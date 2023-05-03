// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>,
// Ryan Henry <rhenry@cs.uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of version 2 of the GNU General Public License as
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// There is a copy of the GNU General Public License in the COPYING file
// packaged with this plugin; if you cannot find it, write to the Free
// Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "percyclient.h"
#include "nullclient.h"
#include "agclient.h"
#include "itclient.h"
#include "recursiveclient.h"

// Simple protected constructor for PercyClient
PercyClient::PercyClient (const PercyClientParams * clientparams, 
	nservers_t num_servers, nservers_t t, PercyStats * stats) :
    clientparams(clientparams),
    num_servers(num_servers),
    t(t),
    goodservers(),
    next_request_identifier(1),
    requested_blocks(),
    blocks_per_query(),
    statuses(),
    undecoded_blocks(),
    decoded_blocks(),
    stats(stats),
    request_id_to_stats_batch_number()
{}

PercyClient::~PercyClient () {}

// A factory method used to get a PercyClient for the given mode
PercyClient * PercyClient::make_client (const PercyClientParams * clientparams, 
	nservers_t num_servers, nservers_t t, sid_t * sids,
	PercyStats * stats)
{   
    PercyClient * retptr = NULL;
    if (clientparams->is_null()) {
	return new NullClient(clientparams, num_servers);
    }
    switch (clientparams->percy_params()->get_mode()) {
    case MODE_ZZ_P:
	retptr = new PercyClient_ZZ_p(clientparams, num_servers, t, sids, 
		stats);
	break;
    case MODE_CHOR:
	retptr = new PercyClient_Chor(clientparams, num_servers, stats);
	break;
    case MODE_GF28:
	retptr = new PercyClient_GF2E<GF28_Element>(clientparams, num_servers, 
		t, sids, stats);
	break;
    case MODE_GF216:
	retptr = new PercyClient_GF2E<GF216_Element>(clientparams, num_servers,
		t, sids, stats);
	break;
    case MODE_AG:
	retptr = new PercyAGClient(clientparams, stats);
	break;
    case MODE_RECURSIVE_AG:
	retptr = new RecursiveClient(
		dynamic_cast<const RecursiveClientParams*>(clientparams), 1, 0, 
		sids, stats);
	break;
    case MODE_HYBRID:
    /*
	retptr = new PercyHybridClient(
		static_cast<const HybridClientParams*>(params), num_servers, 
		t, sids, stats);
    */
	break;
    default:
	break;
    }
    return retptr;
}

// Encode a request for the given block numbers (0-based).
// Returns a positive request identifier on success and returns 0 on
// failure.
nqueries_t PercyClient::encode_request (vector<dbsize_t> block_numbers, 
	nqueries_t querybsize)
{
    // Add to requested_blocks
    nqueries_t request_identifier = next_request_identifier++;
    requested_blocks[request_identifier] = block_numbers;
    statuses[request_identifier] = ENCODE;

    if (querybsize > 1) {
	nqueries_t max_qbs = num_servers - t - clientparams->percy_params()->tau();
	if (querybsize > max_qbs) {
	    std::cerr << "qbs is too large.  Using the maximum qbs = " << max_qbs << "\n";
	    querybsize = max_qbs;
	}
	switch (clientparams->percy_params()->get_mode()) {
	case MODE_ZZ_P:
	case MODE_GF28:
	case MODE_GF216:
	    blocks_per_query[request_identifier] = querybsize;
	    break;
	default:
	    std::cerr << "Cannot have qbs > 1 for this mode.  Using qbs=1.\n";
	    blocks_per_query[request_identifier] = 1;
	    querybsize = 1;
	    break;
	}
    } else {
	blocks_per_query[request_identifier] = 1;
    }

    // Get encode start time
    if (stats) {
	nqueries_t batch_number = stats->start_query_batch(block_numbers.size(), 
		querybsize);
	request_id_to_stats_batch_number[request_identifier] = batch_number;
	stats->encode_start(batch_number);
    }

    // Run implementation
    encode_request_impl(request_identifier);
    // Get encode finish time
    if (stats) {
	stats->encode_done(request_id_to_stats_batch_number[request_identifier]);
    }

    statuses[request_identifier] = SEND;
    return request_identifier;
}

// Send the request for the given identifier to the
// servers connected with the ostreams in the given vector.
dbsize_t PercyClient::send_request(nqueries_t request_identifier, 
	std::vector<ostream*> &osvec, bool send_num_queries) 
{
    // Check for good identifier
    RequestStatuses::iterator status_iter;
    status_iter = statuses.find(request_identifier);
    if (status_iter == statuses.end()) {
	std::cerr << "Invalid request identifier\n";
	return 0;
    }
    if (status_iter->second != SEND) {
	std::cerr << "Specified request not at send stage\n";
	return 0;
    }

    // Get send start time
    if (stats) {
	stats->client_to_server_start(request_id_to_stats_batch_number[request_identifier]);
    }

    // Run implementation
    dbsize_t bytes_sent = send_request_impl(request_identifier, osvec,
	    send_num_queries);

    // Get send done time
    if (stats) {
	stats->client_to_server_done(request_id_to_stats_batch_number[request_identifier], bytes_sent);
    }

    // Flush streams
    for (nservers_t i = 0; i < osvec.size(); ++i) {
	osvec[i]->flush();
    }

    status_iter->second = RECEIVE;
    return bytes_sent;
}

// Receive the servers' replies for the next request.
// Return the number of bytes received
dbsize_t PercyClient::receive_replies(nqueries_t request_identifier,
	std::vector<istream*> &isvec)
{
    // Check for good identifier
    RequestStatuses::iterator status_iter;
    status_iter = statuses.find(request_identifier);
    if (status_iter == statuses.end()) {
	std::cerr << "Invalid request identifier\n";
	return 0;
    }
    if (status_iter->second != RECEIVE) {
	std::cerr << "Specified request not at receive stage\n";
	return 0;
    }

    // Get receive start time
    if (stats) {
	// Wait until we see a character
	isvec[0]->peek();
	stats->server_to_client_start(request_id_to_stats_batch_number[request_identifier]);
    }

    // Run implementation
    dbsize_t bytes_received = receive_replies_impl(request_identifier, isvec);
    
    // Get receive done time
    if (stats) {
	stats->server_to_client_done(request_id_to_stats_batch_number[request_identifier], bytes_received);
    }

    undecoded_blocks[request_identifier] = requested_blocks[request_identifier];
    decoded_blocks[request_identifier] =
	    vector<PercyBlockResults>(requested_blocks[request_identifier].size());
    status_iter->second = DECODE;
    return bytes_received;
}

// Process the servers' replies for all undecoded replies.  The
// successfully decoded blocks are put in decoded_blocks and the undecoded
// block numbers will remain in undecoded_blocks.
nqueries_t PercyClient::process_replies(nservers_t h)
{
#ifdef VERBOSE_PERCYCLIENT
    std::cerr << "IN PROCESS_REPLIES\n";
#endif
    vector<nqueries_t> request_identifiers_decode;
    RequestStatuses::iterator status_iter;
    for (status_iter = statuses.begin(); status_iter != statuses.end(); ++status_iter) {
	if (status_iter->second == DECODE) {
	    request_identifiers_decode.push_back(status_iter->first);
	}
    }

    // Get decode start time
    if (stats) {
	for (nqueries_t r = 0; r < request_identifiers_decode.size(); ++r) {
	    stats->decode_start(request_id_to_stats_batch_number[request_identifiers_decode[r]]);
	}
    }

    // Run implementation
    vector<vector<PercyResult> > results;
    nqueries_t num_unsuccessful = process_replies_impl(h, results);

    // Process results
    nqueries_t i = 0;
    BlockNumbers::iterator ub_iter;
    for (ub_iter = undecoded_blocks.begin(); ub_iter != undecoded_blocks.end();
	    /* increment within */) {
	nqueries_t req_id = ub_iter->first;
#ifdef VERBOSE_PERCYCLIENT
	std::cerr << "req_id = " << req_id << endl;
#endif
	vector<dbsize_t>& undecoded = ub_iter->second;
#ifdef VERBOSE_PERCYCLIENT
	std::cerr << "undecoded = " << undecoded << endl;
#endif
	vector<dbsize_t>& blocks = requested_blocks[req_id];
#ifdef VERBOSE_PERCYCLIENT
	std::cerr << "blocks = " << blocks << endl;
#endif
	vector<PercyBlockResults>& decoded = decoded_blocks[req_id];
#ifdef VERBOSE_PERCYCLIENT
	std::cerr << "decoded(before) = " << decoded << endl;
#endif
	nqueries_t u = 0;
	for (nqueries_t q = 0; q < blocks.size(); ++q) {
	    if (u < undecoded.size() && undecoded[u] != blocks[q]) {
		continue;
	    }
	    if (results.size() > i && results[i].size() > 0) {
		decoded[q].block_number = blocks[q];
		decoded[q].results = results[i];
		undecoded.erase(undecoded.begin() + u);
	    } else {
		++u;
	    }
	    ++i;
	}
	// Advance the iterator *before* we potentially delete what
	// it's pointing to
	++ub_iter;
	if (undecoded.size() == 0) {
	    undecoded_blocks.erase(req_id);
	    statuses[req_id] = DONE;
	} else {
	    statuses[req_id] = INCOMPLETE;
	}
#ifdef VERBOSE_PERCYCLIENT
	std::cerr << "decoded(after) = " << decoded << endl;
#endif
    }

    // Get decode done time
    if (stats) {
	for (nqueries_t r = 0; r < request_identifiers_decode.size(); ++r) {
	    nqueries_t req_id = request_identifiers_decode[r];
	    if (undecoded_blocks.find(req_id) == undecoded_blocks.end()) {
		stats->decode_done(request_id_to_stats_batch_number[req_id], 0);
	    } else {
		stats->decode_done(request_id_to_stats_batch_number[req_id], undecoded_blocks[req_id].size());
	    }
	    stats->finish_query_batch(request_id_to_stats_batch_number[req_id]);
	}
    }

#ifdef VERBOSE_PERCYCLIENT
    std::cerr << "LEAVING PROCESS_REPLIES\n";
#endif
    return num_unsuccessful;
}

// If the request is finished (fully and successfully processed by
// process_replies), return true.  Otherwise, return false.
// results will be filled with the fully decoded queries.
bool PercyClient::get_result (nqueries_t request_identifier,
	vector<PercyBlockResults>& results)
{
    // Check for good identifier
    RequestStatuses::iterator status_iter;
    status_iter = statuses.find(request_identifier);
    if (status_iter == statuses.end()) {
	std::cerr << "Invalid request identifier\n";
	return false;
    }

    RequestStatus status = status_iter->second;
    if (status == DONE) {
	results = decoded_blocks[request_identifier];
	requested_blocks.erase(request_identifier);
	blocks_per_query.erase(request_identifier);
	statuses.erase(request_identifier);
	undecoded_blocks.erase(request_identifier);
	decoded_blocks.erase(request_identifier);
	request_id_to_stats_batch_number.erase(request_identifier);
	return true;
    } else if (status == INCOMPLETE) {
	vector<PercyBlockResults>& stored_results = decoded_blocks[request_identifier];
	nqueries_t num_queries = stored_results.size();
	results.resize(num_queries);
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    if (!(stored_results[q].results.empty())) {
		results.push_back(stored_results[q]);
	    }
	}
	return false;
    } else {
	std::cerr << "Server replies have not been processed\n";
	return false;
    }
    return true;
}

// Do all of the above in one shot.
// results will be set to be the results of decoding the blocks
//     requested with the current call.
bool PercyClient::fetch_blocks(nqueries_t& request_identifier,
	vector<dbsize_t> block_numbers, vector<ostream*> &osvec, 
	vector<istream*> &isvec, vector<PercyBlockResults> &results, nqueries_t querybsize)
{
    request_identifier = encode_request(block_numbers, querybsize);
    send_request(request_identifier, osvec);
    receive_replies(request_identifier, isvec);
    nservers_t k = goodservers.size();
    if(k == 0) {
	std::cerr << "ERROR: No good replies received, returning" << endl;
	return false;
    }
    // Calculate what h to use
    nservers_t h = (nservers_t)(floor(sqrt((t+clientparams->percy_params()->tau()+querybsize-1)*k)))+1;
    const char *envh = getenv("PIRC_H");
    if (envh) {
        nservers_t override_h = atoi(envh);
        if (override_h > 0) {
            h = override_h;
        }
    }
    process_replies(h);
    return get_result(request_identifier, results);
}

