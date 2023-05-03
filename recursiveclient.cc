// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>
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
#include "recursiveclient.h"
#include "streams.h"
#include "percyio.h"

RecursiveClient::RecursiveClient (const RecursiveClientParams * clientparams, 
	nservers_t num_servers, nservers_t t, sid_t * sids, PercyStats * stats)
:
    PercyClient(clientparams, num_servers, t, stats),
    params(static_cast<const RecursiveParams*>(clientparams->percy_params()))
{
    // Make iteration clients
    nqueries_t depth = params->depth();
    for (nqueries_t d = 0; d < depth; ++d) {
	PercyClient * client_d = PercyClient::make_client(
		clientparams->get_iteration(d), num_servers, t, sids);
	iteration_clients.push_back(client_d);
    }
}

RecursiveClient::~RecursiveClient ()
{
    nqueries_t depth = params->depth();
    for (nqueries_t d = 0; d < depth; ++d) {
	if (iteration_clients[d]) {
	    delete iteration_clients[d];
	}
    }
}

void RecursiveClient::encode_request_impl (nqueries_t request_identifier)
{
    const vector<dbsize_t>& block_numbers = get_block_numbers(request_identifier);
    nqueries_t num_queries = block_numbers.size();
    nqueries_t depth = params->depth();
    std::vector<std::vector<nqueries_t> >& this_req_ids =
	    req_ids[request_identifier];
    this_req_ids = std::vector<std::vector<nqueries_t> >(num_queries,
	    std::vector<nqueries_t>(depth, 0));

    for (nqueries_t qind = 0; qind < num_queries; ++qind) {
	dbsize_t block_number = block_numbers[qind];
	std::vector<dbsize_t> iter_indices = params->iteration_indices(block_number);
	for (nqueries_t d = 0; d < depth; ++d) {
	    nqueries_t req_id = iteration_clients[d]->encode_request(
		    std::vector<dbsize_t>(1, iter_indices[d]));
	    this_req_ids[qind][d] = req_id;
	}
    }
}

dbsize_t RecursiveClient::send_request_impl (nqueries_t request_identifier, 
	vector<ostream*> &osvec, bool send_num_queries)
{
    if (osvec.size() != num_servers) {
	std::cerr << "Not enough servers to sent to\n";
	return 0;
    }

    std::vector<std::vector<nqueries_t> >& this_req_ids =
	    req_ids[request_identifier];
    nqueries_t num_queries = this_req_ids.size();
    nqueries_t depth = params->depth();
    dbsize_t retval = 0;

    if (send_num_queries) {
	for (nservers_t i = 0; i < osvec.size(); ++i) {
	    percy_write_le_uint16(*(osvec[i]), num_queries);
	}
    }

    for (nqueries_t qind = 0; qind < num_queries; ++qind) {
	for (nqueries_t d = 0; d < depth; ++d) {
	    retval += iteration_clients[d]->send_request(this_req_ids[qind][d],
		    osvec, false);
	}
    }

    return retval;
}

dbsize_t RecursiveClient::receive_replies_impl (
	nqueries_t request_identifier, vector<istream*> &isvec)
{
    if (isvec.size() != num_servers) {
	std::cerr << "Not enough servers to receive from\n";
	return 0;
    }

    std::vector<std::vector<nqueries_t> >& this_req_ids =
	    req_ids[request_identifier];
    nqueries_t num_queries = this_req_ids.size();
    nqueries_t depth = params->depth();
    dbsize_t retval = 0;

    for (nqueries_t qind = 0; qind < num_queries; ++qind) {
	retval += iteration_clients[depth-1]->receive_replies(
		this_req_ids[qind][depth-1], isvec);
    }
    goodservers = iteration_clients[depth-1]->goodservers;

    unprocessed[request_identifier] = std::vector<nqueries_t>(num_queries, depth);
    return retval;
}

nqueries_t RecursiveClient::process_replies_impl (nservers_t h,
	vector<vector<PercyResult> >& results)
{
    nqueries_t depth = params->depth();
    for (nqueries_t d = depth - 1; d > 0; --d) {
	// Iteration params
	const PercyParams * iter_params = params->get_iteration(d);
	dbsize_t iter_block_size = iter_params->block_size();

	// Process replies for iteration
	iteration_clients[d]->process_replies(h);

	// for each unprocessed request
	std::map<nqueries_t, std::vector<nqueries_t> >::iterator uiter;
	for (uiter = unprocessed.begin(); uiter != unprocessed.end(); ++uiter) {
	    nqueries_t req_id = uiter->first;
	    std::vector<nqueries_t>& this_unprocessed = uiter->second;
	    nqueries_t num_queries = this_unprocessed.size();
	    // for each query in the request
	    for (nqueries_t q = 0; q < num_queries; ++q) {
		if (d < this_unprocessed[q]) {
		    std::vector<PercyBlockResults> results_id_q;
		    if (iteration_clients[d]->get_result(req_ids[req_id][q][d], results_id_q)) {
			// Successful process
			// Next iteration reads
			MemoryStreamBuf * intermediate_bufs = new MemoryStreamBuf[num_servers];
			std::vector<PercyResult>::iterator priter;
			for (priter = results_id_q[0].results.begin(); priter != results_id_q[0].results.end(); ++priter) {
			    std::vector<nservers_t>::iterator Giter;
			    for (Giter = priter->G.begin(); Giter != priter->G.end(); ++Giter) {
				intermediate_bufs[*Giter].add_inbuffer((char*)(priter->sigma.c_str()),
					iter_block_size);
			    }
			}
			std::vector<std::istream*> intermediate_isvec;
			for (nservers_t s = 0; s < num_servers; ++s) {
			    intermediate_isvec.push_back(new std::istream(&(intermediate_bufs[s])));
			}
			iteration_clients[d-1]->receive_replies(req_ids[req_id][q][d-1], intermediate_isvec);
			goodservers = iteration_clients[d-1]->goodservers;
			for (nservers_t s = 0; s < num_servers; ++s) {
			    delete intermediate_isvec[s];
			}
			this_unprocessed[q] = d;
			delete[] intermediate_bufs;
		    }
		}
	    }
	}
    }

    // Last iteration (d = 0)
    iteration_clients[0]->process_replies(h);

    nqueries_t ret = 0;
    // for each unprocessed request
    std::map<nqueries_t, std::vector<nqueries_t> >::iterator uiter;
    for (uiter = unprocessed.begin(); uiter != unprocessed.end(); ++uiter) {
	nqueries_t req_id = uiter->first;
	std::vector<nqueries_t>& this_unprocessed = uiter->second;
	nqueries_t num_queries = this_unprocessed.size();
	// for each query in the request
	nqueries_t u = 0;
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    std::vector<PercyBlockResults> results_id_q;
	    if (0 < this_unprocessed[q-u] &&
		    iteration_clients[0]->get_result(req_ids[req_id][q-u][0], results_id_q)) {
		// Successful process
		// Add to results
		results.push_back(results_id_q[0].results);
		req_ids[req_id].erase(req_ids[req_id].begin() + q - u);
		unprocessed[req_id].erase(unprocessed[req_id].begin() + q - u);
		++u;
	    } else {
		results.push_back(vector<PercyResult>());
		++ret;
	    }
	}
    }

    return ret;
}

