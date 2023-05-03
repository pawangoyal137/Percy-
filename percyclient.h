// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>,
// Ryan Henry <rhenry@cs.uwaterloo.ca>
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

#ifndef __PERCYCLIENT_H__
#define __PERCYCLIENT_H__

#include <vector>
#include <set>
#include <iostream>
#include <string.h>
#include "percyresult.h"
#include "percyparams.h"
#include "percystats.h"
#include "percyio.h"

/// An abstract base class for a PIR client.
class PercyClient {
public:
    /// Destructor
    virtual ~PercyClient ();

    /// Factory method called to get a client object for the given parameters
    /// @param clientparams Parameters for the client.
    /// @param num_servers  The number of servers used
    /// @param t	    The privacy level.  I.e. the maximum number of 
    ///			    servers that can collude and the queries will remain
    ///			    private.
    /// @param sids	    An array of IDs for the servers that will be used.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///			    collected if NULL.
    /// @return		    New dynamically-allocated PercyClient pointer.
    static PercyClient * make_client (const PercyClientParams * clientparams, 
	    nservers_t num_servers, nservers_t t, sid_t * sids = NULL,
	    PercyStats * stats = NULL);

    /// Encode a request for the given block numbers.
    /// @param block_numbers	The requested block indices (0-based).
    /// @param querybsize	The number of blocks to request in one query.
    /// @return			On success, a positive request ID, zero 
    ///				otherwise.
    nqueries_t encode_request (vector<dbsize_t> block_numbers, 
	    nqueries_t querybsize = 1);

    /// Send the request for the given request ID.
    /// @param request_identifier   Request ID.
    /// @param osvec		    Streams for output to the servers.
    /// @param send_num_queries	    If true, send the servers the number of
    ///				    queries prior to sending the queries.
    /// @return			    The total bytes send (to all servers 
    ///				    combined)
    dbsize_t send_request(nqueries_t request_identifier, 
	    std::vector<ostream*> &osvec, bool send_num_queries = true);

    /// Receive the servers' replies for a given request ID.
    /// @param request_identifier   Request ID.
    /// @param isvec		    Streams for input from the servers.
    /// @return			    The total bytes received (from all servers 
    ///				    combined)
    dbsize_t receive_replies(nqueries_t request_identifier,
	    std::vector<istream*> &isvec);

    /// Process the servers' replies for all undecoded replies.  The
    /// successfully decoded blocks are put in decoded_blocks and the undecoded
    /// block numbers will remain in undecoded_blocks.
    /// @param h    The number of honest servers.
    /// @return	    The number of queries that we not successfully decoded.
    nqueries_t process_replies(nservers_t h);

    // If the request is finished (fully and successfully processed by
    // process_replies), return true and put results in the results vector.
    // The results are then removed from the server (no longer stored in this
    // class).  Otherwise, return false and nothing is done to results.
    /// Get the result for a given request ID.
    /// @param request_identifier   Request ID.
    /// @param results		    The results will be put into this vector.
    /// @return			    Return true if the request's response has
    ///				    been received from the servers and has been
    ///				    successfully decoded, false otherwise.
    bool get_result (nqueries_t request_identifier,
	    vector<PercyBlockResults> & results);

    /// Do encode_request(), send_request(), receive_replies(), 
    /// process_replies(), and get_result() all in one shot.
    /// request_identifier will be set to the request_identifier for this
    /// request for use if not all queries are decoded successfully.  In
    /// this case, process another request and call get_result with this
    /// identifier.
    /// @param request_identifier   Will be set to the ID of this request.
    /// @param block_numbers	The requested block indices (0-based).
    /// @param osvec		    Streams for output to the servers.
    /// @param isvec		    Streams for input from the servers.
    /// @param results		    Will be filled with the results of this
    ///				    request.
    /// @param querybsize	The number of blocks to request in one query.
    /// @return True if all queries were successful, false otherwise.
    bool fetch_blocks(nqueries_t& request_identifier, 
	    vector<dbsize_t> block_numbers, vector<ostream*> &osvec, 
	    vector<istream*> &isvec, vector<PercyBlockResults> &results, 
	    nqueries_t querybsize = 1);

protected:
    friend class RecursiveClient;

    /// Get the block numbers for a given request ID.
    /// @param request_identifier   Request ID.
    /// @return			    The blocks requested.
    const vector<dbsize_t>& get_block_numbers (nqueries_t request_identifier) {
	return requested_blocks[request_identifier];
    }

    /// Get the number of blocks to request in a single query for a given
    /// request ID.
    /// @param request_identifier   Request ID.
    /// @return			    The number of blocks to request in a single
    ///				    query.
    nqueries_t get_qbs (nqueries_t request_identifier) {
	return blocks_per_query[request_identifier];
    }

    /// Implementation for encode_request().
    /// @param request_identifier   Request ID.
    virtual void encode_request_impl (nqueries_t request_identifier) = 0;

    /// Implementation for send_request().
    /// @param request_identifier   Request ID.
    /// @param osvec		    Streams for output to the servers
    /// @param send_num_queries	    If true, send the servers the number of
    ///				    queries prior to sending the queries.
    /// @return			    The total bytes send (to all servers 
    ///				    combined)
    virtual dbsize_t send_request_impl (nqueries_t request_identifier,
	    vector<ostream*> &osvec, bool send_num_queries = true) = 0;

    /// Implementation for encode_request().
    /// @param request_identifier   Request ID.
    /// @param isvec		    Streams for input from the servers.
    /// @return			    The total bytes received (from all servers 
    ///				    combined)
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec) = 0;

    /// Implementation for encode_request().
    /// @param h	The number of honest servers.
    /// @param results	Will be filled with the results for all requests for
    ///			which a reply has been received received.
    /// @return		The number of queries that we not successfully decoded.
    virtual nqueries_t process_replies_impl (nservers_t h, 
	    vector<vector<PercyResult> >& results) = 0;

    /// Constructor.  Can only be called by derived classes.  Use make_client()
    /// to create a client object.
    /// @param clientparams Parameters for the client.
    /// @param num_servers  The number of servers used
    /// @param t	    The privacy level.  I.e. the maximum number of 
    ///			    servers that can collude and the queries will remain
    ///			    private.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///			    collected if NULL.
    PercyClient (const PercyClientParams * clientparams, nservers_t num_servers,
	    nservers_t t, PercyStats * stats = NULL);

    /// Parameters for the client.
    const PercyClientParams * clientparams;
    /// The number of servers.
    nservers_t num_servers;
    /// The privacy level.  I.e. the maximum number of servers that can collude 
    /// and the queries will remain private.
    nservers_t t;
    /// The server indices of servers that have not acted up.
    vector<nservers_t> goodservers;
    /// Randomize the requests for IT-PIR clients.
    static const bool randomize = true;

private:
    typedef map<nqueries_t, vector<dbsize_t> > BlockNumbers;
    typedef map<nqueries_t, vector<PercyBlockResults> > BlockResults;
    enum RequestStatus { ENCODE, SEND, RECEIVE, DECODE, INCOMPLETE, DONE };
    typedef map<nqueries_t, RequestStatus> RequestStatuses;

    nqueries_t next_request_identifier;
    BlockNumbers requested_blocks;
    std::map<nqueries_t, nqueries_t> blocks_per_query;
    RequestStatuses statuses;
    BlockNumbers undecoded_blocks;
    BlockResults decoded_blocks;

    PercyStats * stats;
    map<nqueries_t, nqueries_t> request_id_to_stats_batch_number;
};

inline ostream& operator<<(ostream &out, vector<dbsize_t> v) {
    out << "<";
    for (vector<dbsize_t>::const_iterator i=v.begin();
	    i != v.end(); ++i) {
	out << " " << *i;
    }
    out << " >";
    return out;
}

#endif
