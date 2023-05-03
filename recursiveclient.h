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

#ifndef __HYBRIDCLIENT_H__
#define __HYBRIDCLIENT_H__

#include <set>
#include "percyclient.h"
#include "recursiveparams.h"

NTL_CLIENT

/// A PIR client for recursive PIR protocols.

class RecursiveClient : public PercyClient {
public:
    /// Constructor.
    /// @param clientparams Parameters for the client.
    /// @param num_servers  The number of servers used
    /// @param t	    The privacy level.  I.e. the maximum number of 
    ///			    servers that can collude and the queries will remain
    ///			    private.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///			    collected if NULL.
    RecursiveClient (const RecursiveClientParams * clientparams, 
	    nservers_t num_servers, nservers_t t, sid_t * sids = NULL, 
	    PercyStats * stats = NULL);

    /// Destructor.
    ~RecursiveClient ();

private:
    /// Parameters for the protocol.
    const RecursiveParams * params;

    // Virtual members as described in PercyClient class
    virtual void encode_request_impl (nqueries_t request_identifier);
    virtual dbsize_t send_request_impl (nqueries_t request_identifier, 
	    vector<ostream*> &osvec, bool send_num_queries = true);
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec);
    virtual nqueries_t process_replies_impl (nservers_t h,
	    vector<vector<PercyResult> >& results);

    // Clients used for each iteration.
    std::vector<PercyClient*> iteration_clients;

    // inner request identifiers (by req_id, then query index, then iteration
    // index)
    std::map<nqueries_t, std::vector<std::vector<nqueries_t> > > req_ids;

    // For each request_id and query index, the depth at which the reply was
    // fully processed.
    std::map<nqueries_t, std::vector<nqueries_t> > unprocessed;
};

#endif
