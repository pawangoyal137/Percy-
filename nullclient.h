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

#ifndef __NULLCLIENT_H__
#define __NULLCLIENT_H__

#include <NTL/ZZ.h>
#include "percyclient.h"
#include "percyparams.h"

NTL_CLIENT

/// A PIR client that does nothing except send random data to the servers and
/// receive the responses.

class NullClient : public PercyClient {
public:
    /// Constructor.  Can only be called by derived classes.  Use make_client()
    /// to create a client object.
    /// @param params	    Parameters for the client.
    /// @param num_servers  The number of servers used
    NullClient (const PercyClientParams * params, nservers_t num_servers);

    /// Destructor.
    virtual ~NullClient ();

private:
    const PercyParams * params;

    // Virtual members as described in PercyClient class
    virtual void encode_request_impl (nqueries_t request_identifier);
    virtual dbsize_t send_request_impl (nqueries_t request_identifier, 
	    vector<ostream*> &osvec, bool send_num_queries = true);
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec);
    virtual nqueries_t process_replies_impl (nservers_t h,
	    vector<vector<PercyResult> >& results);

    unsigned char * randbuf;
    nqueries_t num_to_process;
};

#endif
