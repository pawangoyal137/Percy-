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

#ifndef __RECURSIVESERVER_H__
#define __RECURSIVESERVER_H__

#include <iostream>
#include "datastore.h"
#include "recursiveparams.h"
#include "percyserver.h"

NTL_CLIENT

/// A PIR server for recursive PIR protocols.

class RecursiveServer : public PercyServer {
public:
    /// Constructor.
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    RecursiveServer (DataStore * datastore, 
	    const RecursiveServerParams * serverparams, PercyStats * stats = NULL);

    /// Destructor.
    virtual ~RecursiveServer ();

    virtual bool handle_request (std::istream &is, std::ostream &os,
	    std::vector<std::iostream*> workers = std::vector<std::iostream*>());

private:
    virtual bool handle_request_impl (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses);

    virtual bool handle_request_distributed (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses, 
	    std::vector<std::iostream*> workers);

    virtual bool handle_request_threaded (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses);

    virtual void combine_results (unsigned char * result, 
	    std::vector<unsigned char*> worker_results);

    const RecursiveServerParams * serverparams;
    const RecursiveParams * params;
    bool dist_first_only;
    std::vector<DataStore*> iteration_datastores;
    std::vector<PercyServer*> iteration_servers;

    bool handle_request_worker (std::istream& is, std::ostream&os);
};

#endif
