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

#ifndef __AGSERVER_H__
#define __AGSERVER_H__

#include <iostream>
#include "datastore.h"
#include "agparams.h"
#include "percyserver.h"

NTL_CLIENT

/// A PIR server for the CPIR protocol by Aguilar Melchor and Gaborit (2007).  
/// This protocol was introduced in 
/// <a href="https://eprint.iacr.org/2007/446">A Lattice-Based Computationally-Efficient Private Information Retrieval Protocol</a> 
/// and revisited in 
/// <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?tp=&arnumber=4622593">HighSpeed Private Information Retrieval Computation on GPU</a>.

class PercyAGServer : public PercyServer {
public:
    /// Constructor.
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    PercyAGServer (DataStore * datastore, const PercyServerParams * params,
	    PercyStats * stats = NULL);
    /// Destructor.
    virtual ~PercyAGServer ();

private:
    virtual bool handle_request_impl (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses);

    virtual void combine_results (unsigned char * result, 
	    std::vector<unsigned char*> worker_results);

    const AGParams * params;
    const AG_Element p;

    bool handle_request_16 (unsigned char * request, unsigned char * response,
	    const unsigned char * database);
    bool handle_request_20 (unsigned char * request, unsigned char * response,
	    const unsigned char * database);
};

#endif
