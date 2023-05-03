// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>,
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

#ifndef __PERCYSERVER_H__
#define __PERCYSERVER_H__

#include <iostream>
#include "datastore.h"
#include "percyparams.h"
#include "percystats.h"
#include "streams.h"

/**
 * An abstract base class for a PIR server.
 */
class PercyServer {
public:
    /// Destructor
    virtual ~PercyServer();

    /// Factory method called to get a server object for the given parameters
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    /// @return		    New dynamically-allocated PercyServer pointer.
    static PercyServer * make_server (DataStore * datastore, 
	    const PercyServerParams * params, PercyStats * stats = NULL);

    /// Tell the server to be Byzantine
    void be_byzantine ();

    /// Special strassen max level, when set to this value, the optimal
    /// strategy is used
    static const nqueries_t STRASSEN_OPTIMAL = -1;

    /// Set the strassen max level
    /// @param depth	The strassen max level to use.
    void set_strassen_max_depth(nqueries_t depth);

    /// Get the strassen max level
    /// @return	    The strassen max level being used.
    nqueries_t get_strassen_max_depth();

    /// Handle a request.
    /// @param is	Input stream from the client.
    /// @param os	Output stream to the client.
    /// @param workers	When the server is a master server, a list of 
    ///			input/output streams to/from the workers.
    /// @return		Returns false if there is an error, true otherwise.
    virtual bool handle_request (std::istream &is, std::ostream &os,
	    std::vector<std::iostream*> workers = std::vector<std::iostream*>());

protected:
    friend class RecursiveServer;
    friend void * thread_work (void * arg);

    /// Implementation of handle_request().  Must be implemented in derived 
    /// classes.
    /// @param requests	    A set of query requests, one for each query.
    /// @param responses    Buffers that must be filled with the query responses.
    /// @return		    Returns false is there is an error, true otherwise.
    virtual bool handle_request_impl (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses) = 0;

    /// Implementation of handle_request() used for master servers.  Must be
    /// implemented in derived classes.
    /// @param requests	    A set of query requests, one for each query.
    /// @param responses    Buffers that must be filled with the query responses.
    /// @param workers	    When the server is a master server, a list of 
    ///			    input/output streams to/from the workers.
    /// @return		    Returns false is there is an error, true otherwise.
    virtual bool handle_request_distributed (std::vector<unsigned char*> requests,
	    std::vector<unsigned char*> responses, 
	    std::vector<std::iostream*> workers);

    /// Implementation of handle_request() to be used when using multithreading.  
    /// Must be implemented in derived classes.
    /// @param requests	    A set of query requests, one for each query.
    /// @param responses    Buffers that must be filled with the query responses.
    /// @return		    Returns false is there is an error, true otherwise.
    virtual bool handle_request_threaded (std::vector<unsigned char*> requests,
	    std::vector<unsigned char*> responses);

    /// Combines the results of threads/workers when partitioning database
    /// records (DIST_SPLIT_RECORDS).
    /// @param result		Buffer to be filled with combined result.
    /// @param worker_results	The results from all threads/workers to be 
    ///				combined.
    virtual void combine_results (unsigned char * result, 
	    std::vector<unsigned char*> worker_results) = 0;

    /// Constructor.  Can only be called by derived classes.
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    PercyServer (DataStore * datastore, const PercyServerParams * serverparams,
	    PercyStats * stats = NULL);

    /// Whether or not the server is Byzantine.
    bool byzantine;
    /// The database used by the server.
    DataStore * datastore;
    /// The parameters for the server.
    const PercyServerParams * serverparams;
    /// Statistics collection object.
    PercyStats * stats;
    /// Maximum depth allowed when using Strassen's matrix multiplication.
    nqueries_t strassen_max_depth;
    /// The strassen level reached in computation.
    nqueries_t strassen_level_reached;
    /// When using multithreading, the server objects for the threads.
    std::vector<PercyServer*> subservers;

    // Lightweight matrix implementation
    template <typename GF2E_Element> class Matrix;
    template <typename GF2E_Element> class SubMatrix;
    template <typename GF2E_Element> class Col;
    template <typename GF2E_Element> class Row;
    template <typename GF2E_Element> class Elem;

};

std::vector<nqueries_t> split_queries (nqueries_t num_queries, 
	nservers_t num_workers);

#include "gf2e_matrix.h"

#endif
