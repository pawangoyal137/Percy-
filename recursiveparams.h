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

#ifndef __RECURSIVEPARAMS_H__
#define __RECURSIVEPARAMS_H__

#include "percyparams.h"

class RecursiveParams : public PercyParams {
public:
    RecursiveParams(dbsize_t num_blocks, dbsize_t block_size, 
	    dbsize_t word_size, PercyMode mode,
	    std::vector<const PercyParams*> iterations, nservers_t tau = 0);

    virtual ~RecursiveParams () {}

    // Find the indices needed at each iteration for a given record index
    std::vector<dbsize_t> iteration_indices (dbsize_t index) const;

    // Accessors
    std::vector<const PercyParams*> get_iterations () const { 
	return iterations;
    }
    const PercyParams * get_iteration (nqueries_t index) const {
	return (index < iterations.size() ? iterations[index] : NULL);
    }
    nqueries_t depth () const {
	return iterations.size();
    }

    // Return the size of the request/response
    dbsize_t request_size (nqueries_t num_queries = 1) const;
    dbsize_t response_size (nqueries_t num_queries = 1) const;

    // Prints the mode-specfic paramaters.  Meant to be overloaded by
    // mode-specific classes
    virtual void print_mode_specific (std::ostream& os, char sep = ',') const {}

    // For use in distributed computation.
    virtual std::vector<const PercyParams*> create_worker_params (
	    std::vector<Dimension> worker_dims) const;

    virtual bool is_recursive () const { return true; }

protected:
    // Write the parameters to a stream to check compatibility
    virtual void write (std::ostream &os) const;

    // Read the parameters from a stream (as written by write()) and check that
    // they are compatible with these parameters.
    virtual bool check_compatible (std::istream &is) const;

    // Members
    std::vector<const PercyParams*> iterations;
};


class RecursiveClientParams : public PercyClientParams {
public:
    RecursiveClientParams (const RecursiveParams * rparams, 
	    nservers_t num_servers = 1, bool is_null = false);

    virtual ~RecursiveClientParams ();

    // Accessors
    const RecursiveParams * recursive_params () const { return rparams; }

    std::vector<const PercyClientParams*> get_iterations () const { 
	return iterations;
    }
    const PercyClientParams * get_iteration (nqueries_t index) const {
	return (index < iterations.size() ? iterations[index] : NULL);
    }
    nqueries_t depth () const {
	return iterations.size();
    }

protected:
    // Recursive params
    const RecursiveParams * rparams;

    // Members
    std::vector<const PercyClientParams*> iterations;
};


class RecursiveServerParams : public PercyServerParams {
public:
    RecursiveServerParams (const RecursiveParams * rparams, nservers_t sid,
	    bool is_worker = false, bool be_byzantine = false);

    // Distributed computation constructor
/*
    RecursiveServerParams (const RecursiveParams * rparams, nservers_t sid,
	    nservers_t num_workers, DistType dtype = DIST_TYPE_NONE,
	    DistSplit dsplit = DIST_SPLIT_RECORDS,
	    std::vector<nservers_t> worker_sids = std::vector<nservers_t>(),
	    bool first_only = false, bool be_byzantine = false);
*/

    RecursiveServerParams (const RecursiveParams * rparams, nservers_t sid,
	    nservers_t num_threads = 0, DistSplit tsplit = DIST_SPLIT_RECORDS,
	    nservers_t num_workers = 0, DistSplit wsplit = DIST_SPLIT_RECORDS,
	    std::vector<nservers_t> worker_sids = std::vector<nservers_t>(),
	    bool fork = false, bool first_only = false,
	    bool be_byzantine = false, bool is_worker = false);

    virtual ~RecursiveServerParams ();

    // Accessors
    const RecursiveParams * recursive_params () const { return rparams; }
    bool dist_first_only () const { return first_only; }
    bool is_worker () const { return _is_worker; }

    std::vector<const PercyServerParams*> get_iterations () const { 
	return iterations;
    }
    const PercyServerParams * get_iteration (nqueries_t index) const {
	return (index < iterations.size() ? iterations[index] : NULL);
    }
    nqueries_t depth () const {
	return iterations.size();
    }

    std::vector<const RecursiveParams*> get_worker_rparams () const {
	return worker_rparams;
    }

    virtual void print_distributed (std::ostream& os) const;

protected:
    // Recursive params
    const RecursiveParams * rparams;

    // Members
    std::vector<const PercyServerParams*> iterations;
    bool first_only;
    bool _is_worker;
    std::vector<const RecursiveParams*> worker_rparams;
};

#endif
