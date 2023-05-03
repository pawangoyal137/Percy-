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

#include <stdio.h>
#include "recursiveparams.h"
#include "percyio.h"

RecursiveParams::RecursiveParams(dbsize_t num_blocks, dbsize_t block_size,
	dbsize_t word_size, PercyMode mode, 
	std::vector<const PercyParams*> iterations, nservers_t tau)
:
    PercyParams(num_blocks, block_size, word_size, mode, tau),
    iterations(iterations)
{}

std::vector<dbsize_t> RecursiveParams::iteration_indices (dbsize_t index) const
{
    std::vector<dbsize_t> result;
    for (nqueries_t d = 0; d < depth(); ++d) {
	dbsize_t vbs_d = iterations[d]->virtual_block_size();
	result.push_back(index);
	index %= vbs_d;
    }
    return result;
}

dbsize_t RecursiveParams::request_size (nqueries_t num_queries) const
{
    dbsize_t total = 0;
    for (nqueries_t d = 0; d < depth(); ++d) {
	total += iterations[d]->request_size(num_queries);
    }
    return total;
}

dbsize_t RecursiveParams::response_size (nqueries_t num_queries) const
{
    return iterations.back()->response_size(num_queries);
}

void RecursiveParams::write (std::ostream &os) const
{
    // Output parent class
    PercyParams::write(os);

    // Output iteration params
    for (nqueries_t d = 0; d < depth(); ++d) {
	if (iterations[d] == NULL) continue;
	percy_write_le_uint16(os, d);
	iterations[d]->write(os);
    }

    os.flush();
}

bool RecursiveParams::check_compatible (std::istream &is) const
{
    // Check parent class
    if (!(PercyParams::check_compatible(is))) return false;

    // Check AG parameters
    for (nqueries_t d = 0; d < depth(); ++d) {
	if (iterations[d] == NULL) continue;
	nqueries_t other_d;
	percy_read_le_uint16(is, other_d);
	if (is.eof() || other_d != d) return false;
	if (!(iterations[d]->check_compatible(is))) return false;
    }

    return true;
}

std::vector<const PercyParams*> RecursiveParams::create_worker_params (
	std::vector<Dimension> worker_dims) const
{
    return std::vector<const PercyParams*>(worker_dims.size(), NULL);
}


RecursiveClientParams::RecursiveClientParams (const RecursiveParams * rparams,
	nservers_t num_servers, bool is_null)
:
    PercyClientParams(rparams, num_servers, is_null),
    rparams(rparams)
{
    for (nqueries_t d = 0; d < rparams->depth(); ++d) {
	iterations.push_back(new PercyClientParams(rparams->get_iteration(d),
		num_servers));
    }
}

RecursiveClientParams::~RecursiveClientParams ()
{
    for (nqueries_t d = 0; d < depth(); ++d) {
	delete iterations[d];
    }
}


RecursiveServerParams::RecursiveServerParams (const RecursiveParams * rparams,
	nservers_t sid, bool is_worker, bool be_byzantine)
:
    PercyServerParams(rparams, sid, be_byzantine),
    rparams(rparams),
    iterations(),
    first_only(false),
    _is_worker(is_worker),
    worker_rparams()
{
    for (nqueries_t d = 0; d < rparams->depth(); ++d) {
	const PercyParams * iter_params = rparams->get_iteration(d);
	if (iter_params) {
	    iterations.push_back(new PercyServerParams(iter_params,
		    sid, be_byzantine));
	} else {
	    iterations.push_back(NULL);
	}
    }
}

RecursiveServerParams::RecursiveServerParams (const RecursiveParams * rparams, 
	nservers_t sid, nservers_t num_threads, DistSplit tsplit, 
	nservers_t num_workers, DistSplit wsplit, 
	std::vector<nservers_t> worker_sids, bool fork, bool first_only, 
	bool be_byzantine, bool is_worker)
:
    PercyServerParams(rparams, sid, be_byzantine),
    rparams(rparams),
    iterations(),
    first_only(first_only),
    _is_worker(is_worker)
{
    // Set distributed values
    _num_threads = num_threads;
    this->tsplit = tsplit;
    is_forked = fork;
    _num_workers = num_workers;
    this->wsplit = wsplit;
    
    // Create iterations
    for (nqueries_t d = 0; d < rparams->depth(); ++d) {
	const PercyParams * iter_params = rparams->get_iteration(d);
	if (iter_params == NULL) {
	    iterations.push_back(NULL);
#ifdef VERBOSE_RECURSIVE
	    std::cerr << "Iteration " << d << ": NONE\n";
#endif
	    continue;
	} else if (first_only && d > 0) {
	    iterations.push_back(new PercyServerParams(iter_params, sid, 
		    num_threads, tsplit, 0, wsplit, worker_sids, fork, 
		    be_byzantine));
	} else {
	    iterations.push_back(new PercyServerParams(iter_params, sid, 
		    num_threads, tsplit, num_workers, wsplit, worker_sids, fork, 
		    be_byzantine));
	}
#ifdef VERBOSE_RECURSIVE
	std::cerr << "Iteration " << d << ": ";
	iterations.back()->print(std::cerr);
	std::cerr << "\n";
#endif
    }

    // Set up worker params
    if (num_workers > 0) {
	for (nqueries_t i = 0; i < _num_workers; ++i) {
	    const PercyParams * worker_iter_0 = iterations[0]->get_worker_params(i);
	    if (!worker_iter_0) break;
	    if (first_only) {
		worker_params.push_back(worker_iter_0);
		worker_serverparams.push_back(iterations[0]->get_worker_serverparams(i));
	    } else {
		std::vector<const PercyParams*> worker_iterations;
		if (depth() == 0) break;
		worker_iterations.push_back(worker_iter_0);
		dbsize_t worker_num_blocks = worker_iter_0->num_blocks();
		dbsize_t worker_block_size = worker_iter_0->block_size();
		for (nqueries_t d = 1; d < depth(); ++d) {
		    if (iterations[d]) {
			worker_iterations.push_back(iterations[d]->get_worker_params(i));
		    } else {
			worker_iterations.push_back(NULL);
		    }
		}
		const RecursiveParams * wrparams = new RecursiveParams(worker_num_blocks,
			worker_block_size, rparams->word_size(), rparams->get_mode(),
			worker_iterations, rparams->tau());
		worker_params.push_back(wrparams);
		worker_rparams.push_back(wrparams);
		worker_serverparams.push_back(new RecursiveServerParams(wrparams,
			(i < worker_sids.size() ? worker_sids[i] : sid),
			_num_threads, tsplit, 0, wsplit, std::vector<nservers_t>(),
			fork, false, be_byzantine, (wsplit == DIST_SPLIT_RECORDS)));
	    }
	}
    } else {
	for (nqueries_t i = 0; i < _num_threads; ++i) {
	    std::vector<const PercyParams*> thread_iterations;
	    const PercyParams * thread_iter_0 = iterations[0]->get_thread_params(i);
	    if (!thread_iter_0) break;
	    thread_iterations.push_back(thread_iter_0);
	    dbsize_t thread_num_blocks = thread_iter_0->num_blocks();
	    dbsize_t thread_block_size = thread_iter_0->block_size();
	    for (nqueries_t d = 1; d < depth(); ++d) {
		if (iterations[d]) {
		    thread_iterations.push_back(iterations[d]->get_thread_params(i));
		} else {
		    thread_iterations.push_back(NULL);
		}
	    }
	    const RecursiveParams * trparams = new RecursiveParams(thread_num_blocks,
		    thread_block_size, rparams->word_size(), rparams->get_mode(),
		    thread_iterations, rparams->tau());
	    thread_params.push_back(trparams);
	    thread_serverparams.push_back(new RecursiveServerParams(trparams,
		    sid, false, be_byzantine));
	}
    }
}

RecursiveServerParams::~RecursiveServerParams ()
{
    for (nqueries_t d = 0; d < depth(); ++d) {
	if (iterations[d]) delete iterations[d];
    }
    if (first_only) {
	worker_params.clear();
	worker_serverparams.clear();
    }
}

void RecursiveServerParams::print_distributed (std::ostream& os) const
{
    if (_is_worker) {
	os << "WORKER:";
	PercyServerParams::print_distributed(os);
    } else if (first_only && _num_workers > 0) {
	os << "first:(DISTRIBUTED;" << _num_workers << ";" << wsplit << ");";
	os << "rest:";
	if (_num_threads > 0) {
	    os << "(" << (is_forked ? "FORKED" : "THREADED") << ";" <<
		    _num_threads << ";" << tsplit << ")";
	} else {
	    os << "None";
	}
    } else {
	PercyServerParams::print_distributed(os);
    }
}

