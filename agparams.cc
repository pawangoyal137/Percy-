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

#include <cmath>
#include <stdio.h>
#include "agparams.h"
#include "percytypes.h"
#include "percyio.h"

AGParams::AGParams(dbsize_t num_blocks, dbsize_t block_size, dbsize_t N, 
	dbsize_t word_size, dbsize_t virtual_block_size)
:
    PercyParams(num_blocks, block_size, word_size, MODE_AG, 0, 
	    virtual_block_size),
    _N(N),
    _block_rows((block_size - 1) / (word_size * N / 8) + 1)
{
    // Initialize p and q
    _q = 1;
    _q <<= 2 * _word_size - 1;
    switch (_word_size) {
    case 16:
	// _p = 281474976710597; // 2^48-59
	_p = 1;
	_p <<= 48;
	_p -= 59;
	break;
    case 20:
	// _p = 1152921504606846883; // 2^60-93
	_p = 1;
	_p <<= 60;
	_p -= 93;
	break;
    case 24:
	// p = 4722366482869645213603 = 2^72-93
	_p = 1;
	_p <<= 72;
	_p -= 93;
	break;
    default:
	break;
    }
}

dbsize_t AGParams::request_size (nqueries_t num_queries) const
{
    return (6 * _word_size * _N * _N / 8) * _num_virtual_blocks * num_queries;
}

dbsize_t AGParams::response_size (nqueries_t num_queries) const
{
    return (6 * _word_size * _N / 8) * _block_rows * num_queries * 
	    _virtual_block_size;
}

void AGParams::print_mode_specific (std::ostream& os) const
{
    os << "N=" << _N;
}

void AGParams::write (std::ostream &os) const
{
    // Output parent class
    PercyParams::write(os);

    // Output AG parameters
    PERCY_WRITE_LE_DBSIZE(os, _N);
    PERCY_WRITE_LE_AGELT(os, _p);
    PERCY_WRITE_LE_AGELT(os, _q);

    os.flush();
}

bool AGParams::check_compatible (std::istream &is) const
{
    // Check parent class
    if (!(PercyParams::check_compatible(is))) return false;

    // Check AG parameters
    dbsize_t readN;
    PERCY_READ_LE_DBSIZE(is, readN);
    AG_Element readp, readq;
    PERCY_READ_LE_AGELT(is, readp);
    PERCY_READ_LE_AGELT(is, readq);
    if (readN != _N ||
	readp != _p ||
	readq != _q) return false;

    return true;
}

std::vector<const PercyParams*> AGParams::create_worker_params (
	std::vector<Dimension> worker_dims) const
{
    std::vector<const PercyParams*> worker_params;
    std::vector<Dimension>::iterator it;
    for (it = worker_dims.begin(); it != worker_dims.end(); ++it) {
	worker_params.push_back(new AGParams(it->first, it->second, 
		_N, _word_size, _virtual_block_size));
    }
    return worker_params;
}


nqueries_t optimal_depth (dbsize_t num_blocks, dbsize_t block_size, dbsize_t N)
{
    nqueries_t depth = 1;
    dbbits_t cost = depth * 6 * N * N * pow(num_blocks, (double)1 / depth) + 
	    pow(6, depth) * block_size;
    while (true) {
	dbbits_t nextcost = (depth + 1) * 6 * N * N * pow(num_blocks, (double)1 / (depth + 1)) + 
		pow(6, depth + 1) * block_size;
	if (nextcost > cost) {
	    return depth;
	}
	depth++;
	cost = nextcost;
    }
}

RecursiveAGParams::RecursiveAGParams (dbsize_t num_blocks, dbsize_t block_size,
	nqueries_t depth, dbsize_t N, dbsize_t word_size)
:
    RecursiveParams(num_blocks, block_size, word_size, MODE_RECURSIVE_AG, 
	    std::vector<const PercyParams*>()),
    _N(N)
{
    // If depth not specified, get optimal depth.
    if (depth == 0) {
	depth = optimal_depth(num_blocks, block_size, N);
    }

    // Make iteration params
    dbsize_t depth_root = pow(num_blocks, (double)1 / depth);
    dbsize_t bound = pow(depth_root, depth);
    dbsize_t virtual_block_size = num_blocks;
    dbsize_t iteration_block_size = block_size;
    nqueries_t d;
    for (d = 0; num_blocks > bound; ++d) {
	dbsize_t iteration_num_blocks = virtual_block_size;
	virtual_block_size = (virtual_block_size - 1) / (depth_root + 1) + 1;
	iterations.push_back(new AGParams(iteration_num_blocks,
		iteration_block_size, _N, _word_size, virtual_block_size));
        bound = bound / depth_root * (depth_root + 1);
	iteration_block_size = iterations.back()->response_size() / 
		virtual_block_size;
    }
    for (; d < depth; ++d) {
	dbsize_t iteration_num_blocks = virtual_block_size;
	virtual_block_size = (virtual_block_size - 1) / depth_root + 1;
	iterations.push_back(new AGParams(iteration_num_blocks,
		iteration_block_size, _N, _word_size, virtual_block_size));
	iteration_block_size = iterations.back()->response_size() / 
		virtual_block_size;
    }
}

RecursiveAGParams::~RecursiveAGParams ()
{
    for (nqueries_t d = 0; d < depth(); ++d) {
	delete iterations[d];
    }
}

void RecursiveAGParams::print_mode_specific (std::ostream& os) const
{
    os << "N=" << _N << ";"
	<< "depth=" << depth();
}

