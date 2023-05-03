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
#include "hybridparams.h"
#include "percyio.h"

// Choose the optimal value for depth
// If it_num_blocks == 0, assumes optimal value
nqueries_t optimal_depth (dbsize_t num_blocks, dbsize_t block_size,
	dbsize_t it_word_size, dbsize_t ag_N, dbsize_t ag_word_size,
	dbsize_t it_num_blocks = 0)
{
    dbbits_t cost = num_blocks * it_word_size + block_size;
    nqueries_t depth = 2;
    while (true) {
	dbbits_t nextcost;
	if (it_num_blocks == 0) {
	    nextcost = depth * pow(6 * ag_N * ag_N * ag_word_size, 
		    (double)(depth - 1) / depth) *
		    pow(it_word_size * num_blocks, (double)1 / depth)
		    + pow(6, depth - 1) * block_size;
	} else {
	    nextcost = it_num_blocks * it_word_size + 
		    6 * (depth - 1) * ag_N * ag_N * ag_word_size * 
		    pow(num_blocks / it_num_blocks, (double)1 / (depth - 1)) + 
		    pow(6, depth - 1) * block_size;
	}
	if (nextcost > cost) {
	    break;
	}
	cost = nextcost;
	depth += 1;
    }
    return depth - 1;
}

// Choose the optimal value for it_num_blocks
dbsize_t optimal_it_num_blocks (dbsize_t num_blocks, dbsize_t it_word_size,
	nqueries_t depth, dbsize_t ag_N, dbsize_t ag_word_size)
{
    if (depth > 1) {
	return (dbsize_t)pow(6 * ag_word_size * ag_N * ag_N * 
		pow(num_blocks, (double)1 / (depth - 1)) / it_word_size,
		(double)(depth - 1) / depth);
    } else {
	return num_blocks;
    }
}

// For Non-ZZ_p IT modes
HybridParams::HybridParams (dbsize_t num_blocks, dbsize_t block_size,
	PercyMode it_mode, nqueries_t depth, nservers_t tau, 
	dbsize_t it_num_blocks, dbsize_t ag_N, dbsize_t c_word_size)
:
    RecursiveParams(num_blocks, block_size, 8, MODE_HYBRID,
	    std::vector<const PercyParams*>(), tau),
    it_mode(it_mode)
{
    dbsize_t it_word_size = 8;
    switch (it_mode) {
    case MODE_CHOR:
	it_word_size = 1;
	break;
    case MODE_ZZ_P:
    case MODE_GF28:
	it_word_size = 8;
	break;
    case MODE_GF216:
	it_word_size = 16;
	break;
    default:
	std::cerr << "Invalid iteration 1 mode: " << it_mode << "\n";
	break;
    }
    if (depth == 0) {
	depth = optimal_depth(num_blocks, block_size, it_word_size, ag_N,
		c_word_size, it_num_blocks);
    }
    if (it_num_blocks == 0) {
	it_num_blocks = optimal_it_num_blocks(num_blocks, it_word_size, 
		depth, ag_N, c_word_size);
    }

    // First iteration params
    dbsize_t virtual_block_size = (_num_blocks - 1) / it_num_blocks + 1;
    switch (it_mode) {
    case MODE_CHOR:
	iterations.push_back(new ChorParams(num_blocks, block_size, 
		virtual_block_size));
	break;
    case MODE_ZZ_P:
	iterations.push_back(new ZZ_pParams(num_blocks, block_size, 
		it_word_size, to_ZZ("257"), tau, NULL, false, 
		virtual_block_size));
	break;
    case MODE_GF28:
    case MODE_GF216:
	iterations.push_back(new GF2EParams(num_blocks, block_size, 
		it_word_size, tau, virtual_block_size));
	break;
    default:
	std::cerr << "Invalid iteration 1 mode: " << it_mode << "\n";
	break;
    }

    // Rest of iterations
    if (depth > 1) {
	dbsize_t vbs_0 = virtual_block_size;
	dbsize_t depth_root = pow(virtual_block_size, (double)1 / (depth-1));
	dbsize_t bound = pow(depth_root, depth-1);
	dbsize_t iteration_block_size = block_size;
	nqueries_t d;
	for (d = 1; vbs_0 > bound; ++d) {
	    dbsize_t iteration_num_blocks = virtual_block_size;
	    virtual_block_size = (virtual_block_size - 1) / (depth_root + 1) 
		    + 1;
	    iterations.push_back(new AGParams(iteration_num_blocks,
		    iteration_block_size, ag_N, c_word_size,
		    virtual_block_size));
	    bound = bound / depth_root * (depth_root + 1);
	    iteration_block_size = iterations.back()->response_size () /
		    virtual_block_size;
	}
	for (; d < depth; ++d) {
	    dbsize_t iteration_num_blocks = virtual_block_size;
	    virtual_block_size = (virtual_block_size - 1) / depth_root + 1;
	    iterations.push_back(new AGParams(iteration_num_blocks,
		    iteration_block_size, ag_N, c_word_size, 
		    virtual_block_size));
	    iteration_block_size = iterations.back()->response_size () /
		    virtual_block_size;
	}
    }
}

// For ZZ_p as the IT mode
HybridParams::HybridParams (dbsize_t num_blocks, dbsize_t block_size,
	dbsize_t it_word_size, ZZ it_modulus, nqueries_t depth, nservers_t tau, 
	dbsize_t it_num_blocks, dbsize_t ag_N, dbsize_t c_word_size)
:
    RecursiveParams(num_blocks, block_size, 8, MODE_HYBRID,
	    std::vector<const PercyParams*>(), tau),
    it_mode(MODE_ZZ_P)
{
    if (depth == 0) {
	depth = optimal_depth(num_blocks, block_size, it_word_size, ag_N,
		c_word_size, it_num_blocks);
    }
    if (it_num_blocks == 0) {
	it_num_blocks = optimal_it_num_blocks(num_blocks, it_word_size, 
		depth, ag_N, c_word_size);
    }

    // First iteration params
    dbsize_t virtual_block_size = (num_blocks - 1) / it_num_blocks + 1;
    ZZ_pParams * iter1_params = new ZZ_pParams(num_blocks, block_size, 
	    it_word_size, it_modulus, tau, NULL, false, virtual_block_size);
    iterations.push_back(iter1_params);

    // Rest of iterations
    if (depth > 1) {
	dbsize_t vbs_0 = virtual_block_size;
	dbsize_t depth_root = pow(virtual_block_size, (double)1 / (depth-1));
	dbsize_t bound = pow(depth_root, depth-1);
	dbsize_t iteration_block_size = iter1_params->response_size() / virtual_block_size;
	std::cerr << "iteration_block_size = " << iteration_block_size << "\n";
	nqueries_t d;
	for (d = 1; vbs_0 > bound; ++d) {
	    dbsize_t iteration_num_blocks = virtual_block_size;
	    virtual_block_size = (virtual_block_size - 1) / (depth_root + 1) + 1;
	    iterations.push_back(new AGParams(iteration_num_blocks,
		    iteration_block_size, ag_N, c_word_size, 
		    virtual_block_size));
	    bound = bound / depth_root * (depth_root + 1);
	    iteration_block_size = iterations.back()->response_size () /
		    virtual_block_size;
	}
	for (; d < depth; ++d) {
	    dbsize_t iteration_num_blocks = virtual_block_size;
	    virtual_block_size = (virtual_block_size - 1) / depth_root + 1;
	    iterations.push_back(new AGParams(iteration_num_blocks,
		    iteration_block_size, ag_N, c_word_size, 
		    virtual_block_size));
	    iteration_block_size = iterations.back()->response_size () /
		    virtual_block_size;
	}
    }
}

HybridParams::~HybridParams ()
{
    for (nqueries_t d = 0; d < depth(); ++d) {
	delete iterations[d];
    }
}

void HybridParams::print_mode_specific (std::ostream& os) const
{
    os << "depth=" << depth() << ";"
	<< "it_mode=" << it_mode;
}

