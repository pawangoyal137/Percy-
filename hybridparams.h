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

#ifndef __HYBRIDPARAMS_H__
#define __HYBRIDPARAMS_H__

#include "recursiveparams.h"
#include "itparams.h"
#include "agparams.h"

class HybridParams : public RecursiveParams {
public:
    // For Non-ZZ_p IT modes
    HybridParams (dbsize_t num_blocks, dbsize_t block_size,
	    PercyMode it_mode, nqueries_t depth = 0, nservers_t tau = 0,
	    dbsize_t it_num_blocks = 0, dbsize_t ag_N = 50, 
	    dbsize_t c_word_size = 20);

    // For ZZ_p as the IT mode
    HybridParams (dbsize_t num_blocks, dbsize_t block_size,
	    dbsize_t it_word_size, ZZ it_modulus, nqueries_t depth = 0, 
	    nservers_t tau = 0, dbsize_t it_num_blocks = 0, 
	    dbsize_t ag_N = 50, dbsize_t c_word_size = 20);

    virtual ~HybridParams ();

    virtual dbsize_t server_block_size () const { 
	return iterations[0]->server_block_size();
    }

    // Prints the mode-specfic paramaters.  Meant to be overloaded by
    // mode-specific classes
    virtual void print_mode_specific (std::ostream& os) const;

protected:
    PercyMode it_mode;
};

#endif
