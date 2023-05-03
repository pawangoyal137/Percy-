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

#ifndef __AGPARAMS_H__
#define __AGPARAMS_H__

#include "percyparams.h"
#include "recursiveparams.h"

#ifdef NEED_UINT128
#include <NTL/ZZ.h>

NTL_CLIENT

typedef ZZ AG_Element;

#define AGELT_TO_ULONG(x) (to_ulong(x))
#define AGELT_TO_UINT64(x) ( \
    (uint64_t)(to_ulong(x)&0xFFFFFFFFUL) + \
    ((uint64_t)(to_ulong(x>>32)&0xFFFFFFFFUL)<<32) )
inline AG_Element UINT64_TO_AGELT(uint64_t x) {
    unsigned long lo = ((unsigned long)x)&0xFFFFFFFFUL;
    unsigned long hi = (x>>32)&0xFFFFFFFFUL;
    AG_Element r = to_ZZ(hi);
    AG_Element loz = to_ZZ(lo);
    r <<= 32;
    r += loz;
    return r;
}
inline void PERCY_WRITE_LE_AGELT(std::ostream &os, const AG_Element &elt) {
    unsigned char buf[16];
    BytesFromZZ(buf, elt, 16);
    os.write((const char *)buf, 16);
}
inline void PERCY_READ_LE_AGELT(std::istream &is, AG_Element &elt) {
    unsigned char buf[16];
    is.read((char *)buf, 16);
    ZZFromBytes(elt, buf, 16);
}
#define ULONG_TO_AGELT(x) (to_ZZ(x))

#else
typedef __uint128_t AG_Element;

#define AGELT_TO_ULONG(x) ((unsigned long)(x))
#define AGELT_TO_UINT64(x) ((uint64_t)(x))
#define UINT64_TO_AGELT(x) ((AG_Element)(x))
#define ULONG_TO_AGELT(x) ((AG_Element)(x))
inline void swap(AG_Element &a, AG_Element &b) {
    AG_Element tmp = a;
    a = b;
    b = tmp;
}
#define PERCY_WRITE_LE_AGELT(os, elt) (os).write((char*)(&(elt)), 16)
#define PERCY_READ_LE_AGELT(is, elt) (is).read((char*)(&(elt)), 16)
#endif

class AGParams : public PercyParams {
public:
    AGParams(dbsize_t num_blocks, dbsize_t block_size, dbsize_t N = 50, 
	    dbsize_t word_size = 20, dbsize_t virtual_block_size = 1);

    virtual ~AGParams () {}

    // Accessors
    dbsize_t N () const { return _N; }
    AG_Element p () const { return _p; }
    AG_Element q () const { return _q; }
    dbsize_t block_rows () const { return _block_rows; }

/*
    // Helpful methods
    std::vector<dbsize_t> iteration_num_subsets () const { return iter_subsets; }
    dbsize_t iteration_num_subsets (nqueries_t i) const { return iter_subsets[i]; }
    dbsize_t iteration_subset_size (nqueries_t i) const { return iter_subset_sizes[i]; }
*/

    // Return the size of the request/response
    virtual dbsize_t request_size (nqueries_t num_queries = 1) const;
    virtual dbsize_t response_size (nqueries_t num_queries = 1) const;

    // Prints the mode-specfic paramaters.  Meant to be overloaded by
    // mode-specific classes
    virtual void print_mode_specific (std::ostream& os) const;

    // For use in distributed computation
    virtual std::vector<const PercyParams*> create_worker_params (
	    std::vector<Dimension> worker_dims) const;

protected:
    // Write the parameters to a stream to check compatibility
    virtual void write (std::ostream &os) const;

    // Read the parameters from a stream (as written by write()) and check that
    // they are compatible with these parameters.
    virtual bool check_compatible (std::istream &is) const;

    dbsize_t _N;
    AG_Element _p;
    AG_Element _q;
    dbsize_t _block_rows;
};


class RecursiveAGParams : public RecursiveParams {
public:
    RecursiveAGParams (dbsize_t num_blocks, dbsize_t block_size, 
	    nqueries_t depth = 0, dbsize_t N = 50, dbsize_t word_size = 20);

    virtual ~RecursiveAGParams ();

    // Accessors
    dbsize_t N () const { return _N; }

    // Prints the mode-specfic paramaters.  Meant to be overloaded by
    // mode-specific classes
    virtual void print_mode_specific (std::ostream& os) const;

protected:
    dbsize_t _N;

private:
    void init_iterations (std::vector<dbsize_t> iteration_num_blocks);
};

#endif
