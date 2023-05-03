// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>
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

#ifndef __ITSERVER_IMPL_H__
#define __ITSERVER_IMPL_H__

#include <iostream>
#include <string.h>
#include "percyparams.h"
#include "gf2e.h"
#include "xor.h"
#include "gf2e_matrix.h"

////// GF2E Server //////

template <typename GF2E_Element>
PercyServer_GF2E<GF2E_Element>::PercyServer_GF2E (DataStore * datastore, 
	const PercyServerParams * params, PercyStats * stats)
:
    PercyServer(datastore, params, stats),
    params(static_cast<const GF2EParams*>(params->percy_params()))
{}

template <typename GF2E_Element>
PercyServer_GF2E<GF2E_Element>::~PercyServer_GF2E ()
{}

template <typename GF2E_Element>
inline void compute_outputvec_single(
        const GF2E_Element *data, const GF2E_Element *inputvec, 
	GF2E_Element *outputvec, dbsize_t num_blocks, dbsize_t words_per_block,
	dbsize_t virtual_block_size = 1);

template <>
inline void compute_outputvec_single<GF28_Element>(
        const GF28_Element *data, const GF28_Element *inputvec, 
	GF28_Element *outputvec, dbsize_t num_blocks, dbsize_t words_per_block,
	dbsize_t virtual_block_size)
{
    const GF28_Element *block = data;
    dbsize_t num_virtual_blocks = (num_blocks - 1) / virtual_block_size + 1;
    dbsize_t words_per_virtual_block = words_per_block * virtual_block_size;
    for (unsigned int j = 0; j < num_virtual_blocks; ++j) {
        const GF28_Element *multrow = GF28_mult_table[inputvec[j]];
	const GF28_Element *blockc = block;
	if (j == num_virtual_blocks-1) {
	    words_per_virtual_block = words_per_block *
		    ((num_blocks - 1) % virtual_block_size + 1);
	}
	GF28_Element *oc = outputvec;
	GF28_Element *oc_end = oc + (words_per_virtual_block & ~7);
	while (oc < oc_end) {
	    uint64_t accum = (uint64_t) multrow[*(blockc++)];
	    accum |= (uint64_t) multrow[*(blockc++)] << 8;
	    accum |= (uint64_t) multrow[*(blockc++)] << 16;
	    accum |= (uint64_t) multrow[*(blockc++)] << 24;
	    accum |= (uint64_t) multrow[*(blockc++)] << 32;
	    accum |= (uint64_t) multrow[*(blockc++)] << 40;
	    accum |= (uint64_t) multrow[*(blockc++)] << 48;
	    accum |= (uint64_t) multrow[*(blockc++)] << 56;
	    *((uint64_t *) oc) ^= accum;
	    oc+=8;
	}
	for (unsigned int c = 0; c < (words_per_virtual_block & 7); ++c) {
	    *(oc++) ^= multrow[*(blockc++)];
	}
	block += words_per_virtual_block;
    }
}

template <>
inline void compute_outputvec_single<GF216_Element>(
        const GF216_Element *data, const GF216_Element *inputvec, 
	GF216_Element *outputvec, dbsize_t num_blocks, dbsize_t words_per_block,
	dbsize_t virtual_block_size)
{
    const GF216_Element *block = data;
    dbsize_t num_virtual_blocks = (num_blocks - 1) / virtual_block_size + 1;
    dbsize_t words_per_virtual_block = words_per_block * virtual_block_size;
    for (dbsize_t j = 0; j < num_virtual_blocks; ++j) {
        GF216_Element inpv_j = inputvec[j];
	if (j == num_virtual_blocks-1) {
	    words_per_virtual_block = words_per_block *
		    ((num_blocks - 1) % virtual_block_size + 1);
	}
        if (inpv_j == 0) {
	    block += words_per_virtual_block;
	    continue;
	}
	GF216_Element log_j = GF216_log_table[inpv_j];
	const GF216_Element *start = GF216_exp_table + log_j;
	const GF216_Element *blockc = block;
	GF216_Element *oc = outputvec;
	GF216_Element *oc_end = oc + (words_per_virtual_block & ~3);
	GF216_Element block_c;
	while(oc < oc_end) {
	    uint64_t accum = 0;
	    block_c = *(blockc++);
	    if (block_c != 0) {
		GF216_Element log_c = GF216_log_table[block_c];
		accum |= (uint64_t) start[log_c];
	    }
	    block_c = *(blockc++);
	    if (block_c != 0) {
		GF216_Element log_c = GF216_log_table[block_c];
		accum |= (uint64_t) start[log_c] << 16;
	    }
	    block_c = *(blockc++);
	    if (block_c != 0) {
		GF216_Element log_c = GF216_log_table[block_c];
		accum |= (uint64_t) start[log_c] << 32;
	    }
	    block_c = *(blockc++);
	    if (block_c != 0) {
		GF216_Element log_c = GF216_log_table[block_c];
		accum |= (uint64_t) start[log_c] << 48;
	    }
	    *((uint64_t *) oc) ^= accum;
	    oc+=4;
	}
	for (dbsize_t c = 0; c < (words_per_virtual_block & 3); ++c, ++oc) {
	    block_c = *(blockc++);
	    if (block_c != 0) {
		GF216_Element log_c = GF216_log_table[block_c];
		*oc ^= start[log_c];
	    }
	}
	block += words_per_virtual_block;
    }
}

template <typename GF2E_Element>
inline void compute_outputvec_multi(
        const GF2E_Element *data, const std::vector<GF2E_Element*> inputvecs, 
	std::vector<GF2E_Element*> outputvecs, nqueries_t num_queries, 
	dbsize_t num_blocks, dbsize_t words_per_block, 
	dbsize_t virtual_block_size = 1);

template <>
inline void compute_outputvec_multi<GF28_Element>(
        const GF28_Element *data, const std::vector<GF28_Element*> inputvecs, 
	std::vector<GF28_Element*> outputvecs, nqueries_t num_queries, 
	dbsize_t num_blocks, dbsize_t words_per_block, 
	dbsize_t virtual_block_size)
{
    dbsize_t num_virtual_blocks = (num_blocks - 1) / virtual_block_size + 1;
    if (num_queries == 2) {
	const GF28_Element *block = data;
	const GF28_Element *inp1 = inputvecs[0];
	const GF28_Element *inp2 = inputvecs[1];
	dbsize_t words_per_virtual_block = words_per_block * virtual_block_size;
	for (unsigned int j = 0; j < num_virtual_blocks; ++j) {
	    const GF28_Element *multrow1 = GF28_mult_table[*(inp1++)];
	    const GF28_Element *multrow2 = GF28_mult_table[*(inp2++)];
	    GF28_Element *oc1 = outputvecs[0];
	    GF28_Element *oc2 = outputvecs[1];
	    if (j == num_virtual_blocks-1) {
		words_per_virtual_block = words_per_block *
			((num_blocks - 1) % virtual_block_size + 1);
	    }
	    const GF28_Element *blockc = block;
	    GF28_Element *oc1_end = oc1 + (words_per_virtual_block & ~7);
	    while (oc1 < oc1_end) {
		GF28_Element v1 = *(blockc++);
		uint64_t accum1 = (uint64_t) multrow1[v1];
		uint64_t accum2 = (uint64_t) multrow2[v1];
		GF28_Element v2 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v2] << 8;
		accum2 |= (uint64_t) multrow2[v2] << 8;
		GF28_Element v3 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v3] << 16;
		accum2 |= (uint64_t) multrow2[v3] << 16;
		GF28_Element v4 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v4] << 24;
		accum2 |= (uint64_t) multrow2[v4] << 24;
		GF28_Element v5 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v5] << 32;
		accum2 |= (uint64_t) multrow2[v5] << 32;
		GF28_Element v6 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v6] << 40;
		accum2 |= (uint64_t) multrow2[v6] << 40;
		GF28_Element v7 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v7] << 48;
		accum2 |= (uint64_t) multrow2[v7] << 48;
		GF28_Element v8 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v8] << 56;
		accum2 |= (uint64_t) multrow2[v8] << 56;
		*((uint64_t *) oc1) ^= accum1;
		*((uint64_t *) oc2) ^= accum2;
		oc1+=8;
		oc2+=8;
	    }
	    for (unsigned int c = 0; c < (words_per_virtual_block & 7); ++c) {
		GF28_Element v = *(blockc++);
		*(oc1++) ^= multrow1[v];
		*(oc2++) ^= multrow2[v];
	    }
	    block += words_per_virtual_block;
	}
    } else if (num_queries == 3) {
	const GF28_Element *block = data;
	const GF28_Element *inp1 = inputvecs[0];
	const GF28_Element *inp2 = inputvecs[1];
	const GF28_Element *inp3 = inputvecs[2];
	dbsize_t words_per_virtual_block = words_per_block * virtual_block_size;
	for (unsigned int j = 0; j < num_virtual_blocks; ++j) {
	    const GF28_Element *multrow1 = GF28_mult_table[*(inp1++)];
	    const GF28_Element *multrow2 = GF28_mult_table[*(inp2++)];
	    const GF28_Element *multrow3 = GF28_mult_table[*(inp3++)];
	    const GF28_Element *blockc = block;
	    GF28_Element *oc1 = outputvecs[0];
	    GF28_Element *oc2 = outputvecs[1];
	    GF28_Element *oc3 = outputvecs[2];
	    if (j == num_virtual_blocks-1) {
		words_per_virtual_block = words_per_block *
			((num_blocks - 1) % virtual_block_size + 1);
	    }
	    GF28_Element *oc1_end = oc1 + (words_per_virtual_block & ~7);
	    while (oc1 < oc1_end) {
		GF28_Element v1 = *(blockc++);
		uint64_t accum1 = (uint64_t) multrow1[v1];
		uint64_t accum2 = (uint64_t) multrow2[v1];
		uint64_t accum3 = (uint64_t) multrow3[v1];
		GF28_Element v2 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v2] << 8;
		accum2 |= (uint64_t) multrow2[v2] << 8;
		accum3 |= (uint64_t) multrow3[v2] << 8;
		GF28_Element v3 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v3] << 16;
		accum2 |= (uint64_t) multrow2[v3] << 16;
		accum3 |= (uint64_t) multrow3[v3] << 16;
		GF28_Element v4 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v4] << 24;
		accum2 |= (uint64_t) multrow2[v4] << 24;
		accum3 |= (uint64_t) multrow3[v4] << 24;
		GF28_Element v5 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v5] << 32;
		accum2 |= (uint64_t) multrow2[v5] << 32;
		accum3 |= (uint64_t) multrow3[v5] << 32;
		GF28_Element v6 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v6] << 40;
		accum2 |= (uint64_t) multrow2[v6] << 40;
		accum3 |= (uint64_t) multrow3[v6] << 40;
		GF28_Element v7 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v7] << 48;
		accum2 |= (uint64_t) multrow2[v7] << 48;
		accum3 |= (uint64_t) multrow3[v7] << 48;
		GF28_Element v8 = *(blockc++);
		accum1 |= (uint64_t) multrow1[v8] << 56;
		accum2 |= (uint64_t) multrow2[v8] << 56;
		accum3 |= (uint64_t) multrow3[v8] << 56;
		*((uint64_t *) oc1) ^= accum1;
		*((uint64_t *) oc2) ^= accum2;
		*((uint64_t *) oc3) ^= accum3;
		oc1+=8;
		oc2+=8;
		oc3+=8;
	    }
	    for (unsigned int c = 0; c < (words_per_virtual_block & 7); ++c) {
		GF28_Element v = *(blockc++);
		*(oc1++) ^= multrow1[v];
		*(oc2++) ^= multrow2[v];
		*(oc3++) ^= multrow3[v];
	    }
	    block += words_per_virtual_block;
	}
    } else {
	const GF28_Element *block = data;
	const GF28_Element **inpv = new const GF28_Element*[num_queries];
	const GF28_Element **multrowv = new const GF28_Element*[num_queries];
	GF28_Element **ocv = new GF28_Element*[num_queries];
	for (unsigned int q=0; q<num_queries; ++q) {
	    inpv[q] = inputvecs[q];
	}
	dbsize_t words_per_virtual_block = words_per_block * virtual_block_size;
	for (unsigned int j = 0; j < num_virtual_blocks; ++j) {
	    for (unsigned int q=0; q<num_queries; ++q) {
		multrowv[q] = GF28_mult_table[*(inpv[q]++)];
	    }
	    const GF28_Element *blockc = block;
	    GF28_Element *ocv0 = outputvecs[0];
	    if (j == num_virtual_blocks-1) {
		words_per_virtual_block = words_per_block *
			((num_blocks - 1) % virtual_block_size + 1);
	    }
	    for (unsigned int q=1; q<num_queries; ++q) {
		ocv[q] = outputvecs[q];
	    }
	    GF28_Element *oc_end = ocv0 + (words_per_virtual_block & ~7);
	    while (ocv0 < oc_end) {
		GF28_Element v1 = *(blockc++);
		GF28_Element v2 = *(blockc++);
		GF28_Element v3 = *(blockc++);
		GF28_Element v4 = *(blockc++);
		GF28_Element v5 = *(blockc++);
		GF28_Element v6 = *(blockc++);
		GF28_Element v7 = *(blockc++);
		GF28_Element v8 = *(blockc++);
		uint64_t accum = (uint64_t) multrowv[0][v1];
		accum |= (uint64_t) multrowv[0][v2] << 8;
		accum |= (uint64_t) multrowv[0][v3] << 16;
		accum |= (uint64_t) multrowv[0][v4] << 24;
		accum |= (uint64_t) multrowv[0][v5] << 32;
		accum |= (uint64_t) multrowv[0][v6] << 40;
		accum |= (uint64_t) multrowv[0][v7] << 48;
		accum |= (uint64_t) multrowv[0][v8] << 56;
		*((uint64_t *) ocv0) ^= accum;
		ocv0 += 8;
		for (unsigned int q=1; q<num_queries; ++q) {
		    uint64_t accum = (uint64_t) multrowv[q][v1];
		    accum |= (uint64_t) multrowv[q][v2] << 8;
		    accum |= (uint64_t) multrowv[q][v3] << 16;
		    accum |= (uint64_t) multrowv[q][v4] << 24;
		    accum |= (uint64_t) multrowv[q][v5] << 32;
		    accum |= (uint64_t) multrowv[q][v6] << 40;
		    accum |= (uint64_t) multrowv[q][v7] << 48;
		    accum |= (uint64_t) multrowv[q][v8] << 56;
		    *((uint64_t *) ocv[q]) ^= accum;
		    ocv[q] += 8;
		}
	    }
	    for (unsigned int c = 0; c < (words_per_virtual_block & 7); ++c) {
		GF28_Element v = *(blockc++);
		*(ocv0++) ^= multrowv[0][v];
		for (unsigned int q=1; q<num_queries; ++q) {
		    *(ocv[q]++) ^= multrowv[q][v];
		}
	    }
	    block += words_per_virtual_block;
	}
	delete[] inpv;
	delete[] multrowv;
	delete[] ocv;
    }
}

template <>
inline void compute_outputvec_multi<GF216_Element>(
        const GF216_Element *data, const std::vector<GF216_Element*> inputvecs, 
	std::vector<GF216_Element*> outputvecs, nqueries_t num_queries, 
	dbsize_t num_blocks, dbsize_t words_per_block, 
	dbsize_t virtual_block_size)
{
    for (nqueries_t q = 0; q < num_queries; q++) {
        compute_outputvec_single<GF216_Element>(data, inputvecs[q],
                outputvecs[q], num_blocks, words_per_block);
    }
}

template<>
inline nqueries_t PercyServer_GF2E<GF28_Element>::optimal_strassen_depth(
	nqueries_t res_rows, dbsize_t inner_dim, dbsize_t res_cols)
{
    nqueries_t depth = 0;
    if( res_rows < 10 ) {
	// Don't do Strassen
	return 0;
    }

    nqueries_t max_number_rows = 3;

    // For larger database don't go so far down yet
    if( (inner_dim >= 32768 && res_rows < 48) ||
	(inner_dim >= 16384 && res_rows < 32)) {
	max_number_rows = 7;
    }

    // The optimal level of 3 for res_rows results from the optimized
    // multiplication routines for 1, 2 and 3 queries
    while(res_rows > max_number_rows && inner_dim > 1 && res_cols > 1) {
	res_rows /= 2;
	inner_dim /= 2;
	res_cols /= 2;

	depth++;
    }

    return depth;
}

template<>
inline nqueries_t PercyServer_GF2E<GF216_Element>::optimal_strassen_depth(
	nqueries_t res_rows, dbsize_t inner_dim, dbsize_t res_cols)
{
    std::cerr << "Warning: optimal levels for GF216 not yet determined"
	<< std::endl;

    nqueries_t depth = 0;
    if( res_rows < 10 ) {
	// Don't do Strassen
	return 0;
    }

    while(res_rows > 4 && inner_dim > 1 && res_cols > 1) {
	res_rows /= 2;
	inner_dim /= 2;
	res_cols /= 2;

	depth++;
    }

    return depth;
}

template <typename GF2E_Element>
GF2E_Element * PercyServer_GF2E<GF2E_Element>::alloc_strassen_memory(
	nqueries_t res_rows, dbsize_t inner_dim, dbsize_t res_cols)
{
    size_t num_elts = 0;

    while(res_rows > 1 && inner_dim > 1 && res_cols > 1) {
	res_rows /= 2;
	inner_dim /= 2;
	res_cols /= 2;

	num_elts += res_rows * inner_dim;
	num_elts += inner_dim * res_cols;
	num_elts += res_rows * res_cols;
    }

    GF2E_Element * memory = new GF2E_Element[num_elts];

    return memory;
}

template <typename GF2E_Element>
void PercyServer_GF2E<GF2E_Element>::free_strassen_memory(
	GF2E_Element * memory)
{
    delete[] memory;
}

template <typename GF2E_Element>
inline void PercyServer_GF2E<GF2E_Element>::strassen_split(
        Matrix<const GF2E_Element> &mat_a,
	Matrix<const GF2E_Element> &mat_b, Matrix<GF2E_Element> &mat_c,
	dbsize_t depth, GF2E_Element * memory)
{
    dbsize_t q, r, s, qrem, rrem, srem;

    dbsize_t res_rows  = mat_a.num_rows;
    dbsize_t inner_dim = mat_a.num_cols;
    dbsize_t res_cols  = mat_b.num_cols;

    mat_c.clear();

    if( res_rows <= 1 || inner_dim <= 1 ||
	    res_cols <= 1 || depth > strassen_real_max_depth) {
	if(res_rows > 1) {
	    // TODO: See if we want to change this.
	    std::vector<GF2E_Element*> inputvec, outputvec;
	    for (nqueries_t q = 0; q < res_rows; ++q) {
		inputvec.push_back(const_cast<GF2E_Element*>(mat_a.data + q * mat_a.num_cols));
		outputvec.push_back(const_cast<GF2E_Element*>(mat_c.data + q * mat_c.num_cols));
	    }
	    compute_outputvec_multi(mat_b.data, inputvec, outputvec,
		    res_rows, inner_dim, res_cols);
	} else {
	    compute_outputvec_single(mat_b.data, mat_a.data, mat_c.data,
		    inner_dim, res_cols);
	}

	return;
    }

    if(depth > strassen_level_reached) {
	strassen_level_reached = depth;
    }

    qrem = res_rows % 2;
    rrem = inner_dim % 2;
    srem = res_cols % 2;

    q = res_rows - qrem;
    r = inner_dim - rrem;
    s = res_cols - srem;

    SubMatrix<const GF2E_Element> a11 =
	SubMatrix<const GF2E_Element> (mat_a, 0, 0, q, r);
    Col<const GF2E_Element> a12 =
	Col<const GF2E_Element> (mat_a, 0, r, q);
    Row<const GF2E_Element> a21 =
	Row<const GF2E_Element> (mat_a, q, 0, r);
    Elem<const GF2E_Element> a22 =
	Elem<const GF2E_Element> (mat_a, q, r);

    SubMatrix<const GF2E_Element> b11 =
	SubMatrix<const GF2E_Element> (mat_b, 0, 0, r, s);
    Col<const GF2E_Element> b12 =
	Col<const GF2E_Element> (mat_b, 0, s, r);
    Row<const GF2E_Element> b21 =
	Row<const GF2E_Element> (mat_b, r, 0, s);
    Elem<const GF2E_Element> b22 =
	Elem<const GF2E_Element> (mat_b, r, s);

    SubMatrix<GF2E_Element> c11 =
	SubMatrix<GF2E_Element> (mat_c, 0, 0, q, s);
    Col<GF2E_Element> c12 =
	Col<GF2E_Element> (mat_c, 0, s, q);
    Row<GF2E_Element> c21 =
	Row<GF2E_Element> (mat_c, q, 0, s);
    Elem<GF2E_Element> c22 =
	Elem<GF2E_Element> (mat_c, q, s);

    // A21 * B11
    if(qrem > 0) {
	c21.add_mult_of(a21, b11);
    }

    // A11 * B12
    if(srem > 0) {
	c12.add_mult_of(a11, b12);
    }

    // A21 * B12
    if(qrem > 0 && srem > 0) {
	c22.add_mult_of(a21, b12);
    }

    if(rrem > 0) {
	// A12 * B21
	c11.add_mult_of(a12, b21);

	// A22 * B21
	if(qrem > 0) {
	    c21.add_mult_of(a22, b21);
	}

	// A12 * B22
	if(srem > 0) {
	    c12.add_mult_of(a12, b22);
	}

	// A22 * B22
	if(qrem > 0 && srem > 0) {
	    c22.add_mult_of(a22, b22);
	}
    }

    strassen_mult(mat_a, mat_b, mat_c, depth, memory);
}


template <typename GF2E_Element>
inline void PercyServer_GF2E<GF2E_Element>::strassen_mult(
        Matrix<const GF2E_Element> &mat_a,
	Matrix<const GF2E_Element> &mat_b, Matrix<GF2E_Element> &mat_c,
	dbsize_t depth, GF2E_Element * memory)
{
    dbsize_t nq = mat_a.num_rows / 2;
    dbsize_t nr = mat_a.num_cols / 2;
    dbsize_t ns = mat_b.num_cols / 2;

    Matrix<const GF2E_Element> matrix_left =
	Matrix<const GF2E_Element>(memory, nq, nr);
    memory += nq * nr;
    Matrix<const GF2E_Element> matrix_right =
	Matrix<const GF2E_Element>(memory, nr, ns);
    memory += nr * ns;
    Matrix<GF2E_Element> matrix_subresult =
	Matrix<GF2E_Element>(memory, nq, ns);
    memory += nq * ns;

    SubMatrix<const GF2E_Element> a11, a12, a21, a22;
    mat_a.strassen_submatrices(a11, a12, a21, a22);

    SubMatrix<const GF2E_Element> b11, b12, b21, b22;
    mat_b.strassen_submatrices(b11, b12, b21, b22);

    SubMatrix<GF2E_Element> c11, c12, c21, c22;
    mat_c.strassen_submatrices(c11, c12, c21, c22);

    // #1
    matrix_left.is_sum_of(a11, a22);
    matrix_right.is_sum_of(b11, b22);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c11.add(matrix_subresult);
    c22.add(matrix_subresult);

    // #2
    matrix_left.is_sum_of(a21, a22);
    matrix_right.copy_from(b11);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c21.add(matrix_subresult);
    c22.add(matrix_subresult);

    // #3
    matrix_left.copy_from(a11);
    matrix_right.is_sum_of(b12, b22);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c12.add(matrix_subresult);
    c22.add(matrix_subresult);

    // #4
    matrix_left.copy_from(a22);
    matrix_right.is_sum_of(b21, b11);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c11.add(matrix_subresult);
    c21.add(matrix_subresult);

    // #5
    matrix_left.is_sum_of(a11, a12);
    matrix_right.copy_from(b22);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c11.add(matrix_subresult);
    c12.add(matrix_subresult);

    // #6
    matrix_left.is_sum_of(a21, a11);
    matrix_right.is_sum_of(b11, b12);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c22.add(matrix_subresult);

    // #7
    matrix_left.is_sum_of(a12, a22);
    matrix_right.is_sum_of(b21, b22);
    strassen_split(matrix_left, matrix_right, matrix_subresult,
	    depth + 1, memory);
    c11.add(matrix_subresult);
}

template <typename GF2E_Element>
nqueries_t PercyServer_GF2E<GF2E_Element>::optimal_strassen_num_queries(nqueries_t num_queries)
{
    int bitlength = 0;
    nqueries_t tmp = num_queries;
    nqueries_t increased;
    double factor;

    while(tmp > 0) {
	bitlength++;
	tmp /= 2;
    }

    nqueries_t ptwo = 2 << (bitlength - 1);
    nqueries_t pthree = 3 << (bitlength - 2);

    if(pthree > num_queries) {
	increased = pthree;
	factor = 1.125;
    } else {
	increased = ptwo;
	factor = 1.05;
    }

    // Only increase size if the cost is not too much
    if (((double) increased) / num_queries < factor) {
	return increased;
    } else {
	return num_queries;
    }
}

template <typename GF2E_Element>
bool PercyServer_GF2E<GF2E_Element>::handle_request_impl (
	std::vector<unsigned char*> requests, 
	std::vector<unsigned char*> responses)
{
    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }

    // Read some values from the params
    dbsize_t words_per_block = params->words_per_block();
    dbsize_t num_blocks = params->num_blocks();
    dbsize_t block_size = params->block_size();
    dbsize_t virtual_block_size = params->virtual_block_size();

    // Strassen is more efficient for certain number of queries,
    // so we test if a resize is necessary.
    dbsize_t optimal_num_queries = optimal_strassen_num_queries(num_queries);

    const GF2E_Element *data = (const GF2E_Element*)(datastore->get_data());

    strassen_level_reached = 0;

    // Compute the output vector and send it back to the client

    //struct timeval ts, te;
    //gettimeofday(&ts, NULL);
    if (num_queries > 1) {
	if (strassen_max_depth == PercyServer::STRASSEN_OPTIMAL) {
	    strassen_real_max_depth =
		optimal_strassen_depth(optimal_num_queries, num_blocks,
			words_per_block);
	} else {
	    strassen_real_max_depth = strassen_max_depth;
	}

	if (strassen_real_max_depth == 0 || virtual_block_size > 1) {
	    std::vector<GF2E_Element*> inputvecs, outputvecs;
	    for (nqueries_t q = 0; q < num_queries; ++q) {
		memset(responses[q], '\0', block_size);
		inputvecs.push_back((GF2E_Element*)(requests[q]));
		outputvecs.push_back((GF2E_Element*)(responses[q]));
	    }
	    compute_outputvec_multi<GF2E_Element>(data, inputvecs, outputvecs,
		    num_queries, num_blocks, words_per_block, 
		    virtual_block_size);
	} else {
	    // TODO: Take a look at this and see if it needs to change
	    GF2E_Element * input = 
		    new GF2E_Element[optimal_num_queries * num_blocks];
	    GF2E_Element * output = 
		    new GF2E_Element[optimal_num_queries * words_per_block];
	    memset(output, '\0', optimal_num_queries * words_per_block * 
		    sizeof(GF2E_Element));
	    for (nqueries_t q = 0; q < num_queries; ++q) {
		memcpy(input + q*num_blocks*sizeof(GF2E_Element), requests[q], 
			num_blocks*sizeof(GF2E_Element));
	    }
	    if (optimal_num_queries > num_queries) {
		memset(input +
			num_queries * num_blocks * sizeof(GF2E_Element),
			'\0', (optimal_num_queries - num_queries) *
			num_blocks * sizeof(GF2E_Element));
	    }

	    GF2E_Element * memory = alloc_strassen_memory(
		    optimal_num_queries, num_blocks, words_per_block);
	    Matrix<const GF2E_Element> mat_input =
		Matrix<const GF2E_Element>(input, optimal_num_queries,
			num_blocks);
	    Matrix<const GF2E_Element> mat_data =
		Matrix<const GF2E_Element>(data, num_blocks,
			words_per_block);
	    Matrix<GF2E_Element> mat_output =
		Matrix<GF2E_Element>(output, optimal_num_queries,
			words_per_block);
	    strassen_split(mat_input, mat_data, mat_output, 1, memory);
	    free_strassen_memory(memory);

	    for (nqueries_t q = 0; q < num_queries; ++q) {
		memcpy(responses[q], output + q*words_per_block*sizeof(GF2E_Element),
			block_size);
	    }
	}
    } else {
	memset(responses[0], '\0', block_size);
	compute_outputvec_single<GF2E_Element>(data, 
		(GF2E_Element*)(requests[0]), (GF2E_Element*)(responses[0]),
		num_blocks, words_per_block, virtual_block_size);
    }

    //gettimeofday(&te, NULL);
    //int td = (te.tv_sec - ts.tv_sec)*1000000 + (te.tv_usec - ts.tv_usec);
    //fprintf(stderr, "%d.%3d msec computation\n", td/1000, td%1000);

    // If the server is Byzantine, give wrong output
    if (byzantine) {
	for (nqueries_t q = 0; q < num_queries; q++) {
	    for (dbsize_t i = 0; i < words_per_block; i++) {
		responses[q][i]++;
	    }
	}
    }

    return true;
}

template <typename GF2E_Element>
void PercyServer_GF2E<GF2E_Element>::combine_results (unsigned char * result, 
	std::vector<unsigned char*> worker_results)
{
    dbsize_t block_size = params->block_size() * params->virtual_block_size();
    for (nservers_t i = 0; i < worker_results.size(); ++i) {
	XOR_equal(result, worker_results[i], block_size);
    }
}

#endif
