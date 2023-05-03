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

#ifndef __ITSERVER_H__
#define __ITSERVER_H__

#include <vec_ZZ_p.h>
#include "percyserver.h"
#include "gf2e.h"
#include "itparams.h"

NTL_CLIENT

/// A PIR server for the IT-PIR protocol by Goldberg (2007) over the integers
/// modulo p.
/// This protocol was introduced in
/// <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=4223220">Improving the Robustness of Private Information Retrieval</a>
/// and revisited in
/// <ul>
/// <li><a href="http://link.springer.com/chapter/10.1007%2F978-3-642-27576-0_13">Revisiting the Computational Practicality of Private Information Retrieval</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2046784">Practical PIR for Electronic Commerce</a></li>
/// <li><a href="https://www.usenix.org/conference/usenixsecurity12/technical-sessions/presentation/devet">Optimally Robust Private Information Retrieval</a></li>
/// <li><a href="http://internetsociety.org/doc/one-block-size-fits-all-pir-and-spir-variable-length-records-multi-block-queries">One (Block) Size Fits All: PIR and SPIR Over Arbitrary-Length Records via Multi-block PIR Queries</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2517840.2517854">Outsourced Private Information Retrieval with Pricing and Access Control</a></li>
/// </ul>

class PercyServer_ZZ_p : public PercyServer {
public:
    /// Constructor.
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    PercyServer_ZZ_p (DataStore * datastore, const PercyServerParams * params,
	    PercyStats * stats = NULL);

    /// Destructor
    virtual ~PercyServer_ZZ_p ();

private:
    virtual bool handle_request_impl (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses);

    virtual void combine_results (unsigned char * result, 
	    std::vector<unsigned char*> worker_results);

    const ZZ_pParams * params;

#ifdef STRASSEN_COUNT_OPERATIONS
    long nr_multiplications;
    long nr_additions;
#endif

    void compute_one(const unsigned char * data, ZZ_p * value, 
	    bool hybrid_protection, dbsize_t num_blocks, 
	    dbsize_t words_per_block, dbsize_t bytes_per_word,
	    nqueries_t num_queries, const vec_ZZ_p *inputvector, dbsize_t c);
    void compute_all(const unsigned char * data, vec_ZZ_p * responses,
            bool hybrid_protection, dbsize_t num_blocks,
	    dbsize_t words_per_block, dbsize_t bytes_per_word, 
	    dbsize_t last_block_words, nqueries_t num_queries, 
	    const vec_ZZ_p * inputvector);
    void compute_mult_strassen(const unsigned char * data,
	    vec_ZZ_p * responses, dbsize_t num_blocks,
	    dbsize_t words_per_block, dbsize_t bytes_per_word,
	    dbsize_t last_block_words, nqueries_t num_queries, 
	    const vec_ZZ_p * inputvector);
    void matrix_mult_naive(const vec_ZZ_p * matrix_a,
	    dbsize_t num_rows_a, const vec_ZZ_p * matrix_b,
	    dbsize_t num_cols_b, vec_ZZ_p * matrix_res,
	    dbsize_t row_offset = 0, dbsize_t col_offset = 0,
	    dbsize_t inner_offset = 0, dbsize_t inner_dim = 0);
    nqueries_t optimal_strassen_depth( nqueries_t num_queries,
	    dbsize_t num_blocks, dbsize_t words_per_block,
	    dbsize_t bytes_per_word);
    void strassen_split(const vec_ZZ_p * matrix_a, dbsize_t num_rows_a,
	    const vec_ZZ_p * matrix_b, dbsize_t num_cols_b,
	    vec_ZZ_p * matrix_res, nqueries_t depth);
    void strassen_mult(const vec_ZZ_p * matrix_a, dbsize_t num_rows_a,
	    const vec_ZZ_p * matrix_b, dbsize_t num_cols_b,
	    vec_ZZ_p * matrix_res, nqueries_t depth);
};

/// A PIR server for the IT-PIR protocol by Chor et al.\ (1995).
/// This protocol was introduced in
/// <a href="http://dl.acm.org/citation.cfm?id=796270">Private Information Retrieval</a>.

class PercyServer_Chor : public PercyServer {
public:
    /// Constructor.
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    PercyServer_Chor (DataStore * datastore, const PercyServerParams * params,
	    PercyStats * stats = NULL);

    /// Destructor.
    virtual ~PercyServer_Chor ();

private:
    virtual bool handle_request_impl (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses);

    virtual void combine_results (unsigned char * result, 
	    std::vector<unsigned char*> worker_results);

    const ChorParams * params;
};

/// A PIR server for the IT-PIR protocol by Goldberg (2007) over GF(2^E).
/// Supported fields are GF(2^8) and GF(2^16).
/// This protocol was introduced in
/// <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=4223220">Improving the Robustness of Private Information Retrieval</a>
/// and revisited in
/// <ul>
/// <li><a href="http://link.springer.com/chapter/10.1007%2F978-3-642-27576-0_13">Revisiting the Computational Practicality of Private Information Retrieval</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2046784">Practical PIR for Electronic Commerce</a></li>
/// <li><a href="https://www.usenix.org/conference/usenixsecurity12/technical-sessions/presentation/devet">Optimally Robust Private Information Retrieval</a></li>
/// <li><a href="http://internetsociety.org/doc/one-block-size-fits-all-pir-and-spir-variable-length-records-multi-block-queries">One (Block) Size Fits All: PIR and SPIR Over Arbitrary-Length Records via Multi-block PIR Queries</a></li>
/// <li><a href="http://dl.acm.org/citation.cfm?id=2517840.2517854">Outsourced Private Information Retrieval with Pricing and Access Control</a></li>
/// </ul>

template <typename GF2E_Element>
class PercyServer_GF2E : public PercyServer {
public:
    /// Constructor.
    /// @param datastore    Database the server will use.
    /// @param params	    Parameters for the server.
    /// @param stats	    Statistics collection object.  No statistics will be
    ///                     collected if NULL.
    PercyServer_GF2E (DataStore * datastore, const PercyServerParams * params,
	    PercyStats * stats = NULL);

    /// Destructor.
    virtual ~PercyServer_GF2E ();

private:
    inline void strassen_mult(Matrix<const GF2E_Element> &mat_a,
	    Matrix<const GF2E_Element> &mat_b,
	    Matrix<GF2E_Element> &mat_c,
	    dbsize_t depth, GF2E_Element * memory);

    inline void strassen_split(Matrix<const GF2E_Element> &mat_a,
	    Matrix<const GF2E_Element> &mat_b,
	    Matrix<GF2E_Element> &mat_c,
	    dbsize_t depth, GF2E_Element * memory);

    GF2E_Element * alloc_strassen_memory(nqueries_t res_rows,
	    dbsize_t inner_dim, dbsize_t res_cols);
    void free_strassen_memory(GF2E_Element * memory);

    nqueries_t optimal_strassen_num_queries(nqueries_t num_queries);
    inline nqueries_t optimal_strassen_depth(nqueries_t res_rows,
	    dbsize_t inner_dim, dbsize_t res_cols);

    nqueries_t strassen_real_max_depth;

    virtual bool handle_request_impl (std::vector<unsigned char*> requests, 
	    std::vector<unsigned char*> responses);

    virtual void combine_results (unsigned char * result, 
	    std::vector<unsigned char*> worker_results);

    const GF2EParams * params;
};

#include "itserver_impl.h"

#endif
