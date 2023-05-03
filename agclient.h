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

#ifndef __AGCLIENT_H__
#define __AGCLIENT_H__

#include <NTL/mat_ZZ_p.h>
#include "percyclient.h"
#include "agparams.h"
#include "percystats.h"

NTL_CLIENT

#ifdef NEED_UINT128
#define AG_NTL
#endif

struct AGDecodeInfo {
#ifdef AG_NTL
    mat_ZZ_p A_inv_B;
    vec_ZZ_p Delta_vec_inv;
#else
    AG_Element * A_inv_B;
    AG_Element * Delta_vec_inv;
#endif

    std::vector<dbsize_t> perm;
};


/// A PIR client for the CPIR protocol by Aguilar Melchor and Gaborit (2007).  
/// This protocol was introduced in 
/// <a href="https://eprint.iacr.org/2007/446">A Lattice-Based Computationally-Efficient Private Information Retrieval Protocol</a> 
/// and revisited in 
/// <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?tp=&arnumber=4622593">HighSpeed Private Information Retrieval Computation on GPU</a>.

class PercyAGClient : public PercyClient {
public:
    /// Constructor.
    /// @param params	Parameters for the client.
    /// @param stats	Statistics collection object.  No statistics will be
    ///			collected if NULL.
    PercyAGClient (const PercyClientParams * params, PercyStats * stats = NULL);

    /// Destructor.
    virtual ~PercyAGClient ();

private:
    /// Parameters for the protocol.
    const AGParams * params;

    // Virtual members as described in PercyClient class
    virtual void encode_request_impl (nqueries_t request_identifier);
    virtual dbsize_t send_request_impl (nqueries_t request_identifier, 
	    vector<ostream*> &osvec, bool send_num_queries = true);
    virtual dbsize_t receive_replies_impl (nqueries_t request_identifier,
	    vector<istream*> &isvec);
    virtual nqueries_t process_replies_impl (nservers_t h,
	    vector<vector<PercyResult> >& results);

#ifndef AG_NTL
    // Generate random elements of ZZ_p for the first n elements of buffer.  If
    // nonzero is set to true, then the random values will be non-zero.
    // NOTE: 
    //  - Assumes that buffer is large enough to hold n elements.
    void random_ag_elements (AG_Element * buffer, dbsize_t n, 
	    bool nonzero = false);

    void soft_noise_matrix (AG_Element * M, dbsize_t size);
    void hard_noise_matrix (AG_Element * M, dbsize_t size);

    // Compute result = A*B where A is rows_A-by-cols_A and B is cols_A-by-cols_B
    // NOTE: 
    //  - Assumes that result is large enough to hold rows_A*cols_B elements.
    //  - Also assumes cols_A*(p-1)^2 fits in an AG_Element (that is, modulo p
    //    operations do not need to be done until the end).
    void mult_matrices (AG_Element * results, AG_Element * A, AG_Element * B, 
	    dbsize_t rows_A, dbsize_t cols_A, dbsize_t cols_B);

    // Find the inverse of n modulo p and put it in inv.  Return true if such and
    // inverse exists and false otherwise.
    bool inv (AG_Element& inv, const AG_Element& n);

    // Do row operations on M to discover if M is invertible.  Return true if M is
    // invertible and false othersize.
    // If aux and/or aux2 are specified, then the row operations will also be
    // performed on them.  If M is invertible, then aux becomes M^(-1)*aux and aux2
    // becomes M^(-1)*aux2.
    // NOTE:
    //  - Assumes that all of M, aux and aux2 are size-by-size matrices.
    //  - All row operations are in-place and change the contents of M, aux and
    //    aux2.  Use copies of these matrices if you will need to keep the
    //    original.
    //  - Assumes that size*(p-1)^2 fits in an AG_Element
    //  - If M is discovered to be singular, the algorithm short-circuits and the
    //    values in M, aux and aux2 are not guaranteed.
    bool is_invertible (AG_Element * M, dbsize_t size, AG_Element * aux = NULL, 
	    AG_Element * aux2 = NULL);
#endif
    
    AG_Element p;
    AG_Element q;
    ZZ q_ntl;
    ZZ p_ntl;

    
    map<nqueries_t, vector<unsigned char *> > queries;
    map<nqueries_t, vector<AGDecodeInfo> > sentinfo;
#ifdef AG_NTL
    map<nqueries_t, map<nservers_t, vec_vec_ZZ_p> > stored_answers;
#else
    map<nqueries_t, map<nservers_t, vector<AG_Element *> > > stored_answers;
#endif
};

#endif
