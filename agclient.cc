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

#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <utility>
#include <sstream>
#include <fstream>
#include "percytypes.h"
#include "agclient.h"
#include "percyio.h"

void print_element (std::ostream& os, AG_Element n);
void print_matrix (std::ostream& os, AG_Element * M, dbsize_t rows, 
	dbsize_t cols);

PercyAGClient::PercyAGClient (const PercyClientParams * params, 
	PercyStats * stats)
:
    PercyClient(params, 1, 0, stats),
    params(static_cast<const AGParams*>(params->percy_params()))
{
#ifdef NEED_UINT128
    q_ntl = this->params->q();
    p_ntl = this->params->p();
#else
    q = this->params->q();
    ZZFromBytes(q_ntl, (unsigned char*)(&q), sizeof(AG_Element));
    p = this->params->p();
    ZZFromBytes(p_ntl, (unsigned char*)(&p), sizeof(AG_Element));
#endif
}

PercyAGClient::~PercyAGClient ()
{
    map<nqueries_t, vector<unsigned char *> >::iterator iter;
    for (iter = queries.begin(); iter != queries.end(); ++iter) {
	for (nqueries_t i = 0; i < iter->second.size(); ++i) {
	    delete[] iter->second[i];
	}
    }
#ifndef AG_NTL
    map<nqueries_t, map<nservers_t, vector<AG_Element*> > >::iterator iter2;
    for (iter2 = stored_answers.begin(); iter2 != stored_answers.end(); ++iter2) {
	map<nservers_t, vector<AG_Element*> >::iterator iterinner;
	for (iterinner = iter2->second.begin(); iterinner != iter2->second.end(); ++iterinner) {
	    vector<AG_Element*>& qanswer = iterinner->second;
	    for (nservers_t j = 0; j < qanswer.size(); ++j) {
		delete[] qanswer[j];
	    }
	}
    }
    map<nqueries_t, vector<AGDecodeInfo> >::iterator iter3;
    for (iter3 = sentinfo.begin(); iter3 != sentinfo.end(); ++iter3) {
	for (nqueries_t i = 0; i < iter3->second.size(); ++i) {
	    AGDecodeInfo& info = iter3->second[i];
	    if (info.A_inv_B) {
		delete[] info.A_inv_B;
	    }
	    if (info.Delta_vec_inv) {
		delete[] info.Delta_vec_inv;
	    }
	}
    }
#endif
}

typedef std::pair<dbsize_t, dbsize_t> bufferloc_t;

bufferloc_t add_to_buffer (unsigned char * buffer, bufferloc_t location, ZZ_p k,
	dbsize_t num_bits)
{
    ZZ n = rep(k);
    n &= power2_ZZ(num_bits) - 1;
    n <<= location.second;
    dbsize_t bytes_needed = (num_bits - 1) / 8 + 2;
    unsigned char * buf = new unsigned char[bytes_needed];
    BytesFromZZ(buf, n, bytes_needed);
    for (dbsize_t i = 0; i < bytes_needed; ++i) {
	if (buf[i] != 0x0) {
	    buffer[location.first + i] |= buf[i];
	}
    }
    delete[] buf;
    return make_pair(location.first + ((location.second + num_bits) / 8),
	    (location.second + num_bits) % 8);
}

bufferloc_t get_from_buffer (unsigned char * buffer, bufferloc_t location, ZZ_p &k,
	dbsize_t num_bits)
{
    ZZ n;
    dbsize_t bytes_needed = (num_bits - 1) / 8 + 2;
    ZZFromBytes(n, buffer + location.first, bytes_needed);
    n >>= location.second;
    n &= power2_ZZ(num_bits) - 1;
    k = to_ZZ_p(n);
    return make_pair(location.first + ((location.second + num_bits) / 8),
	    (location.second + num_bits) % 8);
}

#ifdef AG_NTL
void random_vector (vec_ZZ_p &V)
{
    for (dbsize_t i = 0; i < (dbsize_t)V.length(); ++i) {
	random(V[i]);
	//V[i] = to_ZZ_p(RandomBnd(ZZ_p::modulus()));
    }
}

void random_matrix (mat_ZZ_p &M)
{
    for (dbsize_t i = 0; i < (dbsize_t)M.NumRows(); ++i) {
	random_vector(M[i]);
    }
}

void soft_noise_matrix (mat_ZZ_p &M)
{
    for (dbsize_t i = 0; i < (dbsize_t)M.NumRows(); ++i) {
	vec_ZZ_p& Mi = M[i];
	for (dbsize_t j = 0; j < (dbsize_t)M.NumCols(); ++j) {
	    if (RandomBnd(2) == 0) {
		Mi[j] = -1;
	    } else {
		Mi[j] = 1;
	    }
	}
    }
}

void hard_noise_matrix (mat_ZZ_p &M, ZZ_p q)
{
    soft_noise_matrix(M);
    for (dbsize_t i = 0; i < (dbsize_t)M.NumRows(); ++i) {
	M(i+1, i+1) = q;
    }
}

#else

bufferloc_t add_to_buffer (unsigned char * buffer, bufferloc_t location, 
	AG_Element n, dbsize_t num_bits)
{
    n <<= location.second;
    dbsize_t bytes_needed = (num_bits - 1) / 8 + 2;
    unsigned char * buf = (unsigned char*)(&n);
    for (dbsize_t i = 0; i < bytes_needed; ++i) {
	if (buf[i] != 0x0) {
	    buffer[location.first + i] |= buf[i];
	}
    }
    return make_pair(location.first + ((location.second + num_bits) / 8),
	    (location.second + num_bits) % 8);
}

bufferloc_t get_from_buffer (unsigned char * buffer, bufferloc_t location, 
	AG_Element &n, dbsize_t num_bits)
{
    dbsize_t num_bits_copy = num_bits;
    dbsize_t bytes_needed = (num_bits - 1) / 8 + 1;
    num_bits += location.second;
    unsigned char * buf = (unsigned char*)(&n);
    for (dbsize_t i = 0; i < bytes_needed; ++i) {
	buf[i] |= buffer[location.first + i];
	if (num_bits < 8) {
	    buf[i] &= (unsigned char)((1 << num_bits) - 1);
	    break;
	} else {
	    num_bits -= 8;
	}
    }
    n >>= location.second;
    return make_pair(location.first + ((location.second + num_bits_copy) / 8),
	    (location.second + num_bits_copy) % 8);
}

// Generate random elements of ZZ_p for the first n elements of buffer.  If
// nonzero is set to true, then the random values will be non-zero.
// NOTE: Assumes that buffer is large enough to hold n elements.
void PercyAGClient::random_ag_elements (AG_Element * buffer, dbsize_t n, 
	bool nonzero)
{
    dbsize_t size = sizeof(AG_Element);
    for (dbsize_t i = 0; i < n; ++i) {
	BytesFromZZ((unsigned char*)(buffer + i), RandomBnd(p_ntl), size);
	while (nonzero && buffer[i] == 0) {
	    BytesFromZZ((unsigned char*)(buffer + i), RandomBnd(p_ntl), size);
	}
    }
}

void PercyAGClient::soft_noise_matrix (AG_Element * M, dbsize_t size)
{
    for (dbsize_t i = 0; i < size; ++i) {
	AG_Element * Mi = M + i * size;
	for (dbsize_t j = 0; j < size; ++j) {
	    if (RandomBnd(2) == 0) {
		Mi[j] = p-1;
	    } else {
		Mi[j] = 1;
	    }
	}
    }
}

void PercyAGClient::hard_noise_matrix (AG_Element * M, dbsize_t size)
{
    soft_noise_matrix(M, size);
    for (dbsize_t i = 0; i < size; ++i) {
	M[i*size+i] = q;
    }
}

// Compute result = A*B where A is rows_A-by-cols_A and B is cols_A-by-cols_B
// NOTE: Assumes that result is large enough to hold rows_A*cols_B elements.
//       Also assumes cols_A*(p-1)^2 fist in an AG_Element (that is, modulo p
//       operations do not need to be done until the end).
void PercyAGClient::mult_matrices (AG_Element * results, AG_Element * A, 
	AG_Element * B, dbsize_t rows_A, dbsize_t cols_A, dbsize_t cols_B)
{
    memset(results, 0, sizeof(AG_Element) * rows_A * cols_B);
    AG_Element * a, * b, * r;
    a = A;
    for (dbsize_t i = 0; i < rows_A; ++i) {
	b = B;
	for (dbsize_t j = 0; j < cols_A; ++j) {
	    r = results + i * cols_B;
	    for (dbsize_t k = 0; k < cols_B; ++k) {
		*r += *a * *b;
		++r;
		++b;
	    }
	    ++a;
	}
    }
    for (dbsize_t i = 0; i < rows_A * cols_B; ++i) {
	results[i] %= p;
    }
}

// Find the inverse of n modulo p and put it in inv.  Return true if such and
// inverse exists and false otherwise.
bool PercyAGClient::inv (AG_Element& inv, const AG_Element& n)
{
    if (n == 0 || n >= p) {
	return false;
    }
    AG_Element r2 = p;
    AG_Element r1 = n;
    AG_Element x2 = ULONG_TO_AGELT(1);
    AG_Element x1 = ULONG_TO_AGELT(0);
    AG_Element y2 = ULONG_TO_AGELT(0);
    AG_Element y1 = ULONG_TO_AGELT(1);
    while (r1 > 1) {
	AG_Element q = r2;
	q /= r1;
	AG_Element r = r2;
	r %= r1;
	AG_Element x = p;
	x -= q;
	AG_Element y = x;
	x *= x1;
	x += x2;
	x %= p;
	y *= y1;
	y += y2;
	y %= p;
	r2 = r1;
	r1 = r;
	x2 = x1;
	x1 = x;
	y2 = y1;
	y1 = y;
    }
    inv = 0;
    inv = y1;
    return true;
}

// Do row operations on M to discover if M is invertible.  Return true if M is
// invertible and false othersize.
// If aux and/or aux2 are specified, then the row operations will also be
// performed on them.  If M is invertible, then aux becomes M^(-1)*aux and aux2
// becomes M^(-1)*aux2.
// NOTE:
//  - All row operations are in-place and change the contents of M, aux and
//    aux2.  Use copies of these matrices if you will need to keep the
//    original.
//  - Assumes that size*(p-1)^2 fits in an AG_Element
//  - If M is discovered to be singular, the algorithm short-circuits and the
//    values in M, aux and aux2 are not guaranteed.
bool PercyAGClient::is_invertible (AG_Element * M, dbsize_t size, 
	AG_Element * aux, AG_Element * aux2)
{
    // For each column i, reduce to have only 0s, except 1 at position i
    for (dbsize_t i = 0; i < size; ++i) {
	dbsize_t swaprow = i;
	// Find row for pivot
	while (true) {
	    if (swaprow == size) {
		return false;
	    }
	    if ((M[size*swaprow+i] % p) != 0) {
		break;
	    }
	    ++swaprow;
	}

	// Swap rows if necessary
	if (swaprow != i) {
	    for (dbsize_t j = i; j < size; ++j) {
		swap(M[swaprow*size+j], M[i*size+j]);
	    }
	    if (aux) {
		for (dbsize_t j = 0; j < size; ++j) {
		    swap(aux[swaprow*size+j], aux[i*size+j]);
		}
	    }
	    if (aux2) {
		for (dbsize_t j = 0; j < size; ++j) {
		    swap(aux2[swaprow*size+j], aux2[i*size+j]);
		}
	    }
	}

	AG_Element * row_i = M + i * size;
	AG_Element * aux_row_i = NULL;
	AG_Element * aux2_row_i = NULL;

	// Get inverse of pivot element and multiple pivot row by this
	AG_Element pivot_inv = ULONG_TO_AGELT(1);
	inv(pivot_inv, row_i[i] % p);

	for (dbsize_t j = i; j < size; ++j) {
	    row_i[j] %= p;
	    row_i[j] *= pivot_inv;
	    row_i[j] %= p;
	}
	if (aux) {
	    aux_row_i = aux + i * size;
	    for (dbsize_t j = 0; j < size; ++j) {
		aux_row_i[j] %= p;
		aux_row_i[j] *= pivot_inv;
		aux_row_i[j] %= p;
	    }
	}
	if (aux2) {
	    aux2_row_i = aux2 + i * size;
	    for (dbsize_t j = 0; j < size; ++j) {
		aux2_row_i[j] %= p;
		aux2_row_i[j] *= pivot_inv;
		aux2_row_i[j] %= p;
	    }
	}
	// NOTE: pivot row is now modulo p

	// For each (non-pivot) row, make the i-th element zero.
	for (dbsize_t j = 0; j < size; ++j) {
	    // Don't do for j-th row
	    if (j == i) {
		continue;
	    }
	    AG_Element * row_j = M + j * size;
	    // Not needed if M_j,i == 0 mod p
	    if ((row_j[i] % p) == 0) {
		continue;
	    }
	    AG_Element mult = p - (row_j[i] % p);
	    // Add mult * row_i to row_j
	    for (dbsize_t k = i; k < size; ++k) {
		row_j[k] += mult * row_i[k];
	    }
	    if (aux) {
		AG_Element * aux_row_j = aux + j * size;
		for (dbsize_t k = 0; k < size; ++k) {
		    aux_row_j[k] += mult * aux_row_i[k];
		}
	    }
	    if (aux2) {
		AG_Element * aux2_row_j = aux2 + j * size;
		for (dbsize_t k = 0; k < size; ++k) {
		    aux2_row_j[k] += mult * aux2_row_i[k];
		}
	    }
	}
    }

    // Do mod p on all elements
    for (dbsize_t i = 0; i < size * size; ++i) {
	M[i] %= p;
    }
    if (aux) {
	for (dbsize_t i = 0; i < size * size; ++i) {
	    aux[i] %= p;
	}
    }
    if (aux2) {
	for (dbsize_t i = 0; i < size * size; ++i) {
	    aux2[i] %= p;
	}
    }

    return true;
}

void print_element (std::ostream& os, AG_Element n)
{
    std::vector<int> rev_digits;
    if (n == 0) {
	rev_digits.push_back(0);
    } else {
	while (n > 0) {
	    rev_digits.push_back((int)(n % 10));
	    n /= 10;
	}
    }
    while (!(rev_digits.empty())) {
	os << rev_digits.back();
	rev_digits.pop_back();
    }
}

void print_matrix (std::ostream& os, AG_Element * M, dbsize_t rows, 
	dbsize_t cols)
{
    os << "[";
    dbsize_t k = 0;
    for (dbsize_t i = 0; i < rows; ++i) {
	os << "[";
	for (dbsize_t j = 0; j < cols; ++j, ++k) {
	    if (j > 0) {
		os << " ";
	    }
	    print_element(os, M[k]);
	}
	os << "]\n";
    }
    os << "]\n";
}

void to_mat_ZZ_p (mat_ZZ_p& dest, AG_Element * src, dbsize_t rows, dbsize_t cols)
{
    dest.SetDims(rows, cols);
    dbsize_t k = 0;
    for (dbsize_t i = 0; i < rows; ++i) {
	for (dbsize_t j = 0; j < cols; ++j, ++k) {
	    dest(i+1, j+1) = to_ZZ_p(ZZFromBytes((unsigned char*)(src+k), 
		    sizeof(AG_Element)));
	}
    }
}

void to_ag_matrix (AG_Element * dest, const mat_ZZ_p& src, dbsize_t rows, dbsize_t cols)
{
    dbsize_t k = 0;
    for (dbsize_t i = 0; i < rows; ++i) {
	for (dbsize_t j = 0; j < cols; ++j, ++k) {
	    BytesFromZZ((unsigned char*)(dest+k), rep(src[i][j]), sizeof(AG_Element));
	}
    }
}

#endif

void PercyAGClient::encode_request_impl (nqueries_t request_identifier)
{
    const std::vector<dbsize_t>& block_numbers = get_block_numbers(request_identifier);
    nqueries_t num_queries = block_numbers.size();
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t virtual_block_size = params->virtual_block_size();
    dbsize_t word_size = params->word_size();
    dbsize_t N = params->N();

#ifdef AG_NTL
    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();
    ZZ_p::init(p_ntl);

    // Define and size matrices
    mat_ZZ_p A, A_inv, B, Delta, P, P_inv, D, M1_prime, M2_prime;
    vec_ZZ_p Delta_vec, Delta_vec_inv;
    A.SetDims(N, N);
    B.SetDims(N, N);
    Delta.SetDims(N, N);
    Delta_vec.SetLength(N);
    Delta_vec_inv.SetLength(N);
    P.SetDims(N, N);
    P_inv.SetDims(N, N);
    D.SetDims(N, N);
    M1_prime.SetDims(N, N);
    M2_prime.SetDims(N, N);
#else
    // Define and size matrices
    AG_Element * A = new AG_Element[N*N];
    memset(A, 0, N*N*sizeof(AG_Element));
    AG_Element * B = new AG_Element[N*N];
    memset(B, 0, N*N*sizeof(AG_Element));
    AG_Element * A_copy = new AG_Element[N*N];
    memset(A_copy, 0, N*N*sizeof(AG_Element));
    AG_Element * P_inv = new AG_Element[N*N];
    memset(P_inv, 0, N*N*sizeof(AG_Element));
    AG_Element * PA = new AG_Element[N*N];
    memset(PA, 0, N*N*sizeof(AG_Element));
    AG_Element * PB = new AG_Element[N*N];
    memset(PB, 0, N*N*sizeof(AG_Element));
    AG_Element * D = new AG_Element[N*N];
    memset(D, 0, N*N*sizeof(AG_Element));

    ZZ_p::init(p_ntl);
#endif

    queries[request_identifier] = vector<unsigned char *>();
    sentinfo[request_identifier] = std::vector<AGDecodeInfo>();

    for (nqueries_t qind = 0; qind < num_queries; ++qind) {
	/*
	// Get the indices for each iterator
	std::vector<dbsize_t> iter_indices;
	dbsize_t r = block_numbers[qind];
	for (nqueries_t d = 0; d < depth; ++d) {
	    iter_indices.push_back(r / params->iteration_subset_size(d));
	    r = r % params->iteration_subset_size(d);
	}

	// Add to sentinfo vector
	sentinfo[request_identifier].push_back(vector<AGDecodeInfo>());
	std::vector<AGDecodeInfo>& query_sentinfo = sentinfo[request_identifier].back();
	*/

	//for (nqueries_t d = 0; d < depth; ++d) {
	    // Get random invertible A and random B;
#ifdef AG_NTL
	    random_matrix(B);
	    ZZ_p det = ZZ_p::zero();
	    while (IsZero(det)){
		random_matrix(A);
		inv(det, A_inv, A);
	    }
#else
	    AG_Element * A_inv_B = new AG_Element[N*N];
	    memset(A_inv_B, 0, N*N*sizeof(AG_Element));
	    random_ag_elements(B, N*N);
	    bool keep = false;
	    while (!keep) {
		random_ag_elements(A, N*N);
		memcpy(A_copy, A, N*N*sizeof(AG_Element));
		memcpy(A_inv_B, B, N*N*sizeof(AG_Element));
		keep = is_invertible(A_copy, N, A_inv_B);
	    }
#endif

	    // Get random invertible diagonal matrix Delta
#ifdef AG_NTL
	    for (dbsize_t i = 0; i < N; ++i) {
		ZZ_p current = ZZ_p::zero();
		while (IsZero(current)) {
		    random(current);
		}
		Delta(i+1, i+1) = current;
		Delta_vec[i] = current;
	    }
#else
	    AG_Element * Delta_vec = new AG_Element[N];
	    AG_Element * Delta_vec_inv = new AG_Element[N];
	    memset(Delta_vec, 0, N*sizeof(AG_Element));
	    memset(Delta_vec_inv, 0, N*sizeof(AG_Element));
	    random_ag_elements(Delta_vec, N, true);
#endif
	    for (dbsize_t i = 0; i < N; ++i) {
		inv(Delta_vec_inv[i], Delta_vec[i]);
	    }

	    // Create the random column permutation
	    std::vector<dbsize_t> col_inds;
	    for (dbsize_t i = 0; i < 2*N; ++i) {
		col_inds.push_back(i);
	    }
	    std::vector<dbsize_t> perm;
	    for (dbsize_t i = 0; i < 2*N; ++i) {
		dbsize_t rnd = RandomBnd(2*N - i);
		dbsize_t col = col_inds[rnd];
		perm.push_back(col);
		col_inds.erase(col_inds.begin() + rnd);
	    }

	    // Get block number
	    dbsize_t i0 = block_numbers[qind] / virtual_block_size;

	    for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
		// Get random invertible P
#ifdef AG_NTL
		det = ZZ_p::zero();
		while (IsZero(det)) {
		    random_matrix(P_inv);
		    determinant(det, P);
		    inv(det, P, P_inv);
		}
#else
		keep = false;
		while (!keep) {
		    random_ag_elements(P_inv, N*N);
		    memcpy(PA, A, N*N*sizeof(AG_Element));
		    memcpy(PB, B, N*N*sizeof(AG_Element));
		    keep = is_invertible(P_inv, N, PA, PB);
		}
#endif

		// Generate noise matrix
		if (i == i0) {
#ifdef AG_NTL
		    hard_noise_matrix(D, to_ZZ_p(q_ntl));
#else
		    hard_noise_matrix(D, N);
#endif
		} else {
#ifdef AG_NTL
		    soft_noise_matrix(D);
#else
		    soft_noise_matrix(D, N);
#endif
		}

#ifdef AG_NTL
		// Compute M' = [ M1'=P*A | M2'=P*B+D*Delta ]
		M1_prime = P * A;
		M2_prime = (P * B) + (D * Delta);
#else
		// Comput PB <- PB + D*Delta
		for (dbsize_t k = 0; k < N*N; ++k) {
		    D[k] *= Delta_vec[k % N];
		    PB[k] += D[k];
		    PB[k] %= p;
		}
#endif

		// Permute the columns of M' by perm to get M and add it to the
		// send_buffer
		unsigned char * send_buffer = new unsigned char[6 * N * N * word_size / 8];
		memset(send_buffer, 0, 6*N*N*word_size/8);
		bufferloc_t location;
		for (dbsize_t k = 0; k < N; ++k) {
#ifdef AG_NTL
		    vec_ZZ_p &row1 = M1_prime[k];
		    vec_ZZ_p &row2 = M2_prime[k];
#else
		    AG_Element * row1 = PA + k * N;
		    AG_Element * row2 = PB + k * N;
#endif
		    for (dbsize_t j = 0; j < 2*N; ++j) {
			dbsize_t col = perm[j];
			if (col < N) {
			    location = add_to_buffer(send_buffer, location,
				    row1[col], word_size * 3);
			} else {
			    location = add_to_buffer(send_buffer, location,
				    row2[col-N], word_size * 3);
			}
		    }
		}
		queries[request_identifier].push_back(send_buffer);
	    }

#ifndef AG_NTL
	    delete[] Delta_vec;
#endif

	    // Add info needed for decoding to sent_* vectors
	    AGDecodeInfo decodeinfo;
#ifdef AG_NTL
	    decodeinfo.A_inv_B = A_inv * B;
#else
	    decodeinfo.A_inv_B = A_inv_B;
#endif
	    decodeinfo.Delta_vec_inv = Delta_vec_inv;
	    decodeinfo.perm = perm;
	    sentinfo[request_identifier].push_back(decodeinfo);
	//}
    }

#ifdef AG_NTL
    savectx.restore();
#else
    delete[] A;
    delete[] B;
    delete[] A_copy;
    delete[] P_inv;
    delete[] PA;
    delete[] PB;
    delete[] D;
#endif
}

dbsize_t PercyAGClient::send_request_impl (nqueries_t request_identifier,
	std::vector<ostream*> &osvec, bool send_num_queries)
{
    nservers_t num_servers = osvec.size();
    dbsize_t word_size = params->word_size();
    dbsize_t N = params->N();

    // Send the number of queries
    nqueries_t num_queries = get_block_numbers(request_identifier).size();
    if (send_num_queries) {
	for (nservers_t k = 0; k < num_servers; ++k) {
	    percy_write_le_uint16(*(osvec[k]), num_queries);
	}
    }

    // Send it!
    vector<unsigned char *>& this_queries = queries[request_identifier];
    dbsize_t num_matrices = this_queries.size();
    dbsize_t matrix_size = 6 * N * N * word_size / 8;
    for (dbsize_t i = 0; i < num_matrices; ++i) {
	for (nservers_t k = 0; k < num_servers; ++k) {
	    osvec[k]->write((char*)(this_queries[i]), matrix_size);
	}
	delete[] this_queries[i];
    }

    // Flush 'em
    for (nservers_t k = 0; k < num_servers; ++k) {
	osvec[k]->flush();
    }

    queries.erase(request_identifier);
    return num_servers * num_matrices * matrix_size;
}

dbsize_t PercyAGClient::receive_replies_impl (nqueries_t request_identifier,
	std::vector<istream*> &isvec)
{
#ifdef AG_NTL
    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();
    ZZ_p::init(p_ntl);
#endif

    // Get parameters
    nqueries_t num_queries = get_block_numbers(request_identifier).size();
    dbsize_t block_size = params->block_size();
    dbsize_t word_size = params->word_size();
    dbsize_t N = params->N();
    dbsize_t block_row_size = N * word_size / 8;
    dbsize_t block_num_rows = (block_size - 1) / block_row_size + 1;
    dbsize_t readsize = 6 * block_num_rows * block_row_size;
    dbsize_t size = 2 * block_num_rows * N;

    nservers_t num_servers = isvec.size();
    goodservers.clear();

    // Read answers
    unsigned char * receive_buffer = new unsigned char[readsize];
#ifdef AG_NTL
    map<nservers_t, vec_vec_ZZ_p>& answers = stored_answers[request_identifier];
#else
    map<nservers_t, vector<AG_Element*> >& answers = stored_answers[request_identifier];
#endif
    for (nservers_t k = 0; k < num_servers; ++k) {
	// Add answer matrix
#ifdef AG_NTL
	vec_vec_ZZ_p& sanswers = answers[k];
	sanswers.SetLength(num_queries);
#else
	vector<AG_Element*>& sanswers = answers[k];
#endif
	bool isgood = true;
	for (nqueries_t qind = 0; qind < num_queries; ++qind) {
#ifdef AG_NTL
	    vec_ZZ_p &V = sanswers[qind];
	    V.SetLength(size);
#else
	    AG_Element * V = new AG_Element[size];
	    memset(V, 0, size * sizeof(AG_Element));
	    sanswers.push_back(V);
#endif
	
	    // Read answer matrix
	    isvec[k]->read((char*)receive_buffer, readsize);
	    if (isvec[k]->eof() || isvec[k]->gcount() < 1) {
		std::cerr << "Error reading response from server " << k << "\n";
		isgood = false;
		break;
	    }
	    bufferloc_t location;
	    for (dbsize_t i = 0; i < size; ++i) {
		location = get_from_buffer(receive_buffer, location, V[i],
			word_size * 3);
	    }
	}
	if (isgood) {
	    goodservers.push_back(k);
	} else {
#ifndef AG_NTL
	    for (nqueries_t i = 0; i < sanswers.size(); ++i) {
		delete[] sanswers[i];
	    }
#endif
	    answers.erase(k);
	}
    }
    delete[] receive_buffer;

#ifdef AG_NTL
    savectx.restore();
#endif
    return goodservers.size() * num_queries * readsize;
}

nqueries_t PercyAGClient::process_replies_impl (nservers_t h, 
	std::vector<std::vector<PercyResult> > &results)
{
    // Get parameters
    dbsize_t block_size = params->block_size();
    dbsize_t word_size = params->word_size();
    dbsize_t N = params->N();

#ifdef AG_NTL
    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();
    ZZ_p::init(p_ntl);

    // Word mask
    ZZ mask = power2_ZZ(word_size) - 1;
#else
    AG_Element mask = ULONG_TO_AGELT(1);
    mask <<= word_size;
    mask -= 1;

    AG_Element q_inv;
    inv(q_inv, q);
#endif

    // Get matrix height
    dbsize_t height = 0;
    if (word_size == 16) {
	height = (block_size - 1) / (2 * N) + 1;
    } else if (word_size == 20) {
	height = (block_size - 1) / (5 * N / 2) + 1;
    } else if (word_size == 24) {
	height = (block_size - 1) / (3 * N) + 1;
    } else {
	return stored_answers.size();
    }

#ifdef AG_NTL
    map<nqueries_t, map<nservers_t, vec_vec_ZZ_p> >::iterator iter;
#else
    map<nqueries_t, map<nservers_t, vector<AG_Element*> > >::iterator iter;
#endif
    for (iter = stored_answers.begin(); iter != stored_answers.end(); ++iter) {
	nqueries_t request_identifier = iter->first;
	nqueries_t num_queries = get_block_numbers(request_identifier).size();
#ifdef AG_NTL
	map<nservers_t, vec_vec_ZZ_p>& answers = iter->second;
#else
	map<nservers_t, vector<AG_Element*> >& answers = iter->second;
#endif
	nservers_t num_good_servers = goodservers.size();
	std::vector<AGDecodeInfo>& recvinfo = sentinfo[request_identifier];
	for (dbsize_t qind = 0; qind < num_queries; ++qind) {
	    results.push_back(vector<PercyResult>());
	    vector<PercyResult>& qresult = results.back();

	    // Get stored information
#ifdef AG_NTL
	    mat_ZZ_p &A_inv_B = recvinfo[qind].A_inv_B;
#else
	    AG_Element * A_inv_B = recvinfo[qind].A_inv_B;
#endif
	    std::vector<dbsize_t> &perm = recvinfo[qind].perm;

#ifdef AG_NTL
	    vec_ZZ_p &Delta_inv_vec = recvinfo[qind].Delta_vec_inv;
#else
	    AG_Element * Delta_inv_vec = recvinfo[qind].Delta_vec_inv;
#endif

	    for (nservers_t gserv = 0; gserv < num_good_servers; ++gserv) {
		nservers_t serv = goodservers[gserv];
#ifdef AG_NTL
		vec_ZZ_p &sanswers = answers[serv][qind];
#else
		AG_Element * sanswers = answers[serv][qind];
#endif

		// Get U and D
#ifdef AG_NTL
		mat_ZZ_p U, D;
		U.SetDims(height, N);
		D.SetDims(height, N);
#else
		AG_Element * U = new AG_Element[height * N];
		AG_Element * E = new AG_Element[height * N];
#endif

		dbsize_t ind = 0;
		for (dbsize_t i = 0; i < height; ++i) {
#ifdef AG_NTL
		    vec_ZZ_p &rowU = U[i];
		    vec_ZZ_p &rowD = D[i];
#else
		    AG_Element * rowU = U + i * N;
		    AG_Element * rowD = E + i * N;
#endif
		    for (dbsize_t j = 0; j < 2*N; ++j, ++ind) {
			dbsize_t col = perm[j];
			if (col < N) {
			    rowU[col] = sanswers[ind];
			} else {
			    rowD[col-N] = sanswers[ind];
			}
		    }
		}

		// Find scramble noise
#ifdef AG_NTL
		mat_ZZ_p E = D - (U * A_inv_B);
#else
		AG_Element * U_A_inv_B = new AG_Element[height * N];
		mult_matrices(U_A_inv_B, U, A_inv_B, height, N, N);
		for (dbsize_t i = 0; i < height * N; ++i) {
		    E[i] += p - U_A_inv_B[i];
		}
		delete[] U_A_inv_B;
#endif

		// Unscramble noise
		for (dbsize_t i = 0; i < height; ++i) {
#ifdef AG_NTL
		    vec_ZZ_p &row = E[i];
#else
		    AG_Element * row = E + i * N;
#endif
		    for (dbsize_t j = 0; j < N; ++j) {
			row[j] *= Delta_inv_vec[j];
#ifndef AG_NTL
			row[j] %= p;
#endif
		    }
		}

		// Get record (as ZZ vector)
#ifdef AG_NTL
		vec_ZZ_p sresult;
		sresult.SetLength(height * N);
#else
		AG_Element * sresult = new AG_Element[height * N];
#endif
		ind = 0;
		for (dbsize_t i = 0; i < height; ++i) {
		    for (dbsize_t j = 0; j < N; ++j, ++ind) {
#ifdef AG_NTL
			ZZ cell (rep(E(i+1, j+1)));
			ZZ epsilon = cell % q_ntl;
			if (!(epsilon < q_ntl/2)) {
			    epsilon -= q_ntl;
			}
			sresult[ind] = to_ZZ_p(((cell - epsilon) / q_ntl) & mask);
#else
			AG_Element& cell = E[ind];
			AG_Element epsilon = cell;
			epsilon %= q;
			sresult[ind] = cell;
			sresult[ind] >>= (2 * word_size - 1);
			if (epsilon > (q >> 1)) {
			    sresult[ind] += 1;
			}
			sresult[ind] &= mask;
#endif
		    }
		}
#ifndef AG_NTL
		delete[] U;
		delete[] E;
		delete[] sanswers;
#endif

		// Convert to bytes
		unsigned char * result_char = NULL;
		if (word_size == 16) {
		    result_char = new unsigned char[height * 2 * N];
		    memset(result_char, 0, height * 2 * N);
		} else if (word_size == 20) {
		    result_char = new unsigned char[height * 5 * N / 2];
		    memset(result_char, 0, height * 5 * N / 2);
		} else if (word_size == 24) {
		    result_char = new unsigned char[height * 3 * N];
		    memset(result_char, 0, height * 3 * N);
		}
		bufferloc_t location;
		for (dbsize_t i = 0; i < height * N; ++i) {
		    location = add_to_buffer(result_char, location, sresult[i],
			    word_size);
		}

		// Create result object
		qresult.push_back(PercyResult(std::vector<nservers_t>(1, serv), 
			std::string((char*)result_char, block_size)));

		delete[] result_char;
#ifndef AG_NTL
		delete[] sresult;
#endif
	    }

#ifndef AG_NTL
	    delete[] A_inv_B;
	    delete[] Delta_inv_vec;
#endif
	}
    }

    // Remove blocks from received_blocks
    stored_answers.clear();
    sentinfo.clear();
#ifdef AG_NTL
    savectx.restore();
#endif
    return 0;
}

