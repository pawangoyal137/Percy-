// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>
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

#ifndef __ITCLIENT_IMPL_H__
#define __ITCLIENT_IMPL_H__

#include <iostream>
#include <algorithm>
#include <utility>
#include "gf2e.h"
#include "itparams.h"
#include "percyresult.h"
#include "rsdecoder.h"

template <typename GF2E_Element>
inline void setCoeffs(GF2X &GF2X_P);

template <>
inline void setCoeffs<GF28_Element>(GF2X &GF2X_P) {
    SetCoeff(GF2X_P, 8, 1);
    SetCoeff(GF2X_P, 4, 1);
    SetCoeff(GF2X_P, 3, 1);
    SetCoeff(GF2X_P, 1, 1);
    SetCoeff(GF2X_P, 0, 1);
}

template <>
inline void setCoeffs<GF216_Element>(GF2X &GF2X_P) {
    SetCoeff(GF2X_P, 16, 1);
    SetCoeff(GF2X_P, 14, 1);
    SetCoeff(GF2X_P, 12, 1);
    SetCoeff(GF2X_P, 7, 1);
    SetCoeff(GF2X_P, 6, 1);
    SetCoeff(GF2X_P, 4, 1);
    SetCoeff(GF2X_P, 2, 1);
    SetCoeff(GF2X_P, 1, 1);
    SetCoeff(GF2X_P, 0, 1);
}

// Constructor and destructor
template <typename GF2E_Element>
PercyClient_GF2E<GF2E_Element>::PercyClient_GF2E (const PercyClientParams * params,
	nservers_t num_servers, nservers_t t, sid_t * sids, PercyStats * stats) :
    PercyClient(params, num_servers, t, stats),
    params(static_cast<const GF2EParams*>(params->percy_params())),
    indices(),
    interp_indices(),
    vec_interp_indices(),
    vecs_interp_indices(),
    randmults(),
    answers(),
    unfinished_results(),
    decoded()
{
    // Initialize the GF2E modulus
    GF2X MODULUS_P;
    setCoeffs<GF2E_Element>(MODULUS_P);

    GF2E::init(MODULUS_P);
    GF2X::HexOutput = 1;

    choose_indices(sids);
    nqueries_t max_qbs = num_servers - 1 - this->params->tau();
    choose_interp_indices(max_qbs);
}

template <typename GF2E_Element>
PercyClient_GF2E<GF2E_Element>::~PercyClient_GF2E ()
{
    if (indices != NULL) {
	delete[] indices;
    }
    if (interp_indices != NULL) {
	delete[] interp_indices;
    }
    while (answers.size() > 0) {
	if (answers.back() != NULL) {
	    delete[] answers.back();
	}
	answers.pop_back();
    }
    typename map<nqueries_t, GF2E_Element *>::iterator siter;
    for (siter = stored_shares.begin(); siter != stored_shares.end(); ++siter) {
	delete[] siter->second;
    }
    stored_shares.clear();
}


// Generate t-private (t+1)-of-l shares of a given secret value.
template <typename GF2E_Element>
static void genshares_GF2E(nservers_t t, nservers_t l,
    const GF2E_Element *indices, GF2E_Element *values, GF2E_Element secret)
{
    // Pick a random polynomial of degree t with the right constant term
    GF2E_Element * coeffs = new GF2E_Element[t+1];
    coeffs[0] = secret;
    for (nservers_t i=1;i<=t;++i) {
        coeffs[i] = RandomBits_ulong(8*sizeof(GF2E_Element));
    }

    // Evaluate the polynomial at each of the indices
    for (nservers_t i=0; i<l; ++i) {
        values[i] = evalpoly_GF2E<GF2E_Element>(coeffs, t, indices[i]);
    }
    delete[] coeffs;
}

#ifdef VERBOSE_ITCLIENT
template <typename GF2E_Element>
static void showpoly(const char *label, const GF2E_Element *poly,
			nqueries_t numcoeffs)
{
    if (label) {
	fprintf(stderr, "%s", label);
    }
    fprintf(stderr, "[ ");
    for (nqueries_t i = 0; i < numcoeffs; ++i) {
	fprintf(stderr, "%02x ", poly[i]);
    }
    fprintf(stderr, "]\n");
}
#endif

// Generate t-private (t+1)-of-l shares of a given secret value.
template <typename GF2E_Element>
static void gensharesmulti_GF2E(nservers_t t, nservers_t l,
    const GF2E_Element *indices, GF2E_Element *values,
    const GF2E_Element* secret, const GF2E_Element *interp_indices,
    nqueries_t qbs)
{
    // Special case for t=0 so we don't have to worry about polynomials
    // with degree < 0
    if (t==0) {
	for (nservers_t i=0; i<l; ++i) {
	    values[i] = interpolate_GF2E(interp_indices, secret, qbs,
					    indices[i]);
	}
	return;
    }

    // Pick a random polynomial of degree t-1
    GF2E_Element * randpoly = new GF2E_Element[t];
    for (nservers_t i=0; i<t; ++i) {
        randpoly[i] = RandomBits_ulong(8*sizeof(GF2E_Element));
    }
    
    GF2E_Element * buildroots = new GF2E_Element[qbs+1];
    polyfromroots_GF2E(buildroots, interp_indices, qbs);
    GF2E_Element * final = new GF2E_Element[t+qbs];
    multpoly_GF2E(final, randpoly, t-1, buildroots, qbs);

    // Evaluate the polynomial at each of the indices
    for (nservers_t i=0; i<l; ++i) {
        values[i] = evalpoly_GF2E<GF2E_Element>(final, t+qbs-1, indices[i]) ^
			interpolate_GF2E(interp_indices, secret, qbs,
					    indices[i]);
    }
    delete[] randpoly;
    delete[] buildroots;
    delete[] final;
}

template <typename GF2E_Element>
void PercyClient_GF2E<GF2E_Element>::choose_indices(sid_t *sids) {
    indices = new GF2E_Element[num_servers];
    for (nservers_t j=0; j<num_servers; ++j) {
        if (params->tau()) {
            // Use the indices provided
            indices[j] = (GF2E_Element)sids[j];
        } else {
            // Use random indices
            GF2E_Element r;
            bool ok = false;
            do {
                r = RandomLen_long(8*sizeof(GF2E_Element));
                if (r != 0) {
                    ok = true;
                    for (nservers_t k=0;k<j;++k) {
                        if (indices[k] == r) {
                            ok = false;
                        }
                    }
                }
            } while (!ok);
            indices[j] = r;
        }
    }

    // Create NTL version
    indices_ntl.SetLength(num_servers);
    for (nservers_t ix=0; ix<num_servers; ++ix) {
	conv(indices_ntl[ix],
		GF2XFromBytes((unsigned char *)(indices + ix),
		    sizeof(GF2E_Element)));
    }
}


template <typename GF2E_Element>
void PercyClient_GF2E<GF2E_Element>::choose_interp_indices(nqueries_t qbs) {
    vec_interp_indices.SetLength(qbs);
    interp_indices = new GF2E_Element[qbs];
    GF2E_Element candidate = 0x00;
    nqueries_t i = 0;
    while (i < qbs){
	bool in_server_indices = false;
        for (nservers_t s=0; s<num_servers; s++){
	    if (indices[s] == candidate){
		in_server_indices = true;
	    }
	}	
	if (!in_server_indices) {
            interp_indices[i] = candidate;
	    conv(vec_interp_indices[i],
		    GF2XFromBytes((unsigned char *)(interp_indices + i),
			    sizeof(GF2E_Element)));
	    ++i;
   	}
	--candidate;
    }
}

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
template <typename GF2E_Element>
void PercyClient_GF2E<GF2E_Element>::encode_request_impl (
	nqueries_t request_identifier)
{
    const vector<dbsize_t>& block_numbers = get_block_numbers(request_identifier);
    nqueries_t qbs = get_qbs(request_identifier);
    nqueries_t num_block_requests = block_numbers.size();
    nqueries_t num_queries = (num_block_requests+qbs-1)/qbs;
    // Generate random multiples (!= 0)
    if (randomize) {
	randmults[request_identifier] = vector<GF2E_Element>();
	vector<GF2E_Element>& request_randmults = randmults[request_identifier];
        for (nqueries_t i = 0; i < num_queries * num_servers; ++i) {
	    GF2E_Element r = 0;
	    while (r == 0) {
		r = RandomLen_long(8*sizeof(GF2E_Element));
	    }
	    request_randmults.push_back(r);
        }
    }

    // Construct the shares of the e_{index} vector
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t virtual_block_size = params->virtual_block_size();
    GF2E_Element * shares = new GF2E_Element[num_queries * num_virtual_blocks * num_servers];
    nqueries_t current_query = 0;
    GF2E_Element* secrets = new GF2E_Element[qbs];
    for (nqueries_t q=0; q<num_block_requests; q += qbs) {
	nqueries_t tmpblock = (num_block_requests - q) > qbs ?
				qbs : (num_block_requests - q);
	vec_GF2E interps = vec_interp_indices;
	interps.SetLength(tmpblock);
	vecs_interp_indices.push_back(interps);
        for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
            if (tmpblock == 1) {
                dbsize_t virtual_block_number = block_numbers[q] / virtual_block_size;
                genshares_GF2E<GF2E_Element>(t, num_servers, indices,
                        shares + (current_query * num_virtual_blocks + i) * num_servers, 
                        i == virtual_block_number);
            } else {
                // create GF2E_Element
                for (nqueries_t k = 0; k < tmpblock; k++) {
                    dbsize_t virtual_block_number = block_numbers[q+k] / virtual_block_size;
                    if (virtual_block_number == i) {
                        secrets[k] = (GF2E_Element)(1);
                    } else {
                        secrets[k] = (GF2E_Element)(0);
                    }
                }
                gensharesmulti_GF2E<GF2E_Element>(t, num_servers, indices,
		      shares + (current_query * num_virtual_blocks + i) * num_servers,
		      secrets, interp_indices, tmpblock);
            }
        }
        current_query += 1;
    }
    delete[] secrets;

    // Multiply shares by random multiples
    if (randomize) {
	vector<GF2E_Element>& request_randmults = randmults[request_identifier];
        for (nservers_t p = 0; p < num_servers; ++p) {
            for (nqueries_t i = 0; i < num_queries; ++i) {
                for (dbsize_t j = 0; j < num_virtual_blocks; ++j) {
		    dbsize_t index = (i * num_virtual_blocks + j) * num_servers
			    + p;
		    shares[index] = 
                        multiply_GF2E<GF2E_Element>(shares[index],
                                request_randmults[p*num_queries + i]);
                }
            }
        }
    }

    // Stored shares
    stored_shares[request_identifier] = shares;
}

template <typename GF2E_Element>
dbsize_t PercyClient_GF2E<GF2E_Element>::send_request_impl (
	nqueries_t request_identifier, vector<ostream*> &osvec,
	bool send_num_queries)
{
    // Get some parameters and stored values
    GF2E_Element * shares = stored_shares[request_identifier];
    nqueries_t qbs = get_qbs(request_identifier);
    nqueries_t num_block_requests = get_block_numbers(request_identifier).size();
    nqueries_t num_queries = (num_block_requests+qbs-1)/qbs;
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();

    // Check number of streams
    if (num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request_GF2E.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << num_servers
            << ".\n";
        return 0;
    }

    // Send the params and query to each server
    for (nservers_t j = 0; j < num_servers; ++j) {
	if (send_num_queries) {
	    percy_write_le_uint16(*osvec[j], num_queries);
	}
        for (nqueries_t q=0; q<num_queries; ++q) {
            for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
		dbsize_t index = (q * num_virtual_blocks + i) * num_servers + j;
                char *shareptr = (char *)&(shares[index]);
                osvec[j]->write(shareptr, sizeof(GF2E_Element));
            }
        }
        osvec[j]->flush();
    }

    // Clean up
    delete[] shares;
    stored_shares.erase(request_identifier);

    return num_servers * num_queries * num_virtual_blocks * sizeof(GF2E_Element);
}

// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
template <typename GF2E_Element>
dbsize_t PercyClient_GF2E<GF2E_Element>::receive_replies_impl (
	nqueries_t request_identifier, vector<istream*> &isvec)
{
    dbsize_t words_per_block = params->words_per_block();
    nqueries_t qbs = get_qbs(request_identifier);
    nqueries_t num_block_requests = get_block_numbers(request_identifier).size();
    nqueries_t num_queries = (num_block_requests+qbs-1)/qbs;
    nqueries_t q;

    // The vector of servers that have responded properly
    goodservers.clear();

    // The responses from the servers
    nqueries_t prev = answers.size();
    for (nqueries_t q = 0; q < num_queries; ++q) {
	answers.push_back(new GF2E_Element[words_per_block * num_servers]);
    }

    // Read the replies
    for (nservers_t j = 0; j < num_servers; ++j) {
	bool isgood = true;
	for (nqueries_t q = prev; q < prev + num_queries; ++q) {
	    for (dbsize_t i = 0; isgood && i < words_per_block; ++i) {
		isvec[j]->read((char *)(answers[q] + i * num_servers + j),
			sizeof(GF2E_Element));
		if (isvec[j]->eof()) {
		    std::cerr << "Server " << j+1 << " did not send complete reply.\n";
		    std::cerr << "Marking server " << j+1 << " as bad.\n";
		    isgood = false;
		    break;
		}
		if ((dbsize_t)(isvec[j]->gcount()) < 1) {
		    // Mark this server as bad
		    std::cerr << "Marking server " << j+1 << " as bad.\n";
		    isgood = false;
		    break;
		}
	    }
	}
	if (isgood) {
	    goodservers.push_back(j);
	}
    }

    // Remove random multiple
    if (randomize) {
	vector<GF2E_Element>& request_randmults = randmults[request_identifier];
        for (nservers_t j = 0; j < num_servers; ++j) {
            for (nqueries_t q = 0; q < num_queries; ++q) {
                GF2E_Element randmult_inv =
                        inverse_GF2E<GF2E_Element>(request_randmults[j*num_queries + q]);
                for (dbsize_t i = 0; i < words_per_block; ++i) {
                    answers[q+prev][i*num_servers + j] = multiply_GF2E<GF2E_Element>(
			    answers[q+prev][i*num_servers + j], randmult_inv);
                }
            }
        }
	randmults.erase(request_identifier);
    }

    // Add to unfinished_results and decoded
    for (q=0; q<num_queries; ++q) {
	unfinished_results.push_back(
		vector<DecoderResult<GF2E> >(1, DecoderResult<GF2E>(goodservers, vector<map<dbsize_t, GF2E> >())));
	decoded.push_back(std::set<dbsize_t>());
    }

    return num_servers * num_queries * words_per_block * sizeof(GF2E_Element);
}

// Construct the numpoints Lagrange coefficents to interpolate at the
// point alpha from the source points
// indices[goodservers[firstpoint .. firstpoint+numpoints-1]]
template <typename GF2E_Element>
void PercyClient_GF2E<GF2E_Element>::construct_lagrange_coeffs(GF2E_Element *coeffs,
	GF2E_Element alpha, nservers_t firstpoint, nservers_t numpoints)
{
    for (nservers_t i=0; i < numpoints; ++i) {
        GF2E_Element numer = 1, denom = 1;
        for (nservers_t j=0; j < numpoints; ++j) {
            if (j==i) continue;
            GF2E_Element numerdiff = indices[goodservers[firstpoint+j]] ^
	    		             alpha;
            GF2E_Element denomdiff = indices[goodservers[firstpoint+j]] ^
				     indices[goodservers[firstpoint+i]];
            numer = multiply_GF2E<GF2E_Element>(numer, numerdiff);
            denom = multiply_GF2E<GF2E_Element>(denom, denomdiff);
        }
        coeffs[i] = multiply_GF2E<GF2E_Element>(numer,
			    inverse_GF2E<GF2E_Element>(denom));
    }
}

// Do Lagrange interpolation of source values
// answers[goodservers[firstpoint..firstpoint+numpoints-1]] with the
// coefficients coeffs[0..t]
template <typename GF2E_Element>
inline GF2E_Element PercyClient_GF2E<GF2E_Element>::interpolate(
	const GF2E_Element *word_answers, const GF2E_Element *coeffs, 
	nservers_t firstpoint, nservers_t numpoints)
{
    GF2E_Element res = 0;

    for (nservers_t i = 0; i < numpoints; ++i) {
	res ^= multiply_GF2E<GF2E_Element>(coeffs[i],
		word_answers[goodservers[firstpoint+i]]);
    }

    return res;
}

template <typename GF2E_Element>
bool PercyClient_GF2E<GF2E_Element>::try_fast_recover(nservers_t h,
	vector<vector<PercyResult> > &results)
{
    dbsize_t words_per_block = params->words_per_block();
    nqueries_t num_queries = vecs_interp_indices.size();
    nservers_t k = goodservers.size();

    if (k <= t) return false;

    results.clear();
    for (nqueries_t q = 0; q < num_queries; ++q) {

	// blocks per query; the polynomial is of degree t+bpq-1
	nqueries_t bpq = vecs_interp_indices[q].length();

	// bpq == 0 is the old case of just retrieving one block,
	// interpolated at 0
	const bool justzero = (bpq == 0);
	if (justzero) {
	    bpq = 1;
	}
	GF2E_Element *resblock = new GF2E_Element[words_per_block * bpq];

	// Construct the Lagrange coefficients for the interp indices as
	// well as the target indices[goodservers[0..(k-t-bpq-1)]] from
	// the source indices[goodservers[k-t-bpq .. k-1]]

	// result_coeffs[a] is the array of t+bpq coeffs to interpolate
	//                  interp_indices[a]
	GF2E_Element *result_coeffs = new GF2E_Element[bpq*(t+bpq)];
	for (nservers_t a = 0; a < bpq; ++a) {
	    construct_lagrange_coeffs(result_coeffs + a*(t+bpq),
					justzero ? 0 : interp_indices[a],
					k-t-bpq, t+bpq);
	}

	// check_coeffs[a] is the array of t+bpq coeffs to interpolate
	//                 indices[goodservers[a]]
	GF2E_Element *check_coeffs = new GF2E_Element[(k-t-bpq)*(t+bpq)];
	for (nservers_t a = 0; a < k-t-bpq; ++a) {
	    construct_lagrange_coeffs(check_coeffs + a*(t+bpq),
					indices[goodservers[a]],
					k-t-bpq, t+bpq);
	}

	for (dbsize_t word = 0; word < words_per_block; ++word) {
	    // Check if all the servers agree
	    bool match = true;
	    for (nservers_t a = 0; a < k-t-bpq; ++a) {
		GF2E_Element expectedans =
			answers[q][word*num_servers + goodservers[a]];
		GF2E_Element actualans = interpolate(
			answers[q] + word*num_servers,
			check_coeffs + a*(t+bpq), k-t-bpq, t+bpq);

		if (expectedans != actualans) {
		    match = false;
		    break;
		}
	    }
	    if (!match) {
		delete[] check_coeffs;
		delete[] result_coeffs;
		delete[] resblock;
		results.clear();
		return false;
	    }
	    for (nservers_t a = 0; a < bpq; ++a) {
		resblock[a*words_per_block+word] = interpolate(
			answers[q] + word*num_servers, 
			result_coeffs + a*(t+bpq), k-t-bpq, t+bpq);
	    }
	}
	for (nservers_t a = 0; a < bpq; ++a) {
	    PercyResult thisresult(goodservers,
		string((char *)(resblock + a*words_per_block), 
		    words_per_block*sizeof(GF2E_Element)));
	    results.push_back(vector<PercyResult>(1, thisresult));
	}

	delete[] result_coeffs;
	delete[] check_coeffs;
	delete[] resblock;
    }

    return true;
}

//Process the received replies and return the decoded results as a
//PercyResults object.
template <typename GF2E_Element>
nqueries_t PercyClient_GF2E<GF2E_Element>::process_replies_impl (nservers_t h,
	vector<vector<PercyResult> >& results)
{
    dbsize_t words_per_block = params->words_per_block();
    nservers_t tau = params->tau();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return false;
    }

    if (try_fast_recover(h, results)) {
	// Clean up
	while (answers.size() > 0) {
	    if (answers.back() != NULL) {
		delete[] answers.back();
	    }
	    answers.pop_back();
	}
	answers_ntl.clear();
	unfinished_results.clear();
	decoded.clear();
	vecs_interp_indices.clear();
	return 0;
    }

    // Create NTL answers not yet created
    for (nqueries_t q = answers_ntl.size(); q < answers.size(); ++q) {
	answers_ntl.push_back(vector<vec_GF2E>(words_per_block));
	for (dbsize_t c = 0; c < words_per_block; ++c) {
	    answers_ntl[q][c].SetLength(num_servers);
	    for (nservers_t j = 0; j < num_servers; ++j) {
		conv(answers_ntl[q][c][j],
			GF2XFromBytes((unsigned char *)(answers[q] +
			    c * num_servers + j), sizeof(GF2E_Element)));
	    }
	}
    }

    RSDecoder_GF2E decoder;

    // Call the decoder's Recover method
    bool res = decoder.Recover(sizeof(GF2E_Element), t+tau, h, goodservers,
	    answers_ntl, indices_ntl, unfinished_results, decoded,
	    vecs_interp_indices);

    // Convert the results to PercyBlockResults
    // Remove queries that are completely decoded
    GF2E_Element *sigma = new GF2E_Element[words_per_block];
    for (nqueries_t q = 0; q < answers.size(); ++q) {
	for (dbsize_t i = 0; i < unfinished_results[q].size(); ++i) {
	    ssize_t indices_len = vecs_interp_indices[q].length();
	    const int num_subqueries = indices_len > 0 ? indices_len : 1;
	    for (int sub_query = 0; sub_query < num_subqueries; sub_query++) {
		results.push_back(vector<PercyResult>());
		vector<PercyResult> &block_results = results.back();
		if (decoded[q].size() == words_per_block) {
		    for (dbsize_t j = 0; j < words_per_block; ++j) {
			BytesFromGF2X((unsigned char *)(sigma + j),
				      rep(unfinished_results[q][i].recovered[sub_query][j]),
				      sizeof(GF2E_Element));
		    }
		    block_results.push_back(PercyResult(unfinished_results[q][i].G,
			string((char *)sigma, words_per_block * sizeof(GF2E_Element))));
		}
	    }
	}

        if (decoded[q].size() == words_per_block) {
            // Remove query
            if (answers[q] != NULL) {
                delete[] answers[q];
            }
            answers.erase(answers.begin() + q);
            answers_ntl.erase(answers_ntl.begin() + q);
            unfinished_results.erase(unfinished_results.begin() + q);
            decoded.erase(decoded.begin() + q);
            vecs_interp_indices.erase(vecs_interp_indices.begin() + q);
            --q;
        }
    }
    delete[] sigma;

    (void)res;
    return answers.size();
}

#endif
