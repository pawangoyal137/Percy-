// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>,
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

#include <vec_vec_ZZ_p.h>
#include <ZZ_pX.h>
#include "itclient.h"
#include "rsdecoder.h"
#include "gf2e.h"

#ifdef SPIR_SUPPORT
#include "spirclient.h"
#endif

// Some helper methods

void XOR_equal(unsigned char *dst, const unsigned char *src, unsigned const int len) {
    for(unsigned int i=0; i<len; i++) {
        *(dst++) ^= *(src++);
    }
}


// Generate t-private (t+1)-of-l shares of a given secret value.
// If a non-NULL pointer *polyp is provided, the random polynomial used is 
// copied to that address.

static void genshares(nservers_t t, nservers_t l,
	const vec_ZZ_p &indices, vec_ZZ_p &values, dbsize_t secret,
	ZZ_pX *polyp = NULL)
{
      // Pick a random polynomial of degree t
      ZZ_pX randpoly = random_ZZ_pX(t+1);
      // Set the constant term to the secret
      SetCoeff(randpoly, 0, secret);
  
      // Evaluate the polynomial at each of the indices
      for (nservers_t i=0; i<l; ++i) {
        eval(values[i], randpoly, indices[i]);
      }
 
     // Store this polynomial if requested by the caller
     if (polyp != NULL) {
  	*polyp = randpoly;
     } 

}

// Generate t-private shares, degree is t+qbs-1
static void gensharesmulti(nservers_t t, nservers_t l,
        const vec_ZZ_p &indices, vec_ZZ_p &values, const vec_ZZ_p& secret,
	const vec_ZZ_p &index, ZZ_pX *polyp = NULL)
{
    // Pick a random polynomial of degree t-1
    ZZ_pX randpoly = random_ZZ_pX(t);
    // Modification of the polynomial
    ZZ_pX addpoly = BuildFromRoots(index); 
    ZZ_pX modpoly = randpoly*addpoly;
    ZZ_pX lagrange = interpolate(index, secret);
    ZZ_pX final = modpoly+lagrange;

    // Evaluate the polynomial at each of the indices
    for (nservers_t i=0; i<l; ++i) {
        eval(values[i], final, indices[i]);
    }

    // Store this polynomial if requested by the caller
    if (polyp != NULL) {
        *polyp = final;
    }
}


// ZZ_P CLASS

// Constructor and destructor
PercyClient_ZZ_p::PercyClient_ZZ_p (const PercyClientParams * params, 
	nservers_t num_servers, nservers_t t, sid_t * sids, PercyStats * stats) :
    PercyClient(params, num_servers, t, stats),
    params(static_cast<const ZZ_pParams*>(params->percy_params())),
    indices(),
    interp_indices(),
    vecs_interp_indices(),
    randmults(),
    answers(),
    unfinished_results(),
    decoded()
{
    choose_indices(sids);
    nqueries_t max_qbs = num_servers - t - this->params->tau();
    choose_interp_indices(max_qbs);
}

PercyClient_ZZ_p::~PercyClient_ZZ_p ()
{
    map<nqueries_t, vec_vec_ZZ_p *>::iterator qiter;
    for (qiter = queries.begin(); qiter != queries.end(); ++qiter) {
	delete[] qiter->second;
    }
    queries.clear();
#ifdef SPIR_SUPPORT
    map<nqueries_t, SPIRClientQuery *>::iterator siter;
    for (siter = spir_query_info.begin(); siter != spir_query_info.end(); ++siter) {
	delete[] siter->second;
    }
    spir_query_info.clear();
#endif
}

// TODO
void PercyClient_ZZ_p::choose_indices (sid_t * sids) 
{
    indices.SetLength(num_servers);
    for (nservers_t j = 0; j < num_servers; ++j) {
        if (params->tau() > 0 || params->spir()) {
            // Use the constant indices 1, 2, ..., num_servers.
            indices[j] = to_ZZ_p(sids[j]);
        } else {
            // Use random indices
            ZZ_p r;
            bool ok = false;
            do {
                r = random_ZZ_p();
                if (r != 0) {
                    ok = true;
                    for (nservers_t k = 0; k < j; ++k) {
                        if (indices[k] == r) {
                            ok = false;
                        }   
                    }   
                }   
            } while (!ok);
            indices[j] = r;
        }
    }
}

void PercyClient_ZZ_p::choose_interp_indices (nqueries_t qbs)
{
    interp_indices.SetLength(qbs);
    ZZ_p candidate = ZZ_p::zero();
    nqueries_t i = 0;
    while (i < qbs) {
	bool in_server_indices = false;
	for (nservers_t s=0; s<num_servers; ++s) {
	    if (indices[s] == candidate) {
		in_server_indices = true;
	    }
	}
	if (!in_server_indices) {
	    interp_indices[i] = candidate;
	    ++i;
	}
	--candidate;
    }
}

// Send a request for the given block number (0-based) to the
// servers connected with the ostreams in the given vector.
void PercyClient_ZZ_p::encode_request_impl (nqueries_t request_identifier)
{
    const vector<dbsize_t>& block_numbers = get_block_numbers(request_identifier);
    nqueries_t qbs = get_qbs(request_identifier);
    nqueries_t num_block_requests = block_numbers.size();
    nqueries_t num_queries = (num_block_requests+qbs-1)/qbs;
    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();
    params->mod_modulus();

    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t virtual_block_size = params->virtual_block_size();

#ifdef SPIR_SUPPORT
    SPIRClientQuery * spir_queries = new SPIRClientQuery[num_queries];
    spir_query_info[request_identifier] = spir_queries;
#endif

    // Generate random multiples (!= 0)
    if (randomize && !params->spir()) {
	randmults[request_identifier] = vec_ZZ_p();
	vec_ZZ_p& request_randmults = randmults[request_identifier];
	request_randmults.SetLength(num_servers * num_queries);
	for (nservers_t j = 0; j < num_servers * num_queries; ++j) {
                ZZ_p r;
                while (IsZero(r)) {
                    r = random_ZZ_p();
                }
                request_randmults[j] = r;
        }
    }

    // Construct the shares of the e_{index} vector
    vec_vec_ZZ_p * shares = new vec_vec_ZZ_p[num_queries];
    vec_ZZ_pX * polyvec = new vec_ZZ_pX[num_queries];
    nqueries_t current_query = 0;
    for (nqueries_t q=0; q<num_block_requests; q += qbs) {
	nqueries_t tmpblock = (num_block_requests - q) > qbs ?
				qbs : (num_block_requests - q);
	vec_ZZ_p interps = interp_indices;
	interps.SetLength(tmpblock);
	vecs_interp_indices.push_back(interps);
	shares[current_query].SetLength(num_virtual_blocks);
	polyvec[current_query].SetLength(num_virtual_blocks);
        // we need to make sure that we have 000, 000 100 010 000 001 for fetching block 2, 3, 5
	for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
	    params->mod_modulussq();
	    shares[current_query][i].SetLength(num_servers);
	    params->mod_modulus();
	    if (tmpblock == 1) {
		dbsize_t virtual_block_number = block_numbers[q] / virtual_block_size;
	        genshares(t, num_servers, indices, shares[current_query][i],
			i == virtual_block_number, &(polyvec[current_query][i]));
	    } else {
	        // create vec_ZZ_p
	        vec_ZZ_p secrets;
	        secrets.SetLength(tmpblock);
	        for (nqueries_t k = 0; k < tmpblock; k++) {
		    dbsize_t virtual_block_number = block_numbers[q+k] / virtual_block_size;
	            if (virtual_block_number == i) {
		      secrets[k] = to_ZZ_p(1);
	   	    } else {
		      secrets[k] = to_ZZ_p(0);
		    } 
	        }
	        gensharesmulti(t, num_servers, indices, 
			shares[current_query][i], secrets, interps, 
			&(polyvec[current_query][i]));
	    }
	}
        
#ifdef SPIR_SUPPORT
        if (params->spir()) {
	    spir_queries[q].init_parameters(params, num_servers, t,
		    virtual_block_number, indices, shares[q], polyvec[q]);
        }
#endif
	
	current_query += 1;
    }

    // Multiply shares by random multiples
    if (randomize && !params->spir()) {
	vec_ZZ_p& request_randmults = randmults[request_identifier];
	for (nservers_t j = 0; j < num_servers; ++j) {
	    for (nqueries_t q = 0; q < num_queries; ++q) {
                for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
                    shares[q][i][j] *= request_randmults[j*num_queries+q];
                }
            }
        }
    }

    // Optionally encrypt the shares
    if (params->hybrid()) {
        for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
            for (nservers_t j = 0; j < num_servers; ++j) {
                for (nqueries_t q=0; q<num_queries; ++q) {
                    params->mod_modulus();
                    ZZ share = rep(shares[q][i][j]);
                    params->mod_modulussq();
                    shares[q][i][j] = params->encrypt(share);
                }
            }
        }
        // If we're encrypting, leave modulus^2 as the active modulus
        params->mod_modulussq();
    }

    // Store queries
    queries[request_identifier] = shares;

    // Clean up
    delete[] polyvec;
    savectx.restore();
}

dbsize_t PercyClient_ZZ_p::send_request_impl (nqueries_t request_identifier,
        std::vector<ostream*> &osvec, bool send_num_queries)
{
    if (num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << num_servers
            << ".\n";
        return -1;
    }

    // Save the current ZZ_p context
    ZZ_pContext savectx;
    savectx.save();
    params->mod_modulus();

    // Get stored info
    vec_vec_ZZ_p * shares = queries[request_identifier];
#ifdef SPIR_SUPPORT
    SPIRClientQuery * spir_queries = spir_query_info[request_identifier];
#endif
    nqueries_t num_block_requests = get_block_numbers(request_identifier).size();
    nqueries_t qbs = get_qbs(request_identifier);
    nqueries_t num_queries = (num_block_requests+qbs-1)/qbs;
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();

    // Send the query to each server
    unsigned int modulus_bytes = params->hybrid() ? params->modulussq_bytes() :
	    params->modulus_bytes();
    for (nservers_t j = 0; j < num_servers; ++j) {
	if (send_num_queries) {
	    percy_write_le_uint16(*osvec[j], num_queries);
	}

        for (nqueries_t q=0; q<num_queries; ++q) {
#ifdef SPIR_SUPPORT
            if (params->spir()) {
		spir_queries[q].send_to_server(params, j, *(osvec[j]));
            }
#endif

	    unsigned char * bytes = new unsigned char[modulus_bytes];
            for (dbsize_t i = 0; i < num_virtual_blocks; ++i) {
                BytesFromZZ(bytes, rep(shares[q][i][j]), modulus_bytes);
                osvec[j]->write((char *)bytes, modulus_bytes);
            }
	    delete[] bytes;
        }
        osvec[j]->flush();
    }

#ifdef SPIR_SUPPORT
    delete[] spir_queries;
    spir_query_info.erase(request_identifier);
#endif
    delete[] shares;
    queries.erase(request_identifier);

    savectx.restore();
    return num_servers * num_queries * num_virtual_blocks * modulus_bytes;
}

// Receive the server's replies, and return a number of servers that
// gave complete (but not necessarily correct) replies.
dbsize_t PercyClient_ZZ_p::receive_replies_impl (nqueries_t request_identifier,
	std::vector<istream*> &isvec)
{
    dbsize_t words_per_block = params->words_per_block();
    nqueries_t num_block_requests = get_block_numbers(request_identifier).size();
    nqueries_t qbs = get_qbs(request_identifier);
    nqueries_t num_queries = (num_block_requests+qbs-1)/qbs;
    nqueries_t q;

    // Choose the right modulus
    unsigned int modulus_bytes;
    if (params->hybrid()) {
        modulus_bytes = params->modulussq_bytes();
        params->mod_modulussq();
    } else {
        modulus_bytes = params->modulus_bytes();
        params->mod_modulus();
    }

    // The vector of servers that have responded properly
    goodservers.clear();

    // Add to the appropriate vectors
    nqueries_t prev = answers.size();
    for (q=0; q<num_queries; ++q) {
	answers.push_back(vector<vec_ZZ_p>(words_per_block));
        for (dbsize_t c = 0; c < words_per_block; ++c) {
            answers.back()[c].SetLength(num_servers);
        }
    }
    // Read the replies
    unsigned char * bytes = new unsigned char[modulus_bytes];
    for (nservers_t j = 0; j < num_servers; ++j) {
        bool isgood = true;
	for (q=0; q<num_queries; ++q) {
	    for (dbsize_t i = 0; isgood && i < words_per_block; ++i) {
                isvec[j]->read((char *)bytes, modulus_bytes);
                if (isvec[j]->eof()) {
                    std::cerr << "Server " << j+1 << " did not send complete reply.\n";
                    std::cerr << "Marking server " << j+1 << " as bad.\n";
                    isgood = false;
                    break;
                }
                if ((dbsize_t)(isvec[j]->gcount()) < modulus_bytes) {
                    // Mark this server as bad
                    std::cerr << "Marking server " << j+1 << " as bad.\n";
                    isgood = false;
                    break;
                }
                ZZ ans;
                ZZFromBytes(ans, bytes, modulus_bytes);
                answers[q+prev][i][j] = to_ZZ_p(ans);
            }
        }
        if (isgood) {
            goodservers.push_back(j);
        }
    }
    delete[] bytes;

    // Add to unfinished_results and decoded
    for (q=0; q<num_queries; ++q) {
	unfinished_results.push_back(
		vector<DecoderResult<ZZ_p> >(1, DecoderResult<ZZ_p>(goodservers, vector<map<dbsize_t, ZZ_p> >())));
	decoded.push_back(std::set<dbsize_t>());
    }

    // Optionally decrypt the answers
    if (params->hybrid()) {
        params->mod_modulussq();
        for (dbsize_t i = 0; i < words_per_block; ++i) {
            for (nservers_t j = 0; j < num_servers; ++j) {
                for (q=0; q<num_queries; ++q) {
                    ZZ_p dec = params->decrypt(answers[q+prev][i][j]);
                    clear(answers[q+prev][i][j]);
                    answers[q+prev][i][j] = dec;
                }
            }
        }
    }

    // Remove random multiple
    if (randomize && !params->spir()) {
	vec_ZZ_p& request_randmults = randmults[request_identifier];
        for (nservers_t j = 0; j < num_servers; ++j) {
            for (q = 0; q < num_queries; ++q) {
                ZZ_p randmult_inv = inv(request_randmults[j*num_queries+q]);
                for (dbsize_t i = 0; i < words_per_block; ++i) {
                    answers[q+prev][i][j] *= randmult_inv;
                }
            }
        }
	randmults.erase(request_identifier);
    }

    // Now we're mod modulus for sure
    params->mod_modulus();

    return num_servers * num_queries * words_per_block * modulus_bytes;
}

nqueries_t PercyClient_ZZ_p::process_replies_impl (nservers_t h, 
	vector<vector<PercyResult> > &results)
{
    dbsize_t words_per_block = params->words_per_block();
    dbsize_t bytes_per_word = (params->word_size() - 1) / 8 + 1;
    nservers_t tau = params->tau();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return false;
    }

    RSDecoder_ZZ_p decoder(params->get_p1(), params->get_p2());
    // Call the decoder's Recover method
    bool res = decoder.Recover(bytes_per_word, t+tau, h, goodservers, answers, 
	    indices, unfinished_results, decoded, vecs_interp_indices);
    
    // Convert the results to PercyBlockResults
    // Remove queries that are completely decoded
    unsigned char * sigma = new unsigned char[words_per_block * bytes_per_word];
    for (nqueries_t q = 0; q < answers.size(); ++q) {
	for (dbsize_t i = 0; i < unfinished_results[q].size(); ++i) {
	    ssize_t indices_len = vecs_interp_indices[q].length();
	    const int num_subqueries = indices_len > 0 ? indices_len : 1;
	    for (int sub_query = 0; sub_query < num_subqueries; sub_query++) {
		results.push_back(vector<PercyResult>());
		vector<PercyResult> &curr_block_res = results.back();
		if (decoded[q].size() == answers[q].size()) {
		    for (dbsize_t j = 0; j < words_per_block; ++j) {
			BytesFromZZ(sigma+(j*bytes_per_word),
				rep(unfinished_results[q][i].recovered[sub_query][j]), 
				bytes_per_word);
		    }
		    curr_block_res.push_back(PercyResult(unfinished_results[q][i].G,
			string((char *)sigma, words_per_block * bytes_per_word)));
		}
	    }
	}

	if (decoded[q].size() == answers[q].size()) {
	    // Remove query
	    answers.erase(answers.begin() + q);
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


// CHOR CLASS

// Constructor and destructor
PercyClient_Chor::PercyClient_Chor (const PercyClientParams * params, 
	nservers_t num_servers, PercyStats * stats) :
    PercyClient(params, num_servers, num_servers-1, stats),
    params(static_cast<const ChorParams*>(params->percy_params())),
    answers()
{}

PercyClient_Chor::~PercyClient_Chor ()
{
    map<nqueries_t, unsigned char *>::iterator qiter;
    for (qiter = queries.begin(); qiter != queries.end(); ++qiter) {
	delete[] qiter->second;
    }
    queries.clear();
    while (answers.size() > 0) {
	if (answers.back() != NULL) {
	    delete[] answers.back();
	}
	answers.pop_back();
    }
}


void PercyClient_Chor::encode_request_impl (nqueries_t request_identifier)
{
    const vector<dbsize_t>& block_numbers = get_block_numbers(request_identifier);
    nqueries_t num_queries = block_numbers.size();

    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t num_bytes = (num_virtual_blocks - 1) / 8 + 1;

    unsigned char * shares = new unsigned char [num_servers*num_queries*num_bytes];
    memset(shares, 0, num_servers * num_queries * num_bytes);
    for (nqueries_t q = 0; q < num_queries; q++) {
        // Generate random bytestrings for each server but the last
        nservers_t last_server = num_servers - 1;
        for (nservers_t j = 0; j < last_server; j++) {
	    BytesFromZZ(shares + (j * num_queries + q) * num_bytes,
		    RandomBits_ZZ(num_virtual_blocks), num_bytes);
        }
        // Create a bytestring out of the block_numbers vector, storing it in the location of the 
        // final server's bytestring
        memset(shares + (last_server*num_queries+q)*num_bytes, '\0', num_bytes);
        unsigned int block_byte = (unsigned int) block_numbers[q] / 8;
        unsigned char block_pos = (unsigned char) block_numbers[q] % 8;
        unsigned char byte_val = 1 << block_pos;
        shares[(last_server*num_queries + q)*num_bytes + block_byte] = byte_val;
        // Compute the final bytestring (such that the XOR of all the 
	// server bytestrings is the bytestring 
        // representing the query vector)
        for (nservers_t i = 0; i < last_server; i++) {
            XOR_equal(
                    shares + (last_server*num_queries + q)*num_bytes,
                    shares + (i*num_queries + q)*num_bytes, num_bytes);
        }
    }
    queries[request_identifier] = shares;
}

dbsize_t PercyClient_Chor::send_request_impl (nqueries_t request_identifier,
        std::vector<ostream*> &osvec, bool send_num_queries)
{
    if (num_servers != osvec.size()) {
        std::cerr << "Incorrect iostream vector size passed to "
            "send_request_Chor.\n";
        std::cerr << "Was " << osvec.size() << ", should be " << num_servers
            << ".\n";
        return 0;
    }

    // Send the query to each server
    nqueries_t num_queries = get_block_numbers(request_identifier).size();
    unsigned char * shares = queries[request_identifier];
    dbsize_t num_bytes = (params->num_blocks() - 1) / 8 + 1;
    for (nservers_t j = 0; j < num_servers; ++j) {
	if (send_num_queries) {
	    percy_write_le_uint16(*osvec[j], num_queries);
	}

        osvec[j]->write((char *)(shares + j*num_queries*num_bytes), num_queries*num_bytes);
        osvec[j]->flush();
    }
    delete[] shares;
    queries.erase(request_identifier);

    return num_queries * num_bytes;
}

dbsize_t PercyClient_Chor::receive_replies_impl (
	nqueries_t request_identifier, std::vector<istream*> &isvec)
{
    //dbsize_t words_per_block = params->words_per_block();
    dbsize_t block_size = params->block_size();
    nqueries_t num_queries = get_block_numbers(request_identifier).size();

    // The vector of servers that have responded properly
    goodservers.clear();

    // Allocate space for answers
    nqueries_t previous_queries = answers.size();
    for (nqueries_t q = 0; q < num_queries; ++q) {
	answers.push_back(new unsigned char[num_servers * block_size]);
    }

    // Read the replies
    for (nservers_t j = 0; j < num_servers; ++j) {
        bool isgood = true;
        for (nqueries_t q=0; q<num_queries; ++q) {
            for (dbsize_t i = 0; isgood && i < block_size; ++i) {
                isvec[j]->read((char *)(answers[q+previous_queries] 
			+ j*block_size + i), 1);
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

    return num_servers * num_queries * block_size;
}

nqueries_t PercyClient_Chor::process_replies_impl (nservers_t h,
	vector<vector<PercyResult> > &results)
{
    //dbsize_t words_per_block = params->words_per_block();
    dbsize_t block_size = params->block_size();
    nqueries_t num_queries = answers.size();

    // Check that results is empty
    if (!(results.empty())) {
	std::cerr << "The results vector must be empty\n";
	return num_queries;
    }

    if (goodservers.size() < num_servers) {
        std::cerr << "Not enough servers responded." << std::endl;
	return num_queries;
    }

    unsigned char *result = new unsigned char[block_size];
    for (nqueries_t q = 0; q < num_queries; q++) {
        // Recover the query block by XORing all the replies together
        memset(result, '\0', block_size);
        for (nservers_t j = 0; j < num_servers; j++) {
            XOR_equal(result, answers[q] + j*block_size, block_size);
        }

        vector<nservers_t> empty;
        PercyResult thisresult (empty, string((char *)result, block_size));
	results.push_back(vector<PercyResult>(1, thisresult));

	delete[] answers[q];
    }
    delete[] result;
    answers.clear();

    return 0;
}

