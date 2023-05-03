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

#ifndef __ITPARAMS_H__
#define __ITPARAMS_H__

#include <iostream>
#include <NTL/ZZ_p.h>
#include <vector>
#ifdef SPIR_SUPPORT
#include <PolyCommitCommon.h>
#endif
#include "percytypes.h"
#include "percyparams.h"

NTL_CLIENT

class ZZ_pParams : public PercyParams {
public:
    // Contructor for ZZ_p
    ZZ_pParams (dbsize_t num_blocks, dbsize_t block_size, dbsize_t word_size,
	    ZZ modulus, nservers_t tau = 0, char * pcparams_file = NULL, 
	    bool do_spir = false, dbsize_t virtual_block_size = 1);

    // Generate a new public/private key pair of the given keysize
    ZZ_pParams (dbsize_t num_blocks, dbsize_t block_size, dbsize_t word_size, 
	    unsigned long modulus_bits, nservers_t tau = 0, 
	    dbsize_t virtual_block_size = 1);

    // Use the given factors to generate a public/private key pair
    ZZ_pParams (dbsize_t num_blocks, dbsize_t block_size, dbsize_t word_size, 
	    ZZ p, ZZ q, nservers_t tau = 0, dbsize_t virtual_block_size = 1);

    virtual ~ZZ_pParams ();

    unsigned long modulus_bytes() const {
	return NumBytes(modulus);
    }
    unsigned long modulussq_bytes() const {
	return NumBytes(modulus*modulus);
    }
    void mod_modulus() const {
	modctx.restore();
    }
    void mod_modulussq() const {
	modsqctx.restore();
    }
    const ZZ &get_modulus() const {
	return modulus;
    }
    const ZZ_p &get_g() const {
	return g;
    }
    bool hybrid() const {
	return hybrid_protection;
    }
    bool spir() const {
	return do_spir;
    }
    bool modulus_match(ZZ testmod) const {
	return modulus == testmod;
    }
    char * get_pcparams_filename() const {
	return pcparams_filename;
    }	
#ifdef SPIR_SUPPORT
    PolyCommitParams * get_pcparamsp() const {
        return pcparamsp;
    }
#endif
    virtual dbsize_t server_block_size () const { 
	return (_tau ? modulus_bytes() * words_per_block() : _block_size);
    }
    // Encrypt the given plaintext.  The current ZZ_p context must be
    // modsqctx.
    ZZ_p encrypt(ZZ plaintext) const {
	ZZ_p r;
	random(r);
	return power(g, plaintext) * power(r, modulus);
    }
    // Decrypt the given ciphertext.  This routine will change the
    // current ZZ_p context to modctx.
    ZZ_p decrypt(ZZ_p ciphertext) const {
	modsqctx.restore();
	ZZ Lval = rep(power(ciphertext, lambda) - 1) / modulus;
	modctx.restore();
	ZZ_p ret = to_ZZ_p(Lval) * mu;
	return ret;
    }		
    ZZ get_p1() const { return p1; }
    ZZ get_p2() const { return p2; }

    // Return the size of the request/response
    virtual dbsize_t request_size (nqueries_t num_queries = 1) const;
    virtual dbsize_t response_size (nqueries_t num_queries = 1) const;

    // For use in distributed computation
    virtual std::vector<const PercyParams*> create_worker_params (
	    std::vector<Dimension> worker_dims) const;

protected:
    // Write the parameters to a stream to check compatibility
    virtual void write (std::ostream &os) const;

    // Read the parameters from a stream (as written by write()) and check that
    // they are compatible with these parameters.
    virtual bool check_compatible (std::istream &is) const;

    void create_ZZ_pContexts();
    bool hybrid_protection;
    bool do_spir;
    ZZ_pContext modctx, modsqctx;
    // Paillier public key
    ZZ modulus;
    ZZ_p g;   // mod modulus^2
    char * pcparams_filename;
#ifdef SPIR_SUPPORT
    PolyCommitParams * pcparamsp;
#endif

    void init_hybrid(ZZ p, ZZ q);
    // Paillier private key
    ZZ lambda;
    ZZ_p mu;  // mod modulus
    ZZ p1, p2;
};


class GF2EParams : public PercyParams {
public:
    GF2EParams (dbsize_t num_blocks, dbsize_t block_size, 
	    dbsize_t word_size = 8, nservers_t tau = 0, 
	    dbsize_t virtual_block_size = 1);

    virtual ~GF2EParams () {}

    // Return the size of the request/response
    virtual dbsize_t request_size (nqueries_t num_queries = 1) const;
    virtual dbsize_t response_size (nqueries_t num_queries = 1) const;

    // For use in distributed computation
    virtual std::vector<const PercyParams*> create_worker_params (
	    std::vector<Dimension> worker_dims) const;
};


class ChorParams : public PercyParams {
public:
    ChorParams (dbsize_t num_blocks, dbsize_t block_size, 
	    dbsize_t virtual_block_size = 1);

    virtual ~ChorParams () {}

    // Return the size of the request/response
    virtual dbsize_t request_size (nqueries_t num_queries = 1) const;
    virtual dbsize_t response_size (nqueries_t num_queries = 1) const;

    // For use in distributed computation
    virtual std::vector<const PercyParams*> create_worker_params (
	    std::vector<Dimension> worker_dims) const;
};

#endif
