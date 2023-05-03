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

#include <string.h>
#include <fstream>
#include "itparams.h"
#include "percyio.h"

#ifdef SPIR_SUPPORT
// Open the specified file, read in the PolyCommit parameters contained 
// within, create a new Params object from the parameters, and return it
PolyCommitParams * pcparams_init(const char * filename) {
    if (!filename) {
        return NULL;
    }
    ifstream ifile(filename);
    if(!ifile.is_open()) {
        std::cerr << "Error: Cannot open params file." << endl;
        exit(1);
    }
    PolyCommitParams * pcparamsp = new PolyCommitParams();
    ifile >> *pcparamsp;
    ifile.close();

    return pcparamsp;
}
#endif

ZZ_pParams::ZZ_pParams (dbsize_t num_blocks, dbsize_t block_size, 
	dbsize_t word_size, ZZ modulus, nservers_t tau, char * pcparams_file, 
	bool do_spir, dbsize_t virtual_block_size)
:
    PercyParams(num_blocks, block_size, word_size, MODE_ZZ_P, tau,
	    virtual_block_size),
    hybrid_protection(false),
    do_spir(do_spir),
    modctx(),
    modsqctx(),
    modulus(modulus),
    g(),
    pcparams_filename(pcparams_file)
#ifdef SPIR_SUPPORT
    , pcparamsp(pcparams_init(pcparams_file))
#endif
{
    create_ZZ_pContexts();
}

ZZ_pParams::~ZZ_pParams ()
{
#ifdef SPIR_SUPPORT
    if (pcparamsp) {
	delete pcparamsp;
    }
#endif
}

ZZ_pParams::ZZ_pParams (dbsize_t num_blocks, dbsize_t block_size, 
	dbsize_t word_size, ZZ p, ZZ q, nservers_t tau, 
	dbsize_t virtual_block_size)
:
    PercyParams(num_blocks, block_size, word_size, MODE_ZZ_P, tau,
	    virtual_block_size),
    hybrid_protection(true),
    do_spir(false),
    modctx(),
    modsqctx(),
    modulus(p * q),
    g(),
    pcparams_filename(NULL)
{
    create_ZZ_pContexts();
    init_hybrid(p, q);
}

ZZ_pParams::ZZ_pParams (dbsize_t num_blocks, dbsize_t block_size, 
	dbsize_t word_size, unsigned long modulus_bits, nservers_t tau, 
	dbsize_t virtual_block_size)
:
    PercyParams(num_blocks, block_size, word_size, MODE_ZZ_P, tau,
	    virtual_block_size),
    hybrid_protection(true),
    do_spir(false),
    modctx(),
    modsqctx(),
    modulus(),
    g(),
    pcparams_filename(NULL)
{
    // Pick the sizes for the primes
    unsigned long qsize = modulus_bits / 2;
    unsigned long psize = modulus_bits - qsize;

    // Generate random primes of the appropriate size.  We ensure the
    // top two bits are set so that their product is of the right
    // bitlength.
    ZZ pbase, qbase, p, q;
    RandomBits(pbase, psize);
    if (psize >= 1) SetBit(pbase, psize-1);
    if (psize >= 2) SetBit(pbase, psize-2);
    NextPrime(p, pbase);
    RandomBits(qbase, qsize);
    if (qsize >= 1) SetBit(qbase, qsize-1);
    if (qsize >= 2) SetBit(qbase, qsize-2);
    NextPrime(q, qbase);

    modulus = p * q;
    create_ZZ_pContexts();

    init_hybrid(p, q);
}

void ZZ_pParams::init_hybrid(ZZ p, ZZ q)
{
    p1 = p;
    p2 = q;

    // Generate the Paillier public and private parts
    ZZ pm1, qm1;
    pm1 = p - 1;
    qm1 = q - 1;
    this->lambda = pm1 * qm1 / GCD(pm1, qm1);

    ZZ_pContext savectx;
    savectx.save();
    mod_modulussq();
    random(this->g);
    ZZ muinv = rep(power(this->g, this->lambda) - 1) / modulus;
    mod_modulus();
    this->mu = inv(to_ZZ_p(muinv));
    savectx.restore();
}

void ZZ_pParams::create_ZZ_pContexts()
{
    // Create the ZZ_pContexts
    ZZ_pContext modctx(modulus);
    ZZ_pContext modsqctx(modulus * modulus);
    this->modctx = modctx;
    this->modsqctx = modsqctx;
}

dbsize_t ZZ_pParams::request_size (nqueries_t num_queries) const
{
    // NOTE: this does not include SPIR
    return num_queries * _num_virtual_blocks * modulus_bytes();
}

dbsize_t ZZ_pParams::response_size (nqueries_t num_queries) const
{
    return num_queries * _words_per_block * _virtual_block_size * 
	    modulus_bytes();
}

void ZZ_pParams::write (ostream& os) const
{
    // Output parent class
    PercyParams::write(os);

    // Output the version number and flags
    unsigned char c = (hybrid_protection ? 1 : 0) | 
			(do_spir ? 2 : 0);
    os.write((char *)(&c), 1);

    // Output the modulus
    percy_write_ZZ(os, modulus);

    // Output g, if appropriate
    if (hybrid_protection) {
	percy_write_ZZ(os, rep(g));
    }

    os.flush();
}

bool ZZ_pParams::check_compatible (istream &is) const 
{
    // Test parent class
    if (!(PercyParams::check_compatible(is))) return false;

    // Input the version number and flags
    unsigned char c;
    is.read((char *)(&c), 1);
    if (((c & 1) == 1) != hybrid_protection ||
	((c & 2) == 2) != do_spir) return false;

    // Input the modulus
    ZZ value;
    percy_read_ZZ(is, value);
    if (!modulus_match(value)) return false;

    // Input g, if appropriate
    if (hybrid_protection) {
	percy_read_ZZ(is, value);
	if (rep(g) != value) return false;
    }

    return true;
}

std::vector<const PercyParams*> ZZ_pParams::create_worker_params (
	std::vector<Dimension> worker_dims) const
{
    std::vector<const PercyParams*> worker_params;
    std::vector<Dimension>::iterator it;
    for (it = worker_dims.begin(); it != worker_dims.end(); ++it) {
	worker_params.push_back(new ZZ_pParams(it->first, it->second, 
		_word_size, modulus, _tau, pcparams_filename, do_spir, 
		_virtual_block_size));
    }
    return worker_params;
}


GF2EParams::GF2EParams (dbsize_t num_blocks, dbsize_t block_size,
	dbsize_t word_size, nservers_t tau, dbsize_t virtual_block_size)
:
    PercyParams(num_blocks, block_size, word_size, MODE_GF28, tau,
	    virtual_block_size)
{
    switch (word_size) {
    case 8:
	break;
    case 16:
	mode = MODE_GF216;
	break;
    default:
	std::cerr << "Warning: not a valid GF2E word size: " << word_size
		<< ".  Will use word size of 8.\n";
	_word_size = 8;
	break;
    }
}

dbsize_t GF2EParams::request_size (nqueries_t num_queries) const
{
    return num_queries * _num_virtual_blocks * _word_size / 8;
}

dbsize_t GF2EParams::response_size (nqueries_t num_queries) const
{
    return num_queries * _block_size * _virtual_block_size;
}

std::vector<const PercyParams*> GF2EParams::create_worker_params (
	std::vector<Dimension> worker_dims) const
{
    std::vector<const PercyParams*> worker_params;
    std::vector<Dimension>::iterator it;
    for (it = worker_dims.begin(); it != worker_dims.end(); ++it) {
	worker_params.push_back(new GF2EParams(it->first, it->second, 
		_word_size, _tau, _virtual_block_size));
    }
    return worker_params;
}


ChorParams::ChorParams (dbsize_t num_blocks, dbsize_t block_size, 
	dbsize_t virtual_block_size)
:
    PercyParams(num_blocks, block_size, 1, MODE_CHOR, 0, virtual_block_size)
{}

dbsize_t ChorParams::request_size (nqueries_t num_queries) const
{
    return num_queries * ((_num_virtual_blocks - 1) / 8 + 1);
}

dbsize_t ChorParams::response_size (nqueries_t num_queries) const
{
    return num_queries * _block_size * _virtual_block_size;
}

std::vector<const PercyParams*> ChorParams::create_worker_params (
	std::vector<Dimension> worker_dims) const
{
    std::vector<const PercyParams*> worker_params;
    std::vector<Dimension>::iterator it;
    for (it = worker_dims.begin(); it != worker_dims.end(); ++it) {
	worker_params.push_back(new ChorParams(it->first, it->second, 
		_virtual_block_size));
    }
    return worker_params;
}

