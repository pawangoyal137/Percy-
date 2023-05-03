// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Ryan Henry <rhenry@cs.uwaterloo.ca>
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
#include "percyparams.h"
#include "percyio.h"

const char * PercyModeStrings[] = {
    "NONE",
    "ZZ_P",
    "GF28",
    "GF216",
    "CHOR",
    "AG",
    "RECURSIVE_AG",
    "HYBRID"
};

std::ostream& operator<< (std::ostream& os, const PercyMode mode) {
    return os << PercyModeStrings[mode];
}

/*
const char * DistTypeStrings[] = {
    "NONE",
    "THREADED",
    "MASTER",
    "FORKED"
};

std::ostream& operator<< (std::ostream& os, const DistType disttype) {
    return os << DistTypeStrings[disttype];
}
*/

const char * DistSplitStrings[] = {
    "QUERIES",
    "RECORDS",
    "RECORD_BYTES"
};

std::ostream& operator<< (std::ostream& os, const DistSplit distsplit) {
    return os << DistSplitStrings[distsplit];
}

PercyParams::PercyParams (dbsize_t num_blocks, dbsize_t block_size, 
	dbsize_t word_size, PercyMode mode, nservers_t tau, 
	dbsize_t virtual_block_size)
:
    version(),
    _num_blocks(num_blocks),
    _block_size(block_size),
    _word_size(word_size),
    _words_per_block((8 * block_size - 1) / word_size + 1),
    mode(mode),
    _tau(tau),
    _num_virtual_blocks((num_blocks - 1) / virtual_block_size + 1),
    _virtual_block_size(virtual_block_size)
{
    int vbuf[3];
    sscanf(PERCY_VERSION, "%d.%d.%d", vbuf, vbuf+1, vbuf+2);
    version[0] = (unsigned char)vbuf[0];
    version[1] = (unsigned char)vbuf[1];
    version[2] = (unsigned char)vbuf[2];

#ifdef VERBOSE_PARAMS
    fprintf(stderr, "num_blocks = %d\n", _num_blocks);
    fprintf(stderr, "block_size = %d\n", _block_size);
    fprintf(stderr, "word_size = %d\n", _word_size);
    fprintf(stderr, "words_per_block = %d\n", _words_per_block);
    fprintf(stderr, "mode = %s\n", PercyModeStrings[mode]);
    fprintf(stderr, "tau = %d\n", _tau);
#endif
}

void PercyParams::print (std::ostream& os) const
{
    os << mode << ","
	<< _num_blocks << ","
	<< _block_size << ","
	<< _word_size << ","
	<< _tau << ","
	<< _virtual_block_size << ",";
    print_mode_specific(os);
    os << ",";
}

void PercyParams::write (ostream& os) const
{
    // Output the version number and mode
    unsigned char minibuf[5];
    minibuf[0] = version[0];
    minibuf[1] = version[1];
    minibuf[2] = version[2];
    minibuf[3] = mode;
    minibuf[4] = (_tau ? 1 : 0);
    os.write((char *)minibuf, 5);

    
    // Output the database parameters
    PERCY_WRITE_LE_DBSIZE(os, _num_blocks);
    PERCY_WRITE_LE_DBSIZE(os, _block_size);
    PERCY_WRITE_LE_DBSIZE(os, _word_size);
    PERCY_WRITE_LE_DBSIZE(os, _virtual_block_size);

    os.flush();
}

bool PercyParams::check_compatible (istream &is) const
{
    // Input the version number and mode
    unsigned char minibuf[5];
    is.read((char *)minibuf, 5);
    if ((version[0] != minibuf[0]) ||
	    (version[1] != minibuf[1]) ||
	    (version[2] != minibuf[2]) ||
	    (mode != (PercyMode)minibuf[3]) ||
	    ((_tau ? 1 : 0) != minibuf[4])) return false;

    // Input the database parameters
    dbsize_t value;
    PERCY_READ_LE_DBSIZE(is, value);
    if (value != _num_blocks) return false;
    PERCY_READ_LE_DBSIZE(is, value);
    if (value != _block_size) return false;
    PERCY_READ_LE_DBSIZE(is, value);
    if (value != _word_size) return false;
    PERCY_READ_LE_DBSIZE(is, value);
    if (value != _virtual_block_size) return false;

    return true;
}


PercyClientParams::PercyClientParams (const PercyParams * params, 
	nservers_t num_servers, bool is_null)
:
    params(params),
    _num_servers(num_servers),
    null(is_null)
{}

void PercyClientParams::send (std::ostream &os, nservers_t sid) const
{
    // Output the magic header
    os.write("PIRC", 4);

    // Send sid
    percy_write_le_uint16(os, sid);

    // Send rest
    params->write(os);
}

bool PercyClientParams::is_compatible (std::istream &is, nservers_t sid) const
{
    // Input the magic header
    unsigned char minibuf[4];
    is.read((char *)minibuf, 4);
    if (memcmp(minibuf, "PIRS", 4)) {
	std::cerr << "Did not find expected PercyParams server header.\n";
	return false;
    }

    // Get the sid
    nservers_t readsid;
    percy_read_le_uint16(is, readsid);
    if (readsid != sid) {
	std::cerr << "sids did not match: local=" << sid << ", remote=" <<
		readsid << "\n";
	return false;
    }

    // Get rest of info
    return params->check_compatible(is);
}

void PercyClientParams::print (std::ostream& os) const
{
    params->print(os);
    os << _num_servers << ",";
}


PercyServerParams::PercyServerParams (const PercyParams * params, nservers_t sid, 
	bool be_byzantine)
:
    params(params),
    sid(sid),
    be_byzantine(be_byzantine),
    _num_threads(0),
    tsplit(DIST_SPLIT_RECORDS),
    is_forked(false),
    _num_workers(0),
    wsplit(DIST_SPLIT_RECORDS),
    worker_params(),
    worker_serverparams()
{}

std::vector<Dimension> get_worker_dimensions (const PercyParams * params,
	nservers_t& num_workers, DistSplit split)
{
    std::vector<Dimension> worker_dims;
    dbsize_t num_blocks = params->num_blocks();
    dbsize_t block_size = params->block_size();
    dbsize_t virtual_block_size = params->virtual_block_size();
    if (params->get_mode() == MODE_CHOR) {
	// Make sure for CHOR, we have worker num_blocks as multiples of 8
	// We will fake a virtual_block_size that is a multiple of 8
	while (virtual_block_size % 8 != 0) {
	    virtual_block_size *= 2;
	}
    }
    dbsize_t num_virtual_blocks = (num_blocks - 1) / virtual_block_size + 1;
    switch (split) {
    case DIST_SPLIT_QUERIES:
	worker_dims.resize(num_workers, Dimension(num_blocks, block_size));
	break;
    case DIST_SPLIT_RECORDS:
	{
	    if (num_virtual_blocks < num_workers) {
		num_workers = num_virtual_blocks;
	    }
	    dbsize_t num_vblocks_div = num_virtual_blocks / num_workers;
	    dbsize_t num_vblocks_mod = num_virtual_blocks % num_workers;
	    worker_dims.resize(num_vblocks_mod, Dimension(virtual_block_size
		    * (num_vblocks_div + 1), block_size));
	    worker_dims.resize(num_workers - 1, Dimension(virtual_block_size
		    * num_vblocks_div, block_size));
	    worker_dims.push_back(Dimension(num_blocks - virtual_block_size
		    * ((num_workers-1) * num_vblocks_div + num_vblocks_mod),
		    block_size));
	}
	break;
    case DIST_SPLIT_RECORD_BYTES:
	{
	    if (block_size < num_workers) {
		num_workers = block_size;
	    }
	    dbsize_t block_size_div = block_size / num_workers;
	    dbsize_t block_size_mod = block_size % num_workers;
	    worker_dims.resize(block_size_mod, Dimension(num_blocks, block_size_div + 1));
	    worker_dims.resize(num_workers, Dimension(num_blocks, block_size_div));
	}
	break;
    }
    return worker_dims;
}

/*
PercyServerParams::PercyServerParams (const PercyParams * params, 
	nservers_t sid, nservers_t num_workers, DistType disttype,
	DistSplit distsplit, std::vector<nservers_t> worker_sids, 
	bool be_byzantine)
:
    params(params),
    sid(sid),
    be_byzantine(be_byzantine),
    dtype(disttype),
    _num_workers(num_workers),
    dsplit(distsplit)
{
    if (dtype != DIST_TYPE_NONE) {
	// Split num_blocks and block_size
	std::vector<Dimension> worker_dims = get_worker_dimensions(params, 
		_num_workers, distsplit);
	worker_params = params->create_worker_params(worker_dims);

	worker_sids.resize(_num_workers, 1);
	for (nservers_t i = 0; i < worker_params.size(); ++i) {
	    if (worker_params[i]) {
		worker_serverparams.push_back(new PercyServerParams(worker_params[i], 
			worker_sids[i], be_byzantine));
#ifdef VERBOSE_DISTRIBUTED
		std::cerr << "Worker " << i << ": ";
		worker_serverparams.back()->print(std::cerr);
		std::cerr << "\n";
#endif
	    } else {
		worker_serverparams.push_back(NULL);
	    }
	}
    }
}
*/

PercyServerParams::PercyServerParams (const PercyParams * params, nservers_t sid,
	nservers_t num_threads, DistSplit tsplit, nservers_t num_workers, 
	DistSplit wsplit, std::vector<nservers_t> worker_sids, bool fork, 
	bool be_byzantine)
:
    params(params),
    sid(sid),
    be_byzantine(be_byzantine),
    _num_threads(num_threads),
    tsplit(tsplit),
    is_forked(fork),
    _num_workers(num_workers),
    wsplit(wsplit),
    worker_params(),
    worker_serverparams()
{
    // Split num_blocks and block_size
    if (_num_workers > 0) {
	std::vector<Dimension> worker_dims;
	worker_dims = get_worker_dimensions(params, _num_workers, wsplit);
#ifdef VERBOSE_DISTRIBUTED
	if (_num_workers < num_workers) {
	    std::cerr << "More workers than ";
	    switch (wsplit) {
	    case DIST_SPLIT_RECORDS:
		std::cerr << "virtual blocks";
		break;
	    case DIST_SPLIT_RECORD_BYTES:
		std::cerr << "block bytes";
		break;
	    default:
		break;
	    }
	    std::cerr << ".  Only using " << _num_workers << " workers.\n";
	}
#endif
	worker_sids.resize(_num_workers, sid);
	worker_params = params->create_worker_params(worker_dims);
	for (nservers_t i = 0; i < worker_params.size(); ++i) {
	    if (worker_params[i]) {
		worker_serverparams.push_back(new PercyServerParams(
			worker_params[i], worker_sids[i], num_threads,
			tsplit, 0, wsplit, vector<nservers_t>(), fork,
			be_byzantine));
#ifdef VERBOSE_DISTRIBUTED
		std::cerr << "Worker " << i << ": ";
		worker_serverparams.back()->print(std::cerr);
		std::cerr << "\n";
#endif
	    } else {
		worker_serverparams.push_back(NULL);
	    }
	}
    } else if (_num_threads > 0) {
	std::vector<Dimension> thread_dims;
	thread_dims = get_worker_dimensions(params, _num_threads, tsplit);
#ifdef VERBOSE_THREADED
	if (_num_threads < num_threads) {
	    std::cerr << "More threads than ";
	    switch (wsplit) {
	    case DIST_SPLIT_RECORDS:
		std::cerr << "virtual blocks";
		break;
	    case DIST_SPLIT_RECORD_BYTES:
		std::cerr << "block bytes";
		break;
	    default:
		break;
	    }
	    std::cerr << ".  Only using " << _num_threads << " threads.\n";
	}
#endif
	worker_sids.resize(_num_threads, sid);
	thread_params = params->create_worker_params(thread_dims);
	for (nservers_t i = 0; i < thread_params.size(); ++i) {
	    if (thread_params[i]) {
		thread_serverparams.push_back(new PercyServerParams(
			thread_params[i], worker_sids[i], be_byzantine));
#ifdef VERBOSE_THREADED
		std::cerr << "Thread " << i << ": ";
		thread_serverparams.back()->print(std::cerr);
		std::cerr << "\n";
#endif
	    } else {
		thread_serverparams.push_back(NULL);
	    }
	}
    }
}

PercyServerParams::~PercyServerParams ()
{
    std::vector<const PercyParams*>::iterator it1;
    for (it1 = worker_params.begin(); it1 != worker_params.end(); ++it1) {
	if (*it1 != NULL && *it1 != params) {
	    delete *it1;
	}
    }
    for (it1 = thread_params.begin(); it1 != thread_params.end(); ++it1) {
	if (*it1 != NULL && *it1 != params) {
	    delete *it1;
	}
    }
    std::vector<const PercyServerParams*>::iterator it2;
    for (it2 = worker_serverparams.begin(); it2 != worker_serverparams.end(); ++it2) {
	if (*it2 != NULL && *it2 != this) {
	    delete *it2;
	}
    }
    for (it2 = thread_serverparams.begin(); it2 != thread_serverparams.end(); ++it2) {
	if (*it2 != NULL && *it2 != this) {
	    delete *it2;
	}
    }
}

void PercyServerParams::send (std::ostream &os, bool to_worker) const
{
    // Output the magic header
    if (to_worker) {
	os.write("PIRM", 4);
    } else {
	os.write("PIRS", 4);
    }

    // Send sid
    percy_write_le_uint16(os, sid);

    // Send rest
    params->write(os);
}

bool PercyServerParams::is_compatible (std::istream &is) const
{
    // Input the magic header
    unsigned char minibuf[4];
    is.read((char *)minibuf, 4);
    if (!memcmp(minibuf, "PIRC", 4) && !memcmp(minibuf, "PIRM", 4)) {
	std::cerr << "Did not find expected PercyParams header.\n";
	return false;
    }

    // Get the sid
    nservers_t readsid;
    percy_read_le_uint16(is, readsid);
    if (readsid != sid) {
	std::cerr << "sids did not match: local=" << sid << ", remote=" <<
		readsid << "\n";
	return false;
    }

    // Get rest of info
    return params->check_compatible(is);
}

void PercyServerParams::print (std::ostream& os) const 
{
    params->print(os);
    os << sid << ","
	<< (be_byzantine ? "byzantine" : "honest") << ",";
    print_distributed(os);
    os << ",";
}

void PercyServerParams::print_distributed (std::ostream& os) const
{
    if (_num_workers > 0) {
	os << "(DISTRIBUTED:" << _num_workers << ";" << wsplit << ")";
    } else if (_num_threads > 0) {
	os << "(" << (is_forked ? "FORKED" : "THREADED") << ";" << _num_threads
		<< ";" << tsplit << ")";
    } else {
	os << "None";
    }
}


