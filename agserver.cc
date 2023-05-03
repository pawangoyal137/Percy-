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

#include <NTL/mat_ZZ_p.h>
#include <cstdio>
#include <string.h>
#include <fstream>
#include "agserver.h"
#include "percyio.h"

#include <sys/time.h>

// This class implements the server-side computation for the computational PIR
// scheme introduced in
//
//     C. Aguilar-Melchor and P. Gaborit. A Lattice-Based 
//     Computationally-Efficient Private Information Retrieval Protocol. In 
//     WEWORC 2007, July 2007.
//
// and revisited in 
//
//     C. Aguilar Melchor, B. Crespin, P. Gaborit, V. Jolivet, and P. Rousseau.
//     HighSpeed Private Information Retrieval Computation on GPU. In
//     SECURWARE'08, pages 263-272, 2008.


// Initialize a server with the given DataStore and params.
PercyAGServer::PercyAGServer (DataStore * datastore, 
	const PercyServerParams * params, PercyStats * stats)
:
    PercyServer(datastore, params, stats),
    params(static_cast<const AGParams*>(params->percy_params())),
    p(this->params->p())
{}

PercyAGServer::~PercyAGServer () {}

// Handle a request.  Returns true if successfully completed, and false
// otherwise.
bool PercyAGServer::handle_request_impl (std::vector<unsigned char*> requests,
	std::vector<unsigned char*> responses)
{
    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }

    // Get parameters
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t word_size = params->word_size();

    // Check that there are not too many records for depth
    dbsize_t max_num_blocks = pow(2, word_size-2);
    if (num_virtual_blocks >= max_num_blocks) {
	std::cerr << "The maximum number of blocks for this word size is " 
		<< max_num_blocks - 1 << ".  You have " << num_virtual_blocks 
		<< "\n";
	return false;
    }

    // Get datastore pointer
    const unsigned char * data = datastore->get_data();

    for (nqueries_t qind = 0; qind < num_queries; ++qind) {
	bool ret;
	switch (word_size) {
	case 16:
	    ret = handle_request_16(requests[qind], responses[qind], data);
	    if (!ret) {
		return false;
	    }
	    break;
	case 20:
	    ret = handle_request_20(requests[qind], responses[qind], data);
	    if (!ret) {
		return false;
	    }
	    break;
	default:
	    std::cerr << "Invalid word size: " << word_size << "\n";
	    return false;
	}
    }

    return true;
}

bool PercyAGServer::handle_request_16 (unsigned char * request, 
	unsigned char * response, const unsigned char * database)
{
    // Get parameters
    dbsize_t num_blocks = params->num_blocks();
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t block_size = params->block_size();
    dbsize_t virtual_block_size = params->virtual_block_size();
    dbsize_t N = params->N();

    uint64_t * query = new uint64_t[2 * N * N];

    dbsize_t num_rows_per_record = params->block_rows();
    dbsize_t last_row_size = block_size % (2 * N);
    dbsize_t num_full_rows = ( last_row_size == 0 ? num_rows_per_record : num_rows_per_record - 1 );
    dbsize_t result_words = 2 * N * num_rows_per_record * virtual_block_size;
    AG_Element * result = new AG_Element[result_words];
#ifndef NEED_UINT128
    memset(result, 0, 16 * result_words);
#endif
    const unsigned char * database_current = database;
    unsigned char * request_location = request;
    dbsize_t block_index = 0;
    uint64_t mask48 = 1;
    mask48 <<= 48;
    mask48 -= 1;
    for (dbsize_t a = 0; a < num_virtual_blocks; ++a) {
	// Read query
	for (dbsize_t i = 0; i < 2*N*N; ++i) {
	    query[i] = (*(uint64_t*)(request_location + 6 * i)) & mask48;
	}
	request_location += 12 * N * N;

	// Apply query to database
	AG_Element * result_current = result;
	for (dbsize_t b = 0; 
		block_index < num_blocks && b < virtual_block_size;
		++b, ++block_index) {
	    uint64_t * query_current;
	    for (dbsize_t i = 0; i < num_full_rows; ++i) {
		query_current = query;
		for (dbsize_t j = 0; j < N; ++j) {
		    uint16_t datacell = (*(uint16_t*)database_current);
		    database_current += 2;
		    for (dbsize_t k = 0; k < 2*N; ++k) {
			result_current[k] += UINT64_TO_AGELT(query_current[k]) * datacell;
		    }
		    query_current += 2*N;
		}
		result_current += 2*N;
	    }
	    if (last_row_size > 0) {
		query_current = query;
		for (dbsize_t bl = last_row_size; bl > 0; ) { // decrement inside
		    uint16_t datacell;
		    if (bl == 1) {
			datacell = (*(uint16_t*)database_current) & 0xFF;
			database_current += 1;
			bl = 0;
		    } else {
			datacell = (*(uint16_t*)database_current);
			database_current += 2;
			bl -= 2;
		    }
		    for (dbsize_t k = 0; k < 2*N; ++k) {
			result_current[k] += UINT64_TO_AGELT(query_current[k]) * datacell;
		    }
		    query_current += 2*N;
		}
		result_current += 2*N;
	    }
	}
    }

    // Apply modulo
    for (dbsize_t i = 0; i < result_words; ++i) {
	result[i] %= p;
    }

    // Convert result to str
    memset(response, 0, 6 * result_words);
    for (dbsize_t i = 0; i < result_words - 1; ++i) {
	*(uint64_t*)(response + 6 * i) |= AGELT_TO_UINT64(result[i]);
    }
    *(uint64_t*)(response + 6 * result_words - 8) |= AGELT_TO_UINT64((result[result_words-1] << 16));

    // Clean up
    delete[] result;
    delete[] query;

    return true;
}

bool PercyAGServer::handle_request_20 (unsigned char * request, 
	unsigned char * response, const unsigned char * database)
{
    // Get parameters
    dbsize_t num_blocks = params->num_blocks();
    dbsize_t num_virtual_blocks = params->num_virtual_blocks();
    dbsize_t block_size = params->block_size();
    dbsize_t virtual_block_size = params->virtual_block_size();
    dbsize_t N = params->N();

    uint64_t * query = new uint64_t[2 * N * N];

    dbsize_t num_rows_per_record = params->block_rows();
    dbsize_t last_row_size = block_size % (5 * N / 2);
    dbsize_t num_full_rows = ( last_row_size == 0 ? num_rows_per_record : num_rows_per_record - 1 );
    dbsize_t result_words = 2 * N * num_rows_per_record * virtual_block_size;
    AG_Element * result = new AG_Element[result_words];
#ifndef NEED_UINT128
    memset(result, 0, 16 * result_words);
#endif
    const unsigned char * database_current = database;
    unsigned char * request_location = request;
    dbsize_t block_index = 0;
    uint64_t mask60 = 1;
    mask60 <<= 60;
    mask60 -= 1;
    for (dbsize_t a = 0; a < num_virtual_blocks; ++a) {
	// Read query
	for (dbsize_t i = 0; i < N*N; ++i) {
	    query[2*i] = (*(uint64_t*)(request_location + 15 * i)) & mask60;
	    query[2*i+1] = (*(uint64_t*)(request_location + 15 * i + 7)) >> 4;
	}
	request_location += 15 * N * N;

	// Apply query to database
	AG_Element * result_current = result;
	for (dbsize_t b = 0; 
		block_index < num_blocks && b < virtual_block_size; 
		++b, ++block_index) {
	    uint64_t * query_current;
	    for (dbsize_t i = 0; i < num_full_rows; ++i) {
		query_current = query;
		for (dbsize_t j = 0; j < N/2; ++j) {
		    uint32_t datacell = (*(uint32_t*)database_current) & 0xFFFFF;
		    database_current += 2;
		    for (dbsize_t k = 0; k < 2*N; ++k) {
			result_current[k] += UINT64_TO_AGELT(query_current[k]) * datacell;
		    }
		    query_current += 2*N;
		    datacell = ((*(uint32_t*)database_current) >> 4) & 0xFFFFF;
		    database_current += 3;
		    for (dbsize_t k = 0; k < 2*N; ++k) {
			result_current[k] += UINT64_TO_AGELT(query_current[k]) * datacell;
		    }
		    query_current += 2*N;
		}
		result_current += 2*N;
	    }
	    if (last_row_size > 0) {
		query_current = query;
		for (dbsize_t bl = last_row_size; bl > 0; ) { // decrement inside
		    uint32_t datacell;
		    if (bl == 1) {
			datacell = (*(uint32_t*)database_current) & 0xFF;
			database_current += 1;
			bl = 0;
		    } else if (bl == 2) {
			datacell = (*(uint32_t*)database_current) & 0xFFFF;
			database_current += 2;
			bl = 0;
		    } else {
			datacell = (*(uint32_t*)database_current) & 0xFFFFF;
			database_current += 2;
			bl -= 2;
		    }
		    for (dbsize_t k = 0; k < 2*N; ++k) {
			result_current[k] += UINT64_TO_AGELT(query_current[k]) * datacell;
		    }
		    query_current += 2*N;
		    if (bl == 0) {
			continue;
		    } else if (bl == 1) {
			datacell = ((*(uint32_t*)database_current) >> 4) & 0xF;
			database_current += 1;
			bl = 0;
		    } else if (bl == 2) {
			datacell = ((*(uint32_t*)database_current) >> 4) & 0xFFF;
			database_current += 2;
			bl = 0;
		    } else {
			datacell = ((*(uint32_t*)database_current) >> 4) & 0xFFFFF;
			database_current += 3;
			bl -= 3;
		    }
		    for (dbsize_t k = 0; k < 2*N; ++k) {
			result_current[k] += UINT64_TO_AGELT(query_current[k]) * datacell;
		    }
		    query_current += 2*N;
		}
		result_current += 2*N;
	    }
	}
    }

    // Apply modulo
    for (dbsize_t i = 0; i < result_words; ++i) {
	result[i] %= p;
    }

    // Convert result to str
    memset(response, 0, 15 * result_words / 2);
    for (dbsize_t i = 0; i < result_words / 2; ++i) {
	*(uint64_t*)(response + 15 * i) |= AGELT_TO_UINT64(result[2*i]);
	*(uint64_t*)(response + 15 * i + 7) |= AGELT_TO_UINT64((result[2*i+1] << 4));
    }

    // Clean up
    delete[] result;
    delete[] query;

    return true;
}

void PercyAGServer::combine_results (unsigned char * result, 
	std::vector<unsigned char*> worker_results)
{
    dbsize_t response_size = params->response_size();
    dbsize_t word_size = params->word_size();
    dbsize_t num_words = response_size * 8 / (3 * word_size);
    uint64_t mask = 1;
    AG_Element * result_ag = new AG_Element[num_words];
#ifndef NEED_UINT128
    memset(result_ag, 0, 16 * num_words);
#endif

    // Add results
    switch (word_size) {
    case 16: {
	std::vector<unsigned char*>::iterator writ;
	mask <<= 48;
	mask -= 1;
	for (writ = worker_results.begin(); writ != worker_results.end(); ++writ) {
	    for (dbsize_t w = 0; w < num_words-1; ++w) {
		result_ag[w] += (*(uint64_t*)(*writ + 6 * w)) & mask;
	    }
	    result_ag[num_words-1] += (*(uint64_t*)(*writ + 6 * num_words - 8) >> 16);
	}
    } break;
    case 20: {
	std::vector<unsigned char*>::iterator writ;
	mask <<= 60;
	mask -= 1;
	for (writ = worker_results.begin(); writ != worker_results.end(); ++writ) {
	    for (dbsize_t w = 0; w < num_words/2; ++w) {
		result_ag[2*w] += (*(uint64_t*)(*writ + 15 * w)) & mask;
		result_ag[2*w+1] += (*(uint64_t*)(*writ + 15 * w + 7)) >> 4;
	    }
	}
    } break;
    default:
	break;
    }

    // Reduce modulo p
    AG_Element p = params->p();
    for (dbsize_t w = 0; w < num_words; ++w) {
	result_ag[w] %= p;
    }

    // Convert back to bytes
    switch (word_size) {
    case 16: {
	for (dbsize_t w = 0; w < num_words-1; ++w) {
	    *(uint64_t*)(result + 6 * w) |= AGELT_TO_UINT64(result_ag[w]);
	}
	*(uint64_t*)(result + 6 * num_words - 8) |= AGELT_TO_UINT64((result_ag[num_words-1] << 16));
    } break;
    case 20: {
	for (dbsize_t w = 0; w < num_words/2; ++w) {
	    *(uint64_t*)(result + 15 * w) |= AGELT_TO_UINT64(result_ag[2*w]);
	    *(uint64_t*)(result + 15 * w + 7) |= AGELT_TO_UINT64((result_ag[2*w+1] << 4));
	}
    } break;
    default:
	break;
    }

    delete[] result_ag;
}

