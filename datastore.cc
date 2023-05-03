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

#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fstream>
#include "datastore.h"
#include "percyio.h"
#include "streams.h"
#include "itparams.h"
#include "recursiveparams.h"

#define FILENAME_SIZE 120

DataStore::DataStore (unsigned char * database, 
	const PercyServerParams * serverparams)
:
    database(database),
    serverparams(serverparams),
    params(serverparams->percy_params()),
    block_size(params->server_block_size()),
    num_blocks(params->num_blocks()),
    database_size(params->num_blocks() * params->server_block_size()),
    num_workers(0)
{
    // Set up any sub datastores for distributed computation
    DistSplit split = DIST_SPLIT_RECORDS;
    std::vector<const PercyServerParams*> wserverparams;
    std::vector<const PercyParams*> wparams;
    if (serverparams->is_distributed()) {
	split = serverparams->worker_split();
	wserverparams = serverparams->get_all_worker_serverparams();
	wparams = serverparams->get_all_worker_params();
    } else if (serverparams->is_threaded()) {
	split = serverparams->thread_split();
	wserverparams = serverparams->get_all_thread_serverparams();
	wparams = serverparams->get_all_thread_params();
    }
    num_workers = wparams.size();

    switch (split) {
    case DIST_SPLIT_QUERIES:
	for (nservers_t i = 0; i < num_workers; ++i) {
	    subdb_offsets.push_back(0);
	    subdatastores.push_back(new DataStore(database, wserverparams[i]));
	}
	break;
    case DIST_SPLIT_RECORDS:
	{
	    dbsize_t offset = 0;
	    for (nservers_t i = 0; i < num_workers; ++i) {
		subdb_offsets.push_back(offset);
		subdatastores.push_back(
			new DataStore(database + offset, wserverparams[i]));
		dbsize_t worker_num_blocks = wparams[i]->num_blocks();
		offset += worker_num_blocks * block_size;
	    }
	}
	break;
    case DIST_SPLIT_RECORD_BYTES:
	for (nservers_t i = 0; i < num_workers; ++i) {
	    subdb_offsets.push_back(0);
	    subdatastores.push_back(NULL);
	}
	break;
    }
}

DataStore::~DataStore()
{
    std::vector<DataStore*>::iterator it;
    for (it = subdatastores.begin(); it != subdatastores.end(); ++it) {
	if (*it != NULL && *it != this) {
	    delete *it;
	}
    }
}

void DataStore::set_database (unsigned char * data)
{
    database = data;

    // Update worker datastores
    for (nservers_t i = 0; i < num_workers; ++i) {
	if (subdatastores[i] == this) continue;
	subdatastores[i]->set_database(data + subdb_offsets[i]);
    }
}


FileDataStore::FileDataStore(const char *filename, 
	const PercyServerParams * serverparams, dboffset_t offset)
:
    DataStore(0, serverparams)
{
    map_offset = 0;

    if (params->tau()) {
	// Process the PIRD header
	ifstream dbhead(filename);

	// We'll turn these into exceptions later
	char headbuf[6];
	dbhead.read(headbuf, 6);
	if (memcmp(headbuf, "PIRD\x01\x00", 6)) {
	    std::cerr << "Split database not in PIRD format.\n";
	    exit(1);
	}

	PercyMode taumode = params->get_mode();
	const PercyParams * tauparams = params;
	if (taumode == MODE_HYBRID) {
	    tauparams = dynamic_cast<const RecursiveParams*>(params)->get_iteration(0);
	    taumode = tauparams->get_mode();
	}
	switch (taumode) {
	case MODE_ZZ_P: {
		char c;
		dbhead.read(&c, 1);
		if (c != 'z') {
		    std::cerr << "Incorrent mode for split database.\n";
		    exit(1);
		}
		ZZ dbmodulus;
		percy_read_ZZ(dbhead, dbmodulus);
		if (!dynamic_cast<const ZZ_pParams*>(tauparams)->modulus_match(dbmodulus)) {
		    std::cerr << "Incorrect modulus for split database.\n";
		    exit(1);
		}
	    } break;
	case MODE_GF28:
	case MODE_GF216:
		char c;
		dbhead.read(&c, 1);
		if (c != 'g') {
		    std::cerr << "Incorrent mode for split database.\n";
		    exit(1);
		}
		dbsize_t dbws;
		PERCY_READ_LE_DBSIZE(dbhead, dbws);
		if (dbws != tauparams->word_size()) {
		    std::cerr << "Incorrect word size for split database.\n";
		    exit(1);
		}
	    break;
	default:
	    std::cerr << "Invalid mode for tau-independence: " << taumode << "\n";
	    exit(1);
	}

	map_offset = dbhead.tellg();
	dbhead.close();
    }

    // Open the file so we can mmap it
    dbfd = open(filename, O_RDONLY);

    if (dbfd < 0) {
	fprintf(stderr, "Could not open database %s\n", filename);
	exit(1);
    }
    struct stat st;
    fstat(dbfd, &st);

#if 0
    std::cerr << "database = " << filename << "\n";
    std::cerr << "st.st_size = " << st.st_size << "\n";
    std::cerr << "database_size = " << database_size << "\n";
    std::cerr << "offset = " << offset << "\n";
    std::cerr << "map_offset = " << map_offset << "\n";
#endif
    if ((dboffset_t)(st.st_size) < (dboffset_t)database_size + offset + 
	    map_offset) {
        fprintf(stderr, "Database too small!\n");
        exit(1);
    }

    long pagesize = sysconf(_SC_PAGE_SIZE);
    dboffset_t offset_pages = offset / pagesize;
    map_offset += offset % pagesize;

    mapptr = (unsigned char *)MMAP(NULL, database_size + map_offset, PROT_READ,
	    MAP_SHARED, dbfd, offset_pages * pagesize);
    if (mapptr == MAP_FAILED) {
	perror("mmap");
	exit(1);
    }
    set_database(mapptr + map_offset);
}

FileDataStore::~FileDataStore()
{
    munmap(mapptr, database_size + map_offset);
    close(dbfd);
}

