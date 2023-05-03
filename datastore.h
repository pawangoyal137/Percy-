// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
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

#ifndef __DATASTORE_H__
#define __DATASTORE_H__
#include <ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <set>
#include <map>
#include <vector>
#include "percyparams.h"

NTL_CLIENT

/// A simple database object.  Stores a database that is a contiguous block
/// of memory.

class DataStore {
private:
    unsigned char *database;

protected:
    /// Parameters for the server.
    const PercyServerParams * serverparams;
    /// Parameters for the protocol.
    const PercyParams * params;
    /// Size of each database block in bytes.
    dbsize_t block_size;
    /// Number of blocks in the database.
    dbsize_t num_blocks;
    /// Total size of the database in bytes.
    dbsize_t database_size;

    // Distributed computation stuff
    /// Number of threads/workers.
    nservers_t num_workers;
    /// The distance from the beginning of the database that each 
    /// thread/worker's portion of the database begins.
    std::vector<dbsize_t> subdb_offsets;
    /// The database objects for the threads/workers.
    std::vector<DataStore*> subdatastores;

public:
    /// Constructor.
    /// @param database	A pointer to the database in memory.
    /// @param params	The server parameters.
    DataStore (unsigned char * database, const PercyServerParams * params);

    /// Destructor
    virtual ~DataStore ();

    /// Get a pointer to the database in memory.
    /// @return	The database pointer.
    const unsigned char *get_data() const {return database;}

    /// Change the database pointer.
    /// @param data The location of the new database in memory.
    void set_database (unsigned char * data = NULL);

    /// Get the database objects for the threads/workers.
    /// @return	The database objects for the workers.
    vector<DataStore*> get_worker_datastores () const { 
	return subdatastores;
    }
    /// Get the database object for a thread/worker.
    /// @param worker_index The index of a worker.
    /// @return		    The database object for the worker.
    DataStore * get_worker_datastore (nservers_t worker_index) const {
	return (worker_index < num_workers ? subdatastores[worker_index] : NULL);
    }
};


/// A database that is backed by one contiguous file.

class FileDataStore : public DataStore {
private:
    int dbfd;
    unsigned char *mapptr;
    dboffset_t map_offset;

public:
    /// Constructor.
    /// @param filename	The name of the database file.
    /// @param params	The server parameters.
    /// @param offset	The distance from the beginning of the file to the start
    ///                 of the database.  (Default: 0)
    FileDataStore(const char *filename, const PercyServerParams * params,
	    dboffset_t offset = 0);

    /// Destructor.
    virtual ~FileDataStore();
};

#endif
