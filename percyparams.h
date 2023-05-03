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

#ifndef __PERCYPARAMS_H__
#define __PERCYPARAMS_H__

#include <iostream>
#include <vector>
#include <utility>
#include "version.h"
#include "percytypes.h"

/** @file
Defines the basic structure of protocol parameters (PercyParams), client
parameters (PercyClientParams), and server parameters (PercyServerParams).
*/

/******************************************************************************
  ENUMERATED TYPES AND TYPEDEFS
******************************************************************************/

/// A PIR protocol.
enum PercyMode {
    MODE_NONE = 0,
    MODE_ZZ_P,		///< Goldberg (2007) IT-PIR over the integers modulo p.
    MODE_GF28,		///< Goldberg (2007) IT-PIR over GF(2^8).
    MODE_GF216,		///< Goldberg (2007) IT-PIR over GF(2^16).
    MODE_CHOR,		///< Chor et al. (1995) IT-PIR.
    MODE_AG,		///< Aguilar Melchor and Gaborit (2007) CPIR.
    MODE_RECURSIVE_AG,	///< Aguilar Melchor and Gaborit (2007) CPIR with recursion.
    MODE_HYBRID		///< Devet and Goldberg (2014) Hybrid PIR.
};

/// Strings associated with PercyMode values.
extern const char * PercyModeStrings[];

/// Prints a PercyMode string to a stream.
extern std::ostream& operator<< (std::ostream &os, PercyMode mode);

/// The method used to partition work between threads/workers.
enum DistSplit {
    DIST_SPLIT_QUERIES,	    ///< Partition the queries.
    DIST_SPLIT_RECORDS,	    ///< Partition the database records.
    DIST_SPLIT_RECORD_BYTES ///< Partition the bytes of the records.
};

/// Strings associated with DistSplit values.
extern const char * DistSplitStrings[];

/// Prints a DistSplit string to a stream.
extern std::ostream& operator<< (std::ostream &os, DistSplit mode);

/// A 2-dimensional coordinate.
typedef std::pair<dbsize_t, dbsize_t> Dimension;


/******************************************************************************
  CLASS DEFINITIONS
******************************************************************************/

/// An abstract base class for a protocol's parameters.

class PercyParams {
public:
    /// Constructor.
    /// @param num_blocks	    Number of blocks in the database.
    /// @param block_size	    Size of the blocks in the database in bytes.
    /// @param word_size	    Word size used to split the blocks.
    /// @param mode		    PIR protocol to use.
    /// @param tau		    Databases are tau-independent.  No coalition
    ///				    of up to tau servers can determine the 
    ///				    contents of the database.
    /// @param virtual_block_size   Number of actual blocks contained in each
    ///				    virtual block when being used as one 
    ///				    iteration of a recursive protocol.
    PercyParams (dbsize_t num_blocks, dbsize_t block_size, dbsize_t word_size,
	    PercyMode mode, nservers_t tau = 0, 
	    dbsize_t virtual_block_size = 1);

    /// Destructor
    virtual ~PercyParams () {}

    // Access parameters
    /// Get the number of blocks in the database.
    dbsize_t num_blocks () const { return _num_blocks; }
    /// Get the size of each block in the database in bytes.
    dbsize_t block_size () const { return _block_size; }
    /// Get the size of each block in the database that the datastore will 
    /// actually use.  This may be different than block_size() when 
    /// tau-independence is used.
    virtual dbsize_t server_block_size () const { return _block_size; }
    /// Get the word size used to split blocks.
    dbsize_t word_size () const { return _word_size; }
    /// Get the number of words per database block.
    dbsize_t words_per_block () const { return _words_per_block; }
    /// Get the protocol being used.
    PercyMode get_mode () const { return mode; }
    /// Get the level of tau-independence.  A value of zero indicates that any
    /// server can determine the contents of the database.
    nservers_t tau () const { return _tau; }
    /// Get the number of virtual blocks when being used as one iteration of a
    /// recursive protocol.
    dbsize_t num_virtual_blocks () const { return _num_virtual_blocks; }
    /// Get the number of actual blocks in a virtual block when being used as 
    /// one iteration of a recursive protocol.
    dbsize_t virtual_block_size () const { return _virtual_block_size; }

    /// Get the size of a client to server request.
    /// @param num_queries  The number of queries in the request.  (Default: 1)
    virtual dbsize_t request_size (nqueries_t num_queries = 1) const = 0;
    /// Get the size of a server to client response.
    /// @param num_queries  The number of queries in the request.  (Default: 1)
    virtual dbsize_t response_size (nqueries_t num_queries = 1) const = 0;

    /// Prints the parameters in CSV form.  The fields are:
    /// <pre>mode,num_blocks,block_size,word_size,tau,virtual_block_size,mode_specific,</pre>
    /// where mode_specific is the result of calling print_mode_specific()
    /// @param os   Stream to print to.
    void print (std::ostream& os) const;

    /// Print mode-specific parameters.
    /// Meant to be overloaded by subclasses to print more details about that
    /// mode.  Is used by print() to print the mode_specific column.
    /// @param os   Stream to print to.
    virtual void print_mode_specific (std::ostream& os) const {}

    /// Write the parameters to a stream to check compatibility.
    /// @param os   Stream to write to.
    virtual void write (std::ostream &os) const;

    /// Read the parameters from a stream (as written by write()) and check that
    /// they are compatible with these parameters.
    /// @param is   Stream to read from.
    /// @return	    Return true if the read parameters are compatitble with 
    ///		    these parameters; false otherwise.
    virtual bool check_compatible (std::istream &is) const;

    /// Create protocol parameters for threads/workers.
    /// @params worker_dims The database dimensions for each thread/worker.
    /// @return		    The protocol parameters for each thread/worker.
    virtual std::vector<const PercyParams*> create_worker_params (
	    std::vector<Dimension> worker_dims) const = 0;

    /// Check if the protocol is recursive.
    virtual bool is_recursive () const { return false; }

protected:
    /// Percy++ version.
    unsigned char version[3];
    /// Number of database blocks.
    dbsize_t _num_blocks;
    /// Size of database blocks in bytes.
    dbsize_t _block_size;
    /// Word size used to split blocks.
    dbsize_t _word_size;
    /// Number of words per database block.
    dbsize_t _words_per_block;
    /// Protocol being used.
    PercyMode mode;
    /// Level of tau-independence.
    nservers_t _tau;
    /// Number of virtual blocks when part of a recursive protocol.
    dbsize_t _num_virtual_blocks;
    /// Number of actual blocks per virtual block when part of a recursive
    /// protocol.
    dbsize_t _virtual_block_size;
};


/// Client parameters.

class PercyClientParams {
public:
    /// Constructor.
    /// @param params	    The protocol parameters.
    /// @param num_servers  The number of servers being queried.
    /// @param is_null	    If true, the client does nothing except send random
    ///			    data to the server(s) and read the response(s).
    PercyClientParams (const PercyParams * params, nservers_t num_servers = 1, 
	    bool is_null = false);

    /// Destructor
    virtual ~PercyClientParams () {};

    /// Send parameters and sid to server to test compatibility.
    /// @param os   Stream to write to.
    /// @param sid  Server ID of server being sent params.
    void send (std::ostream &os, nservers_t sid) const;

    /// Test compatibility with a server.  Read parameters from a server and
    /// check compatibility with the server.
    /// @param is   Stream to read from.
    /// @param sid  Server ID of server begin checked.
    /// @return	    Returns true if server is compatible with this client; false
    ///		    otherwise.
    bool is_compatible (std::istream &is, nservers_t sid) const;

    // Access parameters
    /// Get the number of servers being queried.
    nservers_t num_servers () const { return _num_servers; }

    /// Prints the parameters in CSV form.  Prints the protocol parameters and
    /// add the following client fields:
    /// <pre>num_servers,</pre>
    /// @param os   Stream to print to.
    void print (std::ostream& os) const;

    /// Get the protocol parameters.
    const PercyParams * percy_params () const { return params; }

    /// Check if the client is a null client.  If true, the client does nothing
    /// except send random data to the server(s) and read the response(s).
    bool is_null () const { return null; }

protected:
    /// The protocol parameters.
    const PercyParams * params;

    /// The number of servers being queried.
    nservers_t _num_servers;
    /// Whether or not the client a null client.
    bool null;
};


/// Information needed to connect to a server.
struct serverinfo {
    nservers_t sid; ///< Server ID
    char * addr;    ///< Address
    uint16_t port;  ///< Port
};

/// Server parameters.
class PercyServerParams {
public:
    /// Constructor
    /// @param params	    The protocol parameters.
    /// @param sid	    The server ID.
    /// @param be_byzantine If true, do not send correct responses to the client.
    PercyServerParams (const PercyParams * params, nservers_t sid,
	    bool be_byzantine = false);

    /// Constructor for parallel server computation.
    /// @param params	    The protocol parameters.
    /// @param sid	    The server ID.
    /// @param num_threads  The number of threads to distribute computation over.
    /// @param tsplit	    The method of partitioning work between threads.
    /// @param num_workers  The number of workers to distributed computation over.
    /// @param wsplit	    The method of partitioning work between workers.
    /// @param worker_sids  The server IDs of workers.
    /// @param fork	    If true, use process forking instead of threading.
    /// @param be_byzantine If true, do not send correct responses to the client.
    PercyServerParams (const PercyParams * params, nservers_t sid,
	    nservers_t num_threads = 0, DistSplit tsplit = DIST_SPLIT_RECORDS,
	    nservers_t num_workers = 0, DistSplit wsplit = DIST_SPLIT_RECORDS,
	    std::vector<nservers_t> worker_sids = std::vector<nservers_t>(),
	    bool fork = false, bool be_byzantine = false);

    /// Destructor.
    virtual ~PercyServerParams ();

    /// Send parameters and sid to client to test compatibility.
    /// @param os	    Stream to write to.
    /// @param to_worker    If true, the parameters are being sent to a worker.
    void send (std::ostream &os, bool to_worker = false) const;

    /// Test compatibility with client.  Read parameters and a server ID from a
    /// client and test compatibility.
    /// @param is   Stream to read from.
    /// @return	    Returns true if the client is compatible with this server;
    ///		    false otherwise.
    bool is_compatible (std::istream &is) const;

    /// Check if the server is recursive.
    bool is_recursive () const { return params->is_recursive(); }

    // Access parameters
    /// Get the server ID.
    nservers_t get_sid () const { return sid; }
    /// Check if the server is Byzantine (responses incorrectly/maliciously).
    bool is_byzantine () const { return be_byzantine; }

    /// Check if the server is using multithreading.
    bool is_threaded () const {
	return (_num_workers == 0 && _num_threads > 0);
    }
    /// Get the number of threads being used.
    nservers_t num_threads () const { return _num_threads; }
    /// Get the method of partitioning work amoung threads.
    DistSplit thread_split () const { return tsplit; }
    /// Check if the server is using forked processes instead of threading.
    bool use_forked_threads () const { return is_forked; }

    /// Check if the server is using distributed computation.
    bool is_distributed () const { 
	return _num_workers > 0; 
    }
    /// Get the number of worker processes being used.
    nservers_t num_workers () const { return _num_workers; }
    /// Get the method of partitioning work amoung workers.
    DistSplit worker_split () const { return wsplit; }

    /// Get the protocol parameters for a worker.
    /// @param wid  Index of the worker.
    const PercyParams * get_worker_params (nservers_t wid) const {
	return (wid < worker_params.size() ? worker_params[wid] : NULL);
    }
    /// Get the protocol parameters for all workers.
    std::vector<const PercyParams*> get_all_worker_params () const {
	return worker_params;
    }

    /// Get the server parameters for a worker.
    /// param wid   Index of the worker.
    const PercyServerParams * get_worker_serverparams (nservers_t wid) const {
	return (wid < worker_serverparams.size() ? worker_serverparams[wid] : NULL);
    }
    /// Get the server parameters for all workers.
    std::vector<const PercyServerParams*> get_all_worker_serverparams () const {
	return worker_serverparams;
    }

    /// Get the protocol parameters for a thread.
    /// @param wid  Index of the thread.
    const PercyParams * get_thread_params (nservers_t wid) const {
	return (wid < thread_params.size() ? thread_params[wid] : NULL);
    }
    /// Get the protocol parameters for all threads.
    std::vector<const PercyParams*> get_all_thread_params () const {
	return thread_params;
    }

    /// Get the server parameters for a thread.
    /// @param wid  Index of the thread.
    const PercyServerParams * get_thread_serverparams (nservers_t wid) const {
	return (wid < thread_serverparams.size() ? thread_serverparams[wid] : NULL);
    }
    /// Get the server parameters for all threads.
    std::vector<const PercyServerParams*> get_all_thread_serverparams () const {
	return thread_serverparams;
    }

    /// Prints the parameters in CSV form.  Prints the protocol parameters and
    /// add the following server fields:
    /// <pre>sid,byzantine,distributed_information,</pre>
    /// where distributed_information is the result of print_distributed().
    /// @param os   Stream to print to.
    void print (std::ostream& os) const;

    /// Print the distributed properties of the server
    virtual void print_distributed (std::ostream& os) const;

    /// Get the protocol parameters.
    const PercyParams * percy_params () const { return params; }

protected:
    /// The protocol parameters.
    const PercyParams * params;

    /// Server ID
    nservers_t sid;
    /// Whether or not the server is Byzantine.
    bool be_byzantine;

    // Distributed members
    /// The number of threads
    nservers_t _num_threads;
    /// The method of partitioning work between threads
    DistSplit tsplit;
    /// Whether or not the server is using forked processes instead of threading
    bool is_forked;
    /// The number of workers
    nservers_t _num_workers;
    /// The method of partitioning work between workers
    DistSplit wsplit;

    /// The workers' protocol parameters
    std::vector<const PercyParams*> worker_params;
    /// The workers' server parameters
    std::vector<const PercyServerParams*> worker_serverparams;

    /// The threads' protocol parameters
    std::vector<const PercyParams*> thread_params;
    /// The threads' server parameters
    std::vector<const PercyServerParams*> thread_serverparams;
};

#endif
