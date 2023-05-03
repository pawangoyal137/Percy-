// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>
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

#ifndef __PERCYSTATS_H__
#define __PERCYSTATS_H__

#include <sys/time.h>
#include <vector>
#include <map>
#include "percytypes.h"
#include "percyparams.h"

extern void sub_timevals (timeval& result, const timeval& t1, const timeval& t2);
extern std::ostream& operator<< (std::ostream& os, const timeval& t);

struct QueryBatchStats {
    // Batch Size
    nqueries_t batch_size;
    nqueries_t blocks_per_query;

    // Timing
    timeval start_time;
    timeval encode_start_time;
    timeval encode_done_time;
    timeval ctos_start_time;   // client -> server communication start
    timeval ctos_done_time;    // client -> server communication done
    timeval stoc_start_time;   // server -> client communication start
    timeval stoc_done_time;    // server -> client communication done
    timeval decode_start_time;
    timeval decode_done_time;
    timeval end_time;

    // Communication Amounts
    dbsize_t ctos_bytes; // client -> server bytes sent
    dbsize_t stoc_bytes; // server -> client bytes sent

    // Success/Failure
    bool num_unsuccessful;

    // Maximal strassen level reached
    nqueries_t strassen_level_reached;

    QueryBatchStats(): batch_size(), blocks_per_query(), start_time(), 
	encode_start_time(), encode_done_time(), ctos_start_time(), 
	ctos_done_time(), stoc_start_time(), stoc_done_time(), 
	decode_start_time(), decode_done_time(), end_time(), ctos_bytes(), 
	stoc_bytes(), num_unsuccessful(), strassen_level_reached() {}
};

class PercyStats {
public:
    PercyStats (const char * appendix = NULL);
    PercyStats (std::ostream& os, const char * appendix = NULL);
    virtual ~PercyStats ();

    void clear_batches () { query_batches.clear(); }

    nqueries_t start_query_batch (nqueries_t batch_size, 
	    nqueries_t blocks_per_query = 1);
    bool encode_start (nqueries_t batch_number);
    bool encode_done (nqueries_t batch_number);
    bool client_to_server_start (nqueries_t batch_number);
    bool client_to_server_done (nqueries_t batch_number, dbsize_t bytes_sent);
    bool server_to_client_start (nqueries_t batch_number);
    bool server_to_client_done (nqueries_t batch_number, dbsize_t bytes_sent);
    bool decode_start (nqueries_t batch_number);
    bool decode_done (nqueries_t batch_number, nqueries_t num_unsuccessful = 0);
    bool strassen_level_reached (nqueries_t batch_number,
	nqueries_t strassen_level_reached = 0);
    bool finish_query_batch (nqueries_t batch_number, bool print_head = false);

    bool time_communication_separately () const { return time_communication; }

    virtual void print_header () = 0;

protected:
    nqueries_t next_batch_number;

    // Query information
    std::map<nqueries_t, QueryBatchStats> query_batches;

    // Print functions
    virtual bool print (nqueries_t batch_number, bool print_head = false) = 0;
    bool print_and_delete (nqueries_t batch_number, bool print_head = false);
    void print_all (bool print_head = false);

    std::ostream& os;
    const char * appendix;
    bool time_communication;
};


class PercyClientStats : public PercyStats {
public:
    PercyClientStats (const PercyClientParams * params, 
	    const char * appendix = NULL);
    PercyClientStats (const PercyClientParams * params, std::ostream& os,
	    const char * appendix = NULL);
    virtual ~PercyClientStats ();

    virtual void print_header ();

protected:
    virtual bool print (nqueries_t batch_number, bool print_head = false);

    const PercyClientParams * params;
};


class PercyServerStats : public PercyStats {
public:
    PercyServerStats (const PercyServerParams * params, 
	    const char * appendix = NULL);
    PercyServerStats (const PercyServerParams * params, std::ostream& os,
	    const char * appendix = NULL);
    virtual ~PercyServerStats ();

    virtual void print_header ();

protected:
    virtual bool print (nqueries_t batch_number, bool print_head = false);

    const PercyServerParams * params;
};

#endif
