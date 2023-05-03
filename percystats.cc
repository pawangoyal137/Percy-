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

#include <unistd.h>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include "percystats.h"

// result = t1 - t2
void sub_timevals (timeval& result, const timeval& t1, const timeval& t2)
{
    result.tv_sec = t1.tv_sec - t2.tv_sec;
    result.tv_usec = t1.tv_usec - t2.tv_usec;
    if (result.tv_usec < 0 && result.tv_sec > 0) {
	result.tv_sec--;
	result.tv_usec += 1000000;
    } else if (result.tv_usec > 0 && result.tv_sec < 0) {
	result.tv_sec++;
	result.tv_usec -= 1000000;
    }
}

std::ostream& operator<< (std::ostream& os, const timeval& t)
{
    return os << t.tv_sec << "." << std::setfill('0') << std::setw(6) << t.tv_usec << std::setw(0);
}

PercyStats::PercyStats (const char * appendix)
:
    next_batch_number(1),
    os(std::cerr),
    appendix(appendix)
{
    time_communication = (getenv("PIR_TIME_COMM") != NULL);
}

PercyStats::PercyStats (std::ostream& os, const char * appendix)
:
    next_batch_number(1),
    os(os),
    appendix(appendix)
{
    time_communication = (getenv("PIR_TIME_COMM") != NULL);
}

PercyStats::~PercyStats ()
{}

bool PercyStats::print_and_delete (nqueries_t batch_number, bool print_head)
{
    bool ret = print(batch_number, print_head);
    if (ret) {
	query_batches.erase(batch_number);
    }
    return ret;
}

void PercyStats::print_all (bool print_head)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    for (iter = query_batches.begin(); iter != query_batches.end(); ++iter) {
	print(iter->first, print_head);
	print_head = false;
    }
}

nqueries_t PercyStats::start_query_batch (nqueries_t batch_size,
	nqueries_t blocks_per_query)
{
    nqueries_t start_batch_number = next_batch_number;
    nqueries_t batch_number = next_batch_number++;
    while (query_batches.find(batch_number) != query_batches.end()) {
	if (next_batch_number == start_batch_number) {
	    return 0;
	}
	if (next_batch_number == 0) {
	    next_batch_number = 1;
	}
	batch_number = next_batch_number++;
    }
    query_batches[batch_number] = QueryBatchStats();
    QueryBatchStats& qbstats = query_batches[batch_number];
    qbstats.batch_size = batch_size;
    qbstats.blocks_per_query = blocks_per_query;
    gettimeofday(&qbstats.start_time, NULL);
    return batch_number;
}

bool PercyStats::encode_start (nqueries_t batch_number)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.encode_start_time), NULL);
    return true;
}

bool PercyStats::encode_done (nqueries_t batch_number)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.encode_done_time), NULL);
    return true;
}

bool PercyStats::client_to_server_start (nqueries_t batch_number)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.ctos_start_time), NULL);
    return true;
}

bool PercyStats::client_to_server_done (nqueries_t batch_number, dbsize_t bytes_sent)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.ctos_done_time), NULL);
    iter->second.ctos_bytes = bytes_sent;
    return true;
}

bool PercyStats::server_to_client_start (nqueries_t batch_number)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.stoc_start_time), NULL);
    return true;
}

bool PercyStats::server_to_client_done (nqueries_t batch_number, dbsize_t bytes_sent)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.stoc_done_time), NULL);
    iter->second.stoc_bytes = bytes_sent;
    return true;
}

bool PercyStats::decode_start (nqueries_t batch_number)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.decode_start_time), NULL);
    return true;
}

bool PercyStats::decode_done (nqueries_t batch_number, nqueries_t
	num_unsuccessful)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.decode_done_time), NULL);
    iter->second.num_unsuccessful = num_unsuccessful;
    return true;
}

bool PercyStats::strassen_level_reached (nqueries_t batch_number,
	nqueries_t strassen_level_reached)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    iter->second.strassen_level_reached = strassen_level_reached;
    return true;
}

bool PercyStats::finish_query_batch (nqueries_t batch_number, bool print_head)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (iter == query_batches.end()) {
	return false;
    }
    gettimeofday(&(iter->second.end_time), NULL);
    print_and_delete(batch_number);
    return true;
}


PercyClientStats::PercyClientStats (const PercyClientParams * params, 
	const char * appendix)
:
    PercyStats(appendix),
    params(params)
{}

PercyClientStats::PercyClientStats (const PercyClientParams * params, 
	std::ostream& os, const char * appendix)
:
    PercyStats(os, appendix),
    params(params)
{}

PercyClientStats::~PercyClientStats ()
{}

bool PercyClientStats::print (nqueries_t batch_number, bool print_head)
{
    std::map<nqueries_t, QueryBatchStats>::iterator iter;
    iter = query_batches.find(batch_number);
    if (query_batches.find(batch_number) == query_batches.end()) {
	return false;
    }

    if (print_head) {
	print_header();
    }

    // Client params
    params->print(os);

    // Query properties
    QueryBatchStats& qbstats = query_batches[batch_number];
    os << qbstats.batch_size << ",";
    os << qbstats.blocks_per_query << ",";
    timeval diff;
    sub_timevals(diff, qbstats.encode_done_time, qbstats.encode_start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    sub_timevals(diff, qbstats.ctos_done_time, qbstats.ctos_start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    os << qbstats.ctos_bytes << ",";
    sub_timevals(diff, qbstats.stoc_start_time, qbstats.ctos_done_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    sub_timevals(diff, qbstats.stoc_done_time, qbstats.stoc_start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    os << qbstats.stoc_bytes << ",";
    sub_timevals(diff, qbstats.decode_done_time, qbstats.decode_start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    sub_timevals(diff, qbstats.end_time, qbstats.start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    os << qbstats.num_unsuccessful << ",";

    // Appendix
    if (appendix) {
	os << appendix;
    }

    os << "\n";
    os.flush();

    return true;
}

void PercyClientStats::print_header()
{
    os << "#";
    os << "mode,num_blocks,block_size,word_size,tau,virtual_block_size,mode_specific,";
    os << "num_servers,";
    os << "batch_size,";
    os << "blocks_per_query,";
    os << "encode_time,";
    os << "request_time,";
    os << "request_size,";
    os << "server_time,";
    os << "response_time,";
    os << "response_size,";
    os << "decode_time,";
    os << "total_time,";
    os << "num_unsuccessful,";
    os << "\n";
    os.flush();
}


PercyServerStats::PercyServerStats (const PercyServerParams * params, 
	const char * appendix)
:
    PercyStats(appendix),
    params(params)
{}

PercyServerStats::PercyServerStats (const PercyServerParams * params, 
	std::ostream& os, const char * appendix)
:
    PercyStats(os, appendix),
    params(params)
{}

PercyServerStats::~PercyServerStats ()
{}

bool PercyServerStats::print (nqueries_t batch_number, bool print_head)
{
    if (query_batches.find(batch_number) == query_batches.end()) {
	return false;
    }

    if (print_head) {
	print_header();
    }

    // Server params
    params->print(os);

    // Query properties
    QueryBatchStats& qbstats = query_batches[batch_number];
    os << qbstats.batch_size << ",";
    timeval diff;
    sub_timevals(diff, qbstats.ctos_done_time, qbstats.ctos_start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    os << qbstats.ctos_bytes << ",";
    sub_timevals(diff, qbstats.stoc_start_time, qbstats.ctos_done_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    sub_timevals(diff, qbstats.stoc_done_time, qbstats.stoc_start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    os << qbstats.stoc_bytes << ",";
    sub_timevals(diff, qbstats.end_time, qbstats.start_time);
    os << diff.tv_sec << "." << std::setfill('0') << std::setw(6) << diff.tv_usec << std::setw(0) << ",";
    os << qbstats.num_unsuccessful << ",";
    os << qbstats.strassen_level_reached << ",";

    // Appendix
    if (appendix) {
	os << appendix;
    }

    os << "\n";
    os.flush();

    // Output
    return true;
}

void PercyServerStats::print_header()
{
    os << "#";
    os << "mode,num_blocks,block_size,word_size,tau,virtual_block_size,mode_specific,";
    os << "sid,byzantine,distributed,";
    os << "batch_size,";
    os << "request_time,";
    os << "request_size,";
    os << "server_time,";
    os << "response_time,";
    os << "response_size,";
    os << "total_time,";
    os << "num_unsuccessful,";
    os << "strassen_level,";
    os << "\n";
    os.flush();
}

