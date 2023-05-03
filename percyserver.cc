// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Ryan Henry <rhenry@cs.uwaterloo.ca>,
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <pthread.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#include "percytypes.h"
#include "percyio.h"
#include "percyserver.h"
#include "agserver.h"
#include "itserver.h"
#include "recursiveserver.h"
#include "percystats.h"
#include <time.h>

#define READ_SIZE 1048576

PercyServer::PercyServer (DataStore * datastore, 
	const PercyServerParams * serverparams, PercyStats * stats)
:
    byzantine(false),
    datastore(datastore),
    serverparams(serverparams),
    stats(stats),
    strassen_max_depth(PercyServer::STRASSEN_OPTIMAL)
{
    if (serverparams->is_threaded() && !(serverparams->is_recursive())) {
	nservers_t num_threads = serverparams->num_threads();
	std::vector<const PercyServerParams*> wparams = 
		serverparams->get_all_thread_serverparams();
	std::vector<DataStore*> wdatastores = 
		datastore->get_worker_datastores();
	for (nservers_t i = 0; i < num_threads; ++i) {
	    PercyServer * wserver = make_server(wdatastores[i], wparams[i]);
	    subservers.push_back(wserver);
	}
    }
}

PercyServer::~PercyServer ()
{
    for (nservers_t i = 0; i < subservers.size(); ++i) {
	if (subservers[i]) delete subservers[i];
    }
}

PercyServer * PercyServer::make_server (DataStore * datastore,
	const PercyServerParams * serverparams, PercyStats * stats)
{
    PercyServer * retptr = NULL;
    switch (serverparams->percy_params()->get_mode()) {
    case MODE_ZZ_P:
	retptr = new PercyServer_ZZ_p(datastore, serverparams, stats);
	break;
    case MODE_CHOR:
	retptr = new PercyServer_Chor(datastore, serverparams, stats);
	break;
    case MODE_GF28:
	retptr = new PercyServer_GF2E<GF28_Element>(datastore, serverparams,
		stats);
	break;
    case MODE_GF216:
	retptr = new PercyServer_GF2E<GF216_Element>(datastore, serverparams,
		stats);
	break;
    case MODE_AG:
	retptr = new PercyAGServer(datastore, serverparams, stats);
	break;
    case MODE_RECURSIVE_AG:
    case MODE_HYBRID:
	retptr = new RecursiveServer(datastore, 
		dynamic_cast<const RecursiveServerParams*>(serverparams), 
		stats);
	break;
    default:
	break;
    }
    return retptr;
}

void PercyServer::be_byzantine ()
{
    byzantine = true;
}

void PercyServer::set_strassen_max_depth (nqueries_t depth)
{
    strassen_max_depth = depth;
}

nqueries_t PercyServer::get_strassen_max_depth ()
{
    return strassen_max_depth;
}

// Handle a request.
bool PercyServer::handle_request(istream &is, ostream &os, 
	std::vector<std::iostream*> workers)
{
    // Read number of queries
    nqueries_t num_queries;
    percy_read_le_uint16(is, num_queries);
    if (is.eof() || num_queries == 0) {
	return false;
    }

    // Not distributed
    nqueries_t batch_number = 0;
    if (stats) {
	batch_number = stats->start_query_batch(num_queries);
	stats->client_to_server_start(batch_number);
    }

    // Read request
    const PercyParams * params = serverparams->percy_params();
    dbsize_t query_request_size = params->request_size();
    dbsize_t request_size = num_queries * query_request_size;
    unsigned char * request = new unsigned char[request_size];
    // read in chunks to compensate for bug in socket++
    dbsize_t read_chunks = request_size / READ_SIZE;
    dbsize_t read_leftover = request_size % READ_SIZE;
    for (dbsize_t i = 0; i < read_chunks; ++i) {
	is.read((char*)request + i * READ_SIZE, READ_SIZE);
    }
    if (read_leftover > 0) {
	is.read((char*)request + read_chunks * READ_SIZE, read_leftover);
    }
    if (is.eof()) {
	std::cerr << "Request not large enough.\n";
	if (stats) {
	    stats->decode_done(batch_number, num_queries);
	    stats->finish_query_batch(batch_number);
	}
	delete[] request;
	return false;
    }

    if (stats) {
	stats->client_to_server_done(batch_number, request_size);
    }

    // Run the correct implementation functions
    dbsize_t query_response_size = params->response_size();
    dbsize_t response_size = num_queries * query_response_size;
    unsigned char * response = new unsigned char[response_size];
    memset(response, 0, response_size);

    std::vector<unsigned char*> requests, responses;
    for (nqueries_t q = 0; q < num_queries; ++q) {
	requests.push_back(request + q * query_request_size);
	responses.push_back(response + q * query_response_size);
    }

	clock_t start = clock();
    bool ret;
    if (serverparams->is_distributed()) {
	ret = handle_request_distributed(requests, responses, workers);
    } else if (serverparams->is_threaded()) {
	ret = handle_request_threaded(requests, responses);
    } else {
	ret = handle_request_impl(requests, responses);
    }
	clock_t end = clock();
	float totalTime = (float)(end - start) / (CLOCKS_PER_SEC / 1000);
	std::cerr << "Time to generate response: " << totalTime << " ms \n";

    if (stats) {
	stats->server_to_client_start(batch_number);
    }

    // Send response
    os.write((char*)response, response_size);
    os.flush();

    if (stats) {
	stats->server_to_client_done(batch_number, response_size);
	stats->decode_done(batch_number, ( ret ? 0 : num_queries ));
	stats->strassen_level_reached(batch_number,
		strassen_level_reached);
	stats->finish_query_batch(batch_number);
    }

    // Clean up
    delete[] response;
    delete[] request;

    return ret;
}


struct ThreadParams {
    ThreadParams (PercyServer * subserver, std::vector<unsigned char*> requests,
	    std::vector<unsigned char*> responses)
    :
	subserver(subserver),
	requests(requests),
	responses(responses),
	exitstatus(false)
    {}

    PercyServer * subserver;
    std::vector<unsigned char*> requests;
    std::vector<unsigned char*> responses;
    bool exitstatus;
};

void * thread_work (void * arg)
{
#ifdef VERBOSE_THREADED
    std::cerr << "Starting thread...\n";
#endif

    ThreadParams * tparams = (ThreadParams*)arg;
    PercyServer * subserver = tparams->subserver;
    std::vector<unsigned char*> requests = tparams->requests;
    std::vector<unsigned char*> responses = tparams->responses;

    tparams->exitstatus = subserver->handle_request_impl(requests, responses);

#ifdef VERBOSE_THREADED
    std::cerr << "Ending thread...\n";
#endif
    return NULL;
}

// Close fds 3 and up, except for the one given (pass -1 to close them
// all)
static void close_highfds_except(int exceptfd)
{
    // Find the max fd number
    struct rlimit limit;
    getrlimit(RLIMIT_NOFILE, &limit);
    for (int fd = 3; fd < (int)(limit.rlim_cur); ++fd) {
        if (fd != exceptfd) {
            // There's no ill effect from closing a non-open fd, so just
            // do it
            close(fd);
        }
    }
}

bool PercyServer::handle_request_threaded (std::vector<unsigned char*> requests,
	std::vector<unsigned char*> responses)
{
    nqueries_t num_queries = requests.size();
    if (responses.size() != requests.size()) {
	return false;
    }

    nservers_t num_threads = serverparams->num_threads();
    DistSplit split = serverparams->thread_split();
    std::vector<const PercyParams*> wparams = serverparams->get_all_thread_params();
    const PercyParams * params = serverparams->percy_params();

#ifdef VERBOSE_THREADED
    std::cerr << "NUM_THREADS: " << num_threads << "\n";
    timeval st, et, diff;
    gettimeofday(&st, NULL);
#endif

    // Get threads' num_queries
    std::vector<nqueries_t> thread_num_queries;
    if (split == DIST_SPLIT_QUERIES) {
	if (num_threads > num_queries) {
	    num_threads = num_queries;
#ifdef VERBOSE_THREADED
	    std::cerr << "More threads than queries.  Only using " << 
		    num_threads << " thread(s).\n";
#endif
	}
	thread_num_queries = split_queries(num_queries, num_threads);
    } else {
	thread_num_queries = std::vector<nqueries_t>(num_threads, num_queries);
    }

    // Get request and response buffers
    std::vector<std::vector<unsigned char*> > thread_requests (num_threads);
    std::vector<std::vector<unsigned char*> > thread_responses (num_threads);
    dbsize_t query_response_size = params->response_size();
    for (nqueries_t q = 0; q < num_queries; ++q) {
	memset(responses[q], 0, query_response_size);
    }
    unsigned char * intermediate = NULL;
    switch (split) {
    case DIST_SPLIT_QUERIES: {
	nqueries_t q = 0;
	for (nservers_t i = 0; i < num_threads; ++i) {
	    for (nqueries_t wq = 0; wq < thread_num_queries[i]; ++wq) {
		thread_requests[i].push_back(requests[q+wq]);
		thread_responses[i].push_back(responses[q+wq]);
	    }
	    q += thread_num_queries[i];
	}
    } break;
    case DIST_SPLIT_RECORDS: {
	intermediate = new unsigned char[num_queries * num_threads * query_response_size];
	memset(intermediate, 0, num_queries * num_threads * query_response_size);
	dbsize_t request_offset = 0;
	for (nservers_t i = 0; i < num_threads; ++i) {
	    dbsize_t thread_query_request_size = wparams[i]->request_size();
	    unsigned char * wresponse = intermediate + 
		    i * num_queries * query_response_size;
	    for (nqueries_t j = 0; j < num_queries; ++j) {
		thread_requests[i].push_back(requests[j] + request_offset);
		thread_responses[i].push_back(wresponse + j * query_response_size);
	    }
	    request_offset += thread_query_request_size;
	}
    } break;
    case DIST_SPLIT_RECORD_BYTES: {
	dbsize_t response_offset = 0;
	for (nservers_t i = 0; i < num_threads; ++i) {
	    dbsize_t thread_query_response_size = wparams[i]->response_size();
	    for (nqueries_t j = 0; j < num_queries; ++j) {
		thread_requests[i].push_back(requests[j]);
		thread_responses[i].push_back(responses[j] + response_offset);
	    }
	    response_offset += thread_query_response_size;
	}
    } break;
    }

    // Create ThreadParams objects
    std::vector<ThreadParams*> tparams;
    for (nservers_t i = 0; i < num_threads; ++i) {
	tparams.push_back(new ThreadParams(subservers[i], thread_requests[i],
		thread_responses[i]));
    }

#ifdef VERBOSE_THREADED
    gettimeofday(&et, NULL);
    sub_timevals(diff, et, st);
    std::cerr << "SET UP: " << diff << "\n";
    gettimeofday(&st, NULL);
#endif

    // Do the work and wrap up
    bool ret = true;
    if (serverparams->use_forked_threads()) {
	// Forked threads
	std::vector<int> childfds (num_threads, -1);
	pid_t * pid = new pid_t[num_threads];
	for (nservers_t i = 0; i < num_threads; ++i) {
	    int fds[2];
	    if (pipe(fds)) {
		std::cerr << "Pipe failed for thread " << i << "\n";
		pid[i] = -1;
		continue;
	    }
	    pid[i] = fork();
	    if (pid[i] < 0) {
		std::cerr << "Fork failed for thread " << i << "\n";
		continue;
	    } else if (pid[i] == 0) {
		// Child
		int parentfd = fds[1];
		close(fds[0]);
		// For compatibility with MPI
		close_highfds_except(parentfd);
		thread_work((void*)(tparams[i]));
		if (!(tparams[i]->exitstatus)) {
		    exit(1);
		}
		// Send responses to parent
		dbsize_t thread_query_response_size = wparams[i]->response_size();
		for (nqueries_t q = 0; q < thread_num_queries[i]; ++q) {
		    dbsize_t written = write(parentfd, 
			    (char*)(thread_responses[i][q]),
			    thread_query_response_size);
		    if (written != thread_query_response_size) {
			std::cerr << "Thread " << i << "failed when sending to parent\n";
			exit(1);
		    }
		}
		close(parentfd);
		exit(0);
	    } else {
		// Parent
		childfds[i] = fds[0];
		close(fds[1]);
	    }
	}

	// Get responses from children and do post-processing
	switch (split) {
	case DIST_SPLIT_QUERIES: {
	    for (nservers_t i = 0; i < num_threads; ++i) {
		for (nqueries_t q = 0; q < thread_num_queries[i]; ++q) {
		    dbsize_t read_bytes = read(childfds[i], 
			    (char*)(thread_responses[i][q]), 
			    query_response_size);
		    if (read_bytes != query_response_size) {
			std::cerr << "Did not get entire response from thread " 
				<< i << ".\n";
			ret = false;
			break;
		    }
		}
		if (!ret) break;
	    }
#ifdef VERBOSE_THREADED
    gettimeofday(&et, NULL);
    sub_timevals(diff, et, st);
    std::cerr << "THREADED COMPUTATION: " << diff << "\n";
    gettimeofday(&st, NULL);
#endif
	} break;
	case DIST_SPLIT_RECORDS: {
	    for (nservers_t i = 0; i < num_threads; ++i) {
		for (nqueries_t q = 0; q < num_queries; ++q) {
		    dbsize_t read_bytes = read(childfds[i], 
			    (char*)(thread_responses[i][q]),
			    query_response_size);
		    if (read_bytes != query_response_size) {
			std::cerr << "Did not get entire response from thread " 
				<< i << ".\n";
			ret = false;
			break;
		    }
		}
		if (!ret) break;
	    }
#ifdef VERBOSE_THREADED
    gettimeofday(&et, NULL);
    sub_timevals(diff, et, st);
    std::cerr << "THREADED COMPUTATION: " << diff << "\n";
    gettimeofday(&st, NULL);
#endif
	    if (ret) {
		for (nqueries_t q = 0; q < num_queries; ++q) {
		    std::vector<unsigned char*> query_responses;
		    for (nservers_t i = 0; i < num_threads; ++i) {
			query_responses.push_back(thread_responses[i][q]);
		    }
		    combine_results(responses[q], query_responses);
		}
	    }
	} break;
	case DIST_SPLIT_RECORD_BYTES: {
	    dbsize_t response_offset = 0;
	    for (nservers_t i = 0; i < num_threads; ++i) {
		dbsize_t thread_query_response_size = 
			wparams[i]->response_size();
		for (nqueries_t q = 0; q < num_queries; ++q) {
		    dbsize_t read_bytes = read(childfds[i], 
			    (char*)(thread_responses[i][q]),
			    thread_query_response_size);
		    if (read_bytes != thread_query_response_size) {
			std::cerr << "Did not get entire response from thread " 
				<< i << ".\n";
			ret = false;
			break;
		    }
		}
		if (!ret) break;
		response_offset += thread_query_response_size;
	    }
#ifdef VERBOSE_THREADED
    gettimeofday(&et, NULL);
    sub_timevals(diff, et, st);
    std::cerr << "THREADED COMPUTATION: " << diff << "\n";
    gettimeofday(&st, NULL);
#endif
	} break;
	}

	// Wait for children
	int status;
	for (dbsize_t i = 0; i < num_threads; ++i) {
	    close(childfds[i]);
	    waitpid(pid[i], &status, 0);
	    if (!WIFEXITED(status)) {
		std::cerr << "Thread " << i << " did not exit properly\n";
		ret = false;
		continue;
	    }
	    int exitstatus = WEXITSTATUS(status);
	    if (exitstatus) {
		std::cerr << "Thread " << i << " returned a bad exit status: " << exitstatus << "\n";
		ret = false;
	    }
	}

	delete[] pid;

    } else {
	// POSIX threads
	// Run threads
	pthread_t * threads = new pthread_t[num_threads];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int rc;
	for (nservers_t i = 0; i < num_threads; ++i) {
	    rc = pthread_create(&threads[i], &attr, thread_work, 
		    (void*)tparams[i]);
	    if (rc) {
		std::cerr << "Error creating thread " << i << "\n";
		ret = false;
	    }
	}
	pthread_attr_destroy(&attr);
	for (nservers_t i = 0; i < num_threads; ++i) {
	    void * status;
	    rc = pthread_join(threads[i], &status);
	    if (rc) {
		std::cerr << "Error joining thread " << i << "\n";
		ret = false;
	    }
	    if (!(tparams[i]->exitstatus)) {
		std::cerr << "Thread " << i << " did not finish correctly\n";
		ret = false;
	    }
	}
	delete[] threads;
#ifdef VERBOSE_THREADED
    gettimeofday(&et, NULL);
    sub_timevals(diff, et, st);
    std::cerr << "THREADED COMPUTATION: " << diff << "\n";
    gettimeofday(&st, NULL);
#endif

	// Do necessary post-processing
	if (ret) {
	    if (split == DIST_SPLIT_RECORDS) {
		for (nqueries_t q = 0; q < num_queries; ++q) {
		    std::vector<unsigned char*> query_responses;
		    for (nservers_t i = 0; i < num_threads; ++i) {
			query_responses.push_back(thread_responses[i][q]);
		    }
		    combine_results(responses[q], query_responses);
		}
	    }
	}
    }

    if (intermediate) delete[] intermediate;

#ifdef VERBOSE_THREADED
    gettimeofday(&et, NULL);
    sub_timevals(diff, et, st);
    std::cerr << "POST PROCESSING: " << diff << "\n";
#endif

    return ret;
}

// Do distributed computation
bool PercyServer::handle_request_distributed (std::vector<unsigned char*> requests, 
	std::vector<unsigned char*> responses, std::vector<std::iostream*> workers)
{
    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }
    
    // Get some parameters
    nservers_t num_workers = serverparams->num_workers();
    DistSplit split = serverparams->worker_split();
    std::vector<const PercyParams*> wparams = serverparams->get_all_worker_params();
    const PercyParams * params = serverparams->percy_params();

    // Sanity checks
    if (workers.size() < num_workers) {
	std::cerr << "Not enough workers given\n";
	return false;
    }

    // Add worker's number of queries to query_buffers
    std::vector<nqueries_t> worker_num_queries;
    if (split == DIST_SPLIT_QUERIES) {
	if (num_workers > num_queries) {
	    num_workers = num_queries;
#ifdef VERBOSE_DISTRIBUTED
	    std::cerr << "More workers than queries.  Only using " << 
		    num_workers << " worker(s).\n";
#endif
	}
	worker_num_queries = split_queries(num_queries, num_workers);
    } else {
	worker_num_queries = std::vector<nqueries_t>(num_workers, num_queries);
    }
    for (nservers_t i = 0; i < num_workers; ++i) {
	workers[i]->write((char*)&(worker_num_queries[i]), 2);
    }

    // Send the workers their queries
    dbsize_t query_request_size = params->request_size();
    switch (split) {
    case DIST_SPLIT_QUERIES: {
	nqueries_t q = 0;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    for (nqueries_t j = 0; j < worker_num_queries[i]; ++j) {
		workers[i]->write((char*)(requests[q+j]), query_request_size);
	    }
	    q += worker_num_queries[i];
	    workers[i]->flush();
	}
    } break;
    case DIST_SPLIT_RECORDS: {
	dbsize_t request_offset = 0;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    dbsize_t worker_request_size = wparams[i]->request_size();
	    for (nqueries_t j = 0; j < num_queries; ++j) {
		workers[i]->write((char*)(requests[j]) + request_offset, 
			worker_request_size);
	    }
	    workers[i]->flush();
	    request_offset += worker_request_size;
	}
    } break;
    case DIST_SPLIT_RECORD_BYTES: {
	for (nservers_t i = 0; i < num_workers; ++i) {
	    for (nqueries_t j = 0; j < num_workers; ++j) {
		workers[i]->write((char*)(requests[j]), query_request_size);
	    }
	    workers[i]->flush();
	}
    } break;
    }

    // Receive responses from workers
    dbsize_t query_response_size = params->response_size();
    for (nqueries_t q = 0; q < num_queries; ++q) {
	memset(responses[q], 0, query_response_size);
    }
    switch (split) {
    case DIST_SPLIT_QUERIES: {
	nqueries_t q = 0;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    for (nqueries_t j = 0; j < worker_num_queries[i]; ++j) {
		workers[i]->read((char*)(responses[q+j]), query_response_size);
		if (workers[i]->eof()) {
		    std::cerr << "Worker " << i << " did not send a full response\n";
		    return false;
		}
	    }
	    q += worker_num_queries[i];
	}
    } break;
    case DIST_SPLIT_RECORDS: {
	unsigned char * intermediate = new unsigned char[num_workers * query_response_size];
	std::vector<unsigned char*> query_responses;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    query_responses.push_back(intermediate + i * query_response_size);
	}
	for (nqueries_t q = 0; q < num_queries; ++q) {
	    for (nservers_t i = 0; i < num_workers; ++i) {
		workers[i]->read((char*)(query_responses[i]), query_response_size);
		if (workers[i]->eof()) {
		    std::cerr << "Worker " << i << " did not send a full response\n";
		    delete[] intermediate;
		    return false;
		}
	    }
	    combine_results(responses[q], query_responses);
	}
	delete[] intermediate;
    } break;
    case DIST_SPLIT_RECORD_BYTES: {
	dbsize_t response_offset = 0;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    dbsize_t worker_response_size = wparams[i]->response_size();
	    for (nqueries_t j = 0; j < num_queries; ++j) {
		workers[i]->read((char*)(responses[j]) + response_offset,
			worker_response_size);
		if (workers[i]->eof()) {
		    std::cerr << "Worker " << i << " did not send a full response\n";
		    return false;
		}
	    }
	    response_offset += worker_response_size;
	}
    } break;
    }

    return true;
}


std::vector<nqueries_t> split_queries (nqueries_t num_queries, 
	nservers_t num_workers)
{
    nqueries_t num_queries_div = num_queries / num_workers;
    nqueries_t num_queries_mod = num_queries % num_workers;
    std::vector<nqueries_t> split;
    for (nservers_t i = 0; i < num_queries_mod; ++i) {
	split.push_back(num_queries_div + 1);
    }
    for (nservers_t i = num_queries_mod; i < num_workers; ++i) {
	split.push_back(num_queries_div);
    }
    return split;
}

