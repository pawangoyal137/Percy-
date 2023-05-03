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

#include <vector>
#include <fstream>
#include <cstring>
#include "recursiveserver.h"
#include "streams.h"
#include "percyio.h"

RecursiveServer::RecursiveServer (DataStore * datastore, 
	const RecursiveServerParams * serverparams, PercyStats * stats)
:
    PercyServer(datastore, serverparams, stats),
    serverparams(serverparams),
    params(static_cast<const RecursiveParams*>(serverparams->percy_params())),
    dist_first_only(serverparams->dist_first_only())
{
    // Make iteration servers
    nqueries_t depth = params->depth();
    for (nqueries_t d = 0; d < depth; ++d) {
	const PercyServerParams * serverparams_d = serverparams->get_iteration(d);
	DataStore * datastore_d;
	if (serverparams_d == NULL) {
	    datastore_d = NULL;
	} else if (d == 0) {
	    if (datastore == NULL) {
		datastore_d = NULL;
	    } else {
		datastore_d = new DataStore(const_cast<unsigned char*>(
			datastore->get_data()), serverparams_d);
	    }
	} else {
	    datastore_d = new DataStore(NULL, serverparams_d);
	}
	iteration_datastores.push_back(datastore_d);
	PercyServer * server_d = NULL;
	if (serverparams_d) {
	    server_d = PercyServer::make_server(datastore_d, serverparams_d);
	}
	iteration_servers.push_back(server_d);
    }

    // If distributed/threaded and split by queries we need subservers.
    if (serverparams->is_threaded() && serverparams->thread_split() == DIST_SPLIT_QUERIES) {
	std::vector<const PercyServerParams*> wsparams = 
		serverparams->get_all_thread_serverparams();
	std::vector<DataStore*> wdata = datastore->get_worker_datastores();
	for (nservers_t i = 0; i < serverparams->num_threads(); ++i) {
	    const RecursiveServerParams * childparams = 
		    dynamic_cast<const RecursiveServerParams*>(wsparams[i]);
	    subservers.push_back(new RecursiveServer(wdata[i], childparams));
	}
    }
}

RecursiveServer::~RecursiveServer ()
{
    nqueries_t depth = params->depth();
    for (nqueries_t d = 0; d < depth; ++d) {
	if (iteration_datastores[d]) {
	    delete iteration_datastores[d];
	}
	if (iteration_servers[d]) {
	    delete iteration_servers[d];
	}
    }
}

bool RecursiveServer::handle_request (std::istream &is, std::ostream &os,
	std::vector<std::iostream*> workers)
{
    if (serverparams->is_worker()) {
	return handle_request_worker(is, os);
    }
    return PercyServer::handle_request(is, os, workers);
}

bool RecursiveServer::handle_request_impl (std::vector<unsigned char*> requests, 
	std::vector<unsigned char*> responses)
{
    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }
    nqueries_t depth = params->depth();
    std::vector<const PercyParams*> iparams = params->get_iterations();

    // Get iteration request/response sizes
    std::vector<dbsize_t> iteration_request_sizes;
    std::vector<dbsize_t> iteration_response_sizes;
    for (nqueries_t d = 0; d < depth; ++d) {
	iteration_request_sizes.push_back(iparams[d]->request_size());
	iteration_response_sizes.push_back(iparams[d]->response_size());
    }

    // Allocate memory for intermediate datastores
    std::vector<unsigned char *> iteration_responses;
    for (nqueries_t d = 0; d < depth-1; ++d) {
	unsigned char * response_d = new unsigned char[iteration_response_sizes[d]];
	iteration_responses.push_back(response_d);
	iteration_datastores[d+1]->set_database(response_d);
    }

    // Call handle_request for each query and each iteration
    bool good = true;
    for (nqueries_t q = 0; q < num_queries; ++q) {
	unsigned char * request_location = requests[q];
	for (nqueries_t d = 0; d < depth; ++d) {
	    std::vector<unsigned char*> irequests, iresponses;
	    irequests.push_back(request_location);
	    if (d < depth-1) {
		iresponses.push_back(iteration_responses[d]);
		memset(iteration_responses[d], 0, iteration_response_sizes[d]);
	    } else {
		iresponses.push_back(responses[q]);
	    }
	    good = iteration_servers[d]->handle_request_impl(irequests, 
		    iresponses);
	    if (!good) break;
	    request_location += iteration_request_sizes[d];
	}
	if (!good) break;
    }

    // Clean up
    for (nqueries_t d = 0; d < depth-1; ++d) {
	if (iteration_responses[d]) {
	    delete iteration_responses[d];
	}
    }

    return good;
}

bool RecursiveServer::handle_request_threaded (std::vector<unsigned char*> requests, 
	std::vector<unsigned char*> responses)
{
    if (serverparams->thread_split() == DIST_SPLIT_QUERIES) {
	return PercyServer::handle_request_threaded(requests, responses);
    }

    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }
    nqueries_t depth = params->depth();
    std::vector<const PercyParams*> iparams = params->get_iterations();

    // Get iteration request/response sizes
    std::vector<dbsize_t> iteration_request_sizes;
    std::vector<dbsize_t> iteration_response_sizes;
    for (nqueries_t d = 0; d < depth; ++d) {
	iteration_request_sizes.push_back(iparams[d]->request_size());
	iteration_response_sizes.push_back(iparams[d]->response_size());
    }

    // Allocate memory for intermediate datastores
    std::vector<unsigned char *> iteration_responses;
    for (nqueries_t d = 0; d < depth-1; ++d) {
	unsigned char * response_d = new unsigned char[iteration_response_sizes[d]];
	iteration_responses.push_back(response_d);
	iteration_datastores[d+1]->set_database(response_d);
    }

    // Call handle_request for each query and each iteration
    bool good = true;
    for (nqueries_t q = 0; q < num_queries; ++q) {
	unsigned char * request_location = requests[q];
	for (nqueries_t d = 0; d < depth; ++d) {
	    std::vector<unsigned char*> irequests, iresponses;
	    irequests.push_back(request_location);
	    if (d < depth-1) {
		iresponses.push_back(iteration_responses[d]);
		memset(iteration_responses[d], 0, iteration_response_sizes[d]);
	    } else {
		iresponses.push_back(responses[q]);
	    }
	    good = iteration_servers[d]->handle_request_threaded(irequests, 
		    iresponses);
	    if (!good) break;
	    request_location += iteration_request_sizes[d];
	}
	if (!good) break;
    }

    // Clean up
    for (nqueries_t d = 0; d < depth-1; ++d) {
	if (iteration_responses[d]) {
	    delete iteration_responses[d];
	}
    }

    return good;
}

bool RecursiveServer::handle_request_distributed (std::vector<unsigned char*> requests, 
	std::vector<unsigned char*> responses, std::vector<std::iostream*> workers)
{
    if (serverparams->worker_split() == DIST_SPLIT_QUERIES) {
	return PercyServer::handle_request_distributed(requests, responses, workers);
    }

    nqueries_t num_queries = requests.size();
    if (responses.size() != num_queries) {
	return false;
    }
    nqueries_t depth = params->depth();
    std::vector<const PercyParams*> iparams = params->get_iterations();
    std::vector<const PercyServerParams*> isparams = serverparams->get_iterations();
    std::vector<const RecursiveParams*> wrparams = serverparams->get_worker_rparams();
    nservers_t num_workers = serverparams->num_workers();

    // Get iteration request/response sizes
    std::vector<dbsize_t> iteration_request_sizes;
    std::vector<dbsize_t> iteration_response_sizes;
    for (nqueries_t d = 0; d < depth; ++d) {
	iteration_request_sizes.push_back(iparams[d]->request_size());
	iteration_response_sizes.push_back(iparams[d]->response_size());
    }

    // Allocate memory for intermediate datastores
    std::vector<unsigned char *> iteration_responses;
    for (nqueries_t d = 0; d < depth-1; ++d) {
	unsigned char * response_d = new unsigned char[iteration_response_sizes[d]];
	iteration_responses.push_back(response_d);
	iteration_datastores[d+1]->set_database(response_d);
    }

    // Call handle_request for each query and each iteration
    bool good = true;
    for (nqueries_t q = 0; q < num_queries; ++q) {
	unsigned char * request_location = requests[q];
	for (nqueries_t d = 0; d < depth; ++d) {
	    if (d > 0 && !dist_first_only) {
		dbsize_t offset = 0;
		for (nservers_t i = 0; i < num_workers; ++i) {
		    const PercyParams * worker_iter_params = wrparams[i]->get_iteration(d);
		    if (!worker_iter_params) continue;
		    dbsize_t worker_iter_db_size = worker_iter_params->num_blocks()
			    * worker_iter_params->block_size();
		    workers[i]->write((char*)iteration_responses[d-1] + offset,
			    worker_iter_db_size);
		    offset += worker_iter_db_size;
		}
	    }
	    std::vector<unsigned char*> irequests, iresponses;
	    irequests.push_back(request_location);
	    if (d < depth-1) {
		iresponses.push_back(iteration_responses[d]);
		memset(iteration_responses[d], 0, iteration_response_sizes[d]);
	    } else {
		iresponses.push_back(responses[q]);
	    }
	    if (d > 0 && dist_first_only) {
		if (isparams[d]->is_threaded()) {
		    good = iteration_servers[d]->handle_request_threaded(irequests,
			    iresponses);
		} else {
		    good = iteration_servers[d]->handle_request_impl(irequests,
			    iresponses);
		}
	    } else {
		good = iteration_servers[d]->handle_request_distributed(irequests, 
			iresponses, workers);
	    }
	    if (!good) break;
	    request_location += iteration_request_sizes[d];
	}
	if (!good) break;
    }

    // Clean up
    for (nqueries_t d = 0; d < depth-1; ++d) {
	if (iteration_responses[d]) {
	    delete iteration_responses[d];
	}
    }

    return good;
}

bool RecursiveServer::handle_request_worker (std::istream& is, std::ostream& os)
{
    nqueries_t depth = params->depth();
    std::vector<const PercyParams*> iparams = params->get_iterations();

    // Get iteration database sizes
    std::vector<dbsize_t> iteration_db_sizes;
    for (nqueries_t d = 0; d < depth; ++d) {
	if (iparams[d] == NULL) {
	    iteration_db_sizes.push_back(0);
	} else {
	    iteration_db_sizes.push_back(iparams[d]-> num_blocks() * 
		    iparams[d]->block_size());
	}
    }

    // Allocate memory for intermediate datastores
    std::vector<unsigned char *> iteration_databases;
    iteration_databases.push_back(NULL);
    for (nqueries_t d = 1; d < depth; ++d) {
	if (iteration_db_sizes[d] == 0) {
	    iteration_databases.push_back(NULL);
	} else {
	    unsigned char * database_d = new unsigned char[iteration_db_sizes[d]];
	    iteration_databases.push_back(database_d);
	    iteration_datastores[d]->set_database(database_d);
	}
    }

    // Call handle_request for each query and each iteration
    bool good = true;
    while (good) {
	good = iteration_servers[0]->handle_request(is, os);
	if (!good) break;
	for (nqueries_t d = 1; d < depth; ++d) {
	    if (iteration_servers[d] == NULL) continue;
	    is.read((char*)iteration_databases[d], iteration_db_sizes[d]);
	    if (is.eof()) {
		std::cerr << "Error reading intermediate database at iteration " << d << "\n";
		good = false;
		break;
	    }
	    good = iteration_servers[d]->handle_request(is, os);
	    if (!good) break;
	}
	if (!good) break;
    }

    // Clean up
    for (nqueries_t d = 1; d < depth; ++d) {
	if (iteration_databases[d]) {
	    delete iteration_databases[d];
	}
    }

    return good;
}

/*
    DistSplit dsplit = serverparams->dist_split_method();
    nqueries_t depth = params->depth();
    
    if (dsplit == DIST_SPLIT_QUERIES) {
	return PercyServer::handle_request_distributed(request, response, 
		workers, num_queries);
    }

    // Get iteration request/response sizes
    std::vector<dbsize_t> iteration_request_sizes;
    std::vector<dbsize_t> iteration_response_sizes;
    for (nqueries_t d = 0; d < depth; ++d) {
	iteration_request_sizes.push_back(params->get_iteration(d)->request_size());
	iteration_response_sizes.push_back(params->get_iteration(d)->response_size());
    }

    // Allocate memory for intermediate datastores
    std::vector<unsigned char *> iteration_responses;
    for (nqueries_t d = 0; d < depth-1; ++d) {
	unsigned char * response_d = new unsigned char[iteration_response_sizes[d]];
	iteration_responses.push_back(response_d);
	iteration_datastores[d+1]->set_database(response_d);
    }

    // Call handle_request for each query and each iteration
    unsigned char * request_location = request;
    unsigned char * response_location = response;
    for (nqueries_t q = 0; q < num_queries; ++q) {
	for (nqueries_t d = 0; d < depth - 1; ++d) {
	    if (dist_first_only && d > 0) {
		iteration_servers[d]->handle_request_impl(request_location,
			iteration_responses[d]);
	    } else {
		iteration_servers[d]->handle_request_distributed(
			request_location, iteration_responses[d], workers);
	    }
	    request_location += iteration_request_sizes[d];
	}
	if (dist_first_only && depth > 1) {
	    iteration_servers.back()->handle_request_impl(
		    request_location, response_location);
	} else {
	    iteration_servers.back()->handle_request_distributed(
		    request_location, response_location, workers);
	}
	request_location += iteration_request_sizes.back();
	response_location += iteration_response_sizes.back();
    }

    // Clean up
    for (nqueries_t d = 0; d < depth-1; ++d) {
	if (iteration_responses[d]) {
	    delete iteration_responses[d];
	}
    }
    return true;
}
*/

void RecursiveServer::combine_results (unsigned char * result, 
	std::vector<unsigned char*> worker_results)
{}

