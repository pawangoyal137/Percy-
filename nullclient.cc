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

#include <stdlib.h>
#include <stdio.h>
#include "percytypes.h"
#include "nullclient.h"
#include "percyio.h"

#define WRITESIZE 1048576

NullClient::NullClient (const PercyClientParams * params, 
	nservers_t num_servers)
:
    PercyClient(params, num_servers, 0, NULL),
    params(params->percy_params()),
    randbuf(NULL),
    num_to_process(0)
{
    long bufsize = 2 * WRITESIZE;
    randbuf = new unsigned char [bufsize];
    BytesFromZZ(randbuf, RandomLen_ZZ(8 * bufsize), bufsize);
}

NullClient::~NullClient ()
{
    delete[] randbuf;
}

void NullClient::encode_request_impl (nqueries_t request_identifier)
{}

dbsize_t NullClient::send_request_impl (nqueries_t request_identifier,
	std::vector<ostream*> &osvec, bool send_num_queries)
{
    nqueries_t num_queries = get_block_numbers(request_identifier).size();
    dbsize_t request_size = params->request_size(num_queries);
    dbsize_t write_chunks = request_size / WRITESIZE;
    dbsize_t write_leftover = request_size % WRITESIZE;

    for (nservers_t i = 0; i < num_servers; ++i) {
	if (send_num_queries) {
	    percy_write_le_uint16(*osvec[i], num_queries);
	}
	for (dbsize_t j = 0; j < write_chunks; ++j) {
	    long offset = RandomBnd(WRITESIZE);
	    osvec[i]->write((char*)randbuf + offset, WRITESIZE);
	}
	if (write_leftover > 0) {
	    long offset = RandomBnd(WRITESIZE);
	    osvec[i]->write((char*)randbuf + offset, write_leftover);
	}
	osvec[i]->flush();
    }

    return num_servers * request_size;
}

dbsize_t NullClient::receive_replies_impl (nqueries_t request_identifier,
	std::vector<istream*> &isvec)
{
    nqueries_t num_queries = get_block_numbers(request_identifier).size();
    dbsize_t response_size = params->response_size(num_queries);
    unsigned char * readbuf = new unsigned char[response_size];
    goodservers.clear();
    for (nqueries_t i = 0; i < num_servers; ++i) {
	isvec[i]->read((char*)readbuf, response_size);
	if (!(isvec[i]->eof())) {
	    goodservers.push_back(i);
	}
    }
    num_to_process += get_block_numbers(request_identifier).size();
    return num_servers * response_size;
}

nqueries_t NullClient::process_replies_impl (nservers_t h, 
	std::vector<std::vector<PercyResult> > &results)
{
    dbsize_t block_size = params->block_size();
    for (nqueries_t i = 0; i < num_to_process; ++i) {
	results.push_back(vector<PercyResult>(1, PercyResult(goodservers, 
		std::string(block_size, 0))));
    }
    num_to_process = 0;
    return 0;
}

