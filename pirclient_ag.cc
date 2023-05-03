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
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <socket++/sockinet.h>
#include <vector>
#include <sys/stat.h>
#include "percytypes.h"
#include "percyparams.h"
#include "agparams.h"
#include "agclient.h"
#include "config.h"
#include "percyio.h"

#define PERCY_MAX_CONNECT_ATTEMPTS 5
#define DEFAULT_WORD_SIZE 20
#define DEFAULT_SECURITY_N 50

void usage (const char * bin)
{
    std::cerr << "Usage: " << bin << " NUM_BLOCKS BLOCK_SIZE BLOCK_INDICES SERVERINFO [OPTIONS...]\n\n"
	      << "Required Arguments:\n"
	      << "   NUM_BLOCKS      The number of records in the database.\n"
	      << "   BLOCK_SIZE      The size of each database record in bytes.\n"
	      << "   BLOCK_INDICES   The indices of the records to fetch in the form \n"
	      << "                   \"idx1 ... idx\".\n"
	      << "   SERVERINFO      Server connection information in the form \"SID:ADDRESS:PORT\".\n\n"
	      << "Available Options:\n"
	      << "   -w, --word-size WSIZE       Specify the word size (l0).  WSIZE must be one of\n"
	      << "                               {16, 20}.  (Default: "
		    << DEFAULT_WORD_SIZE << ")\n"
	      << "   -N, --security-N NVAL       Specify the security parameter N.  (Default: " 
		    << DEFAULT_SECURITY_N << ")\n"
	      << "   -r, --recursive             Use the recursive version of the scheme."
	      << "   -d, --depth DEPTH           Specify the recursive depth of the scheme.  Implies --recursive.\n"
	      << "                               (Default: Chosen to optimize communication cost)\n"
	      << "   -l, --log-file FILE         Enable statistics collection.  Output results to FILE\n"
	      << "                               (Default: stderr)\n"
	      << "   -L, --log-append STR        Append the string STR to the end of all log lines.\n"
	      << "   -H, --log-header            Print a header in the log file if there are no other entries.\n"
	      << "       --help                  Display this help message and exit.\n"
	      << "       --version               Output version information and exit.\n";
}

// A method to connect to a server.  Returns a pointer to a socket if successful
// and a NULL pointer if unsuccessful.
// NOTE: The serverinfo class is defined in percyparams.h
iosockinet * connect_to_server (const serverinfo& sinfo) {
    std::cerr << "    Attempting to connect to " << sinfo.addr << ":" << sinfo.port << "...";
    bool connected = false;
    unsigned short attempts = 0;
    iosockinet *socket = new iosockinet(sockbuf::sock_stream);
    while (!connected && (attempts++ < PERCY_MAX_CONNECT_ATTEMPTS))
    {   
        try 
        {   
            (*socket)->connect(sinfo.addr, sinfo.port);
            connected = true;
        }   
        catch (sockerr e)
        {   
            cerr << ".";
            sleep(1);
        }   
    }   
    if (connected) {
        std::cerr << "succeeded!" << std::endl;
        return socket;
    } else {
        std::cerr << "failed!" << std::endl;
        return NULL;
    }
}

struct option long_options[] = {
    {"word-size",	required_argument,  NULL,   'w'},
    {"security-N",	required_argument,  NULL,   'N'},
    {"recursive",	no_argument,	    NULL,   'r'},
    {"depth",		required_argument,  NULL,   'd'},
    {"log-file",	required_argument,  NULL,   'l'},
    {"log-append",	required_argument,  NULL,   'L'},
    {"log-header",	no_argument,	    NULL,   'H'},
    {"help",		no_argument,	    NULL,   'h'},
    {"version",		no_argument,	    NULL,   'v'},
    {"null",		no_argument,	    NULL,   '0'},
    {NULL,		0,                  NULL,   0},
};

const char * short_options = ":w:N:rd:l:L:H0";

int main (int argc, char ** argv)
{
    // Initialize NTL and the random number stream
    ZZ modinit;
    modinit = to_ZZ(257);
    ZZ_p::init(modinit);
    unsigned char randbuf[128];
    ifstream urand("/dev/urandom");
    urand.read((char *)randbuf, sizeof(randbuf));
    urand.close();
    ZZ randzz = ZZFromBytes(randbuf, sizeof(randbuf));
    SetSeed(randzz);

    if (argc == 1) {
	usage(argv[0]);
	return -1;
    }

    // Optional arguments
    dbsize_t word_size = DEFAULT_WORD_SIZE;
    dbsize_t N = DEFAULT_SECURITY_N;

    // Recursive arguments
    bool do_recursive = false;
    nqueries_t depth = 0;

    bool do_logging = false;
    char * logfilename = NULL;
    char * logappend = NULL;
    bool print_header = false;

    bool is_null = false;

    int opt;
    while ((opt = getopt_long(argc, argv, short_options, long_options, 0)) != -1) {
	switch (opt) {
	case 'w':
	    word_size = strtoull(optarg, NULL, 10);
	    break;
	case 'N':
	    N = strtoull(optarg, NULL, 10);
	    if (N % 2 != 0) {
		std::cerr << "Error: N must be even\n";
		return 1;
	    }
	    break;
	case 'r':
	    do_recursive = true;
	    break;
	case 'd':
	    do_recursive = true;
	    depth = strtoul(optarg, NULL, 10);
	    break;
	case 'l':
	    do_logging = true;
	    logfilename = optarg;
	    break;
	case 'L':
	    do_logging = true;
	    logappend = optarg;
	    break;
	case 'H':
	    print_header = true;
	    break;
	case 'v':
	    std::cerr << "Percy++ pirclient_ag version " << PERCY_VERSION << "\n";
	    std::cerr << AUTHOR << "\n";
	    return 0;
	    break;
	case 'h':
	    usage(argv[0]);
	    return 0;
	    break;
	case '0':
	    is_null = true;
	    break;
	default:
	    std::cerr << "Invalid option: " << (char)optopt << "\n\n";
	    usage(argv[0]);
	    return 1;
	    break;
	}
    }

    // Required arguments
    if (argc - optind < 4) {
	usage(argv[0]);
	return 1;
    }
    dbsize_t num_blocks = strtoull(argv[optind++], NULL, 10);
    dbsize_t block_size = strtoull(argv[optind++], NULL, 10);
    char * indicesstr = argv[optind++];
    char * serverinfo = argv[optind++];

    // Check argument values
    switch (word_size) {
    case 16:
    case 20:
	break;
    default:
	std::cerr << "Error: Invalid word size: " << word_size << "\n\n";
	usage(argv[0]);
	return 1;
    }

    // Check indices
    std::vector<dbsize_t> indices;
    std::istringstream indexss(indicesstr);
    while (!indexss.eof()) {
	dbsize_t index;
	indexss >> index;
	if (index >= num_blocks) {
	    std::cerr << "Error: Record index is too high: " << index << "\n";
	    return 1;
	}
	indices.push_back(index);
    }

    // Parse server info
    struct serverinfo sinfo;
    char * sidstr = strtok(serverinfo, ":");
    char * addrstr = strtok(NULL, ":");
    char * portstr = strtok(NULL, ":");
    if (portstr == NULL) {
	std::cerr << "Error: Invalid server information: " << serverinfo << "\n";
	return 1;
    }
    sinfo.addr = addrstr;
    sinfo.sid = strtoul(sidstr, NULL, 10);
    sid_t sid = (sid_t)(sinfo.sid);
    if (!sid) {
	std::cerr << "Error: SID must be positive\n";
	return 1;
    }
    sinfo.port = strtoul(portstr, NULL, 10);
    if (sinfo.port < 1024 || sinfo.port > 65535) {
	std::cerr << "Error: port number must be in range (1024, 65336)\n";
	return 1;
    }

    // Create params
    PercyParams * params = NULL;
    PercyClientParams * clientparams = NULL;
    if (do_recursive) {
	RecursiveAGParams * rparams = new RecursiveAGParams(num_blocks, 
		block_size, depth, N, word_size);
	params = rparams;
	clientparams = new RecursiveClientParams(rparams, 1, is_null);
    } else {
	params = new AGParams(num_blocks, block_size, N, word_size);
	clientparams = new PercyClientParams(params, 1, is_null);
    }
    if (!params || !clientparams) {
	std::cerr << "Error: Error creating params object\n";
	if (clientparams) delete clientparams;
	if (params) delete params;
	return 1;
    }

#ifdef VERBOSE
    clientparams->print(std::cerr);
    std::cerr << "\n";
#endif

    // Connect to server
    iosockinet * socket = connect_to_server(sinfo);
    if (!socket) {
	std::cerr << "Error: Error connecting to server\n";
	delete clientparams;
	delete params;
	return 1;
    }
    std::iostream * sio = socket;
    std::istream& is = *sio;
    std::ostream& os = *sio;
    std::vector<std::ostream*> ostreams (1, sio);
    std::vector<std::istream*> istreams (1, sio);

    // Verify that server's parameters match the client's
    // Send params
    clientparams->send(os, sinfo.sid);
    os.flush();

    // Get response
    unsigned char failure;
    is.read((char*)&failure, 1);
    if (failure || is.eof()) {
	std::cerr << "Error: Server did not accept parameters\n";
	delete socket;
	delete clientparams;
	delete params;
	return 1;
    }

    // Loggind
    PercyStats * stats = NULL;
    ofstream logstream;
    if (do_logging) {
	if (logfilename != NULL) {
	    struct stat buffer;
	    print_header &= (stat (logfilename, &buffer) != 0);
	    logstream.open(logfilename, std::ios_base::app | std::ios_base::out);
	    stats = new PercyClientStats(clientparams, logstream, logappend);
	    if (print_header) {
		stats->print_header();
	    }
	} else {
	    stats = new PercyClientStats(clientparams, logappend);
	}
    }

    // Create the client
    PercyClient * client = PercyClient::make_client(clientparams, 1, 0,
	    &sid, stats);
    if (!client) {
	std::cerr << "Error: Error creating client\n";
	if (stats) {
	    if (logfilename != NULL) {
		logstream.close();
	    }
	    delete stats;
	}
	delete socket;
	delete clientparams;
	delete params;
	return 1;
    }

    // Fetch the blocks
    nqueries_t req_id;
    vector<PercyBlockResults> results;
    bool res = client->fetch_blocks(req_id, indices, ostreams, istreams, 
	    results);
    
    // Output results
    bool ret = true;
    if (!res) {
	std::cerr << "Error: The query failed\n";
	ret = false;
    }
    nqueries_t num_res = results.size();
    for (nqueries_t r=0; ret && r < num_res; ++r) {
	if (results[r].results.empty()) {
	    std::cerr << "PIR query failed.\n";
	    ret = false;
	    break;
	}
	else if (results[r].results.size() > 1) {
	    std::cerr << results[r].results.size() << " possible blocks returned.\n";
	}
	// Output the retrieved block(s)
	vector<PercyResult>::const_iterator resiter;
	for (resiter = results[r].results.begin(); resiter != results[r].results.end(); ++resiter) {
	    std::cout << resiter->sigma;
	}
    }

    // Clean up
    delete client;
    if (stats) {
	if (logfilename != NULL) {
	    logstream.close();
	}
	delete stats;
    }
    delete socket;
    delete clientparams;
    delete params;

    return (ret ? 0 : 1);
}

