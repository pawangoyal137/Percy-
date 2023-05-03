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
#include <sys/wait.h>
#include <sys/types.h>
#include "percytypes.h"
#include "percyparams.h"
#include "agparams.h"
#include "agserver.h"
#include "config.h"
#include "percystats.h"
#include "percyio.h"

#define PERCY_DEFAULT_PORT 31337
#define PERCY_MAX_BIND_ATTEMPTS 10
#define PERCY_MAX_CONNECT_ATTEMPTS 5
#define DEFAULT_WORD_SIZE 20
#define DEFAULT_SECURITY_N 50

// A method to connect to a server.  Returns a pointer to a socket if successful
// and a NULL pointer if unsuccessful.
iosockinet * connect_to_server (const char * addr, const uint16_t port) {
    std::cerr << "    Attempting to connect to " << addr << ":" << port << "...";
    bool connected = false;
    unsigned short attempts = 0;
    iosockinet *socket = new iosockinet(sockbuf::sock_stream);
    while (!connected && (attempts++ < PERCY_MAX_CONNECT_ATTEMPTS))
    {   
        try 
        {   
            (*socket)->connect(addr, port);
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

void usage (const char * bin)
{
    std::cerr << "Usage: " << bin << " DATABSE SID NUM_BLOCKS BLOCK_SIZE [OPTIONS...]\n\n"
	      << "          Start the server on the specified database file DATABASE.\n\n"
	      << "       " << bin << " --master WORKERINFO SID NUM_BLOCKS BLOCK_SIZE [OPTIONS...]\n\n"
	      << "          Start the server as the master of a set of worker servers.\n\n"
	      << "       " << bin << " --worker NUM_WORKERS WORKER_INDEX DATABASE SID NUM_BLOCKS BLOCK_SIZE [OPTIONS...]\n\n"
	      << "          Start the PIR server as a worker number WORKER_INDEX of NUM_WORKERS.\n\n"
	      << "Required Arguments:\n"
	      << "   DATABASE                    The location of the database file.\n"
	      << "   SID                         Server identifier\n"
	      << "   NUM_BLOCKS                  The number of records in the database.\n"
	      << "   BLOCK_SIZE                  The size of each database record in bytes.\n"
	      << "   WORKERINFO                  The SID, address and port of each worker in the form\n"
	      << "                               \"sid1:addr1:port1 sid2:addr2:port2 ...\".\n\n"
	      << "   NUM_WORKERS                 Total number of workers.\n"
	      << "   WORKER_INDEX                The worker index for this worker.\n"
	      << "Available Options:\n"
	      << "   -p, --port PORTNO           Listen for connections on the specified port.\n"
	      << "   -w, --word-size WSIZE       Specify the word size (l0).  WSIZE must be one of\n"
	      << "                               {16, 20}.  (Default: "
		    << DEFAULT_WORD_SIZE << ")\n"
	      << "   -O, --offset OFFSET         The database starts OFFSET bytes from the beginning of the file.\n"
	      << "   -N, --security-N NVAL       Specify the security parameter N.  (Default: " 
		    << DEFAULT_SECURITY_N << ")\n"
	      << "   -r, --recursive             Use the recursive version of the scheme.\n"
	      << "   -d, --depth DEPTH           Specify the recursive depth of the scheme.  Implies --recursive.\n"
	      << "                               (Default: Chosen to optimize communication cost)\n"
	      << "   -T, --num-threads T         Distribute computation over T threads.\n"
	      << "   -F, --forked-threads        Use forked processes instead of threads.\n"
	      << "   -S, --split SPLIT           Specify how the queries are split up between threads.  Supported\n"
	      << "                               types are:\n"
	      << "                               Long form   Short form   Description" << std::endl
	      << "                               records     r            each thread is assigned a subset of database rows" << std::endl
	      << "                               bytes       b            each thread is assigned a section of the blocks" << std::endl
	      << "                               queries     q            each thread is assigned a subset of the queries" << std::endl
	      << "                               (default: records)" << std::endl
	      << "   -D, --dist-split SPLIT      Specify how the queries are split up between workers.  Has the same\n"
	      << "                               possible values as --split.  (default: records)\n"
	      << "   -f, --dist-first-only       Only do distributed computation on the first iteration\n"
	      << "   -l, --log-file FILE         Enable statistics collection.  Output results to FILE\n"
	      << "                               (Default: stderr)\n"
	      << "   -L, --log-append STR        Append the string STR to the end of all log lines.\n"
	      << "   -H, --log-header            Print a header in the log file if there are no other entries.\n"
	      << "   -1, --oneconn               Accept only a single connection.  Do not fork.\n"
	      << "       --help                  Display this help message and exit.\n"
	      << "       --version               Output version information and exit.\n\n";
}

bool bind_to_port (sockinetbuf& sin, uint16_t& port)
{
    sin.reuseaddr(true);
    if (!port) {
        port = PERCY_DEFAULT_PORT;
        retry:
        try {
            sin.bind((unsigned long) INADDR_ANY, port);
        } catch (sockerr) {
            //std::cerr << "Debug: failed to bind to port " << port << "."<< std::endl;
            port++;
            if (port < PERCY_DEFAULT_PORT + PERCY_MAX_BIND_ATTEMPTS) {
                goto retry;
            } else {
                std::cerr << "Error: unable to bind socket to a port. (I tried using ports " << PERCY_DEFAULT_PORT << " through " << (PERCY_DEFAULT_PORT + PERCY_MAX_BIND_ATTEMPTS) << ".)" << std::endl;
		return false;
            }
        }
    } else {
        try {
            sin.bind((unsigned long) INADDR_ANY, port);
        } catch (sockerr) {
            std::cerr << "Error: unable to bind socket on port " << port << "." << std::endl;
	    return false;
        }
    }
    sin.listen();
    std::cerr << "Listening on port " << port << "." << std::endl;
    return true;
}

void handle_requests (PercyServer * server, std::iostream& sio, 
	const PercyServerParams * params, std::vector<std::iostream*> workers)
{
    std::cerr << "Received an incoming connection\n";

    // Read params from client and check against ours.
    unsigned char failure = 0;
    failure = !(params->is_compatible(sio));
    sio.write((char*)&failure, 1);
    sio.flush();

    if (failure) {
	std::cerr << "Error: Client is not compatible with the server\n";
	return;
    }

    // Handle the request(s)
    while (server->handle_request(sio, sio, workers)) {
	std::cerr << "Finished a request\n";
    }
    std::cerr << "Query complete.  Terminating connection.\n";
}

struct option long_options[] = {
    {"master",		no_argument,	    NULL,   'M'},
    {"worker",		no_argument,	    NULL,   'W'},
    {"port",		required_argument,  NULL,   'p'},
    {"word-size",	required_argument,  NULL,   'w'},
    {"offset",		required_argument,  NULL,   'O'},
    {"security-N",	required_argument,  NULL,   'N'},
    {"recursive",	no_argument,	    NULL,   'r'},
    {"depth",		required_argument,  NULL,   'd'},
    {"log-file",	required_argument,  NULL,   'l'},
    {"log-append",	required_argument,  NULL,   'L'},
    {"log-header",	no_argument,	    NULL,   'H'},
    {"oneconn",		no_argument,	    NULL,   '1'},
    {"help",		no_argument,	    NULL,   'h'},
    {"version",		no_argument,	    NULL,   'v'},
    {"num-threads",	required_argument,  NULL,   'T'},
    {"forked-threads",	no_argument,	    NULL,   'F'},
    {"split",		required_argument,  NULL,   'S'},
    {"dist-split",	required_argument,  NULL,   'D'},
    {"dist-first-only",	no_argument,	    NULL,   'f'},
    {NULL,		0,		    NULL,   0}
};

const char * short_options = "MWp:w:O:N:rd:l:L:H1T:FS:D:f";

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
#ifdef NORAND
    SetSeed(ZZ::zero());
#else
    SetSeed(randzz);
#endif

    if (argc == 1) {
	usage(argv[0]);
	return -1;
    }

    // Optional arguments
    dbsize_t word_size = DEFAULT_WORD_SIZE;
    dbsize_t N = DEFAULT_SECURITY_N;
    dboffset_t offset = 0;

    // Recursive arguments
    bool do_recursive = false;
    nqueries_t depth = 0;
    
    uint16_t port = 0;
    bool daemon_mode = true;

    // Distributive computation parameters
    bool is_master = false;
    nservers_t num_threads = 0;
    DistSplit tsplit = DIST_SPLIT_RECORDS;
    nservers_t num_workers = 0;
    DistSplit wsplit = DIST_SPLIT_RECORDS;
    bool is_forked = false;
    bool is_worker = false;
    nservers_t worker_index = 0;
    bool first_only = false;
    bool print_header = false;

    // Logging
    bool do_logging = false;
    char * logfilename = NULL;
    char * logappend = NULL;

    // Threading
    //dbsize_t num_threads = 0;
    //PercyThreadingType ttype = THREADING_ROWS;
    //PercyThreadMethod tmethod = THREAD_METHOD_PTHREAD;

    int opt;
    while ((opt = getopt_long(argc, argv, short_options, long_options, 0)) != -1) {
	switch (opt) {
	case 'M':
	    is_master = true;
	    break;
	case 'W':
	    is_worker = true;
	    break;
	case 'p':
	    port = strtoul(optarg, NULL, 10);
	    break;
	case 'w':
	    word_size = strtoull(optarg, NULL, 10);
	    break;
	case 'O':
	    offset = strtoull(optarg, NULL, 10);
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
	case 'T':
	    num_threads = strtoll(optarg, NULL, 10);
	    if (num_threads < 1) {
		fprintf(stderr, "Must specify at least 1 thread!\n");
		return -1;
	    }
	    break;
	case 'F':
	    is_forked = true;
	    break;
	case 'S':
	    if (!strcmp(optarg, "records") || !strcmp(optarg, "r")) {
		tsplit = DIST_SPLIT_RECORDS;
	    } else if (!strcmp(optarg, "bytes") || !strcmp(optarg, "b")) {
		tsplit = DIST_SPLIT_RECORD_BYTES;
		std::cerr << "Splitting record bytes not yet supported\n";
		return -1;
	    } else if (!strcmp(optarg, "queries") || !strcmp(optarg, "q")) {
		tsplit = DIST_SPLIT_QUERIES;
	    } else {
		std::cerr << "Invalid thread splitting method selected.\n";
		usage(argv[0]);
		return 1;
	    }
	    break;
	case 'D':
	    if (!strcmp(optarg, "records") || !strcmp(optarg, "r")) {
		wsplit = DIST_SPLIT_RECORDS;
	    } else if (!strcmp(optarg, "bytes") || !strcmp(optarg, "b")) {
		wsplit = DIST_SPLIT_RECORD_BYTES;
		std::cerr << "Splitting  up record bytes not yet supported\n";
		return -1;
	    } else if (!strcmp(optarg, "queries") || !strcmp(optarg, "q")) {
		wsplit = DIST_SPLIT_QUERIES;
	    } else {
		std::cerr << "Invalid splitting method selected.\n";
		usage(argv[0]);
		return 1;
	    }
	    break;
	case 'f':
	    first_only = true;
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
	case '1':
	    daemon_mode = false;
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
	default:
	    std::cerr << "Invalid option: " << (char)optopt << "\n\n";
	    usage(argv[0]);
	    return 1;
	    break;
	}
    }

    // Cannot be a master and a worker
    if (is_master && is_worker) {
	fprintf(stderr, "Cannot be a master and a worker");
	usage(argv[0]);
	return 1;
    }

    // first_only only applies to wsplit == records and do_recursive == true
    if (wsplit != DIST_SPLIT_RECORDS || !do_recursive) {
	first_only = false;
    }

    // Required arguments
    if ((argc - optind < 4) || (is_worker && argc - optind < 6)) {
	std::cerr << "Error: Not enough arguments\n\n";
	usage(argv[0]);
	return 1;
    }
    if (is_worker) {
	num_workers = strtoul(argv[optind++], NULL, 10);
	worker_index = strtoul(argv[optind++], NULL, 10);
	// Check if these make sense
	if (worker_index >= num_workers) {
	    fprintf(stderr, "The worker index must be less than the number of workers\n");
	    return -1;
	}
    }
    char * dbfile = NULL;
    char * workerinfo = NULL;
    if (is_master) {
	workerinfo = argv[optind++];
    } else {
	dbfile = argv[optind++];
    }
    nservers_t sid = strtoul(argv[optind++], NULL, 10);
    dbsize_t num_blocks = strtoull(argv[optind++], NULL, 10);
    dbsize_t block_size = strtoull(argv[optind++], NULL, 10);

    // Check argument values
    if (!num_blocks) {
	std::cerr << "Error: The number of blocks must be positive\n";
	return 1;
    }
    if (!block_size) {
	std::cerr << "Error: The block size must be positive\n";
	return 1;
    }
    switch (word_size) {
    case 16:
    case 20:
	break;
    default:
	std::cerr << "Error: Invalid word size: " << word_size << "\n\n";
	usage(argv[0]);
	return 1;
    }

    std::vector<std::iostream*> workers;
    std::vector<nservers_t> worker_sids;
    if (is_master) {
	// Get worker addresses and ports
	std::vector<char*> workerstr;
	char * token = strtok(workerinfo, " ");
	while (token != NULL) {
	    workerstr.push_back(token);
	    token = strtok(NULL, " ");
	}
	num_workers = workerstr.size();
	if (num_workers < 1) {
	    std::cerr << "The number of workers must be at least 1\n";
	    return -1;
	}

	// Create worker iostreams
	bool good = true;
	for (dbsize_t i = 0; i < num_workers; ++i) {
	    char * sidchar = strtok(workerstr[i], ":");
	    char * addrchar = strtok(NULL, ":");
	    char * portchar = strtok(NULL, ":");
	    if (addrchar == NULL || portchar == NULL || portchar == NULL) {
		fprintf(stderr, "Error: worker information was incorrectly formatted.\n");
		good = false;
		break;
	    }
	    nservers_t sidnum = strtoul(sidchar, NULL, 10);
	    if (!sidnum) {
		std::cerr << "Error: SID must be an integer greater than 0.\n";
		good = false;
		break;
	    }
	    worker_sids.push_back(sidnum);
	    unsigned long portnum = strtoul(portchar, NULL, 10);
	    if (portnum < 1024 || portnum > 65535) {
		fprintf(stderr, "Error: port number must be an integer greater than 1024 and less than 65535.\n");
		good = false;
		break;
	    }
	    // connect_to_server function is located in distserver.{h,cc}
	    iosockinet * socket = connect_to_server(addrchar, portnum);
	    if (socket == NULL) {
		std::cerr << "Error: cannot connect to worker " << i << ": " << workerstr[i] << ".\n";
		good = false;
		break;
	    }
	    std::iostream * stream = socket;
	    workers.push_back(stream);
	}
	if (!good) {
	    for (nservers_t i = 0; i < workers.size(); ++i) {
		if (workers[i]) delete workers[i];
	    }
	    return -1;
	}

    } else {
	// Check database file
	struct stat filestatus;
	if (stat(dbfile, &filestatus)) {
	    std::cerr << "Error: Cannot file database file: " << dbfile << "\n";
	    for (nservers_t i = 0; i < workers.size(); ++i) {
		if (workers[i]) delete workers[i];
	    }
	    return 1;
	}
	dbsize_t dbsize = num_blocks * block_size;
	if (dbsize > (dbsize_t)filestatus.st_size) {
	    std::cerr << "Error: The database file is not large enough\n";
	    for (nservers_t i = 0; i < workers.size(); ++i) {
		if (workers[i]) delete workers[i];
	    }
	    return 1;
	}

	// Set all worker sids to my sid, so that we're using the right sid if
	// I'm a worker.
	worker_sids = std::vector<nservers_t>(num_workers, sid);
    }

    // Create params
    PercyParams * params = NULL;
    const PercyServerParams * serverparams = NULL;
    if (do_recursive) {
	RecursiveAGParams * rparams = new RecursiveAGParams(num_blocks,
		block_size, depth, N, word_size);
	params = rparams;
	serverparams = new RecursiveServerParams(rparams, sid, num_threads,
		tsplit, num_workers, wsplit, worker_sids, is_forked, first_only);
    } else {
	params = new AGParams(num_blocks, block_size, N, word_size);
	serverparams = new PercyServerParams(params, sid, num_threads, tsplit,
		num_workers, wsplit, worker_sids, is_forked);
    }
    if (!params || !serverparams) {
	std::cerr << "Error: Error creating params object\n";
	if (serverparams) delete serverparams;
	if (params) delete params;
	for (nservers_t i = 0; i < workers.size(); ++i) {
	    if (workers[i]) delete workers[i];
	}
	return 1;
    }

    // Create worker parameters and send params to workers if necessary
    if (is_master) {
	std::vector<const PercyServerParams*> wparams =
		serverparams->get_all_worker_serverparams();
	for (nservers_t i = 0; i < num_workers; ++i) {
#ifdef VERBOSE
	    std::cerr << "WORKER " << i << ": ";
	    wparams[i]->print(std::cerr);
	    std::cerr << "\n";
#endif
	    wparams[i]->send(*workers[i], true);
	}
	unsigned char failure;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    workers[i]->read((char*)&failure, 1);
	    if (failure) {
		std::cerr << "Worker " << i << " did not accept parameters\n";
		delete serverparams;
		delete params;
		for (nservers_t i = 0; i < workers.size(); ++i) {
		    if (workers[i]) delete workers[i];
		}
		return -1;
	    }
	}
    }

    // Create datastore
    DataStore * datastore = NULL;
    if (!is_master) {
	datastore = new FileDataStore(dbfile, serverparams, offset);
	if (!datastore) {
	    std::cerr << "Error creating datastore\n";
	    delete serverparams;
	    delete params;
	    for (nservers_t i = 0; i < workers.size(); ++i) {
		if (workers[i]) delete workers[i];
	    }
	    return -1;
	}
    }

    const PercyServerParams * parentparams = NULL;
    const DataStore * parentdatastore = NULL;
    if (is_worker) {
	parentparams = serverparams;
	serverparams = parentparams->get_worker_serverparams(worker_index);
	parentdatastore = datastore;
	datastore = parentdatastore->get_worker_datastore(worker_index);
	if (wsplit == DIST_SPLIT_RECORDS && tsplit == DIST_SPLIT_QUERIES) {
	    std::cerr << "Warning: Splitting threads by queries will have no effect.\n";
	}
    }

#ifdef VERBOSE
    serverparams->print(std::cerr);
    std::cerr << "\n";
#endif

    // Create socket
    sockinetbuf sin (sockbuf::sock_stream);
    if (!bind_to_port(sin, port)) {
	std::cerr << "Error: Did not successfully bid to a port\n";
	if (datastore) delete (parentdatastore ? parentdatastore : datastore);
	delete (parentparams ? parentparams : serverparams);
	delete params;
	for (nservers_t i = 0; i < workers.size(); ++i) {
	    if (workers[i]) delete workers[i];
	}
	return 1;
    }

    // Logging
    PercyStats * stats = NULL;
    ofstream logstream;
    if (do_logging) {
	if (logfilename) {
	    struct stat buffer;
	    print_header &= (stat (logfilename, &buffer) != 0);
	    logstream.open(logfilename, std::ios_base::app | std::ios_base::out);
	    stats = new PercyServerStats(serverparams, logstream, logappend);
	    if (print_header) {
		stats->print_header();
	    }
	} else {
	    stats = new PercyServerStats(serverparams, logappend);
	}
    }

    // Create the server
    PercyServer * server = PercyServer::make_server(datastore, serverparams, stats);
    if (!server) {
	std::cerr << "Error: Error creating server\n";
	if (stats) {
	    if (logfilename) {
		logstream.close();
	    }
	    delete stats;
	}
	if (datastore) delete (parentdatastore ? parentdatastore : datastore);
	delete (parentparams ? parentparams : serverparams);
	delete params;
	for (nservers_t i = 0; i < workers.size(); ++i) {
	    if (workers[i]) delete workers[i];
	}
	return 1;
    }

    // Start handling requests
    if (daemon_mode) {
	while (true) {
	    iosockinet sio (sin.accept());
	    pid_t childpid = fork();
	    if (childpid) {
		waitpid(childpid, NULL, 0);
	    } else {
		// spawn a grandchild and commit suicide so that the
		// parent doesn't have to wait()
		pid_t grandchildpid = fork();
		if (grandchildpid) {
		    break; // Will exit loop, clean up and exit
		} else {
		    // Handle request
		    handle_requests(server, sio, serverparams, workers);
		    break; // Will exit loop, clean up and exit
		}
	    }
	}
    } else {
	iosockinet sio (sin.accept());
	handle_requests(server, sio, serverparams, workers);
    }

    // Clean up
    delete server;
    if (stats) {
	if (logfilename) {
	    logstream.close();
	}
	delete stats;
    }
    if (datastore) delete (parentdatastore ? parentdatastore : datastore);
    delete (parentparams ? parentparams : serverparams);
    delete params;
    for (nservers_t i = 0; i < workers.size(); ++i) {
	if (workers[i]) delete workers[i];
    }

    return 0;
}

