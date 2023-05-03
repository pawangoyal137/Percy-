// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>,
// Ryan Henry <rhenry@cs.uwaterloo.ca>,
// Femi Olumofin <fgolumof@cs.uwaterloo.ca>
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

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "datastore.h"
#include "itserver.h"
#include "itparams.h"
#include "config.h"
#include <sys/types.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <socket++/sockinet.h>
#include <unistd.h>
#include <getopt.h>
#include "agserver.h"

#define PERCY_DEFAULT_PORT 31337
#define PERCY_MAX_BIND_ATTEMPTS 10
#define PERCY_MAX_CONNECT_ATTEMPTS 5

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

void print_usage (const char * bin, int status = 0) {
    std::cerr << "Usage: " << bin << " DATABASE SID NUM_BLOCKS BLOCK_SIZE [OPTIONS...]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "          Start the PIR server on the specified database file DATABASE." << std::endl;
    std::cerr << std::endl;
    std::cerr << "       " << bin << " --master WORKERINFO SID NUM_BLOCKS BLOCK_SIZE [OPTIONS...]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "          Start the PIR server as the master of a set of worker servers." << std::endl;
    std::cerr << std::endl;
    std::cerr << "       " << bin << " --worker NUM_WORKERS WORKER_INDEX DATABASE SID NUM_BLOCKS BLOCK_SIZE [OPTIONS...]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "          Start the PIR server as worker number WORKER_INDEX of NUM_WORKERS." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Required Arguments:\n";
    std::cerr << "   DATABASE               The location of the database file.\n";
    std::cerr << "   SID                    The server's identifier\n";
    std::cerr << "   NUM_BLOCKS             The number of records in the database\n";
    std::cerr << "   BLOCK_SIZE             The size (in bytes) of each record.\n";
    std::cerr << "   WORKERINFO             The SID, address and port of each worker in the form\n";
    std::cerr << "                          \"sid1:addr1:port1 sid2:addr2:port2 ...\".\n";
    std::cerr << "   NUM_WORKERS            Total number of workers.\n";
    std::cerr << "   WORKER_INDEX           The worker index for this worker.\n";
    std::cerr << std::endl;
    std::cerr << "Available Options:\n";
    std::cerr << "   NOTE: Mandatory arguments to long options are mandatory for short options too." << std::endl;
   // std::cerr << "   -n DBBYTES             use only the first DBBYTES bytes of database (default: entire file)." << std::endl;
    std::cerr << "   -w WORDSIZE            use a word size of WORDSIZE bits (default: 8)." << std::endl;
   // std::cerr << "   -b BLOCKSIZE           use a block size of BLOCKSIZE bytes (default: sqrt(DBBYTES*WORDSIZE)/8)." << std::endl;
    std::cerr << "   -O, --offset OFFSET    the database starts OFFSET bytes from the beginning of the file.\n";
    std::cerr << "   -t, --tau              specify that database is tau independent." << std::endl;
   // std::cerr << "   -S, --sid SERVERID     use the specified SID." << std::endl;
    std::cerr << "   -p, --port PORTNO      listen for connections on the specified port." << std::endl;
    std::cerr << "   -m, --mode MODE        use the specified mode of operation. Supported modes are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          GF28        g            use fast arithmetic in GF(2^8)" << std::endl;
    std::cerr << "                          GF216       s            use fast arithmetic in GF(2^16)" << std::endl;
    std::cerr << "                          ZZ_P        z            use arithmetic in Z mod p" << std::endl;
    std::cerr << "                          CHOR        c            use Chor et al.'s lightweight protocol" << std::endl;
    std::cerr << "                          (default: ZZ_P)" << std::endl;
#ifdef SPIR_SUPPORT
    std::cerr << "   -Z, --spir PCPARAMS    do symmetric PIR with specified PolyCommit" << std::endl;
#endif
    std::cerr << "                          parameters (a file)." << std::endl;
    std::cerr << "   -h, --hybrid           support hybrid security." << std::endl;
    std::cerr << "   -T, --num-threads T    distribute computation over T threads." << std::endl;
    std::cerr << "   -F, --forked-threads   use forked processes instead of threads." << std::endl;
    std::cerr << "   -S, --split SPLIT      specify how the queries are split up between threads.  Supported" << std::endl;
    std::cerr << "                          types are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          records     r            each thread is assigned a subset of database rows" << std::endl;
    std::cerr << "                          bytes       b            each thread is assigned a section of the blocks" << std::endl;
    std::cerr << "                          queries     q            each thread is assigned a subset of the queries" << std::endl;
    std::cerr << "                          (default: records)" << std::endl;
    std::cerr << "   -D, --dist-split SPLIT specify how the queries are split up between workers.  Has the same" << std::endl;
    std::cerr << "                          possible values as --split.  (default: records)" << std::endl;
    std::cerr << "   -l, --log-file FILE    enable statistics collection.  Output results to FILE" << std::endl;
    std::cerr << "                          (Default: stderr)" << std::endl;
    std::cerr << "   -L, --log-append STR   append the string STR to the end of all log lines." << std::endl;
    std::cerr << "   -H, --log-header       print a header in the log file if there are no other entries.\n";
    std::cerr << "   -z, --byzantine        be byzantine." << std::endl;
    std::cerr << "   -1, --oneconn          accept only a single conncetion; do not fork." << std::endl;
    std::cerr << "   -s, --strassen DEPTH   specify the maximum Strassen multiplication depth. Supported"<< std::endl;
    std::cerr << "                          depths are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          optimal     o            the optimal Strassen depth" << std::endl;
    std::cerr << "                          none        0            no Strassen is used" << std::endl;
    std::cerr << "                                      n            Strassen up to depth n (n > 0)" << std::endl;
    std::cerr << "                          (default: optimal)" << std::endl;
    std::cerr << "       --help             display this help and exit." << std::endl;
    std::cerr << "       --version          output version information and exit." << std::endl;
/* TODO
    std::cerr << std::endl;
    std::cerr << "Query Input/Response Output Options:" << std::endl;
    std::cerr << "   -F, --queries-from-file FILE   Use the file FILE as the data sent from the client." << std::endl;
    std::cerr << "   -G, --responses-to-file FILE   Redirect the responses to FILE." << std::endl;
*/
    std::cerr << std::endl;
    std::cerr << "Report bugs to iang+percy@cs.uwaterloo.ca." << std::endl;

    exit(status);
}

// Get the database size.  Returns 0 if the database does not exist.
dbsize_t database_bytes (const char * database) {
    struct stat filestatus;
    int not_exists = stat(database, &filestatus);
    if (not_exists)
    {
	std::cerr << "Error: cannot find database file " << database << std::endl;
	return 0;
    }
    return filestatus.st_size;
}


// Handle all of the requests on a single connection
static void handle_requests(PercyServer * server, std::istream &is, std::ostream &os,
    const PercyServerParams * serverparams, std::vector<std::iostream*> workers)
{
    std::cerr << "Received an incoming connection." << std::endl;
    // Receive the parameters from the client (and the SID, if necessary).
    unsigned char failure = 0;
    nservers_t sid = serverparams->get_sid();

    // First, read the client's query parameters.
//                std::cerr << "Receiving query parameters from client...";
    bool compatible = serverparams->is_compatible(is);
//                std::cerr << "done" << std::endl;
    // check if the two sets of params are compatible...
    if (!compatible) {
	failure = 1;
	std::cerr << "Client is not compatible with server\n";
    }
    os.write((char*)&failure, 1);
    os.flush();

    // Finally, do the PIR query!
    // With probability $PIRS_FAIL/100, fail completely
    // With probability $PIRS_BYZ/100, be Byzantine
    unsigned long rndval = RandomBnd(100);
    unsigned long failat = 0, byzat = 0, byznum = 0;
    const char *failenv = getenv("PIRS_FAIL");
    if (failenv) failat = atoi(failenv);
    const char *byzenv = getenv("PIRS_BYZ");
    if (byzenv) byzat = atoi(byzenv);
    const char *byznumenv = getenv("PIRS_BYZN");
    if (byznumenv) byznum = atoi(byznumenv);

    if (rndval < failat)
    {
	std::cerr << "["<<sid<<"] Failing.\n";
	return;
    }
    if (serverparams->is_byzantine() || rndval < failat + byzat || sid <= byznum)
    {
	std::cerr << "["<<sid<<"] Going Byzantine.\n";
	server->be_byzantine();
    }

    // Handle the request(s)
    // This loop will run until all queries are read and eof is
    // read.  We then gracefully exit the child process.
    while (server->handle_request(is, os, workers)) {
	std::cerr << "Finished a request\n";
    }
    std::cerr << "Query completed." << std::endl;
}

bool bind_to_port (sockinetbuf& sin, uint16_t& port)
{
    sin.reuseaddr(true);
    if (!port)
    {
        port = PERCY_DEFAULT_PORT;
        retry:
        try
        {
            sin.bind((unsigned long) INADDR_ANY, port);
        }
        catch (sockerr)
        {
            //std::cerr << "Debug: failed to bind to port " << port << "."<< std::endl;
            port++;
            if (port < PERCY_DEFAULT_PORT + PERCY_MAX_BIND_ATTEMPTS)
            {
                goto retry;
            }
            else
            {
                std::cerr << "Error: unable to bind socket to a port. (I tried using ports " << PERCY_DEFAULT_PORT << " through " << (PERCY_DEFAULT_PORT + PERCY_MAX_BIND_ATTEMPTS) << ".)" << std::endl;
		return false;
            }
        }
    }
    else
    {
        try
        {
            sin.bind((unsigned long) INADDR_ANY, port);
        }
        catch (sockerr)
        {
            std::cerr << "Error: unable to bind socket on port " << port << "." << std::endl;
	    return false;
        }
    }
    sin.listen();
    std::cerr << "Listening on port " << port << "." << std::endl;
    return true;
}

// List of long (i.e. --) options:
// {"longoptname", no_argument|required_argument|optional_argument, 0, 'shortoptname'}
struct option longopts[] = {
    {"version",		    no_argument,	NULL, 'v'},
    {"help",		    no_argument,	NULL, 'a'},
    {"master",		    no_argument,	NULL, 'M'},
    {"worker",		    no_argument,	NULL, 'W'},
    {"word-size",	    required_argument,	NULL, 'w'},
    {"offset",		    required_argument,	NULL, 'O'},
    {"tau",		    no_argument,	NULL, 't'},
    {"mode",		    required_argument,  NULL, 'm'},
    {"hybrid",		    no_argument,        NULL, 'h'},
    {"byzantine",	    no_argument,        NULL, 'z'},
    {"oneconn",		    no_argument,        NULL, '1'},
    {"strassen",	    required_argument,  NULL, 's'},
#ifdef SPIR_SUPPORT
    {"spir",		    required_argument,  NULL, 'Z'},
#endif
    {"port",		    required_argument,  NULL, 'p'},
    {"num-threads",	    required_argument,	NULL, 'T'},
    {"forked-threads",	    no_argument,	NULL, 'F'},
    {"split",		    required_argument,  NULL, 'S'},
    {"dist-split",	    required_argument,  NULL, 'D'},
    {"partial-database",    required_argument,  NULL, 'P'},
    {"log-file",	    required_argument,  NULL, 'l'},
    {"log-append",	    required_argument,  NULL, 'L'},
    {"log-header",	    no_argument,	NULL, 'H'},
    {NULL,		    0,                  NULL, 0},
};


int main (int argc, char ** argv)
{
    // Ignore SIGPIPE
    signal(SIGPIPE, SIG_IGN);

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
	print_usage(argv[0], -1);
    }

    // Allow hybrid queries?
    bool do_hybrid = false;
    
    // Mode of operation selected (ZZ_p, GF28, GF216 or Chor)
    PercyMode mode = MODE_ZZ_P;

    // Strassen maximum depth
    bool strassen_set = false;
    int strassen_max_depth = PercyServer::STRASSEN_OPTIMAL;

    // Do symmetric PIR?
    bool do_spir = false;

    // PolyCommit Params file to read for SPIR
    char *pcparams_file = NULL;
    
    // Is the database tau independent?
    bool is_tau = false;

    // Should we be byzantine?
    bool be_byzantine = false;

    uint16_t port = 0;
    
    // Word size
    dbsize_t w = 8;
    dboffset_t offset = 0;

    // Distributive computation parameters
    bool is_master = false;
    nservers_t num_threads = 0;
    DistSplit tsplit = DIST_SPLIT_RECORDS;
    nservers_t num_workers = 0;
    DistSplit wsplit = DIST_SPLIT_RECORDS;
    bool is_forked = false;
    bool is_worker = false;
    nservers_t worker_index = 0;
    bool is_partial_database = false;

    // Logging
    bool do_logging = false;
    char * logfilename = NULL;
    char * logappend = NULL;
    bool print_header = false;

    // Run once or as a daemon
    bool daemon_mode = true;

    // Parse arguments
#ifdef SPIR_SUPPORT
    const char * shortopts = "Mk:w:O:tm:hz1p:l:L:HZ:T:FS:D:Ps:";
#else
    const char * shortopts = "Mk:w:O:tm:hz1p:l:L:HT:FS:D:Ps:";
#endif
    int opt;
    while ((opt = getopt_long(argc, argv, shortopts, longopts, 0)) != -1) {
        switch(opt) {
	    case 'M':
		is_master = true;
		break;
	    case 'W':
		is_worker = true;
		break;
            case 'w':
                w = strtoull(optarg, NULL, 10);
                break;
	    case 'O':
		offset = strtoull(optarg, NULL, 10);
		break;
            case 't':
                is_tau = true;
                break;
            case 'm':
                if(!strcmp(optarg, "ZZ_P") || !strcmp(optarg, "z")) {
                    mode = MODE_ZZ_P;
                }
                else if(!strcmp(optarg, "GF28") || !strcmp(optarg, "g")) {
                    mode = MODE_GF28;
                }
                else if(!strcmp(optarg, "GF216") || !strcmp(optarg, "s")) {
                    mode = MODE_GF216;
                }
                else if(!strcmp(optarg, "CHOR") || !strcmp(optarg, "c")) {
                    mode = MODE_CHOR;
                }
                else {
                    std::cerr << "Unknown mode selected. Valid modes are ZZ_P, GF28, GF216 and CHOR.\n\n";
		    print_usage(argv[0], -1);
                }
                break;
            case 'h':
                do_hybrid = true;
                break;
            case 'z':
                be_byzantine = true;
                break;
#ifdef SPIR_SUPPORT
            case 'Z':
                do_spir = true;
                pcparams_file = optarg;
                break;
#endif
	    case 'p':
		port = strtoul(optarg, NULL, 10);
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
		    print_usage(argv[0], -1);
		}
		break;
	    case 'D':
		if (!strcmp(optarg, "records") || !strcmp(optarg, "r")) {
		    wsplit = DIST_SPLIT_RECORDS;
		} else if (!strcmp(optarg, "bytes") || !strcmp(optarg, "b")) {
		    wsplit = DIST_SPLIT_RECORD_BYTES;
		    std::cerr << "Splitting record bytes not yet supported\n";
		    return -1;
		} else if (!strcmp(optarg, "queries") || !strcmp(optarg, "q")) {
		    wsplit = DIST_SPLIT_QUERIES;
		} else {
		    std::cerr << "Invalid worker splitting method selected.\n";
		    print_usage(argv[0], -1);
		}
		break;
	    case 'P':
		is_partial_database = true;
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
	    case 'a':
		print_usage(argv[0]);
            case 'v':
                std::cerr << "Percy++ pirserver version " << PERCY_VERSION << std::endl;
                std::cerr << AUTHOR << std::endl;
		return 0;
                break;
	    case 's':
                if(!strcmp(optarg, "optimal") || !strcmp(optarg, "o")) {
		    strassen_set = true;
                    strassen_max_depth = PercyServer::STRASSEN_OPTIMAL;
                } else {
		    char * endptr;
		    int depth = (int) strtol(optarg, &endptr, 10);
		    if (endptr == optarg) {
			std::cerr << "Please specify a valid integer for Strassen depth";
			return -1;
		    } else {
			if( depth < -1 ) {
			    std::cerr << "Strassen depth should be non-negative";
			    return -1;
			} else {
			    strassen_set = true;
			    strassen_max_depth = depth;
			}
		    }
		}
		break;
            default:
		std::cerr << "Invalid option: " << opt << "\n";
		return -1;
        }
    }

    // Cannot be a master and a worker
    if (is_master && is_worker) {
	fprintf(stderr, "Cannot be a master and a worker");
	print_usage(argv[0], -1);
    }

    // Make sure enough mandatory arguments are present.
    if ((argc - optind < 4) || (is_worker && argc - optind < 6)) {
	fprintf(stderr, "Not enough arguments\n");
        print_usage(argv[0], -1);
    }    

    // Get required arguments
    if (is_worker) {
	num_workers = strtoul(argv[optind++], NULL, 10);
	worker_index = strtoul(argv[optind++], NULL, 10);
	// Check if these make sense
	if (worker_index >= num_workers) {
	    fprintf(stderr, "The worker index must be less than the number of workers\n");
	    return -1;
	}
    }

    char * database = NULL;
    char * workerinfo = NULL;
    if (is_master) {
	workerinfo = argv[optind++];
    } else {
	database = argv[optind++];
    }

    nservers_t sid = strtoul(argv[optind++], NULL, 10);
    dbsize_t num_blocks = strtoull(argv[optind++], NULL, 10);
    dbsize_t block_size = strtoull(argv[optind++], NULL, 10);

#if 0
    std::cerr << "sid = " << sid << "\n";
    std::cerr << "num_blocks = " << num_blocks << "\n";
    std::cerr << "block_size = " << block_size << "\n";
    std::cerr << "workerinfo = " << workerinfo << "\n";
    std::cerr << "database = " << database << "\n";
    std::cerr << "num_workers = " << num_workers << "\n";
    std::cerr << "worker_index = " << worker_index << "\n";
#endif

    dbbits_t n_bytes = num_blocks * block_size;
    dbsize_t b = block_size * 8;

    // Change threading method from pthread to fork if in ZZ_P
    if (num_threads > 0 and mode == MODE_ZZ_P and !is_forked) {
	fprintf(stderr, "The pthread library is not compatible with ZZ_p.  Forking separate processes instead.\n");
	is_forked = true;
    }

    // Sanity checks
    if (do_hybrid && (mode != MODE_ZZ_P)) {
        fprintf(stderr, "Error: hybrid security can only be used with the integers mod p mode of operation.\n");
	return -1;
    }
#ifdef SPIR_SUPPORT
    if (do_hybrid && do_spir) {
        fprintf(stderr, "Error: cannot use hybrid security with symmetric PIR.\n");
	return -1;
    }
    if (do_spir && mode != MODE_ZZ_P) {
        fprintf(stderr, "Error: symmetric PIR can only be used with the integers mod p mode of operation.\n");
	return -1;
    }
#endif
    if (is_tau && mode == MODE_CHOR) {
        fprintf(stderr, "Error: Chor et al.'s PIR scheme does not support tau independence.\n");
	return -1;
    }
    if (strassen_set && mode != MODE_ZZ_P && mode != MODE_GF28 &&
	    mode != MODE_GF216) {
        fprintf(stderr, "Strassen is only implemented for GF28, GF216 and ZZ_p.\n");
	return -1;
    }

    bool badw = false;
    switch (mode) {
    case MODE_ZZ_P:
	badw = (w % 8 != 0);
	break;
    case MODE_CHOR:
	badw = (w != 1);
	break;
    case MODE_GF28:
	badw = (w != 8);
	break;
    case MODE_GF216:
	badw = (w != 16);
	break;
    default:
	break;
    }
    if (badw) {
	std::cerr << "Invalid word size for mode " << mode << ": " << w << "\n";
	return 1;
    }

    // Choose an appropriate modulus.
    ZZ modulus = to_ZZ("256");
    if (mode == MODE_ZZ_P) {
	switch (w) {
	case 2048:
	    modulus = to_ZZ("51162405833378812589599605953260132300166393994651819099454781579567509212081792013783783759303440508155949594262147212874957344953142209597742684263402581129339826752613431877280173074502314648334418584122460414512816448592261381117519846844295394134225624418756277265452922709245846828145574822031541004633366879073894273715489429502290966133193310966178373909137394353164436844312924586836474134940807305776164928781025210917912257206480517698118422827367766257579221703667784216949825206167241852365543481875593117676222875888924950402025039269210778276794873837063438751454865130720887819939394489366347567251243");
	    break;
	case 1536:
	    modulus = to_ZZ("4065256781338999183533854850423382625119065920051798531476300569026463202897155088318466013703570859212040475097762405522038651420119366364979939687154236065682459920101982590074846996306687236388206057475890613264408059472973401701686869808348910896596468985609043697525749128687318350246421674945679872669881805678484464202726328189280359385791023305618545788872763420795247846720674554774196715770302797683129209164871258189464484019233379849839076263862630987");
	    break;
	case 1024:
	    modulus = to_ZZ("343308946066366926839932845260501528909643718159825813630709694160026342456154871924497152436552679706642965502704642456637620829912957820221098686748075257358288200837461739492534713539606088624083011849535450485951774635526473457667739540374042376629835941950802202870595346459371144019363420985729553740241");
	    break;
	case 256:
	    modulus = to_ZZ("115792089237316195423570985008687907853269984665640564039457584007913129640233");
	    badw = do_hybrid;
	    break;
	case 192:
	    modulus = to_ZZ("6277101735386680763835789423207666416102355444464034513029");
	    badw = do_hybrid;
	    break;
	case 160:
	    // NOTE: p2s is the prime from the PolyCommit params; spir
	    // will break if this value gets changed!
	    //
	    // TODO: read the prime from the PolyCommit params and check
	    //          that it is consistent with w.
	    modulus = to_ZZ("2425980306017163398341728799446792216592523285797");
	    badw = do_hybrid;
	    break;
	case 128:
	    modulus = to_ZZ("340282366920938463463374607431768211507");
	    badw = do_hybrid;
	    break;
	case 96:
	    modulus = to_ZZ("79228162514264337593543950397");
	    badw = do_hybrid;
	    break;
	case 32:
	    modulus = to_ZZ("4294967311");
	    badw = do_hybrid;
	    break;
	case 16:
	    modulus = to_ZZ("65537");
	    badw = do_hybrid;
	    break;
	case 8:
	    modulus = to_ZZ("257");
	    badw = do_hybrid;
	    break;
	default:
	    std::cerr << "Error: No modulus available for w = " << w << "." << std::endl;
	    return -1;
	}

	if (badw) {
	    std::cerr << "Error: No hybrid-compatible modulus available for w = " << w << "." << std::endl;
	    return -1;
	}

#ifdef SPIR_SUPPORT
	if (do_spir && w!=160)
	{
	    fprintf(stderr, "Error: symmetric PIR currently supports only w=160.\n");
	    return -1;
	}
#endif
    }

    std::vector<std::iostream*> workers;
    std::vector<nservers_t> worker_sids;
    dbsize_t dbsize = 0;
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
	    for (unsigned i = 0; i < workers.size(); ++i) 
		if (workers[i]) delete workers[i];
	    return -1;
	}

	// Create worker iostreams
	for (dbsize_t i = 0; i < num_workers; ++i) {
	    char * sidchar = strtok(workerstr[i], ":");
	    char * addrchar = strtok(NULL, ":");
	    char * portchar = strtok(NULL, ":");
	    if (addrchar == NULL || portchar == NULL || portchar == NULL) {
		fprintf(stderr, "Error: worker information was incorrectly formatted.\n");
		for (unsigned i = 0; i < workers.size(); ++i) 
		    if (workers[i]) delete workers[i];
		return -1;
	    }
	    nservers_t sidnum = strtoul(sidchar, NULL, 10);
	    if (!sidnum || sidnum > modulus) {
		std::cerr << "Error: SID must be an integer greater than 0 and less than " << modulus << ".\n";
		for (unsigned i = 0; i < workers.size(); ++i) 
		    if (workers[i]) delete workers[i];
		return -1;
	    }
	    worker_sids.push_back(sidnum);
	    unsigned long portnum = strtoul(portchar, NULL, 10);
	    if (portnum < 1024 || portnum > 65535) {
		fprintf(stderr, "Error: port number must be an integer greater than 1024 and less than 65535.\n");
		for (unsigned i = 0; i < workers.size(); ++i) 
		    if (workers[i]) delete workers[i];
		return -1;
	    }
	    // connect_to_server function is located in distserver.{h,cc}
	    iosockinet * socket = connect_to_server(addrchar, portnum);
	    if (socket == NULL) {
		std::cerr << "Error: cannot connect to worker " << i << ": " << workerstr[i] << ".\n";
		for (unsigned i = 0; i < workers.size(); ++i) 
		    if (workers[i]) delete workers[i];
		return -1;
	    }
	    std::iostream * stream = socket;
	    workers.push_back(stream);
	}

	// database size must be specified
	if (!n_bytes) {
	    fprintf(stderr, "Error: Database size (n) must be specified.\n");
	    for (unsigned i = 0; i < workers.size(); ++i) 
		if (workers[i]) delete workers[i];
	    return -1;
	}

    } else {
	// Make sure the specified database file exists.
	dbsize = database_bytes(database);
	if (dbsize == 0) {
	    fprintf(stderr, "Error: the database must exist and be non-empty.\n");
	    return -1;
	}

	// If no value for "n" is specified, then use a default database
	// size of dbsize. Otherwise, just check that 0<n<=dbsize.
	if (!n_bytes) {
	    n_bytes = dbsize;
	} else if (!(is_worker && is_partial_database) && n_bytes > dbsize) {
	    fprintf(stderr, "Error: n cannot be larger than database file.\n");
	    return -1;
	}

	// Set all worker sids to my sid, so that we're using the right sid if
	// I'm a worker.
	worker_sids = std::vector<nservers_t>(num_workers, sid);
    }

    dbbits_t n = n_bytes * 8;
    if (n_bytes > n) {
        fprintf(stderr, "Error: database file is too large for the current architecture!\n");
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	return -1;
    }

    // If no value for "b" is specified, then use a default block size
    // of \sqrt(n * w) bits.
    if (!b) {
        b = sqrt(n * w);
        if (n != b*b/w)
        {
            fprintf(stderr, "Error: optimal parameter choice is invalid for this database. Please specify a value for both of b and w.\n");
	    for (unsigned i = 0; i < workers.size(); ++i) 
		if (workers[i]) delete workers[i];
	    return -1;
        }
    }

    // Sanity checks for (n,b,w).
    if (n % b != 0 || b % w != 0)
    {
        fprintf(stderr, "Error: b must divide n and w must divide b.\n");
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	return -1;
    }
    
    // Compute the number of blocks, and number of words per block.
    dbsize_t words_per_block = b / w;
    //std::cerr << "Number of blocks: " << num_blocks << std::endl;
    //std::cerr << "Words per block:  " << words_per_block << std::endl;
    //std::cerr << "Bits per block:  " << b << std::endl;
    //std::cerr << "Bits per word:  " << w << std::endl;
    if (num_blocks != words_per_block)
    {
        std::cerr << "Warning: non-optimal choice of blocksize detected." << std::endl;
    }

    // Create the PercyServerParams object.
    PercyParams * params = NULL;
    switch (mode) {
    case MODE_ZZ_P:
	params = new ZZ_pParams(num_blocks, block_size, w, modulus, 
		is_tau, pcparams_file, do_spir);
	break;
    case MODE_CHOR:
	params = new ChorParams(num_blocks, block_size);
	break;
    case MODE_GF28:
    case MODE_GF216:
	params = new GF2EParams(num_blocks, block_size, w, is_tau);
	break;
    default:
	std::cerr << "Invalid mode: " << mode << "\n";
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	return -1;
    }
    if (!params) {
	std::cerr << "Error creating parameters object\n";
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	return -1;
    }
    const PercyServerParams * serverparams = new PercyServerParams(params, sid, 
	    num_threads, tsplit, num_workers, wsplit, worker_sids, is_forked, 
	    be_byzantine);
    if (!serverparams) {
	delete params;
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	return -1;
    }

    // Send params to workers if necessary
    if (is_master) {
	std::vector<const PercyServerParams*> wparams =
		serverparams->get_all_worker_serverparams();
	for (nservers_t i = 0; i < num_workers; ++i) {
	    wparams[i]->send(*workers[i], true);
	}
	unsigned char failure;
	for (nservers_t i = 0; i < num_workers; ++i) {
	    workers[i]->read((char*)&failure, 1);
	    if (failure) {
		std::cerr << "Worker " << i << " did not accept parameters\n";
		delete serverparams;
		delete params;
		for (unsigned i = 0; i < workers.size(); ++i) 
		    if (workers[i]) delete workers[i];
		return -1;
	    }
	}
    }

    // Create datastore
    DataStore * datastore = NULL;
    if (!is_master) {
	if (is_worker && is_partial_database) {
	    const PercyParams * wparams = serverparams->get_worker_params(worker_index);
	    const PercyServerParams * wsparams = serverparams->get_worker_serverparams(worker_index);
	    if (wparams->num_blocks() * wparams->block_size() > dbsize) {
		std::cerr << "Error: The partial database is not large enough!\n";
		delete serverparams;
		delete params;
		for (unsigned i = 0; i < workers.size(); ++i) 
		    if (workers[i]) delete workers[i];
		return -1;
	    }
	    datastore = new FileDataStore(database, wsparams, offset);
	} else {
	    datastore = new FileDataStore(database, serverparams, offset);
	}
	if (!datastore) {
	    std::cerr << "Error creating datastore\n";
	    delete serverparams;
	    delete params;
	    for (unsigned i = 0; i < workers.size(); ++i) 
		if (workers[i]) delete workers[i];
	    return -1;
	}
    }

    const PercyServerParams * parentparams = NULL;
    const DataStore * parentdatastore = NULL;
    if (is_worker) {
	parentparams = serverparams;
	serverparams = parentparams->get_worker_serverparams(worker_index);
	if (!is_partial_database) {
	    parentdatastore = datastore;
	    datastore = parentdatastore->get_worker_datastore(worker_index);
	}
    }

#ifdef VERBOSE
    serverparams->print(std::cerr);
    std::cerr << "\n";
#endif

    // Create a socket for clients to connect to.
    sockinetbuf sin(sockbuf::sock_stream);
    if (!bind_to_port(sin, port)) {
	if (datastore) delete (parentdatastore ? parentdatastore : datastore);
	delete (parentparams ? parentparams : serverparams);
	delete params;
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	fprintf(stderr, "Did not successfully bind to a port\n");
	return -1;
    }

    // Logging
    PercyStats * stats = NULL;
    ofstream logstream;
    if (do_logging) {
	if (logfilename != NULL) {
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
	if (!stats) {
	    if (datastore) delete (parentdatastore ? parentdatastore : datastore);
	    delete (parentparams ? parentparams : serverparams);
	    delete params;
	    for (unsigned i = 0; i < workers.size(); ++i) 
		if (workers[i]) delete workers[i];
	    fprintf(stderr, "Error creating stats object\n");
	    return -1;
	}
    }

    // Create the server
    PercyServer * server = PercyServer::make_server(datastore, serverparams, stats);
    if (server == NULL) {
	if (stats) delete stats;
	if (datastore) delete (parentdatastore ? parentdatastore : datastore);
	delete (parentparams ? parentparams : serverparams);
	delete params;
	for (unsigned i = 0; i < workers.size(); ++i) 
	    if (workers[i]) delete workers[i];
	fprintf(stderr, "Server not created successfully\n");
	return -1;
    }

    server->set_strassen_max_depth(strassen_max_depth);

    if (daemon_mode) {
	// Daemon mode
	while(true) {
	    iosockinet sio(sin.accept());

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
		    handle_requests(server, sio, sio, serverparams, workers);
		    break; // Will exit loop, clean up and exit
		}
	    }
	}
    } else {
	// One connection
	// Get incoming socket connection.
	iosockinet sio(sin.accept());

	// Handle request
	handle_requests(server, sio, sio, serverparams, workers);
    }

    // Clean up
    delete server;
    if (stats) delete stats;
    if (datastore) delete (parentdatastore ? parentdatastore : datastore);
    delete (parentparams ? parentparams : serverparams);
    delete params;
    for (unsigned i = 0; i < workers.size(); ++i) 
	if (workers[i]) delete workers[i];

    return 0;
}

