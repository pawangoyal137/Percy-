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
#include "hybridparams.h"
#include "recursiveclient.h"
#include "config.h"
#include "percyio.h"

#define PERCY_MAX_CONNECT_ATTEMPTS 5
#define DEFAULT_AG_WORD_SIZE 20
#define DEFAULT_AG_SECURITY_N 50
#define DEFAULT_IT_MODE MODE_GF28
#define DEFAULT_IT_WORD_SIZE 32

void usage (const char * bin)
{
    std::cerr << "Usage: " << bin << " NUM_BLOCKS BLOCK_SIZE BLOCK_INDICES ELL T SERVERINFO [OPTIONS...]\n\n"
	      << "Required Arguments:\n"
	      << "   NUM_BLOCKS      The number of records in the database.\n"
	      << "   BLOCK_SIZE      The size of each database record in bytes.\n"
	      << "   BLOCK_INDICES   The indices of the records to fetch in the form \n"
	      << "                   \"idx1 ... idx\".\n"
	      << "   ELL             The number of servers.\n"
	      << "   T               The number of servers that can collude.\n"
	      << "   SERVERINFO      Server connection information in the form\n"
	      << "                   SID_0:ADDRESS_0:PORT_0 SID_1:ADDRESS_1:PORT_1 ... SID_ELL:ADDRESS_ELL:PORT_ELL\".\n\n"
	      << "Available Options:\n"
	      << "   -k, --k K               the number of servers that must response.\n"
	      << "                           (default: number of servers in SERVERINFO)\n"
	      << "   -m, --it-mode MODE      Specify the mode used for the IT-PIR step.  Supported modes are:\n"
	      << "                           Long form   Short form   Description\n"
	      << "                           GF28        g            use fast arithmetic in GF(2^8)\n"
	      << "                           GF216       s            use fast arithmetic in GF(2^16)\n"
	      << "                           ZZ_P        z            use arithmetic in Z mod p\n"
	      << "                           CHOR        c            use Chor et al.'s lightweight protocol\n"
	      << "                           (Default: " << DEFAULT_IT_MODE << ")\n"
	      << "   -Z, --it-word-size WS   Specify the word size for the IT-PIR step when the mode is ZZ_P\n"
	      << "                           (Default: " << DEFAULT_IT_WORD_SIZE << ")\n"
	      << "   -I, --it-num-blocks NB  The number of blocks that the database is split into for\n"
	      << "                           the IT-PIR step.  By default this value is optimized to\n"
	      << "                           minimize communication cost.\n"
	      << "   -t, --tau TAU           tau-independence value of database shares (default: 0).\n"
	      << "   -w, --ag-word-size WS   Specify the word size for the CPIR step (l0).\n"
	      << "                           WS must be one of {16, 20}.  (Default: "
		    << DEFAULT_AG_WORD_SIZE << ")\n"
	      << "   -N, --security-N NVAL   Specify the security parameter N.  (Default: " 
		    << DEFAULT_AG_SECURITY_N << ")\n"
	      << "   -d, --depth DEPTH       Specify the recursive depth of the scheme.\n"
	      << "                           By default, this value will be optimized to minimize\n"
	      << "                           communication cost.\n"
	      << "   -l, --log-file FILE     Enable statistics collection.  Output results to FILE\n"
	      << "                           (Default: stderr)\n"
	      << "   -L, --log-append STR    Append the string STR to the end of all log lines.\n"
	      << "   -H, --log-header        Print a header in the log file if there are no other entries.\n"
	      << "       --help              Display this help message and exit.\n"
	      << "       --version           Output version information and exit.\n";
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
    {"k",		required_argument,  NULL,   'k'},
    {"it-mode",		required_argument,  NULL,   'm'},
    {"it-word-size",	required_argument,  NULL,   'Z'},
    {"it-num-blocks",	required_argument,  NULL,   'I'},
    {"tau",		required_argument,  NULL,   't'},
    {"ag-word-size",	required_argument,  NULL,   'w'},
    {"security-N",	required_argument,  NULL,   'N'},
    {"depth",		required_argument,  NULL,   'd'},
    {"log-file",	required_argument,  NULL,   'l'},
    {"log-append",	required_argument,  NULL,   'L'},
    {"log-header",	no_argument,	    NULL,   'H'},
    {"help",		no_argument,        NULL,   'h'},
    {"version",		no_argument,	    NULL,   'v'},
    {NULL,		0,                  NULL,   0},
};

const char * short_options = "k:m:Z:I:t:w:N:d:l:L:H";

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

    // Servers needed
    nservers_t k = 0;
   
    // Optional arguments
    nqueries_t depth = 0;

    PercyMode it_mode = DEFAULT_IT_MODE;
    dbsize_t it_word_size = DEFAULT_IT_WORD_SIZE;
    dbsize_t it_num_blocks = 0;

    nservers_t tau = 0;
    // TODO: implement this stuff
//    char * pcparams_file = NULL;
//    bool do_spir = false;

    dbsize_t ag_word_size = DEFAULT_AG_WORD_SIZE;
    dbsize_t ag_N = DEFAULT_AG_SECURITY_N;

    bool do_logging = false;
    char * logfilename = NULL;
    char * logappend = NULL;
    bool print_header = false;

    int opt;
    while ((opt = getopt_long(argc, argv, short_options, long_options, 0)) != -1) {
	switch (opt) {
	case 'k':
	    k = strtoul(optarg, NULL, 10);
	    if (k == 0) {
		std::cerr << "k must be positive\n";
		exit(1);
	    }
	    break;
	case 'm':
	    if(!strcmp(optarg, "ZZ_P") || !strcmp(optarg, "z")) {
		it_mode = MODE_ZZ_P;
	    }
	    else if(!strcmp(optarg, "GF28") || !strcmp(optarg, "g")) {
		it_mode = MODE_GF28;
	    }
	    else if(!strcmp(optarg, "GF216") || !strcmp(optarg, "s")) {
		it_mode = MODE_GF216;
	    }
	    else if(!strcmp(optarg, "CHOR") || !strcmp(optarg, "c")) {
		it_mode = MODE_CHOR;
	    }
	    else {
		std::cerr << "Unknown mode selected. Valid modes are ZZ_P, GF28, GF216 and CHOR." << std::endl;
		return 1;
	    }
	    break;
	case 'Z':
	    it_word_size = strtoull(optarg, NULL, 10);
	    break;
	case 'I':
	    it_num_blocks = strtoull(optarg, NULL, 10);
	    break;
	case 't':
	    tau = strtoul(optarg, NULL, 10);
	    break;
	case 'w':
	    ag_word_size = strtoull(optarg, NULL, 10);
	    break;
	case 'N':
	    ag_N = strtoull(optarg, NULL, 10);
	    if (ag_N % 2 != 0) {
		std::cerr << "Error: N must be even\n";
		return 1;
	    }
	    break;
	case 'd':
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
	default:
	    std::cerr << "Invalid option: " << (char)optopt << "\n\n";
	    usage(argv[0]);
	    return 1;
	    break;
	}
    }

    // Required arguments
    if (argc - optind < 6) {
	usage(argv[0]);
	return 1;
    }
    dbsize_t num_blocks = strtoull(argv[optind++], NULL, 10);
    dbsize_t block_size = strtoull(argv[optind++], NULL, 10);
    char * indicesstr = argv[optind++];
    nservers_t ell = strtoul(argv[optind++], NULL, 10);
    nservers_t t = strtoul(argv[optind++], NULL, 10);
    char * sinfostr = argv[optind++];

    // Check argument values
    switch (ag_word_size) {
    case 16:
    case 20:
	break;
    default:
	std::cerr << "Error: Invalid AG word size: " << ag_word_size << "\n\n";
	usage(argv[0]);
	return 1;
    }

    // Set it_word_size (if applicable) and it_modulus
    ZZ it_modulus;
    switch (it_mode) {
    case MODE_GF28:
	it_word_size = 8;
	it_modulus = to_ZZ("256");
	break;
    case MODE_GF216:
	it_word_size = 16;
	it_modulus = to_ZZ("65536");
	break;
    case MODE_CHOR:
	it_word_size = 1;
	it_modulus = to_ZZ("256");
	break;
    case MODE_ZZ_P:
	if (it_word_size == 2048) {
	    it_modulus = to_ZZ("51162405833378812589599605953260132300166393994651819099454781579567509212081792013783783759303440508155949594262147212874957344953142209597742684263402581129339826752613431877280173074502314648334418584122460414512816448592261381117519846844295394134225624418756277265452922709245846828145574822031541004633366879073894273715489429502290966133193310966178373909137394353164436844312924586836474134940807305776164928781025210917912257206480517698118422827367766257579221703667784216949825206167241852365543481875593117676222875888924950402025039269210778276794873837063438751454865130720887819939394489366347567251243");
	} else if (it_word_size == 1536) {
	    it_modulus = to_ZZ("4065256781338999183533854850423382625119065920051798531476300569026463202897155088318466013703570859212040475097762405522038651420119366364979939687154236065682459920101982590074846996306687236388206057475890613264408059472973401701686869808348910896596468985609043697525749128687318350246421674945679872669881805678484464202726328189280359385791023305618545788872763420795247846720674554774196715770302797683129209164871258189464484019233379849839076263862630987");
	} else if (it_word_size == 1024) {
	    it_modulus = to_ZZ("343308946066366926839932845260501528909643718159825813630709694160026342456154871924497152436552679706642965502704642456637620829912957820221098686748075257358288200837461739492534713539606088624083011849535450485951774635526473457667739540374042376629835941950802202870595346459371144019363420985729553740241");
	} else if (it_word_size == 256) {
	    it_modulus = to_ZZ("115792089237316195423570985008687907853269984665640564039457584007913129640233");
	} else if (it_word_size == 196) {
	    it_modulus = to_ZZ("6277101735386680763835789423207666416102355444464034513029");
	} else if (it_word_size == 160) {
	    it_modulus = to_ZZ("2425980306017163398341728799446792216592523285797");
	} else if (it_word_size == 128) {
	    it_modulus = to_ZZ("340282366920938463463374607431768211507");
	} else if (it_word_size == 96) {
	    it_modulus = to_ZZ("79228162514264337593543950397");
	} else if (it_word_size == 32) {
	    it_modulus = to_ZZ("4294967311");
	} else if (it_word_size == 16) {
	    it_modulus = to_ZZ("65537");
	} else if (it_word_size == 8) {
	    it_modulus = to_ZZ("257");
	} else {
	    std::cerr << "Invalid word size for ZZ_P: " << it_word_size << "\n";
	    return 1;
	}
	break;
    default:
	return 1;
	break;
    }

    if (tau && it_mode == MODE_CHOR) {
	std::cerr << "Error: Chor et al.'s PIR scheme does not support tau independence." << std::endl;
	exit(1);
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
    vector<char*> sstr;
    sid_t * server_indices = new sid_t[ell];
    char * token = strtok(sinfostr, " ");
    while (token != NULL) {
	sstr.push_back(token);
	token = strtok(NULL, " ");
    }
    if (sstr.size() != ell) {
	std::cerr << "Error: in SERVERINFO, expected " << ell << " servers but found " 
		<< sstr.size() << ".\n";
	exit(1);
    }
    vector<serverinfo> sinfos;
    for (nservers_t i = 0; i < ell; ++i) {
	struct serverinfo sinfo;
	char * sidstr = strtok(sstr[i], ":");
	char * addrstr = strtok(NULL, ":");
	char * portstr = strtok(NULL, ":");
	if (portstr == NULL) {
	    std::cerr << "Error: Invalid server information: " << sstr[i] << "\n";
	    exit(1);
	}
	sinfo.addr = addrstr;
	sinfo.sid = strtoul(sidstr, NULL, 10);
	if (!sinfo.sid) {
	    std::cerr << "Error: SID must be positive\n";
	    exit(1);
	}
	server_indices[i] = (sid_t)(sinfo.sid);
	sinfo.port = strtoul(portstr, NULL, 10);
	if (sinfo.port < 1024 || sinfo.port > 65535) {
	    std::cerr << "Error: port number must be in range (1024, 65336)\n";
	    exit(1);
	}
	sinfos.push_back(sinfo);
    }

    // Connect to server
    vector<iosockinet*> sockets;
    vector<iostream*> streams;
    vector<ostream*> osvec;
    vector<istream*> isvec;
    for (nservers_t i = 0; i < ell; ++i) {
	iosockinet * socket = connect_to_server(sinfos[i]); // Temp change
	if (!socket) continue;
	sockets.push_back(socket);
	streams.push_back(socket);
	osvec.push_back(socket);
	isvec.push_back(socket);
    }

    if (k == 0) {
	k = ell;
    }

    if (streams.size() < k)
    {
	std::cerr << "Error: fewer than k=" << k << " out of ell=" << ell << " servers are online." << std::endl;
	for (nservers_t i = 0; i < streams.size(); ++i) {
	    if (streams[i]) delete streams[i];
	}
	delete[] server_indices;
	exit(1);
    }

    // Sanity checks for (ell,t,k,tau).
    if (t+tau >= k)
    {
	std::cerr << "Error: t+tau must be less than k." << std::endl;
	for (nservers_t i = 0; i < streams.size(); ++i) {
	    if (streams[i]) delete streams[i];
	}
	delete[] server_indices;
	exit(1);
    }
    if (it_mode == MODE_CHOR && t != k-1) {
        std::cerr << "Error: Chor requires that t=k-1." << std::endl;
	for (nservers_t i = 0; i < streams.size(); ++i) {
	    if (streams[i]) delete streams[i];
	}
	delete[] server_indices;
        exit(1);
    }

    // Create params
    RecursiveParams * params = NULL;
    if (it_mode == MODE_ZZ_P) {
	params = new HybridParams (num_blocks, block_size, it_word_size,
		it_modulus, depth, tau, it_num_blocks, ag_N, ag_word_size);
    } else {
	params = new HybridParams (num_blocks, block_size, it_mode,
		depth, tau, it_num_blocks, ag_N, ag_word_size);
    }
    RecursiveClientParams * clientparams = 
	    new RecursiveClientParams(params, streams.size());
    if (!params || !clientparams) {
	std::cerr << "Error creating parameters object\n";
	exit(1);
    }

#ifdef VERBOSE
    clientparams->print(std::cerr);
    std::cerr << "\n";
#endif

    // Verify that server's parameters match the client's
    for (nservers_t i = 0; i < ell; ++i) {
	// Send params
	clientparams->send(*osvec[i], sinfos[i].sid);
	osvec[i]->flush();

	// Get response
	unsigned char failure;
	isvec[i]->read((char*)&failure, 1);
	if (failure || isvec[i]->eof()) {
	    std::cerr << "Error: Server " << i << " (" << sstr[i] << ") did not accept parameters\n";
	    exit(1);
	}
    }

    // Logging
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
    RecursiveClient client (clientparams, ell, t, server_indices, stats);

    // Fetch the blocks
    nqueries_t req_id;
    vector<PercyBlockResults> results;
    bool res = client.fetch_blocks(req_id, indices, osvec, isvec, results);

    // Output results
    int ret = 0;
    if (!res) {
	std::cerr << "Error: The query failed\n";
	ret = 1;
    }
    nqueries_t num_res = results.size();
    for (nqueries_t r=0; ret == 0 && r < num_res; ++r) {
	if (results[r].results.empty()) {
	    std::cerr << "PIR query failed.\n";
	    ret = 1;
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
    if (stats) {
	if (logfilename != NULL) {
	    logstream.close();
	}
	delete stats;
    }
    delete clientparams;
    delete params;
    for (nservers_t i = 0; i < ell; ++i) {
	delete sockets[i];
    }
    delete[] server_indices;

    return ret;
}

