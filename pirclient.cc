// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Paul Hendry <pshendry@uwaterloo.ca>,
// Ann Yang <y242yang@uwaterloo.ca>,
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

#include <signal.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ZZ_p.h>
#include <socket++/sockinet.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <getopt.h>
#include <iterator>
#include <sys/stat.h>
#include "itclient.h"
#include "config.h"

#define PERCY_MAX_CONNECT_ATTEMPTS 5

// Let E = w * (m(h-t-1)-k+h+1)
// where w is the bitlength of a word and m is the number of resent queries.
// Then this value is highest that E can be.  (NOTE: this is chosen so that
// with this many resent queries, TK1 fails approx 1 time in 1,000,000,000,000)
#define QUERY_LIMIT_EXP 40

void print_usage_exit(char *bin, int status = -1)
{
//    std::cerr << "Usage: " << bin << " [OPTIONS...] r s w ell k t \"idx1 idx2 ... idxq\" \"sid1:addr1:port1 ... sidell:addrell:portell\"" << std::endl;
    std::cerr << "Usage: " << bin << " NUM_BLOCKS BLOCK_SIZE SERVERINFO PRIVACY_LEVEL BLOCK_INDICES [OPTIONS...]\n";
    std::cerr << std::endl;
    std::cerr << "          Query the specified PIR servers at the specified indices.\n";
    std::cerr << std::endl;
    std::cerr << "Required Arguments:\n";
    std::cerr << "   NUM_BLOCKS           The number of records in the database\n";
    std::cerr << "   BLOCK_SIZE           The size (in bytes) of each record.\n";
    std::cerr << "   SERVERINFO           The SID, address and port of each server in the form\n";
    std::cerr << "                        \"sid1:addr1:port1 sid2:addr2:port2 ...\".\n";
    std::cerr << "   PRIVACY_LEVEL        The maximum number of servers that can collude\n";
    std::cerr << "   BLOCK_INDICES        The indices of the records to fetch (0-based)\n";
/*
    std::cerr << "   r     number of block s." << std::endl;
    std::cerr << "   s     words per block." << std::endl;
    std::cerr << "   w     word size (in bits)." << std::endl;
    std::cerr << "   ell   number of servers." << std::endl;
    std::cerr << "   k     number of servers that need to respond." << std::endl;
    std::cerr << "   t     number of servers that can collude." << std::endl;
    std::cerr << "   idxi  indices of blocks to fetch (0-based)." << std::endl;
    std::cerr << "   sidi:addri:porti  the SID, address and port number for each of the ell servers." << std::endl;
*/
    std::cerr << std::endl;
    std::cerr << "Available Options:\n";
    std::cerr << "   NOTE: Mandatory or optional arguments to long options are also mandatory or optional" << std::endl;
    std::cerr << "         for any corresponding short options too." << std::endl;
    std::cerr << "   -k, --k K              the number of servers that must respond.\n";
    std::cerr << "                          (default: number of servers in SERVERINFO)\n";
    std::cerr << "   -w, --word-size WS     use a word size of WS bits (default: 8)." << std::endl;
    std::cerr << "   -t, --tau TAU          tau-independence value of database shares (default: 0)." << std::endl;
    std::cerr << "   -m, --mode MODE        use the specified mode of operation. Supported modes are:" << std::endl;
    std::cerr << "                          Long form   Short form   Description" << std::endl;
    std::cerr << "                          GF28        g            use fast arithmetic in GF(2^8)" << std::endl;
    std::cerr << "                          GF216       s            use fast arithmetic in GF(2^16)" << std::endl;
    std::cerr << "                          ZZ_P        z            use arithmetic in Z mod p" << std::endl;
    std::cerr << "                          CHOR        c            use Chor et al.'s lightweight protocol" << std::endl;
    std::cerr << "                          (default: ZZ_P)" << std::endl;
    std::cerr << "   -Q, --batch-query QBS  request QBS blocks in each query." << std::endl;
#ifdef SPIR_SUPPORT
    std::cerr << "   -Z, --spir PCPARAMS    do symmetric PIR with specified PolyCommit" << std::endl;
    std::cerr << "                          parameters (a file)." << std::endl;
#endif
    std::cerr << "   -h, --hybrid [KEY]     enable hybrid security; optionally specify a" << std::endl;
    std::cerr << "                          file containing the keypair." << std::endl;
    std::cerr << "   -l, --log-file FILE    enable statistics collection.  Output results to FILE" << std::endl;
    std::cerr << "                          (Default: stderr)" << std::endl;
    std::cerr << "   -L, --log-append STR   append the string STR to the end of all log lines." << std::endl;
    std::cerr << "   -H, --log-header       print a header in the log file if there are no other entries.\n";
    std::cerr << "       --help             display this help and exit." << std::endl;
    std::cerr << "       --version          output version information and exit." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Report bugs to iang+percy@cs.uwaterloo.ca." << std::endl;
    exit(status);
}

// List of long (i.e. --) options:
// {"longoptname", no_argument|required_argument|optional_argument, 0, 'shortoptname'}
struct option long_options[] = {
    {"version",     no_argument,        NULL, 'v'},
    {"help",        no_argument,        NULL, 'a'},
    {"k",	    required_argument,  NULL, 'k'},
    {"word-size",   required_argument,  NULL, 'w'},
    {"tau",         required_argument,  NULL, 't'},
    {"mode",        required_argument,  NULL, 'm'},
    {"batch-query", required_argument,  NULL, 'Q'},
    // hybrid: optionally specify the public-private key pair.
    {"hybrid",      optional_argument,  NULL, 'h'},
#ifdef SPIR_SUPPORT
    // spir: must specify PolyCommit parameters.
    {"spir",        required_argument,  NULL, 'Z'},
#endif
    {"log-file",    required_argument,  NULL, 'l'},
    {"log-append",  required_argument,  NULL, 'L'},
    {"log-header",  no_argument,	NULL, 'H'},
    {NULL,          0,                  NULL, 0},
};

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


int main(int argc, char **argv)
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
	print_usage_exit(argv[0]);
    }

    // Word size
    dbsize_t w = 8;

    // Servers needed
    nservers_t k = 0;
    
    // Do symmetric PIR?
    bool do_spir = false;

    // PolyCommit Params file to read for SPIR
    char *pcparams_file = NULL;
	
    // Enable hybrid security?
    bool do_hybrid = false;
	
    // Mode of operation selected (ZZ_p, GF28, GF216 or Chor)
    PercyMode mode = MODE_ZZ_P;
   
    // qbs indicates query block size
    nqueries_t qbs = 1;

    nservers_t tau = 0;

    // Logging
    bool do_logging = false;
    char * logfilename = NULL;
    char * logappend = NULL;
    bool print_header = false;
    
    // Read the optional flags.
    int opt;
    //In getopt's short options, ':' indicates the previous argument takes a required
    //argument, '::' indicates an optional one
    while((opt = getopt_long(argc, argv, "k:w:t:m:Q:h::s:l:d:L:H", long_options, 0)) != -1) {
        switch(opt) {
	    case 'k':
		k = strtoul(optarg, NULL, 10);
		if (k == 0) {
		    std::cerr << "k must be positive\n";
		    exit(1);
		}
		break;
            case 'w':
                w = strtoull(optarg, NULL, 10);
                break;
	    case 't':
		tau = strtoul(optarg, NULL, 10);
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
                    std::cerr << "Unknown mode selected. Valid modes are ZZ_P, GF28, GF216 and CHOR." << std::endl;
                    exit(1);
                }
		break;
	    case 'Q':
		qbs = strtoul(optarg, NULL, 10);
		break;
	    case 'h':
                do_hybrid = true;
		if (optarg)
		{
		    std::cerr << "Warning: hybrid with specified keys not yet implementing. Ignoring optional argument." << std::endl;
		}
		// optarg == keyfile
                break;
#ifdef SPIR_SUPPORT
	    case 'Z':
		do_spir = true;
                pcparams_file = optarg;
                break;
#endif
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
                std::cerr << "Percy++ pirclient version " << PERCY_VERSION << std::endl;
                std::cerr << AUTHOR << std::endl;
                exit(0);
                break;
            default:
                print_usage_exit(argv[0], 0);
                break;
        }
    }
	
    // Check for conflicting optional arguments.
    if (do_hybrid && (mode != MODE_ZZ_P)) {
	std::cerr << "Error: hybrid security can only be used with the integers mod p mode "
		"of operation." << std::endl;
	exit(1);
    }
#ifdef SPIR_SUPPORT
    if (do_hybrid && do_spir) {
	std::cerr << "Error: cannot use hybrid security with symmetric PIR." << std::endl;
	exit(1);
    }
    if (do_spir && mode != MODE_ZZ_P) {
	std::cerr << "Error: symmetric PIR can only be used with the integers mod p mode "
		"of operation." << std::endl;
	exit(1);
    }
#endif
    if (tau && mode == MODE_CHOR) {
	std::cerr << "Error: Chor et al.'s PIR scheme does not support tau independence." << std::endl;
	exit(1);
    }
	
    // Make sure enough mandatory arguments are present.
    if(argc - optind < 5)
    {
        print_usage_exit(argv[0]);
    }
	
    dbsize_t num_blocks = strtoull(argv[optind++], NULL, 10);
    dbsize_t block_size = strtoull(argv[optind++], NULL, 10);
    dbsize_t words_per_block = block_size * 8 / w;

    // Sanity checks for (n,b,w).
    bool badw = false;
    switch (mode) {
    case MODE_CHOR:
	badw = ( w != 1 );
	break;
    case MODE_GF28:
	badw = ( w != 8 );
	break;
    case MODE_GF216:
	badw = ( w != 16 );
	break;
    default:
	badw = ( w % 8 != 0 );
	break;
    }
    if (badw) {
	std::cerr << "Error: Invalid word size for mode " << mode << ": " << w << "\n";
	return 1;
    }
    
    // Compute the number of blocks, and number of words per block.
    //std::cerr << "Number of blocks: " << num_blocks << std::endl;
    //std::cerr << "Words per block:  " << words_per_block << std::endl;
    if (num_blocks != words_per_block)
    {
	std::cerr << "Warning: non-optimal choice of blocksize detected." << std::endl;
    }

    // Choose an appropriate modulus.
    ZZ p1, p2;
    if (w == 2048)
    {
	p1 = to_ZZ("208647130951457402363969335056365957472826150618980217460328400485971950387185944410889077723063406198415802830757517777351462262669194793047360775411639408116452523756687066355086195124187048682420529316060567502352699557841412039275095485224490337148164650010000499984813523719988826268799665657866626493329");
	p2 = to_ZZ("245210205383950153265232956846271987008710436579074459102383753214859717124121302267932009072054546430711727811323033561244148876933172687995163379778095734152594201215411509169035373484564340604271927100344464582888777887746436564737355045100633587336239754449508771770564607896955672999950235015535154415867");
    }
    else if (w == 1536)
    {
	p1 = to_ZZ("1762848592595080314705600925431624874456855439794595868418019480189213868063348394981842423875338178991362893712297567682392276281463789141688306484765105096429658863055172316227409205756175078509101834587188923103831602929062176351");
	p2 = to_ZZ("2306072568237159640249655953989533876736033293267891813492402870069702343561490811306173449455816459207943593196801405361355605814646339972518285709494570145269396000374210678514250118174550977925517522811232946347459478425104006037");
    }
    else if (w == 1024)
    {
	p1 = to_ZZ("14710132128541592475387440366744304824352604767753216777226640368050037133836174845369895150342922969891066267019166301546403100960464521216972792406229873");
	p2 = to_ZZ("23338263930359653850870152235447790059566299230626918909126959284529524161146399225204807633841208114717867386116272471253667601589249734409576687328158817");
    }
    else if (mode == MODE_GF28)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("256");
    }
    else if (mode == MODE_GF216)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("65536");
    }
    else if (mode == MODE_CHOR)
    {
	//Important: for Chor we pretend as though a word is 1 byte. This is because many 
	//parts of the code rely on a word being a byte multiple (for example, where bytes_per_word 
	//is used). We set this here since the calculations for the optimal database shape need a 
	//word size of 1 bit.
	//words_per_block /= 8;
	p1 = to_ZZ("1");
	p2 = to_ZZ("256");
    }
    else if (w == 8 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("257");
    }
    else if (w == 16 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("65537");
    }
    else if (w == 32 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("4294967311");
    }
    else if (w == 96 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("79228162514264337593543950397");
    }
    else if (w == 128 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("340282366920938463463374607431768211507");
    }
    else if (w == 160 && !do_hybrid)
    {
	// NOTE: p2s is the prime from the PolyCommit params; spir
	// will break if this value gets changed!
	//
	// TODO: read the prime from the PolyCommit params and check
	// 		 that it is consistent with w.
	p1 = to_ZZ("1");
	p2 = to_ZZ("2425980306017163398341728799446792216592523285797");
    }
    else if (w == 192 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("6277101735386680763835789423207666416102355444464034513029");
    }
    else if (w == 256 && !do_hybrid)
    {
	p1 = to_ZZ("1");
	p2 = to_ZZ("115792089237316195423570985008687907853269984665640564039457584007913129640233");
    }
    else if (do_hybrid)
    {
	std::cerr << "Error: No hybrid-compatible modulus available for w = " << w << "." << std::endl;
	exit(1);
    }
    else
    {
	std::cerr << "Error: No modulus available for w = " << w << "." << std::endl;
	exit(1);
    }
#ifdef SPIR_SUPPORT
    if (do_spir && w!=160)
    {
	std::cerr << "Error: symmetric PIR currently supports only w=160." << std::endl;
	exit(1);
    }
#endif
    ZZ modulus = p1 * p2;
    
    unsigned long modulusbytes = NumBytes(modulus);
    if (modulusbytes <= w/8)
    {
	std::cerr << "Error: w must be at most " << (modulusbytes-1)*8 << " (was " << w << ")." << std::endl;
	exit(1);
    }

    // Parse the server information
    vector<char*> servers;
    char *token = strtok(argv[optind++], " ");
    while (token != NULL) {
	servers.push_back(token);
	token = strtok(NULL, " ");
    }
    nservers_t ell = servers.size();
    sid_t * server_indices = new sid_t[ell];
    vector<serverinfo> sinfos;
    for (nservers_t j = 0; j < ell; j++) {
	char *sids = strtok(servers[j], ":");
	char *addrs = strtok(NULL, ":");
	char *ports = strtok(NULL, ":");
	if (ports == NULL) {
	    std::cerr << "Error: something is wrong with the server information in position " << (j+1) << "." << std::endl;
	    delete[] server_indices;
	    exit(1);
	}
	unsigned long sid = strtoul(sids, NULL, 10);
	server_indices[j] = (sid_t)sid;
	if (!sid || sid > modulus) {
	    std::cerr << "Error: SID must be an integer greater than 0 and less than " << modulus << "." << std::endl;
	    delete[] server_indices;
	    exit(1);
	}
	unsigned long port = strtoul(ports, NULL, 10);
	if (port < 1024 || port > 65535) {
	    std::cerr << "Error: port number must be an integer greater than 1024 and less than 65536." << std::endl;
	    delete[] server_indices;
	    exit(1);
	}
	serverinfo sinfo;
	sinfo.sid = sid;
	sinfo.addr = addrs;
	sinfo.port = port;
	sinfos.push_back(sinfo);
    }
    // Set up an iostream to each of the servers.
    vector<serverinfo> onlinesinfos;
    vector<iosockinet*> serverstreams;
    vector<istream*> istreams;
    vector<ostream*> ostreams;
    for (nservers_t j = 0; j < ell; j++)
    {
	iosockinet *socket = connect_to_server(sinfos[j]);
	std::iostream * stream = socket;
        if ( socket != NULL ) {
	    onlinesinfos.push_back(sinfos[j]);
	    serverstreams.push_back(socket);
	    istreams.push_back(stream);
	    ostreams.push_back(stream);
	}
    }

    if (k == 0) {
	k = ell;
    }
    nservers_t t = strtoul(argv[optind++], NULL, 10);

    if (serverstreams.size() < k)
    {
	std::cerr << "Error: fewer than k=" << k << " out of ell=" << ell << " servers are online." << std::endl;
	for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	    if (serverstreams[i]) delete serverstreams[i];
	}
	delete[] server_indices;
	exit(1);
    }

    // Sanity checks for (ell,t,k,tau).
    if (qbs-1+t+tau >= k)
    {
	std::cerr << "Error: qbs-1+t+tau must be less than k." << std::endl;
	for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	    if (serverstreams[i]) delete serverstreams[i];
	}
	delete[] server_indices;
	exit(1);
    }
    if (mode == MODE_CHOR && t != k-1) {
        std::cerr << "Error: Chor requires that t=k-1." << std::endl;
	for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	    if (serverstreams[i]) delete serverstreams[i];
	}
	delete[] server_indices;
        exit(1);
    }
    if (mode == MODE_CHOR && qbs != 1) {
        std::cerr << "Error: Chor requires that qbs=1." << std::endl;
        exit(1);
    }

    // Get (and perform sanity check on) the query indices.
    vector<dbsize_t> indices;
    istringstream indexss(argv[optind++]);
    while(!indexss.eof())
    {
	dbsize_t index;
	indexss >> index;
	indices.push_back(index);
    }
    for (nqueries_t q = 0; q < indices.size(); q++)
    {
	if (indices[q] >= num_blocks) {
	    std::cerr << "Error: all query indices must be less than number of blocks." << std::endl;
	    for (nservers_t i = 0; i < serverstreams.size(); ++i) {
		if (serverstreams[i]) delete serverstreams[i];
	    }
	    delete[] server_indices;
	    exit(1);
	}
    }
    std::cerr << "Fetching " << indices.size() << " blocks." << std::endl;
    
    // Create the PercyClientParams object.
    PercyParams * params = NULL;
    switch (mode) {
    case MODE_ZZ_P:
	if (do_hybrid) {
	    params = new ZZ_pParams(num_blocks, words_per_block * w / 8,
		    w, p1, p2, tau);
	} else {
	    params = new ZZ_pParams(num_blocks, words_per_block * w / 8,
		    w, modulus, tau, pcparams_file, do_spir);
	}
	break;
    case MODE_CHOR:
	    params = new ChorParams(num_blocks, words_per_block * w / 8);
	break;
    case MODE_GF28:
    case MODE_GF216:
	    params = new GF2EParams(num_blocks, words_per_block * w / 8, w, tau);
	break;
    default:
	std::cerr << "Invalid mode: " << mode << "\n";
	for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	    if (serverstreams[i]) delete serverstreams[i];
	}
	delete[] server_indices;
	exit(1);
	break;
    }
    if (!params) {
	std::cerr << "Error creating parameters object\n";
	for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	    if (serverstreams[i]) delete serverstreams[i];
	}
	delete[] server_indices;
	exit(1);
    }
    PercyClientParams clientparams (params, serverstreams.size());
    

    // Send each server the parameters (and its SID, if necessary).
    for (nservers_t j = 0; j < serverstreams.size(); j++) {
	std::ostream &os = *ostreams[j];
	
	// Send the PercyClientParams.
//		std::cerr << "Sending query parameters to server " << onlinesinfos[j].sid << "...";
	clientparams.send(os, onlinesinfos[j].sid);
//		std::cerr << "done" << std::endl;
    
	std::istream &is = *istreams[j];
	unsigned char failure;
//		std::cerr << "Receiving response from server " << onlinesinfos[j].sid << "...";
	is.read((char*)&failure, 1);
//		std::cerr << "done" << std::endl;
	if (failure) {
	    std::cerr << "Error: " << onlinesinfos[j].addr << ":" << onlinesinfos[j].port << " did not accept parameters." <<  std::endl;
	    delete params;
	    for (nservers_t i = 0; i < serverstreams.size(); ++i) {
		if (serverstreams[i]) delete serverstreams[i];
	    }
	    delete[] server_indices;
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
	    stats = new PercyClientStats(&clientparams, logstream, logappend);
	    if (print_header) {
		stats->print_header();
	    }
	} else {
	    stats = new PercyClientStats(&clientparams, logappend);
	}
    }

    // Finally, do the PIR query!
    int ret = 0;
    PercyClient * client = PercyClient::make_client(&clientparams, ell, t, server_indices, stats);
    delete[] server_indices;
    
    if (client == NULL) {
	std::cerr << "Error: A client was not created\n";
	if (stats) {
	    if (logfilename != NULL) {
		logstream.close();
	    }
	    delete stats;
	}
	delete params;
	for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	    if (serverstreams[i]) delete serverstreams[i];
	}
	exit(1);
    }
    srand(time(NULL));

    // Get the value of h
    nservers_t h = (nservers_t)(floor(sqrt((t+tau)*k)))+1;
    const char *envh = getenv("PIRC_H");
    if (envh) {
        nservers_t override_h = atoi(envh);
        if (override_h > 0) {
            h = override_h;
        }
    }

    vector<PercyBlockResults> results, current, original;
#ifdef VERBOSE_PIRCLIENT
    std::cerr << "INDICES = [";
#endif
    for (nqueries_t q = 0; q < indices.size(); ++q) {
	results.push_back(PercyBlockResults());
	results.back().block_number = indices[q];
#ifdef VERBOSE_PIRCLIENT
	std::cerr << " " << indices[q];
#endif
    }
#ifdef VERBOSE_PIRCLIENT
    std::cerr << " ]\n";
#endif
    std::set<dbsize_t> undecoded (indices.begin(), indices.end());

    // Fetch the blocks
    nqueries_t request_identifier;
    dbsize_t res = client->fetch_blocks(request_identifier, indices, ostreams, 
	    istreams, current, qbs);
#ifdef VERBOSE_PIRCLIENT
    std::cerr << "CURRENT =\n" << current;
#endif

    if (res) {
	results = current;
    } else {
	for (nqueries_t m = 0; ; ++m) {
#ifdef VERBOSE_PIRCLIENT
	    std::cerr << "CURRENT =\n" << current;
	    std::cerr << "ORIGINAL =\n" << original;
#endif
	    vector<PercyBlockResults>::const_iterator iter;
	    // Find decoded queries in original
	    for (iter = original.begin(); iter != original.end(); ++iter) {
		if (!(iter->results.empty())) {
		    if (undecoded.find(iter->block_number) != undecoded.end()) {
			for (nqueries_t q = 0; q < indices.size(); ++q) {
			    if (indices[q] == iter->block_number) {
				results[q].results = iter->results;
			    }
			}
			undecoded.erase(iter->block_number);
		    }
		}
	    }
	    // Find decoded queries in current
	    for (iter = current.begin(); iter != current.end(); ++iter) {
		if (!(iter->results.empty())) {
		    if (undecoded.find(iter->block_number) != undecoded.end()) {
			for (nqueries_t q = 0; q < indices.size(); ++q) {
			    if (indices[q] == iter->block_number) {
				results[q].results = iter->results;
			    }
			}
			undecoded.erase(iter->block_number);
		    }
		}
	    }
#ifdef VERBOSE_PIRCLIENT
	    std::cerr << "RESULTS = " << results;
#endif

	    // Check if done
	    if (undecoded.empty()) {
		break;
	    }

	    // E = w * (m(h-t-1)-k+h+1)
	    // where w is the bitlength of a word and m is the number of
	    // resent queries
	    // But here, the degree of the poly (t in the above formula) is t+tau+qbs-1
	    int E = w * (m * (h - t - params->tau() - qbs) - ell + h + 1);
	    if (E > QUERY_LIMIT_EXP) {
		std::cerr << "Reached the maximum number of re-sent queries.\n";
		std::cerr << "Too few honest servers to recover data!\n";
		delete client;
		if (stats) {
		    if (logfilename != NULL) {
			logstream.close();
		    }
		    delete stats;
		}
		delete params;
		for (nservers_t i = 0; i < serverstreams.size(); ++i) {
		    if (serverstreams[i]) delete serverstreams[i];
		}
		exit(1);
	    }

	    // Choose block to requery
	    vector<dbsize_t> requery_indices;
	    nqueries_t rq_index = rand() % undecoded.size();
	    std::set<dbsize_t>::iterator ud_iter = undecoded.begin();
	    std::advance(ud_iter, rq_index);
	    dbsize_t rq = *ud_iter;
	    requery_indices.push_back(rq);
	    std::cerr << "Requerying for block number '" << rq << "'\n";

	    // Requery
	    current.clear();
	    original.clear();
	    nqueries_t requery_identifier;
	    client->fetch_blocks(requery_identifier, requery_indices, 
		    ostreams, istreams, current);
	    client->get_result(request_identifier, original);
	}
    }

    nqueries_t num_res = results.size();
    for (nqueries_t r=0; r < num_res; ++r) {
	if (results[r].results.empty()) {
	    std::cerr << "PIR query failed.\n";
	    ret = 1;
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
	
    // Tidy up after query.
    delete client;
    if (stats) {
	if (logfilename != NULL) {
	    logstream.close();
	}
	delete stats;
    }
    delete params;
    for (nservers_t i = 0; i < serverstreams.size(); ++i) {
	if (serverstreams[i]) delete serverstreams[i];
    }

    std::cerr << "Client shutting down." << std::endl;
    (void)res;
    return ret;
}
