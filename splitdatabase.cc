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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <NTL/ZZ_pX.h>
#include <ctype.h>
#include <getopt.h>
#include "percyio.h"
#include "percytypes.h"
#include "config.h"
#include "gf2e.h"

NTL_CLIENT

// Split a given database into l pieces with tau-independence: no
// coalition of up to tau of the pieces has any information about the
// contents of the original database.  You must supply a modulus to use
// for the splitting, and that same modulus must be used during the PIR
// protocol to query the pieces.
//
// If the original database name is "dbname", the pieces will be called
// "dbname.1", "dbname.2", etc.

ZZ zz_p_modulus (dbsize_t word_size) {
    if (word_size == 2048) {
	return to_ZZ("51162405833378812589599605953260132300166393994651819099454781579567509212081792013783783759303440508155949594262147212874957344953142209597742684263402581129339826752613431877280173074502314648334418584122460414512816448592261381117519846844295394134225624418756277265452922709245846828145574822031541004633366879073894273715489429502290966133193310966178373909137394353164436844312924586836474134940807305776164928781025210917912257206480517698118422827367766257579221703667784216949825206167241852365543481875593117676222875888924950402025039269210778276794873837063438751454865130720887819939394489366347567251243");
    } else if (word_size == 1536) {
	return to_ZZ("4065256781338999183533854850423382625119065920051798531476300569026463202897155088318466013703570859212040475097762405522038651420119366364979939687154236065682459920101982590074846996306687236388206057475890613264408059472973401701686869808348910896596468985609043697525749128687318350246421674945679872669881805678484464202726328189280359385791023305618545788872763420795247846720674554774196715770302797683129209164871258189464484019233379849839076263862630987");
    } else if (word_size == 1024) {
	return to_ZZ("343308946066366926839932845260501528909643718159825813630709694160026342456154871924497152436552679706642965502704642456637620829912957820221098686748075257358288200837461739492534713539606088624083011849535450485951774635526473457667739540374042376629835941950802202870595346459371144019363420985729553740241");
    } else if (word_size == 256) {
	return to_ZZ("115792089237316195423570985008687907853269984665640564039457584007913129640233");
    } else if (word_size == 196) {
	return to_ZZ("6277101735386680763835789423207666416102355444464034513029");
    } else if (word_size == 160) {
	return to_ZZ("2425980306017163398341728799446792216592523285797");
    } else if (word_size == 128) {
	return to_ZZ("340282366920938463463374607431768211507");
    } else if (word_size == 96) {
	return to_ZZ("79228162514264337593543950397");
    } else if (word_size == 32) {
	return to_ZZ("4294967311");
    } else if (word_size == 16) {
	return to_ZZ("65537");
    } else if (word_size == 8) {
	return to_ZZ("257");
    } else {
	std::cerr << "Invalid word size for ZZ_P: " << word_size << "\n";
	return ZZ::zero();
    }
}

void split_zz_p (std::ifstream& infile, std::ofstream* outfile, nservers_t tau, 
	nservers_t l, ZZ modulus, dbsize_t dbsize)
{
    streamsize outbytes = NumBytes(modulus);
    streamsize inbytes = outbytes-1;
    if (outbytes < 2) {
	std::cerr << "Error: modulus is too small; must be at least 256.\n";
	exit(1);
    }

#if 0
    std::cout << "Converting " << inbytes << "-byte inputs to " << outbytes
	<< "-byte outputs.\n";
#endif

    char c = 'z';
    for(nservers_t i=1;i<=l;++i) {
	outfile[i-1].write(&c, 1);
	percy_write_ZZ(outfile[i-1], modulus);
    }

    unsigned char * inbuf = new unsigned char[inbytes];
    unsigned char * outbuf = new unsigned char[outbytes];

    dbsize_t bytesread = 0;
    while(1) {
	if (dbsize != 0 && bytesread >= dbsize) break;

	// Read the input in chunks of inbytes bytes
	infile.read((char *)inbuf, inbytes);
	bytesread += inbytes;

	if (infile.gcount() == 0) break;

	// Convert the chunk to a ZZ
	ZZ Wz = ZZFromBytes(inbuf, infile.gcount());
	// cout << "reading " << Wz << "\n";

	// Pick a random polynomial of degree tau, and set the constant
	// coefficient to the ZZ read from the input database
	ZZ_pX randpoly = random_ZZ_pX(tau+1);
	SetCoeff(randpoly, 0, to_ZZ_p(Wz));
	// cout << "poly(" << tau << ") = " << randpoly << "\n";

	// For each output file i, output the value of the polynomial
	// evaluated at i.
	for(nservers_t i=1;i<=l;++i) {
	    ZZ_p value;
	    eval(value, randpoly, to_ZZ_p(i));
	    BytesFromZZ(outbuf, rep(value), outbytes);
	    outfile[i-1].write((char *)outbuf, outbytes);
	    // cout << "writing " << i << "/" << tau << ": " << rep(value) << "\n";
	}
    }

    for(nservers_t i=0;i<l;++i) {
	outfile[i].close();
    }
    infile.close();

    delete[] outbuf;
    delete[] inbuf;
}

template <typename GF2E_Element>
void split_gf2e (std::ifstream& infile, std::ofstream* outfile, nservers_t tau, 
	nservers_t l, dbsize_t word_size, dbsize_t dbsize)
{
    char c = 'g';
    for(nservers_t i=1;i<=l;++i) {
	outfile[i-1].write(&c, 1);
	PERCY_WRITE_LE_DBSIZE(outfile[i-1], word_size);
    }

    GF2E_Element * coeffs = new GF2E_Element[tau+1];
    GF2E_Element * indices = new GF2E_Element[l];
    for (nservers_t i = 1; i <= l; ++i) {
	indices[i-1] = (GF2E_Element)i;
    }

    dbsize_t bytesread = 0;
    while (1) {
	if (dbsize != 0 && bytesread >= dbsize) break;

	// Read word
	infile.read((char*)coeffs, sizeof(GF2E_Element));
	bytesread += sizeof(GF2E_Element);

	if (infile.gcount() == 0) break;

	// Pick a random polynomial of degree tau, and set the constant
	// coefficient to the word read from infile
	for (nservers_t i = 0; i < tau; ++i) {
	    coeffs[i+1] = RandomBits_ulong(8*sizeof(GF2E_Element));
	}

	// For each output file i, output the value of the polynomial
	// evaluated at i.
	for(nservers_t i=1;i<=l;++i) {
	    GF2E_Element outword = 
		    evalpoly_GF2E<GF2E_Element>(coeffs, tau, indices[i-1]);
	    outfile[i-1].write((char *)(&outword), sizeof(GF2E_Element));
	    // cout << "writing " << i << "/" << tau << ": " << rep(value) << "\n";
	}
    }

    delete[] coeffs;
    delete[] indices;
}

dbsize_t parse_database_size (const char * cdbsize)
{
    std::stringstream ssdbsize;
    ssdbsize << cdbsize;
    dbsize_t dbsize;
    std::string dbunits;
    ssdbsize >> dbsize;
    ssdbsize >> dbunits;
    if (dbunits == "") {
	// Do nothing
    } else if (dbunits == "kB" || dbunits == "KB") {
	dbsize *= 1000;
    } else if (dbunits == "k" || dbunits == "K" || dbunits == "KiB") {
	dbsize *= 1024;
    } else if (dbunits == "MB") {
	dbsize *= 1000000;
    } else if (dbunits == "M" || dbunits == "MiB") {
	dbsize *= 1048576;
    } else if (dbunits == "GB") {
	dbsize *= 1000000000;
    } else if (dbunits == "G" || dbunits == "GiB") {
	dbsize *= 1073741824;
    } else if (dbunits == "TB") {
	dbsize *= 1000000;
	dbsize *= 1000000;
    } else if (dbunits == "T" || dbunits == "TiB") {
	dbsize *= 1048576;
	dbsize *= 1048576;
    } else {
	std::cerr << "Invalid database size: " << cdbsize << "\n";
	return 0;
    }
    return dbsize;
}

void usage (const char * bin) {
    std::cerr << "Usage: " << bin << " [OPTIONS] database tau l [database_size]\n\n"
	      << "Split the database for tau-independence with l servers.  That is, if the resulting\n"
	      << "databases are given to different PIR servers, then any coalition of at most tau of\n"
	      << "the l servers will not be able to learn the contents of the database.\n\n"
	      << "Available Options:\n"
	      << "   -m, --mode MODE        use the specified mode of operation. Supported modes are:\n"
	      << "                          Long form   Short form   Description\n"
	      << "                          GF28        g            use fast arithmetic in GF(2^8)\n"
	      << "                          GF216       s            use fast arithmetic in GF(2^16)\n"
	      << "                          ZZ_P        z            use arithmetic in Z mod p\n"
	      << "                          (default: ZZ_P)\n"
	      << "   -w, --word-size WS     use a word size of WS bits for ZZ_P.  (default: 1024)\n"
	      << "   -p, --modulus MOD      use a modulus of MOD for ZZ_P.  Will override any value\n"
	      << "                          of --word-size.\n"
	      << "   --help		    display this help and exit.\n"
	      << "   --version		    output version information and exit.\n\n"
	      << "Report bugs to iang+percy@cs.uwaterloo.ca.\n"
	      ;
}

struct option long_options[] = {
    {"mode",		required_argument,	NULL,   'm'},
    {"word-size",	required_argument,	NULL,   'w'},
    {"modulus",		required_argument,	NULL,	'p'},
    {"help",		no_argument,		NULL,   'h'},
    {"version",		no_argument,		NULL,   'v'},
    {NULL,		0,			NULL,   0}
};

const char * short_options = "m:w:p:";

int main(int argc, char **argv)
{
    char mode = 'z';
    dbsize_t word_size = 1024;
    ZZ modulus = ZZ::zero();

    int opt;
    while ((opt = getopt_long(argc, argv, short_options, long_options, 0)) != -1) {
	switch (opt) {
	case 'm':
	    if(!strcmp(optarg, "ZZ_P") || !strcmp(optarg, "z")) {
		mode = 'z';
	    }
	    else if(!strcmp(optarg, "GF28") || !strcmp(optarg, "g")) {
		mode = 'g';
	    }
	    else if(!strcmp(optarg, "GF216") || !strcmp(optarg, "s")) {
		mode = 's';
	    }
	    else {
		std::cerr << "Unknown mode selected.\n\n";
		usage(argv[0]);
		exit(1);
	    }
	    break;
	case 'w':
	    word_size = strtoull(optarg, NULL, 10);
	    break;
	case 'p': {
		stringstream ss (stringstream::in | stringstream::out);
		ss << argv[5];
		ss >> modulus;
	    } break;
	case 'h':
	    usage(argv[0]);
	    exit(0);
	    break;
	case 'v':
	    std::cerr << "Percy++ splitdatabase version " << PERCY_VERSION << "\n";
	    std::cerr << AUTHOR << "\n";
	    exit(0);
	    break;
	default:
	    std::cerr << "Invalid option: " << opt << "\n";
	    exit(1);
	}
    }

    if (argc - optind < 3) {
	usage(argv[0]);
	exit(1);
    }

    char * database = argv[optind++];
    nservers_t tau = strtoul(argv[optind++], NULL, 10);
    nservers_t l = strtoul(argv[optind++], NULL, 10);
    dbsize_t dbsize = 0;
    if (argc - optind > 0) {
	dbsize = parse_database_size(argv[optind++]);
	if (!dbsize) {
	    exit(1);
	}
    }

    if (IsZero(modulus)) {
	modulus = zz_p_modulus(word_size);
	if (IsZero(modulus)) {
	    exit(1);
	}
    }

    // Open the input database
    ifstream infile (database);
    int flen = strlen(database);

    // Create l output databases
    ofstream * outfile = new ofstream[l];
    nservers_t i;
    char * ofile = new char[flen+13];
    for(i=1;i<=l;++i) {
	sprintf(ofile, "%s.%d", database, i);
	outfile[i-1].open(ofile);
	// Write the tau-independence header for each output database
	outfile[i-1].write("PIRD\x01\x00", 6);
    }
    delete[] ofile;

    switch (mode) {
    case 'z': // ZZ_P
	ZZ_p::init(modulus);
	split_zz_p(infile, outfile, tau, l, modulus, dbsize);
	break;
    case 'g': // GF28
	word_size = 8;
	split_gf2e<GF28_Element>(infile, outfile, tau, l, word_size, 
		dbsize);
	break;
    case 's': // GF216
	word_size = 16;
	split_gf2e<GF216_Element>(infile, outfile, tau, l, word_size, 
		dbsize);
	break;
    default:
	std::cerr << "Invalid mode: " << mode << "\n\n";
	usage(argv[0]);
	exit(1);
    }

    delete[] outfile;

    return 0;
}

