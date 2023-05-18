# Getting Percy++ working
Percy++ uses the following libraries; so will need to install them 
before we can use Percy++:

- NTL (http://www.shoup.net/ntl/)
- Socket++ (http://www.linuxhacker.at/socketxx/)
- libgcrypt (http://www.gnu.org/software/libgcrypt/)

Install GMP if not already
```
sudo apt-get install libgmp-dev
```

### Get NTL library
```
wget https://libntl.org/ntl-11.5.1.tar.gz
gunzip ntl-11.5.1.tar.gz 
tar xf ntl-11.5.1.tar
cd ntl-xxx/src
./configure 
make 
make check
sudo make install
```
### Get socketxx
```
git clone https://github.com/pawangoyal137/socketxx
sudo apt-get install libtool
cd socketxx
./autogen
./configure --prefix=/usr --enable-debug
# Switch sys_siglist to strsignal in socket++/config.h and remove sys_errlist in  socket++/local.h socket++/config.h and socket++/socketstream.ccp
make
sudo make install
```

### Get libgcrypt

```
# install libgpg error required for libgcrypt
wget https://www.gnupg.org/ftp/gcrypt/libgpg-error/libgpg-error-1.47.tar.bz2
tar -xjf libgpg-error-1.47.tar.bz2 
cd libgpg-error-1.47
./configure --prefix=/usr
make
make check 
sudo make install
# get the readme
sudo install -v -m644 -D README /usr/share/doc/libgpg-error-1.47/README

# get libgcrypt
wget https://www.gnupg.org/ftp/gcrypt/libgcrypt/libgcrypt-1.10.2.tar.bz2
tar -xjf libgcrypt-1.10.2.tar.bz2
cd libgcrypt-1.10.2/
./configure --prefix=/usr
make
make -C doc html
sudo apt install texinfo
makeinfo --html --no-split -o doc/gcrypt_nochunks.html doc/gcrypt.texi
makeinfo --plaintext       -o doc/gcrypt.txt           doc/gcrypt.texi
make check
sudo make install
```

------------------------------------------------
-----------------
----------------
# Below is the original Readme and can be found at http://downloads.sourceforge.net/percy/percy++-1.0.0.tar.gz

Percy++ / PIR in C++

Ian Goldberg <iang@cs.uwaterloo.ca>
Casey Devet <cjdevet@cs.uwaterloo.ca>
Wouter Lueks <wouter@telox.net>
Ann Yang <y242yang@uwaterloo.ca>
Paul Hendry <pshendry@uwaterloo.ca>
Ryan Henry <rhenry@cs.uwaterloo.ca>

Version 1.0: 2014-10-17


About Percy++
-------------

Percy++ is an implementation of the private information retrieval (PIR)
protocols from the papers:

    Ian Goldberg.  Improving the Robustness of Private Information
    Retrieval.  Proc. of 2007 IEEE Symposium on Security and Privacy
    (Oakland 2007), May 2007.

    Ryan Henry, Femi Olumofin, Ian Goldberg.  Practical PIR for
    Electronic Commerce.  18th ACM Conference on Computer and
    Communications Security, October 2011.

    Casey Devet, Ian Goldberg, and Nadia Heninger.  Optimally Robust
    Private Information Retrieval.  21st USENIX Security Symposium,
    August 2012.

    Carlos Aguilar Melchor and Philippe Gaborit.  A Lattice-Based
    Computationally-Efficient Private Information Retrieval Protocol.
    WEWORC 2007, July 2007.

    Casey Devet and Ian Goldberg. The Best of Both Worlds: Combining
    Information-Theoretic and Computational PIR for Communication
    Efficiency. 14th Privacy Enhancing Technologies Symposium (PETS
    2014), July 2014.

    Wouter Lueks, Ian Goldberg. "Sublinear Scaling for Multi-Client
    Private Information Retrieval". CACR Tech Report 2014-19,
    September 2014.


Briefly, private information retrieval is the task of fetching a block
of data from a database server (or a group of distributed servers)
without the server(s) learning which block it was that you were
interested in.

These protocols provide t-private v-Byzantine-robust tau-independent
k-out-of-l private information retrieval.  This means:

k-out-of-l: there are l distributed database servers, and we only need
            to receive replies from k of them (the rest might be down,
	    overloaded, unreachable, etc.)

t-private:  no coalition of up to t servers receives *any information at
            all* about the block you are interested in

v-Byzantine-robust: up to v of the servers that do reply might give
		    *incorrect* answers; we will want to detect which
		    servers did that, and to determine the correct
		    database block

tau-independent: the database is split between the servers so that no
                 coalition of up to tau of them can determine the
		 contents of the database itself (tau=0 means all the
		 servers just have a complete copy of the database)

All of the above are "information-theoretic"; that is, the protections
hold even if the servers have unlimited computational power, so long as
no more than t server are colluding to determine the client's query.

Any choice of t, v, tau, k and l will work, so long as they satisfy the
following conditions:

- They are all nonnegative integers.
- 0 < t <= t + tau < k <= l
- 0 <= v < k - t - tau - 1

This library also supports "computational" PIR, in which there is a
single server, and it cannot learn the contents of the query, but the
scheme has no robustness.  The security of this scheme relies on certain
lattice-based cryptographic assumptions.

Finally, it support "hybrid" PIR, which combines an
information-theoretic PIR scheme with a computational one in order to
hedge against either the non-collusion assumption or the cryptographic
assumption being violated.

Percy++ is written in C++, using Victor Shoup's NTL library.


Building and Using Percy++
--------------------------

Percy++ uses the following libraries; you will need to install them 
before you can use Percy++:

- NTL (http://www.shoup.net/ntl/)
- Socket++ (http://www.linuxhacker.at/socketxx/)
- libgcrypt (http://www.gnu.org/software/libgcrypt/)

If Symmetric PIR support is desired, then the following additional 
libraries are required:

- PBC (http://crypto.stanford.edu/pbc/)
- PBCWrapper and PolyCommit (see below)

In addition, if you tell NTL to use the GMP library 
when you install it, Percy++'s calculations will be faster.

Percy++ assumes that NTL headers are located in /usr/local/include/NTL,
and that the above libraries are located in /usr/local/lib. If this is
not the case, you should edit the Makefile accordingly.

If you are building on a 32-bit machine, set BIT32_SUPPORT to true in
the Makefile.  (Many parts of Percy++ are considerably faster on 64-bit
machines.)

Once the libraries are installed, running "make" should build the
programs pirclient, pirserver, agclient, agserver, hybridclient,
hybridserver, and splitdatabase.

The *client and *server programs are sample frontends to the Percy++
library, meant primarily as a demonstration.  Usage information can be
obtained by running the commands with no arguments.

The pirclient/pirserver pair implement information-theoretic PIR, as
presented in the above papers [Goldberg 2007], [Henry et al. 2011],
[Devet et al. 2012].  To overwrite the default value of h (the minimum
number of honest servers), set the environment variable PIRC_H to be the
desired value.

The agclient/agserver pair implement computational PIR, as presented in
[Aguilar Melchor and Gaborit 2007].

The hybridclient/hybridserver pair implement a hybrid
information-theoretic/computation PIR, as presented in [Devet and
Goldberg 2014].


The database is kept in a single file called "database".  You can
populate that file however you like; copying some amount of data from
/dev/urandom is fine.

The runtest and runtestset script are powerful wrappers that allow you
to run tests on large numbers of combinations of modes and parameters,
and gather detailed timing information.  For more information about
running these programs, run them with no arguments.

To try the tau-independence property, you will need to split the
database file into l pieces using the splitdatabase command (run it with
no arguments to see its invocation options).  For example, running
"./splitdatabase database 1 5" will split the database into 5 files,
named database.1, database.2, etc., with 1-independence.  That is, no
single database server could determine the contents of the database, but
more than one could.  Note that this is entirely separate from the t
parameter, which controls how many servers can collude without being
able to learn the value of your query (as opposed to the contents of the
database).  Then using the "-t 1" option to runtest (the parameter after
-t is the value of tau) will test tau-independence.

To try symmetric PIR, first ensure that SPIR_SUPPORT is set to 'true' in
the Makefile, and additionally that PBCWrapper and PolyCommit (available
in a separate package from Percy++) are located in the same directory as
the percy diretory. To try SPIR, run pirclient and pirserver with the
additional option '-Z polycommit.params'.  Note that using symmetric PIR
is much more computationally expensive than not using it.

The server can distribute its computation over a group of worker
servers.  This is done by running the workers and then running
pirserver_master with the correct information.  (Use "make
pirserver_master" to build this program.  See the usage message for
pirserver_master for more).

The server can also use threading to distribute its computation.  This
can be done using the -T tag for pirserver. (See the usage message for
pirserver for more details.)

Feel free to send question, bug reports, patches, etc. to Ian Goldberg
<iang+percy@cs.uwaterloo.ca>.


Reed-Solomon Decoding
---------------------

The RSDecoder class is an implementation of Reed-Solomon decoding.  We
have implemented many decoding algorithms, including the algorithm
described in the paper by Devet, Goldberg and Heninger (listed above).
For more information about the API of RSDecoder, see the
doxygen-generated documentation.


Changelog
---------

Version 1.0 (2014-10-17):
    - Added API documentation.  Install doxygen 1.8.x and use 'make
      docs' to create the documentation.
    - Significantly changed the API of the classes.  See the
      documentation for details.
    - Changed to the usage of executables, including pirclient and
      pirserver.  Run the executables with --help flag to see usage.
    - Implemented the computational PIR protocol by Aguilar Melchor and
      Gaborit (2007).  See agclient and agserver executables.
    - Implemented the hybrid PIR protocol by Devet and Goldberg (2014).
      See hybridclient and hybridserver executables.
    - Added multiblock querying for Goldberg's IT-PIR protocol (modes
      GF28, GF216, ZZ_P).  Can now request multiple blocks in a single
      query.  Use --batch-query option of pirclient.
    - Implemented Strassen's fast matrix multiplication algorithm for
      server computation in Goldberg's IT-PIR protocol (modes GF28,
      GF216, ZZ_P), based on Lueks and Goldberg (2014).
    - Added statistics collection to the servers and clients.  See usage
      messages of executables for details.
    - Replaced testclient script with runtest and runtestset scripts.
      Run these scripts with --help to see usage messages.

    Known bugs: SPIR support may not be working correctly in this
    version.  The doxygen documentation is incomplete.

Version 0.9 (2013-06-07):
    - Added support for PIR servers that use multiple worker hosts to do
      their computations faster, and/or use multiple threads or
      processes per host
    - Considerably cleaned up the APIs on the client and server sides
    - "make libs" will now build separate Percy++ client and server
      libraries for linking into your own programs
    - Removed variable-length arrays for compiler portability
    - Improved the speed of queries for multiple blocks in GF(2^8)
    - Added "-1 / --oneconn" option to pirserver to accept a single
      client connection and not fork (useful for debugging)

Version 0.8 (2012-06-29):
    - Added support for Symmetric PIR, fast arithmetic in GF(2^16), 
      and Chor et. al's lightweight protocol
    - Implemented many Reed-Solomon decoding algorithms, including 
      Berlekamp-Welch, Cohn-Heninger Single-Polynomial Decoding, 
      Cohn-Heninger Multi-Polynomial Decoding, a dynamic programming 
      approach and a portfolio algorithm of all of the above.  This 
      allows for successful decoding with a higher number of Byzantine 
      servers.
    - Modified command-line usage of testclient, pirclient and pirserver.
    - Improved testclient; testclient now kills pirserver processes after
      the test is completed.
    - Modified pirserver and pirclient to use sockets for communication;
      pirserver processes are now launched separately from pirclient.

Version 0.7.1 (2007-06-17):

    (Based on patches from Len Sassaman <Len.Sassaman@esat.kuleuven.be>)
    Added support for *BSD stat(1) in testclient, and testclient now
    does additional sanity checks and auto-generates the test database
    if it doesn't exist (or isn't readable).  Added the makefile
    argument "distclean" to clean up extraneous files.  Utilities now
    display the current version number when given the argument
    --version.  When recovering from Byzantine servers and HardRecover
    is invoked, a command-line message is displayed.

Version 0.7 (2007-04-03):
    The Guruswami-Sudan implementation has been changed to a much more
    effecient algorithm.  This saves about 70% of the runtime in the
    presence of Byzantine servers.  Set the environment variable
    PIRC_NAIVE=1 to revert to the old algorithm for comparison.

Version 0.6 (2007-03-14):
    Thanks to M. Jason Hinek <mjhinek@alumni.uwaterloo.ca>, the
    dependency on MuPAD has been removed.  All computations are now done
    natively in C++ using NTL.

Version 0.5 (2007-03-02):
    Initial release

Copyright
---------

Percy++ is covered by the following (GPL) license:

    Percy++
    Copyright 2007,2012,2013,2014
    Ian Goldberg <iang@cs.uwaterloo.ca>,
    Casey Devet <cjdevet@uwaterloo.ca>,
    Wouter Lueks <wouter@telox.net>,
    Ann Yang <y242yang@uwaterloo.ca>,
    Paul Hendry <pshendry@uwaterloo.ca>,
    Ryan Henry <rhenry@cs.uwaterloo.ca>

    This program is free software; you can redistribute it and/or modify
    it under the terms of version 2 of the GNU General Public License as
    published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    There is a copy of the GNU General Public License in the COPYING file
    packaged with this plugin; if you cannot find it, write to the Free
    Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
    02110-1301 USA

The ZZ_pXY.cc and ZZ_pXY.h files are adapted from files in version
5.4 of the NTL package by Victor Shoup <victor@shoup.net>, with
modifications by M. Jason Hinek <mjhinek@alumni.uwaterloo.ca> and
Ian Goldberg <iang@cs.uwaterloo.ca>, and are also under the above GPL
version 2 license.  Portions of pirclient.cc and pirserver.cc are by
Femi Olumofin <fgolumof@cs.uwaterloo.ca>, under the same GPL version 2
license.
