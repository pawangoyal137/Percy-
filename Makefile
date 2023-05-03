## Percy++ Copyright 2007,2012,2013,2014
## Ian Goldberg <iang@cs.uwaterloo.ca>,
## Casey Devet <cjdevet@uwaterloo.ca>,
## Wouter Lueks <wouter@telox.net>,
## Ann Yang <y242yang@uwaterloo.ca>,
## Paul Hendry <pshendry@uwaterloo.ca>,
## Ryan Henry <rhenry@cs.uwaterloo.ca>
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of version 2 of the GNU General Public License as
## published by the Free Software Foundation.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## There is a copy of the GNU General Public License in the COPYING file
## packaged with this plugin; if you cannot find it, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301 USA

# Set to `true` if the build should support SPIR
SPIR_SUPPORT=false

# Set to `true` if building on a 32-bit machine
BIT32_SUPPORT=false

CXXFLAGS=-Wall -g -O2 -pedantic -I/usr/local/include/NTL #-DVERBOSE -DVERBOSE_DISTRIBUTED -DVERBOSE_RECURSIVE -DVERBOSE_THREADED
LDLIBS=-lntl -lgmp -pthread -lsocket++ -L/usr/local/lib -lgcrypt

TESTS=findpolys_test rr_test time_findpolys time_findpolys_gf28 time_findpolys_w8 time_findpolys_w16 time_findpolys_w32 testdistserver hybridtest agtest

TARGETS=pirserver agserver hybridserver pirclient agclient hybridclient splitdatabase

RUNFILES=database database.* out.client out.real

CLIENT_O=percyclient.o itclient.o percyparams.o itparams.o recover.o \
	rsdecoder.o percyio.o FXY.o gf2e.o subset_iter.o subset.o portfolio.o \
	agparams.o agclient.o percystats.o streams.o recursiveparams.o \
	recursiveclient.o hybridparams.o nullclient.o
SERVER_O=percyserver.o itserver.o percyparams.o itparams.o datastore.o \
	percyio.o gf2e.o agparams.o agserver.o percystats.o \
	streams.o recursiveparams.o recursiveserver.o \
	hybridparams.o
SRCS=$(subst .o,.cc,$(CLIENT_O) $(SERVER_O) pirclient.o pirserver.o pirclient_ag.o pirserver_ag.o splitdatabase.o percyio.o)
LIBS=libpercyserver.a libpercyclient.a

VERSION=$(shell ./updateheaders.py --version-number)

ifeq ($(BIT32_SUPPORT),true)
	CXXFLAGS += -DNEED_UINT128
endif

ifeq ($(SPIR_SUPPORT),true)
	CXXFLAGS+= -I../PolyCommit -I../PBCWrapper -DSPIR_SUPPORT
	LDLIBS+= -L../PBCWrapper -lPBC -lpbc
	CLIENT_O+= ../PolyCommit/PolyCommitCommon.o pspir_crypt.o spirclient.o
	SERVER_O+= ../PolyCommit/PolyCommitCommon.o pspir_crypt.o spirserver.o
endif

all: $(TARGETS)

libs: $(LIBS)

libpercyserver.a: $(SERVER_O)
	ar rcs $@ $^

libpercyclient.a: $(CLIENT_O)
	ar rcs $@ $^

tests: $(TESTS)

pirserver: pirserver.o libpercyserver.a
	g++ -o $@ $^ $(LDLIBS)

pirserver_master: pirserver.cc libpercyserver.a
	g++ $(CXXFLAGS) -DDIST_MASTER -o $@ $^ $(LDLIBS)

pirserver_mpi: pirserver.cc mpiserver.o libpercyserver.a
	mpic++ $(CXXFLAGS) -DMPI_DIST_SERVER -o $@ $^ $(LDLIBS)

pirclient: pirclient.o libpercyclient.a
	g++ -o $@ $^ $(LDLIBS)

agclient: pirclient_ag.o libpercyclient.a
	g++ -o $@ $^ $(LDLIBS)

agserver: pirserver_ag.o libpercyserver.a
	g++ -o $@ $^ $(LDLIBS)

hybridclient: pirclient_hybrid.o libpercyclient.a
	g++ -o $@ $^ $(LDLIBS)

hybridserver: pirserver_hybrid.o libpercyserver.a
	g++ -o $@ $^ $(LDLIBS)

splitdatabase: splitdatabase.o percyio.o gf2e.o
	g++ -o $@ $^ $(LDLIBS)

time_findpolys: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -o $@ $^ $(LDLIBS) # w128

time_findpolys_gf28: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ -static $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_GF28 -o $@ $^ $(LDLIBS)

time_findpolys_gf24: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ -static $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_GF24 -o $@ $^ $(LDLIBS)

time_findpolys_w8: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_W8 -o $@ $^ $(LDLIBS)

time_findpolys_w16: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_W16 -o $@ $^ $(LDLIBS)

time_findpolys_w32: rsdecoder.cc FXY.o subset.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTIME_FINDPOLYS -DUSE_W32 -o $@ $^ $(LDLIBS)

rr_test: rsdecoder.cc FXY.o gf2e.o
	g++ $(CXXFLAGS) -DTEST_RR -o $@ $^ $(LDLIBS)

findpolys_test: rsdecoder.cc FXY.o gf2e.o subset_iter.o portfolio.o
	g++ $(CXXFLAGS) -DTEST_FINDPOLYS -o $@ $^ $(LDLIBS)

testdistserver: testdistserver.cc cmdtools.o
	g++ -o $@ $^

agtest: agtest.o libpercyclient.a libpercyserver.a
	g++ -o $@ $^ $(LDLIBS)

hybridtest: hybridtest.o libpercyclient.a libpercyserver.a
	g++ -o $@ $^ $(LDLIBS)

gf2etest: gf2etest.cc gf2e.o
	g++ $(CXXFLAGS) -o $@ $^ $(LDLIBS)

.PHONY: docs

docs: 
	doxygen docs/Doxyfile
	@docs/fixhtml

clean:
	-rm -f *.o

cleandocs:
	-rm -rf docs/html docs/latex

veryclean: clean
	-rm -f $(TARGETS) $(TESTS) $(LIBS)

distclean: veryclean cleandocs
	-rm -f $(RUNFILES)

dist: 
	#@./updateheaders.py
	tar -czf percy++-$(VERSION).tar.gz --transform 's,^,percy++-$(VERSION)/,' -T releaselist

depend:
	makedepend -Y -- $(CXXFLAGS) -- $(SRCS) 2>/dev/null

# DO NOT DELETE

percyclient.o: percyclient.h percyresult.h percytypes.h percyparams.h
percyclient.o: version.h percystats.h percyio.h /usr/local/include/NTL/ZZ.h
percyclient.o: nullclient.h agclient.h agparams.h recursiveparams.h
percyclient.o: itclient.h /usr/local/include/NTL/vec_vec_ZZ_p.h
percyclient.o: /usr/local/include/NTL/vec_GF2E.h gf2e.h itparams.h
percyclient.o: rsdecoder.h FXY.h portfolio.h recover.h rsdecoder_impl.h
percyclient.o: subset.h subset_iter.h itclient_impl.h recursiveclient.h
itclient.o: /usr/local/include/NTL/vec_vec_ZZ_p.h
itclient.o: /usr/local/include/NTL/ZZ_pX.h itclient.h
itclient.o: /usr/local/include/NTL/vec_GF2E.h percyclient.h percyresult.h
itclient.o: percytypes.h percyparams.h version.h percystats.h percyio.h
itclient.o: /usr/local/include/NTL/ZZ.h gf2e.h itparams.h rsdecoder.h FXY.h
itclient.o: portfolio.h recover.h rsdecoder_impl.h subset.h subset_iter.h
itclient.o: itclient_impl.h
percyparams.o: percyparams.h version.h percytypes.h percyio.h
percyparams.o: /usr/local/include/NTL/ZZ.h
itparams.o: itparams.h percytypes.h percyparams.h version.h percyio.h
itparams.o: /usr/local/include/NTL/ZZ.h
recover.o: subset.h percytypes.h subset_iter.h recover.h rsdecoder.h
recover.o: percyresult.h FXY.h portfolio.h rsdecoder_impl.h gf2e.h
rsdecoder.o: rsdecoder.h percyresult.h percytypes.h FXY.h portfolio.h
rsdecoder.o: recover.h gf2e.h rsdecoder_impl.h subset.h subset_iter.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
FXY.o: FXY.h
gf2e.o: gf2e.h
subset_iter.o: subset_iter.h percytypes.h
subset.o: subset.h percytypes.h
portfolio.o: portfolio.h
agparams.o: agparams.h percyparams.h version.h percytypes.h recursiveparams.h
agparams.o: percyio.h /usr/local/include/NTL/ZZ.h
agclient.o: percytypes.h agclient.h percyclient.h percyresult.h percyparams.h
agclient.o: version.h percystats.h percyio.h /usr/local/include/NTL/ZZ.h
agclient.o: agparams.h recursiveparams.h
percystats.o: percystats.h percytypes.h percyparams.h version.h
streams.o: streams.h percytypes.h percyparams.h version.h
recursiveparams.o: recursiveparams.h percyparams.h version.h percytypes.h
recursiveparams.o: percyio.h /usr/local/include/NTL/ZZ.h
recursiveclient.o: recursiveclient.h percyclient.h percyresult.h percytypes.h
recursiveclient.o: percyparams.h version.h percystats.h percyio.h
recursiveclient.o: /usr/local/include/NTL/ZZ.h recursiveparams.h streams.h
hybridparams.o: hybridparams.h recursiveparams.h percyparams.h version.h
hybridparams.o: percytypes.h itparams.h agparams.h percyio.h
hybridparams.o: /usr/local/include/NTL/ZZ.h
nullclient.o: percytypes.h nullclient.h percyclient.h percyresult.h
nullclient.o: percyparams.h version.h percystats.h percyio.h
nullclient.o: /usr/local/include/NTL/ZZ.h
percyserver.o: percytypes.h percyio.h /usr/local/include/NTL/ZZ.h
percyserver.o: percyserver.h datastore.h percyparams.h version.h percystats.h
percyserver.o: streams.h gf2e_matrix.h gf2e.h xor.h itserver.h
percyserver.o: /usr/local/include/NTL/vec_ZZ_p.h itparams.h itserver_impl.h
percyserver.o: gf2e_matrix_impl.h agserver.h agparams.h recursiveparams.h
percyserver.o: recursiveserver.h
itserver.o: /usr/local/include/NTL/vec_ZZ_p.h itserver.h percyserver.h
itserver.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h version.h
itserver.o: percytypes.h percystats.h streams.h gf2e_matrix.h gf2e.h xor.h
itserver.o: gf2e_matrix_impl.h itparams.h itserver_impl.h
percyparams.o: percyparams.h version.h percytypes.h percyio.h
percyparams.o: /usr/local/include/NTL/ZZ.h
itparams.o: itparams.h percytypes.h percyparams.h version.h percyio.h
itparams.o: /usr/local/include/NTL/ZZ.h
datastore.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h version.h
datastore.o: percytypes.h percyio.h streams.h itparams.h recursiveparams.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
gf2e.o: gf2e.h
agparams.o: agparams.h percyparams.h version.h percytypes.h recursiveparams.h
agparams.o: percyio.h /usr/local/include/NTL/ZZ.h
agserver.o: agserver.h datastore.h /usr/local/include/NTL/ZZ.h percyparams.h
agserver.o: version.h percytypes.h agparams.h recursiveparams.h percyserver.h
agserver.o: percystats.h streams.h gf2e_matrix.h gf2e.h xor.h itserver.h
agserver.o: /usr/local/include/NTL/vec_ZZ_p.h itparams.h itserver_impl.h
agserver.o: gf2e_matrix_impl.h percyio.h
percystats.o: percystats.h percytypes.h percyparams.h version.h
streams.o: streams.h percytypes.h percyparams.h version.h
recursiveparams.o: recursiveparams.h percyparams.h version.h percytypes.h
recursiveparams.o: percyio.h /usr/local/include/NTL/ZZ.h
recursiveserver.o: recursiveserver.h datastore.h /usr/local/include/NTL/ZZ.h
recursiveserver.o: percyparams.h version.h percytypes.h recursiveparams.h
recursiveserver.o: percyserver.h percystats.h streams.h gf2e_matrix.h gf2e.h
recursiveserver.o: xor.h itserver.h /usr/local/include/NTL/vec_ZZ_p.h
recursiveserver.o: itparams.h itserver_impl.h gf2e_matrix_impl.h percyio.h
hybridparams.o: hybridparams.h recursiveparams.h percyparams.h version.h
hybridparams.o: percytypes.h itparams.h agparams.h percyio.h
hybridparams.o: /usr/local/include/NTL/ZZ.h
pirclient.o: /usr/local/include/NTL/ZZ_p.h itclient.h
pirclient.o: /usr/local/include/NTL/vec_vec_ZZ_p.h
pirclient.o: /usr/local/include/NTL/vec_GF2E.h percyclient.h percyresult.h
pirclient.o: percytypes.h percyparams.h version.h percystats.h percyio.h
pirclient.o: /usr/local/include/NTL/ZZ.h gf2e.h itparams.h rsdecoder.h FXY.h
pirclient.o: portfolio.h recover.h rsdecoder_impl.h subset.h subset_iter.h
pirclient.o: itclient_impl.h config.h
pirserver.o: datastore.h /usr/local/include/NTL/ZZ.h percyparams.h version.h
pirserver.o: percytypes.h itserver.h /usr/local/include/NTL/vec_ZZ_p.h
pirserver.o: percyserver.h percystats.h streams.h gf2e_matrix.h gf2e.h xor.h
pirserver.o: gf2e_matrix_impl.h itparams.h itserver_impl.h config.h
pirserver.o: agserver.h agparams.h recursiveparams.h
pirclient_ag.o: percytypes.h percyparams.h version.h agparams.h
pirclient_ag.o: recursiveparams.h agclient.h percyclient.h percyresult.h
pirclient_ag.o: percystats.h percyio.h /usr/local/include/NTL/ZZ.h config.h
pirserver_ag.o: percytypes.h percyparams.h version.h agparams.h
pirserver_ag.o: recursiveparams.h agserver.h datastore.h
pirserver_ag.o: /usr/local/include/NTL/ZZ.h percyserver.h percystats.h
pirserver_ag.o: streams.h gf2e_matrix.h gf2e.h xor.h itserver.h
pirserver_ag.o: /usr/local/include/NTL/vec_ZZ_p.h itparams.h itserver_impl.h
pirserver_ag.o: gf2e_matrix_impl.h config.h percyio.h
splitdatabase.o: percyio.h /usr/local/include/NTL/ZZ.h percytypes.h config.h
splitdatabase.o: version.h gf2e.h
percyio.o: percyio.h /usr/local/include/NTL/ZZ.h
