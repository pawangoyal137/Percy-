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

#ifndef __MPICOMM_H__
#define __MPICOMM_H__

#ifdef MPI_DIST_SERVER
#include <mpi.h>
#endif
#include <streambuf>
#include <iostream>
#include <vector>
#include <stdio.h>
#include "percytypes.h"
#include "percyparams.h"

#ifdef MPI_DIST_SERVER
class MPIStreamBuf : public std::streambuf
{
public:
    MPIStreamBuf (int other_rank, MPI_Comm comm = MPI_COMM_WORLD);
    ~MPIStreamBuf ();

protected:
    virtual int underflow ();
    virtual int overflow (int c = EOF);
    virtual std::streamsize xsputn (const char * s, std::streamsize n);

private:
    char * buffer;
    int buffersize;

    int other_rank;
    int my_rank;
    MPI_Comm comm;
};
#endif

struct BufferInfo {
    BufferInfo () : addr(), size() {}
    BufferInfo (char * addr, dbsize_t size) : addr(addr), size(size) {}
    char * addr;
    dbsize_t size;
};

typedef std::vector<BufferInfo> BufferList;

class MemoryStreamBuf : public std::streambuf {
public:
    MemoryStreamBuf ();
    MemoryStreamBuf (BufferList inbuffers, BufferList outbuffers);
    virtual ~MemoryStreamBuf ();

    void add_inbuffer (char * addr, dbsize_t size);
    void add_outbuffer (char * addr, dbsize_t size);

    virtual std::streamsize showmanyc ();
    virtual int underflow ();
    virtual int uflow ();
    virtual int overflow (int c = EOF);

    bool in_eof ();
    bool out_eof ();

private:
    BufferList inbuffers;
    BufferList outbuffers;
    std::streamsize inbufferindex;
    std::streamsize outbufferindex;
    dbsize_t incharindex;
    dbsize_t outcharindex;
};

#endif
