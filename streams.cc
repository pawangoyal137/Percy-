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

#include "streams.h"

#ifdef MPI_DIST_SERVER
// MPI Streaming implementation

MPIStreamBuf::MPIStreamBuf (int other_rank, MPI_Comm comm) :
    buffer(NULL),
    buffersize(0),
    other_rank(other_rank),
    comm(comm)
{
    MPI_Comm_rank(comm, &my_rank);
}

MPIStreamBuf::~MPIStreamBuf ()
{
    if (buffer != NULL) {
	delete[] buffer;
    }
}

int MPIStreamBuf::underflow ()
{
    if (gptr() < egptr()) {
	// The buffer is not all read
	return traits_type::to_int_type(*gptr());	
    }

    // We need to do another MPI receive
    // Get length of next receive.
    MPI_Status s1;
    MPI_Recv(&buffersize, 1, MPI_INT, other_rank, 0, comm, &s1);

    // Set up buffer
    if (buffer != NULL) {
	delete[] buffer;
    }
    buffer = new char[buffersize];

    // Do receive
    MPI_Status s2;
    MPI_Recv(buffer, buffersize, MPI_CHAR, other_rank, 0, comm, &s2);

#ifdef VERBOSE_MPI
    fprintf(stderr, "[MPI] %d <- %d - Received %d character(s)\n", my_rank, other_rank, buffersize);
#endif

    // Change pointers
    setg(buffer, buffer, buffer + buffersize);

    // Return current character
    return traits_type::to_int_type(*gptr());
}

int MPIStreamBuf::overflow (int c)
{
    if (c != EOF) {
	// Send the character
	char cc = traits_type::to_char_type(c);
	if (xsputn(&cc, 1) < 1) {
	    return EOF;
	}
    }
    return c;
}

std::streamsize MPIStreamBuf::xsputn (const char * s, std::streamsize n)
{
    // Send the size
    int size = (int)n;
    MPI_Send(&size, 1, MPI_INT, other_rank, 0, comm);

    // Send the string
    MPI_Send((char *)s, size, MPI_CHAR, other_rank, 0, comm);

#ifdef VERBOSE_MPI
    fprintf(stderr, "[MPI] %d -> %d - Sent %d character(s)\n", my_rank, other_rank, size);
#endif

    return size;
}
#endif


MemoryStreamBuf::MemoryStreamBuf () :
    inbuffers(),
    outbuffers(),
    inbufferindex(0),
    outbufferindex(0),
    incharindex(0),
    outcharindex(0)
{}

MemoryStreamBuf::MemoryStreamBuf (BufferList inbuffers, BufferList outbuffers)
:
    inbuffers(inbuffers),
    outbuffers(outbuffers),
    inbufferindex(0),
    outbufferindex(0),
    incharindex(0),
    outcharindex(0)
{}

MemoryStreamBuf::~MemoryStreamBuf () {}

void MemoryStreamBuf::add_inbuffer (char * addr, dbsize_t size)
{
    BufferInfo bi;
    bi.addr = addr;
    bi.size = size;
    inbuffers.push_back(bi);
}

void MemoryStreamBuf::add_outbuffer (char * addr, dbsize_t size)
{
    BufferInfo bi;
    bi.addr = addr;
    bi.size = size;
    outbuffers.push_back(bi);
}

std::streamsize MemoryStreamBuf::showmanyc ()
{
    std::streamsize total = 0;
    BufferList::iterator iter;
    for (iter = inbuffers.begin(); iter != inbuffers.end(); ++iter) {
	total += ( iter->size < 0 ? 0 : iter->size );
    }
    return total;
}

int MemoryStreamBuf::underflow ()
{
    if (in_eof()) {
	return EOF;
    }
    char * currbuf = inbuffers[inbufferindex].addr;
    return traits_type::to_int_type(currbuf[incharindex]);
}

int MemoryStreamBuf::uflow ()
{
    int retval = underflow();
    if (retval == EOF) {
	return EOF;
    }
    ++incharindex;
    if (incharindex >= inbuffers[inbufferindex].size) {
	++inbufferindex;
	incharindex = 0;
    }
    return retval;
}

int MemoryStreamBuf::overflow (int c)
{
    if (out_eof()) {
	return EOF;
    }
    char * currbuf = outbuffers[outbufferindex].addr;
    currbuf[outcharindex] = traits_type::to_char_type(c);
    ++outcharindex;
    if (outcharindex >= outbuffers[outbufferindex].size) {
	++outbufferindex;
	outcharindex = 0;
    }
    return c;
}

bool MemoryStreamBuf::in_eof ()
{
    return inbufferindex >= (std::streamsize)(inbuffers.size());
}

bool MemoryStreamBuf::out_eof ()
{
    return outbufferindex >= (std::streamsize)(outbuffers.size());
}

