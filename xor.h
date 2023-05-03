// Percy++ Copyright 2007,2012,2013,2014
// Ian Goldberg <iang@cs.uwaterloo.ca>,
// Casey Devet <cjdevet@uwaterloo.ca>,
// Wouter Lueks <wouter@telox.net>
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

#ifndef __XOR_H__
#define __XOR_H__

typedef unsigned long chunk;

inline void XOR_equal(unsigned char *dst, const unsigned char *src, size_t len) {
    chunk *dl = (chunk *)dst;
    const chunk *sl = (const chunk *)src;
    size_t ll = len / (8*sizeof(chunk));
    size_t lr = len % (8*sizeof(chunk));
    for(size_t i=0; i<ll; i++) {
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
        *(dl++) ^= *(sl++);
    }
    if (lr) {
        dst = (unsigned char *)dl;
        src = (const unsigned char *)sl;
        for(size_t i=0; i<lr; i++) {
            *(dst++) ^= *(src++);
        }
    }
}

inline void XOR_memblocks(unsigned char *dst, const unsigned char *src1,
	const unsigned char *src2, size_t len) {
    chunk *dl = (chunk *)dst;
    const chunk *sl1 = (const chunk *)src1;
    const chunk *sl2 = (const chunk *)src2;
    size_t ll = len / (8*sizeof(chunk));
    size_t lr = len % (8*sizeof(chunk));
    for(size_t i=0; i<ll; i++) {
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
        *(dl++) = *(sl1++) ^ *(sl2++);
    }
    if (lr) {
        dst = (unsigned char *)dl;
        src1 = (const unsigned char *)sl1;
        src2 = (const unsigned char *)sl2;
        for(size_t i=0; i<lr; i++) {
            *(dst++) = *(src1++) ^ *(src2++);
        }
    }
}

#endif

