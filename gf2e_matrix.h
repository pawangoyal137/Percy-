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

#ifndef __GF2E_MATRIX_H__
#define __GF2E_MATRIX_H__

#include "gf2e.h"
#include <iostream>
#include <stdint.h>
#include <string.h>
#include "xor.h"
#include "itserver.h"

template <typename GF2E_Element>
class PercyServer::Matrix {
public:
    GF2E_Element *data;
    dbsize_t num_rows;
    dbsize_t num_cols;

    // Constructors
    Matrix () {}
    Matrix(GF2E_Element *data, dbsize_t num_rows, dbsize_t num_cols) :
	data(data),
	num_rows(num_rows),
	num_cols(num_cols) {}

    void strassen_submatrices(SubMatrix<GF2E_Element> &m11,
	    SubMatrix<GF2E_Element> &m12,
	    SubMatrix<GF2E_Element> &m21,
	    SubMatrix<GF2E_Element> &m22);

    void clear();

    void is_sum_of(const SubMatrix<const GF2E_Element> &m1,
	    const SubMatrix<const GF2E_Element> &m2);
    void copy_from(const SubMatrix<const GF2E_Element> &src);

    template <typename GF2E_Element_>
    friend std::ostream &operator<<(std::ostream &os,
	    const Matrix<GF2E_Element_> &m);
};

template <typename GF2E_Element>
std::ostream &operator<<(std::ostream &os,
	const PercyServer::Matrix<GF2E_Element> &m);

template <typename GF2E_Element>
class PercyServer::SubMatrix {
public:
    Matrix<GF2E_Element> matrix;
    dbsize_t row;
    dbsize_t col;
    dbsize_t num_rows;
    dbsize_t num_cols;

    // Constructors
    SubMatrix () {}
    SubMatrix (Matrix<GF2E_Element> matrix, dbsize_t row, dbsize_t col,
	    dbsize_t num_rows, dbsize_t num_cols) :
	matrix(matrix),
	row(row),
	col(col),
	num_rows(num_rows),
	num_cols(num_cols) {}

    GF2E_Element* get(dbsize_t row, dbsize_t col) const;
    GF2E_Element* first() const;

    void add(const Matrix<GF2E_Element> &m);

    inline void add_mult_of(
	    const Col<const GF2E_Element> &col_a,
	    const Row<const GF2E_Element> &row_b);

    template <typename GF2E_Element_>
    friend std::ostream &operator<<(std::ostream &os,
	    const SubMatrix<GF2E_Element_> &m);
};

template <typename GF2E_Element>
std::ostream &operator<<(std::ostream &os,
	const PercyServer::SubMatrix<GF2E_Element> &m);


template <typename GF2E_Element>
class PercyServer::Row : public virtual SubMatrix<GF2E_Element> {
public:
    // Constructor
    Row (Matrix<GF2E_Element> matrix,
	    dbsize_t row, dbsize_t col, dbsize_t num_cols) :
	SubMatrix<GF2E_Element>(matrix, row, col, 1, num_cols) {}

    GF2E_Element* get(dbsize_t row) const;

    inline void add_mult_of(
	    const Row<const GF2E_Element> &row_a,
	    const SubMatrix<const GF2E_Element> &mat_b);
};

template <typename GF2E_Element>
class PercyServer::Col : public virtual SubMatrix<GF2E_Element> {
public:
    // Constructor
    Col (Matrix<GF2E_Element> matrix,
	    dbsize_t row, dbsize_t col, dbsize_t num_rows) :
	SubMatrix<GF2E_Element>(matrix, row, col, num_rows, 1) {}

    GF2E_Element* get(dbsize_t col) const;

    inline void add_mult_of(
	const SubMatrix<const GF2E_Element> &mat_a,
	const Col<const GF2E_Element> &col_b);

    inline void add_mult_of(
	const Col<const GF2E_Element> &col_a,
	const Elem<const GF2E_Element> &elem_b);
};

template <typename GF2E_Element>
class PercyServer::Elem : public Row<GF2E_Element>, public Col<GF2E_Element> {
public:
    // Constructor
    Elem (Matrix<GF2E_Element> matrix,
	    dbsize_t row, dbsize_t col) :
	SubMatrix<GF2E_Element>(matrix, row, col, 1, 1),
	Row<GF2E_Element>(matrix, row, col, 1),
	Col<GF2E_Element>(matrix, row, col, 1) {}

    GF2E_Element* get() const;

    inline void add_mult_of(
	const Row<const GF2E_Element> &row_a,
	const Col<const GF2E_Element> &col_b);

    inline void add_mult_of(
	const Elem<const GF2E_Element> &elem_a,
	const Elem<const GF2E_Element> &elem_b);
};

#include "gf2e_matrix_impl.h"

#endif
