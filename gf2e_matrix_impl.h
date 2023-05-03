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

template <typename GF2E_Element>
inline void PercyServer::Matrix<GF2E_Element>::strassen_submatrices(
	SubMatrix<GF2E_Element> &m11,
	SubMatrix<GF2E_Element> &m12,
	SubMatrix<GF2E_Element> &m21,
	SubMatrix<GF2E_Element> &m22)
{
    dbsize_t nrows = this->num_rows / 2;
    dbsize_t ncols = this->num_cols / 2;

    m11 = SubMatrix<GF2E_Element>(*this,     0,     0, nrows, ncols);
    m12 = SubMatrix<GF2E_Element>(*this,     0, ncols, nrows, ncols);
    m21 = SubMatrix<GF2E_Element>(*this, nrows,     0, nrows, ncols);
    m22 = SubMatrix<GF2E_Element>(*this, nrows, ncols, nrows, ncols);
}

template <typename GF2E_Element>
inline void PercyServer::Matrix<GF2E_Element>::clear()
{
    memset(this->data, 0,
	    this->num_cols * this->num_rows * sizeof(GF2E_Element));
}

template <typename GF2E_Element>
inline void PercyServer::Matrix<GF2E_Element>::is_sum_of(
	const SubMatrix<const GF2E_Element> &sub1,
	const SubMatrix<const GF2E_Element> &sub2)
{
    // WARNING: only works for equal sized submatrices
    GF2E_Element * dst = this->data;
    const GF2E_Element * src1 = sub1.first();
    const GF2E_Element * src2 = sub2.first();

    size_t nc_bytes = sub1.num_cols * sizeof(GF2E_Element);

    for(unsigned int row = 0; row < sub1.num_rows; ++row) {
	XOR_memblocks((unsigned char*) dst,
		(const unsigned char*) src1,
		(const unsigned char*) src2, nc_bytes);

	dst  += sub1.num_cols;
	src1 += sub1.matrix.num_cols;
	src2 += sub2.matrix.num_cols;
    }
}

template <typename GF2E_Element>
inline void PercyServer::Matrix<GF2E_Element>::copy_from(
	const SubMatrix<const GF2E_Element> &other)
{
    const GF2E_Element *src = other.first();
    GF2E_Element *dst = this->data;
    for(unsigned int row = 0; row < other.num_rows; ++row) {
	memcpy((void *) dst, (const void *) src,
		other.num_cols * sizeof(GF2E_Element));
	dst += this->num_cols;
	src += other.matrix.num_cols;
    }
}

template <typename GF2E_Element>
std::ostream &operator<<(std::ostream &os,
	const PercyServer::Matrix<GF2E_Element> &m)
{
    os << "[";
    for(unsigned int r = 0; r < m.num_rows; ++r) {
	for(unsigned int c = 0; c < m.num_cols; ++c) {
	    os << " " << hex <<
		(int) *(m.data + r * m.num_cols + c);
	}
	os << ",";
    }
    os << "]";
    return os;
}

template <typename GF2E_Element>
inline void PercyServer::SubMatrix<GF2E_Element>::add(
	const Matrix<GF2E_Element> &m)
{
    GF2E_Element *dst = this->first();
    GF2E_Element *src = m.data;

    size_t nc_bytes = this->num_cols * sizeof(GF2E_Element);

    for(unsigned int row = 0; row < this->num_rows; ++row) {
	XOR_equal((unsigned char *)dst, (const unsigned char*)src,
		nc_bytes);
	dst += this->matrix.num_cols;
	src += m.num_cols;
    }

}

template <typename GF2E_Element>
std::ostream &operator<<(std::ostream &os,
	const PercyServer::SubMatrix<GF2E_Element> &m)
{
    GF2E_Element *data = m.matrix.data;
    dbsize_t num_cols = m.matrix.num_cols;

    os << "[";
    for(dbsize_t r = m.row; r < m.num_rows + m.row; ++r) {
	for(dbsize_t c = m.col; c < m.num_cols + m.col; ++c) {
	    os << " " << hex << (int) *(data + r * num_cols + c);
	}
	os << ",";
    }
    os << "]" << endl;
    return os;
}

template <typename GF2E_Element>
inline GF2E_Element* PercyServer::SubMatrix<GF2E_Element>::
get(dbsize_t row, dbsize_t col) const
{
    return this->matrix.data +
	this->matrix.num_cols * (this->row + row) + (this->col + col);
}

template <typename GF2E_Element>
inline GF2E_Element*
PercyServer::SubMatrix<GF2E_Element>::first() const
{
    return this->matrix.data +
	this->matrix.num_cols * this->row + this->col;
}

template <typename GF2E_Element>
inline GF2E_Element*
PercyServer::Row<GF2E_Element>::get(dbsize_t col) const
{
    return SubMatrix<GF2E_Element>::get(0, col);
}

template <typename GF2E_Element>
inline GF2E_Element*
PercyServer::Col<GF2E_Element>::get(dbsize_t row) const
{
    return SubMatrix<GF2E_Element>::get(row, 0);
}
template <typename GF2E_Element>
inline GF2E_Element*
PercyServer::Elem<GF2E_Element>::get() const
{
    return SubMatrix<GF2E_Element>::get(0, 0);
}

template <>
inline void PercyServer::Row<GF28_Element>::add_mult_of(
	const Row<const GF28_Element> &row_a,
	const SubMatrix<const GF28_Element> &mat_b)
{
    const GF28_Element *row = row_a.first();
    const GF28_Element *block = mat_b.first();

    GF28_Element *oc_start = this->first();
    GF28_Element *oc_end = this->get(mat_b.num_cols & ~7);

    const GF28_Element *multrow;
    const GF28_Element *blockc;
    GF28_Element *oc;

    for (unsigned int j = 0; j < row_a.num_cols; ++j) {
        multrow = GF28_mult_table[row[j]];
        blockc = block;
	oc = oc_start;

        while (oc < oc_end) {
            uint64_t accum = (uint64_t) multrow[*(blockc++)];
            accum |= (uint64_t) multrow[*(blockc++)] << 8;
            accum |= (uint64_t) multrow[*(blockc++)] << 16;
            accum |= (uint64_t) multrow[*(blockc++)] << 24;
            accum |= (uint64_t) multrow[*(blockc++)] << 32;
            accum |= (uint64_t) multrow[*(blockc++)] << 40;
            accum |= (uint64_t) multrow[*(blockc++)] << 48;
            accum |= (uint64_t) multrow[*(blockc++)] << 56;
            *((uint64_t *) oc) ^= accum;
            oc+=8;
        }
        for (unsigned int c = 0; c < (mat_b.num_cols & 7); ++c) {
            *(oc++) ^= multrow[*(blockc++)];
        }

        block += mat_b.matrix.num_cols;
    }
}

template <>
inline void PercyServer::Row<GF216_Element>::add_mult_of(
	const Row<const GF216_Element> &row_a,
	const SubMatrix<const GF216_Element> &mat_b)
{
    const GF216_Element *row = row_a.first();
    const GF216_Element *block = mat_b.first();

    GF216_Element *oc;
    GF216_Element *oc_start = this->first();
    GF216_Element *oc_end = this->get(mat_b.num_cols & ~3);
    for (dbsize_t j = 0; j < row_a.num_cols; ++j) {
        GF216_Element inpv_j = row[j];
        if (inpv_j != 0) {
            const GF216_Element *blockc = block;
            GF216_Element log_j = GF216_log_table[inpv_j];
            const GF216_Element *start = GF216_exp_table + log_j;
            oc = oc_start;
            GF216_Element block_c;
            while(oc < oc_end) {
                uint64_t accum = 0;
                block_c = *(blockc++);
                if (block_c != 0) {
                    GF216_Element log_c = GF216_log_table[block_c];
                    accum |= (uint64_t) start[log_c];
                }
                block_c = *(blockc++);
                if (block_c != 0) {
                    GF216_Element log_c = GF216_log_table[block_c];
                    accum |= (uint64_t) start[log_c] << 16;
                }
                block_c = *(blockc++);
                if (block_c != 0) {
                    GF216_Element log_c = GF216_log_table[block_c];
                    accum |= (uint64_t) start[log_c] << 32;
                }
                block_c = *(blockc++);
                if (block_c != 0) {
                    GF216_Element log_c = GF216_log_table[block_c];
                    accum |= (uint64_t) start[log_c] << 48;
                }
                *((uint64_t *) oc) ^= accum;
                oc+=4;
            }
            for (dbsize_t c = 0; c < (mat_b.num_cols & 3); ++c, ++oc) {
                block_c = *(blockc++);
                if (block_c != 0) {
                    GF216_Element log_c = GF216_log_table[block_c];
                    *oc ^= start[log_c];
                }
            }
        }
        block += mat_b.matrix.num_cols;
    }
}

template <typename GF2E_Element>
inline void PercyServer::Col<GF2E_Element>::add_mult_of(
	const SubMatrix<const GF2E_Element> &mat_a,
	const Col<const GF2E_Element> &col_b)
{
    dbsize_t r, c;
    for(r = 0; r < mat_a.num_rows; ++r) {
	GF2E_Element *res = get(r);
	for(c = 0; c < mat_a.num_cols; ++c) {
	    *res ^=
		multiply_GF2E(*mat_a.get(r,c), *col_b.get(c));
	}
    }
}

template <typename GF2E_Element>
inline void PercyServer::Elem<GF2E_Element>::add_mult_of(
	const Row<const GF2E_Element> &row_a,
	const Col<const GF2E_Element> &col_b)
{
    GF2E_Element *acc = this->first();

    const GF2E_Element *a_ptr = row_a.first();
    const GF2E_Element *b_ptr = col_b.first();
    const GF2E_Element *a_end = row_a.get(row_a.num_cols);

    while(a_ptr < a_end) {
	*acc ^= multiply_GF2E(*(a_ptr++), *b_ptr);
	b_ptr += col_b.matrix.num_cols;
    }
}

template <>
inline void PercyServer::SubMatrix<GF28_Element>::add_mult_of(
	const Col<const GF28_Element> &col_a,
	const Row<const GF28_Element> &row_b)
{
    const GF28_Element *row_start = row_b.first();

    const GF28_Element *multcol;
    const GF28_Element *row;
    GF28_Element *oc, *oc_end;

    for (dbsize_t r = 0; r < col_a.num_rows; ++r) {
        multcol = GF28_mult_table[*col_a.get(r)];

	oc = this->get(r,0);
	oc_end = this->get(r,row_b.num_cols);
	row = row_start;

	while(oc < oc_end) {
	    *(oc++) ^= multcol[*(row++)];
	}
    }
}

template <>
inline void PercyServer::SubMatrix<GF216_Element>::add_mult_of(
	const Col<const GF216_Element> &col_a,
	const Row<const GF216_Element> &row_b)
{
    const GF216_Element *row;
    const GF216_Element *row_start = row_b.first();
    GF216_Element row_elem;

    GF216_Element col_elem;
    GF216_Element log_c;
    const GF216_Element *start_c;

    GF216_Element *oc, *oc_end;

    for (dbsize_t r = 0; r < col_a.num_rows; ++r) {
	col_elem = *col_a.get(r);
	if(col_elem != 0) {
	    dbsize_t c = 0;
	    oc = this->get(r,0);
	    oc_end = this->get(r,row_b.num_cols);

	    log_c = GF216_log_table[col_elem];
	    start_c = GF216_exp_table + log_c;

	    row = row_start;

	    while(oc < oc_end) {
		row_elem = *(row++);
		if(row_elem != 0) {
		    *oc ^= start_c[GF216_log_table[row_elem]];
		}
		++c;
		++oc;
	    }
	}
    }
}

template <>
inline void PercyServer::Col<GF28_Element>::add_mult_of(
	const Col<const GF28_Element> &col_a,
	const Elem<const GF28_Element> &elem_b)
{
    const GF28_Element *mult = GF28_mult_table[*elem_b.first()];

    const GF28_Element *col = col_a.first();
    GF28_Element *oc, *oc_end;

    oc = this->get(0);
    oc_end = this->get(col_a.num_rows);

    while(oc < oc_end) {
	*oc ^= mult[*col];
	oc  += this->matrix.num_cols;
	col += col_a.matrix.num_cols;
    }
}

template <>
inline void PercyServer::Col<GF216_Element>::add_mult_of(
	const Col<const GF216_Element> &col_a,
	const Elem<const GF216_Element> &elem_b)
{
    GF216_Element b = *elem_b.first();

    if(b != 0) {
	const GF216_Element log_mult = GF216_log_table[b];
	const GF216_Element *start_mult = GF216_exp_table + log_mult;

	const GF216_Element *col = col_a.first();
	GF216_Element *oc, *oc_end;

	oc = this->get(0);
	oc_end = this->get(col_a.num_rows);

	while(oc < oc_end) {
	    if(*col != 0) {
		*oc ^= start_mult[GF216_log_table[*col]];
	    }

	    oc  += this->matrix.num_cols;
	    col += col_a.matrix.num_cols;
	}
    }
}

template <typename GF2E_Element>
inline void PercyServer::Elem<GF2E_Element>::add_mult_of(
	const Elem<const GF2E_Element> &elem_a,
	const Elem<const GF2E_Element> &elem_b)
{
    *this->first() ^= multiply_GF2E(*elem_a.first(), *elem_b.first());
}
