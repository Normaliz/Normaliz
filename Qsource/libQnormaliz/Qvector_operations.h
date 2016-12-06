/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */
//---------------------------------------------------------------------------
#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H
//---------------------------------------------------------------------------

#include <vector>
#include <ostream>
#include <list>

#include <libQnormaliz/libQnormaliz.h>
#include <libQnormaliz/Qinteger.h>
#include <libQnormaliz/Qconvert.h>

namespace libQnormaliz {
using std::vector;

//---------------------------------------------------------------------------
//							Data access
//---------------------------------------------------------------------------

template <typename T>
std::ostream& operator<< (std::ostream& out, const vector<T>& vec) {
    for (size_t i=0; i<vec.size(); ++i) {
        out << vec[i] << " ";
    }
    out << std::endl;
    return out;
}

//---------------------------------------------------------------------------
//					    	Vector operations
//---------------------------------------------------------------------------
template<typename Number>
Number v_scalar_product(const vector<Number>& a,const vector<Number>& b);

//returns the scalar product of the vector a with the end of the vector b
template<typename Number>
Number v_scalar_product_unequal_vectors_end(const vector<Number>& a,const vector<Number>& b);

//returns the addition a + b, vectors must be of equal size
template<typename Number>
vector<Number> v_add(const vector<Number>& a,const vector<Number>& b);
template<typename Number>
vector<Number> v_add_overflow_check(const vector<Number>& a,const vector<Number>& b);
template<typename Number>
void v_add_result(vector<Number>& result, const size_t length, const vector<Number>& a,const vector<Number>& b);



//---------------------------------------------------------------------------
//							Scalar operations
//---------------------------------------------------------------------------

//v = v * scalar
template<typename Number>
void v_scalar_multiplication(vector<Number>& v, const Number& scalar){
    size_t i,size=v.size();
    for (i = 0; i <size; i++) {
        v[i] *= scalar;
    }
}

template<typename Number>
void v_scalar_division(vector<Number>& v, const Number& scalar);

//---------------------------------------------------------------------------
//							   General vector operations
//---------------------------------------------------------------------------

//returns a new vector with the content of a extended by b
template<typename T>
vector<T> v_merge(const vector<T>& a, const T& b);

//returns a new vector with the content of a and b
template<typename T>
vector<T> v_merge(const vector<T>& a, const vector<T>& b);

//returns a new vector with the last size entries of v
template<typename T>
vector<T> v_cut_front(const vector<T>& v, size_t size);

//the input vectors must be ordered of equal size
//if u is different from v by just one element, it returns that element
//else returns 0 (the elements of u and v are >0)
//int v_difference_ordered_fast(const vector<size_t>& u,const vector<size_t>& v);


template<typename Number>
bool compare_last (const vector<Number>& a, const vector<Number>& b)
{
    return a.back() < b.back();
}

//returns a key vector containing the positions of non-zero entrys of v
template<typename Number>
vector<key_t> v_non_zero_pos(const vector<Number>& v);

// counts the number of positive entries
template<typename Number>
size_t v_nr_positive(const vector<Number>& v);

// check whether the vector only contains 0
template<typename Number>
bool v_is_zero(const vector<Number>& v);

template<typename Number>
bool v_is_symmetric(const vector<Number>& v);

template<typename Number>
bool v_is_nonnegative(const vector<Number>& v);

template<typename Number>
Number v_max_abs(const vector<Number>& v){
    Number tmp = 0;
    for (size_t i=0; i<v.size(); i++){
            if (Iabs(v[i])>tmp) tmp=Iabs(v[i]);
    }
    return tmp;
}

//---------------------------------------------------------------------------
//							   bool vector operations
//---------------------------------------------------------------------------

vector<bool> v_bool_andnot(const vector<bool>& a, const vector<bool>& b);

// swaps entry i and j of the vector<bool> v
void v_bool_entry_swap(vector<bool>& v, size_t i, size_t j);

//---------------------------------------------------------------------------
//							  Special
//---------------------------------------------------------------------------

// computes integral simplex containing a rational vector
template<typename Number>
void approx_simplex(const vector<Number>& q, std::list<vector<Number> >& approx,const long k);

vector<key_t> identity_key(size_t n);

//---------------------------------------------------------------------------
//                            Sorting
//---------------------------------------------------------------------------

template <typename T>
void order_by_perm(vector<T>& v, const vector<key_t>& permfix);

} // namespace

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
