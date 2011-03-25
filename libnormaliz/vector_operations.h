/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 */
//---------------------------------------------------------------------------
#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H
//---------------------------------------------------------------------------

#include <vector>

#include "libnormaliz.h"

namespace libnormaliz {
using std::vector;

//---------------------------------------------------------------------------
//							Data access
//---------------------------------------------------------------------------

template <typename T> void v_write(vector<T>& v);        //used for tests
template <typename T> size_t v_read(const vector<T>& v,std::ostream& out=std::cout);  //used for tests, returns size of v

//---------------------------------------------------------------------------
//					    	Vector operations
//---------------------------------------------------------------------------
template<typename Integer>
Integer v_scalar_product(const vector<Integer>& a,const vector<Integer>& b);

//returns the scalar product of the vector a with the end of the vector b
template<typename Integer>
Integer v_scalar_product_unequal_vectors_end(const vector<Integer>& a,const vector<Integer>& b);

//returns the addition a + b, vectors must be of equal size
template<typename Integer>
vector<Integer> v_add(const vector<Integer>& a,const vector<Integer>& b);

//---------------------------------------------------------------------------
//							abs, gcd and lcm
//---------------------------------------------------------------------------

//returns a vector with the absolute value of the elements
template<typename Integer>
vector<Integer> v_abs(const vector<Integer>& v);

template<typename Integer>
Integer v_gcd(const vector<Integer>& v); //returns gcd of the elements of v

template<typename Integer>
Integer v_lcm(const vector<Integer>& v); //returns lcm of the elements of v

template<typename Integer>
vector<Integer> v_make_prime(const vector<Integer>& v); //returns v divided by the gcd of elements

template<typename Integer>
vector<Integer> v_make_prime(const vector<Integer>& v,Integer& g); //returns v divided by the gcd of elements, g saves the gcd


//---------------------------------------------------------------------------
//							Scalar operations
//---------------------------------------------------------------------------

template<typename Integer>
void v_scalar_multiplication(vector<Integer>& v, const Integer& scalar); //v = v * scalar

template<typename Integer>
vector<Integer> v_scalar_multiplication_two(const vector<Integer>& v, const Integer& scalar);
//returns v * scalar

template<typename Integer>
void v_scalar_division(vector<Integer>& v, const Integer& scalar);
//v = v / scalar, all the elements of v must be divisible with the scalar

template<typename Integer>
void v_reduction_modulo(vector<Integer>& v, const Integer& modulo);
//v = v mod modulo

//---------------------------------------------------------------------------
//								Test
//---------------------------------------------------------------------------

template<typename Integer>
bool v_test_scalar_product(const vector<Integer>& a,const vector<Integer>& b, const Integer& result, const long& m);
// test the main computation for arithmetic overflow
// uses multiplication mod m

//---------------------------------------------------------------------------
//							   General vector operations
//---------------------------------------------------------------------------

//returns a new vector with the content of a and b
template<typename T>
vector<T> v_merge(const vector<T>& a,const vector<T>& b);

//returns a new vector with the last size entries of v
template<typename T>
vector<T> v_cut_front(const vector<T>& v, size_t size);

//the input vectors must be ordered of equal size
//if u is different from v by just one element, it returns that element
//else returns 0 (the elements of u and v are >0)
int v_difference_ordered_fast(const vector<size_t>& u,const vector<size_t>& v);

template<typename Integer>
vector<size_t> v_non_zero_pos(vector<Integer> v); //returns a key vector containing the positions of non-zero entrys of v (counting from 1 to v.size())


}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
