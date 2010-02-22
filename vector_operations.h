/*
 * Normaliz 2.2
 * Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim
 * With contributions by Christof Soeger
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

#include <stdlib.h>
#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;


//---------------------------------------------------------------------------
//							Data acces
//---------------------------------------------------------------------------

void v_write(vector<Integer>& v);        //used for tests
int v_read(const vector<Integer>& v);  //used for tests, returns size of v
int v_read(const vector<int>& v);  //used for tests, returns size of v
int v_read(const vector<bool>& v);  //used for tests, returns size of v

//---------------------------------------------------------------------------
//					    	Vector operations
//---------------------------------------------------------------------------

Integer v_scalar_product(const vector<Integer>& a,const vector<Integer>& b);
Integer v_scalar_product_unequal_vectors_end(const vector<Integer>& a,const vector<Integer>& b);
//returns the scalar product of the vector a with the end of the vector b
vector<Integer> v_merge(const vector<Integer>& a,const vector<Integer>& b);
//returns a vector with the content of a and b
vector<Integer> v_add(const vector<Integer>& a,const vector<Integer>& b);
//returns the addition a + b, vectors must be of equal size
vector<int> v_complement(const int& a, const vector<int>& v);
//returns a vector containg all the elements of v, less v[a]

//returns a new vector with the last size entries of v
vector<Integer> v_cut_front(const vector<Integer>& v, int size);

vector<int> v_non_zero_pos(vector<Integer> v); //returns a key vector containing the positions of non-zero entrys of v (counting from 1 to v.size())

//---------------------------------------------------------------------------
//							abs, gcd and lcm
//---------------------------------------------------------------------------

vector<Integer> v_abs(const vector<Integer>& v); //returns a vector with the absolute
//value of the elements
Integer v_gcd(const vector<Integer>& v); //returns gcd of the elements of v
Integer v_lcm(const vector<Integer>& v); //returns lcm of the elements of v
vector<Integer> v_make_prime(const vector<Integer>& v);
//returns v divided by the gcd of elements
vector<Integer> v_make_prime(const vector<Integer>& v,Integer& g);
//returns v divided by the gcd of elements, g saves the  gcd

//---------------------------------------------------------------------------
//							Scalar operations
//---------------------------------------------------------------------------

void v_scalar_multiplication(vector<Integer>& v, const Integer& scalar);
//v = v * scalar
vector<Integer> v_scalar_multiplication_two(const vector<Integer>& v, const Integer& scalar);
//returns v * scalar
void v_scalar_division(vector<Integer>& v, const Integer& scalar);
//v = v / scalar, all the elements of v must be divisible with the scalar
void v_reduction_modulo(vector<Integer>& v, const Integer& modulo);
//v = v mod modulo

//---------------------------------------------------------------------------
//								Comparation
//---------------------------------------------------------------------------

int v_diference_ordered_fast(const vector<int>& u,const vector<int>& v);
//the input vectors must be ordered of equal size
//if u is different from v by just one element, it returns that element
//else returns 0 (the elements of u and v are >0)

//---------------------------------------------------------------------------
//								Test
//---------------------------------------------------------------------------

bool v_test_scalar_product(const vector<Integer>& a,const vector<Integer>& b, const Integer& result, const int& m);
// test the main computation for arithmetic overflow
// uses multiplication mod m

//---------------------------------------------------------------------------
//							   Error msg
//---------------------------------------------------------------------------

void v_error(string s);

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
