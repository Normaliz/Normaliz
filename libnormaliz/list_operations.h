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
#ifndef LIST_OPERATIONS_H
#define LIST_OPERATIONS_H


//---------------------------------------------------------------------------
                  
#include <vector>
#include <list>

#include "libnormaliz.h"
#include "simplex.h"

namespace libnormaliz {

//---------------------------------------------------------------------------
//							Data access
//---------------------------------------------------------------------------

template<typename Integer>
 int l_read(const list< vector<Integer> >& l);  //used for tests, returns size of l
template<typename Integer>
 int l_read(const list< vector<int> >& l);  //used for tests, returns size of l

template<typename Integer>
int l_read(const list <pair<vector<size_t>,Integer> >& l);  //used for tests, returns size of l
 
template<typename Integer>
 int l_read(const list<int>& l);  //used for tests, returns size of l

//---------------------------------------------------------------------------
//						   List operations
//---------------------------------------------------------------------------

template<typename Integer>
 vector<Integer> l_multiplication(const list< vector<Integer> >& l,const vector<Integer>& v);
 //the list shall contain only vectors of size=v.size(). Returns a vector
 //containing all the scalar products  (we see l as as matrix and return l*v).
template<typename Integer>
 list< vector<Integer> > l_list_x_matrix(const list< vector<Integer> >& l,const Matrix<Integer>& M);
 //the list shall contain only vectors of size=M.nr_of_rows(). Returns a list
 //containing the product  (we see l as as matrix and return l*M).
template<typename Integer>
 void  l_cut(list<  vector<Integer> >& l,int size );
 //cuts all the vectors in l to a given size.
template<typename Integer>
 void  l_cut_front(list<  vector<Integer> >& l,int size );
 //cuts all the vectors in l to a given size, maintaining the back

}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
