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

#include <iostream>
#include <string>

#include "libQnormaliz/Qinteger.h"
#include "libQnormaliz/Qvector_operations.h"
#include "libQnormaliz/Qmatrix.h"
// #include "libQnormaliz/Qsimplex.h"
#include "libQnormaliz/Qlist_operations.h"

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> l_multiplication(const list< vector<Number> >& l,const vector<Number>& v){
    int s=l.size();
    vector<Number> p(s);
    typename list< vector<Number> >::const_iterator i;
    s=0;
    for (i =l.begin(); i != l.end(); ++i, ++s) {
        p[s]=v_scalar_product(*i,v);             //maybe we loose time here?
    }
    return p;
}

//---------------------------------------------------------------------------

template<typename Number>
list< vector<Number> > l_list_x_matrix(const list< vector<Number> >& l,const Matrix<Number>& M){
    list< vector<Number> > result;
    vector<Number> p;
    typename list< vector<Number> >::const_iterator i;
    for (i =l.begin(); i != l.end(); i++) {
        p=M.VxM(*i);
        result.push_back(p);
    }
    return result;
}
//---------------------------------------------------------------------------

template<typename Number>
void  l_cut(list<  vector<Number> >& l, int size){
    typename list< vector<Number> >::iterator i;
    for (i =l.begin(); i != l.end(); i++) {
        (*i).resize(size);
    }
}

//---------------------------------------------------------------------------


template<typename Number>
void  l_cut_front(list<  vector<Number> >& l, int size){
    typename list< vector<Number> >::iterator i;
    vector<Number> tmp;
    for (i =l.begin(); i != l.end(); ) {
        tmp=v_cut_front(*i, size);
        i=l.erase(i);  //important to decrease memory consumption
        l.insert(i,tmp);
    }
}

//---------------------------------------------------------------------------

}
