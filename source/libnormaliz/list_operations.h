/*
 * Normaliz
 * Copyright (C) 2007-2019  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
#ifndef LIBNORMALIZ_LIST_OPERATIONS_H
#define LIBNORMALIZ_LIST_OPERATIONS_H

//---------------------------------------------------------------------------

#include <vector>
#include <list>
#include <ostream>

#include "libnormaliz/general.h"

namespace libnormaliz {
using std::list;
using std::vector;

//---------------------------------------------------------------------------
//                          Data access
//---------------------------------------------------------------------------

template <typename T>
std::ostream& operator<<(std::ostream& out, const list<T>& l) {
    for (const auto& i : l) {
        out << i << " ";
    }
    out << std::endl;
    return out;
}

//---------------------------------------------------------------------------
//                         List operations
//---------------------------------------------------------------------------

template <typename Integer>
vector<Integer> l_multiplication(const list<vector<Integer> >& l, const vector<Integer>& v);
// the list shall contain only vectors of size=v.size(). Returns a vector
// containing all the scalar products  (we see l as as matrix and return l*v).
template <typename Integer>
list<vector<Integer> > l_list_x_matrix(const list<vector<Integer> >& l, const Matrix<Integer>& M);
// the list shall contain only vectors of size=M.nr_of_rows(). Returns a list
// containing the product  (we see l as as matrix and return l*M).
template <typename Integer>
void l_cut(list<vector<Integer> >& l, int size);
// cuts all the vectors in l to a given size.

/*
template <typename Integer>
void l_cut_front(list<vector<Integer> >& l, int size);
// cuts all the vectors in l to a given size, maintaining the back
*/

//---------------------------------------------------------------------------

template <typename T>
void random_order(list<T>& LL) {
    vector<typename list<T>::iterator> list_order;
    size_t nrLL = LL.size();
    list_order.reserve(nrLL);
    auto p = LL.begin();
    for (size_t k = 0; k < nrLL; ++k, ++p) {
        list_order.push_back(p);
    }
    for (size_t k = 0; k < 10 * nrLL; ++k) {
        swap(list_order[rand() % nrLL], list_order[rand() % nrLL]);
    }
    list<T> new_order;
    for (size_t k = 0; k < nrLL; ++k) {
        new_order.push_back(*list_order[k]);
    }
    LL.clear();
    LL.splice(LL.begin(), new_order);
}

//---------------------------------------------------------------------------

template <typename T>
void random_order(list<T>& LL, typename list<T>::iterator from, typename list<T>::iterator to) {
    list<T> MM;
    MM.splice(MM.begin(), LL, from, to);
    random_order(MM);
    LL.splice(LL.begin(), MM);
}

//--------------------------------------------------------------------------
// remove all entries that appear exactly twice (or with even multiplicity)
// it must be possible to sorten the list

template <typename T>
void remove_twins(list<T>& L){
    
    L.sort();
    auto S = L.begin(); // remove all subfacets that appear twice
    for(; S != L.end();){
        auto del = S;
        ++del;
        if(del != L.end() && *S == *del){
            S = L.erase(S);
            S = L.erase(S);
        }
        else
            S++;                
    }
}

//--------------------------------------------------------------------------
// remove all entries whose "first" appears twice (or with even multiplicity)
// it must be possible to sorten the list
// L must be a list of pairs

template <typename T>
void remove_twins_in_first(list<T>& L){
    
    L.sort();
    auto S = L.begin(); // remove all subfacets that appear twice
    for(; S != L.end();){
        auto del = S;
        ++del;
        if(del != L.end() && S->first == del->first){
            S = L.erase(S);
            S = L.erase(S);
        }
        else
            S++;                
    }
}


}  // namespace libnormaliz

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
