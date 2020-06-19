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
#include <set>
#include <ostream>

#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h"

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
vector<Integer> l_multiplication(const list<vector<Integer> >& l, const vector<Integer>& v) {
    int s = l.size();
    vector<Integer> p(s);
    s = 0;
    for (const auto& i : l) {
        p[s++] = v_scalar_product(*i, v);  // maybe we loose time here?
    }
    return p;
}

//---------------------------------------------------------------------------

template <typename Integer>
list<vector<Integer> > l_list_x_matrix(const list<vector<Integer> >& l, const Matrix<Integer>& M) {
    list<vector<Integer> > result;
    vector<Integer> p;
    for (const auto& i : l) {
        p = M.VxM(i);
        result.push_back(p);
    }
    return result;
}
//---------------------------------------------------------------------------

template <typename Integer>
void l_cut(list<vector<Integer> >& l, int size) {
    for (auto& i : l) {
        i.resize(size);
    }
}

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

// formerly map_operations.h

template <typename key, typename T>
std::ostream& operator<<(std::ostream& out, const map<key, T>& M) {
    for (const auto& it : M) {
        out << it.first << ": " << it.second << "  ";
    }
    out << std::endl;
    return out;
}

//---------------------------------------------------------------------------

template <typename key>
bool contains(const set<key>& m, const key& k) {
    return (m.find(k) != m.end());
}

//---------------------------------------------------------------------------

template <typename key, typename T>
bool contains(const map<key, T>& m, const key& k) {
    return (m.find(k) != m.end());
}

//---------------------------------------------------------------------------

template <typename key, typename T>
map<key, T> count_in_map(const vector<key>& v) {
    map<key, T> m;
    T size = v.size();
    for (T i = 0; i < size; ++i) {
        m[v[i]]++;
    }
    return m;
}

template <typename key, typename T>
vector<key> to_vector(const map<key, T>& M) {
    vector<key> v;
    for (const auto& it : M) {
        for (T i = 0; i < it.second; i++) {
            v.push_back(it.first);
        }
    }
    return v;
}

}  // namespace libnormaliz

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
