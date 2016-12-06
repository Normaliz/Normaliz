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

#include <algorithm>
#include <sstream>
#include "libQnormaliz/Qinteger.h"
#include "libQnormaliz/Qvector_operations.h"

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

bool try_convert(long& ret, const long long& val) {
    assert(false); return true;
}

bool try_convert(long& ret, const mpq_class& val) {
    assert(false); return true;
}

bool try_convert(long long& ret, const mpq_class& val) {
    assert(false); return true;
}

bool try_convert(mpq_class& ret, const long long& val) {
    assert(false); return true;
}

bool fits_long_range(long long a) {
    return sizeof(long long) == sizeof(long) || (a <= LONG_MAX && a >= LONG_MIN);
}

//---------------------------------------------------------------------------

template <typename Number>
size_t decimal_length(Number a){
    /* size_t l=1;
    if (a<0) {
        a=-a;
        l++;
    }
    while((a/=10)!=0)
        l++;*/

    ostringstream test;
    test << a;
    return test.str().size();
}

//---------------------------------------------------------------------------

template <typename Number>
Number permutations(const size_t& a, const size_t& b){
    unsigned long i;
    Number P=1;
    for (i = a+1; i <= b; i++) {
        P*=i;
    }
    return P;
}

//---------------------------------------------------------------------------

template<typename Number> 
Number permutations_modulo(const size_t& a, const size_t& b, long m) {
    unsigned long i;
    Number P=1;
    for (i = a+1; i <= b; i++) {
        P*=i; P%=m;
    }
    return P;
}




} //end namespace libQnormaliz
