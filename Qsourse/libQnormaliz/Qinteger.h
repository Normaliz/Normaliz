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

#ifndef INTEGER_H_
#define INTEGER_H_

#include <libQnormaliz/Qgeneral.h>

#include <list>
#include <vector>
#include <string>
#include <limits.h>

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//                     Basic functions
//---------------------------------------------------------------------------

// returns the absolute value of a
template<typename Integer> inline Integer Iabs(const Integer& a) {
    return (a>=0) ? (a) : Integer(-a);
}

//---------------------------------------------------------------------------
//                     Conversions and checks
//---------------------------------------------------------------------------

// convert val to ret
// does the conversion and returns false if it fails
bool try_convert(long& ret, const long long& val);
inline bool try_convert(long long& ret, const long& val) {assert(false); return true;}
bool try_convert(long& ret, const mpz_class& val);
bool try_convert(long long& ret, const mpz_class& val);
inline bool try_convert(mpz_class& ret, const long& val) {assert(false); return true;}
bool try_convert(mpz_class& ret, const long long& val);

// template for same typ "conversion"
template<typename Type>
inline bool try_convert(Type& ret, const Type& val) {ret = val; return true;}


bool fits_long_range(long long a);

template<typename Integer>
inline bool using_GMP() {
  return false;
}

template<>
inline bool using_GMP<mpq_class>() {
  return true;
} 
//---------------------------------------------------------------------------

// Should be completely remoced:
template<typename Integer>
inline bool check_range(const Integer& m) {
    return true;
}

//---------------------------------------------------------------------------
//                     Special functions
//---------------------------------------------------------------------------

//return the number of decimals, needed to write the Integer a
template<typename Integer> size_t decimal_length(Integer a);

//returns b!/a!
template<typename Integer> Integer permutations(const size_t& a, const size_t& b);
template<typename Integer> Integer permutations_modulo(const size_t& a, const size_t& b, long m);

//---------------------------------------------------------------------------
//                     String conversion functions
//---------------------------------------------------------------------------

// forward declaration to silence clang error:
// 'operator<<' should be declared prior to the call site or in the global namespace
template <typename T> std::ostream& operator<< (std::ostream& out, const vector<T>& vec);

template<typename Integer> string toString(Integer a) {
    ostringstream ostream;
    ostream << a;
    return ostream.str();
}
template<> inline string toString(mpz_class a) {
    return a.get_str();
}
template<> inline string toString(mpq_class a) {
    return a.get_str();
}

} // end libnormaliz

//---------------------------------------------------------------------------
#endif /* INTEGER_H_ */
//---------------------------------------------------------------------------
