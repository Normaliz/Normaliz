/*
 * Normaliz 2.5
 * Copyright (C) 2007-2010  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#ifndef INTEGER_H
#define INTEGER_H

#include "libnormaliz.h"

// Integer should (may) support:
// Integer abs(Integer); here implemented as Iabs
// Integer min(Integer, Integer); here we use the template min in <algorithm>
// It provides abs, gcd and lcm
//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

template<typename Integer> inline long explicit_cast_to_long(const Integer& a) {
	return a;
}
template<> inline long explicit_cast_to_long<mpz_class> (const mpz_class& a) {
	return a.get_si();
}

//---------------------------------------------------------------------------
//                     Basic functions
//---------------------------------------------------------------------------

// returns the absolute value of a
template<typename Integer> inline Integer Iabs(const Integer& a) {
	return (a>=0) ? (a) : Integer(-a);
}

//returns gcd of a and b, 	if one is 0 returns the nonzero one
template<typename Integer> Integer gcd(const Integer& a, const Integer& b);

//returns lcm of a and b, 	returns 0 if one is 0
template<typename Integer> Integer lcm(const Integer& a, const Integer& b);
//---------------------------------------------------------------------------
//                     Special functions
//---------------------------------------------------------------------------

//return the number of decimals, needed to write the Integer a
template<typename Integer> int decimal_length(Integer a);

//returns b!/a!
template<typename Integer> Integer permutations(const int& a, const int& b);

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
