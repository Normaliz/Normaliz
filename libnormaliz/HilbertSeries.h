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

/*
 * HilbertSeries represents a Hilbert series of a N-graded algebra
 * with non-homogeneous generators.
 * It is represented as a polynomial divided by a product of (1-t^i).
 * The nominator is represented as vector of coefficients, the h-vector
 * h vector repr.: sum of h[i]*t^i
 * and the denominator is represented as a vector of the exponents of (1-t^i)
 * d vector repr.: product of (1-t^i)^d[i] over i > 0
 *
 * The class offers basic operations on the series and a simplification which
 * transforms the series in the simplest presentation of such a form in terms
 * of the degrees of the nominator and denominator.
 *
 * Furthermore this file include operations for the polynomials used.
 */

//---------------------------------------------------------------------------

#ifndef HILBERT_SERIES_H
#define HILBERT_SERIES_H

//---------------------------------------------------------------------------

#include <vector>
#include <ostream>

//#include "general.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using std::vector;
using std::ostream;

template<typename Integer> class HilbertSeries;

// write a readable representation to the stream
template<typename Integer>
ostream& operator<< (ostream& out, const HilbertSeries<Integer>& HS);


template<typename Integer>
class HilbertSeries {

public:
	// Constructor, creates 0/1
	HilbertSeries();
	// Constructor, creates nom/denom, see class description for format
	HilbertSeries(const vector<Integer>& nom, const vector<Integer>& denom);

	// add another HilbertSeries to this
	HilbertSeries<Integer>& operator+=(const HilbertSeries<Integer>& other);

/*	// multiply the series with 1/(1-t^i)
	void multiply_denom(size_t i);
	// add t^i to the nominator
	void add_to_nom(size_t i);
*/

	// simplify, see class description
	void simplify();

	// returns the nominator, repr. as vector of coefficients, the h-vector
	const vector<Integer>& getNominator() const;
	// returns the denominator, repr. as a vector of the exponents of (1-t^i)^e
	const vector<Integer>& getDenominator() const;

private:
	// the nominator, repr. as vector of coefficients, the h-vector
	vector<Integer> nom;
	// the denominator, repr. as a vector of the exponents of (1-t^i)^e
	vector<Integer> denom;

	friend ostream& operator<< <> (ostream& out, const HilbertSeries<Integer>& HS);

};
//class end *****************************************************************


//---------------------------------------------------------------------------
// polynomial operations, for polynomials repr. as vector of coefficients
//---------------------------------------------------------------------------

// a += b
template<typename Integer>
void poly_add_to (vector<Integer>& a, const vector<Integer>& b);

// a -= b
template<typename Integer>
void poly_subt_to (vector<Integer>& a, const vector<Integer>& b);


// a * b
template<typename Integer>
vector<Integer> poly_mult(const vector<Integer>& a, const vector<Integer>& b);

// a *= (1-t^d)^e
template<typename Integer>
void poly_mult_to(vector<Integer>& a, size_t d, size_t e = 1);


// division with remainder, a = q * b + r
template<typename Integer>
void poly_div(vector<Integer>& q, vector<Integer>& r, const vector<Integer>& a, const vector<Integer>&b);


// remove leading zero coefficients, 0 polynomial leads to empty list
template<typename Integer>
void remove_zeros(vector<Integer>& a);


// Returns the n-th cyclotomic polynomial, all smaller are computed and stored.
// The n-th cyclotomic polynomial is the product of (X-w) over all 
// n-th primitive roots of unity w.
template<typename Integer>
vector<Integer> cyclotomicPoly(size_t n);

// returns the coefficient vector of 1-t^i
template<typename Integer>
vector<Integer> coeff_vector(size_t i);

} //end namespace libnormaliz

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------

