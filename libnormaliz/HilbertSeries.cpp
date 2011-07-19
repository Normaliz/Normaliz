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
#include <cassert>
#include <iostream>
#include "HilbertSeries.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using std::cout; using std::endl;

template <typename Class>
ostream& operator<< (ostream& out, const vector<Class>& vec) {
    for (size_t i=0; i<vec.size(); ++i) {
        out << " " << vec[i];
    }
    out << endl;
    return out;
}


// Constructor, creates 0/1
HilbertSeries::HilbertSeries() {
	nom   = vector<long64>(1,0);
	denom = vector<long64>(1,0);
}

// Constructor, creates nom/denom, see class description for format
HilbertSeries::HilbertSeries(const vector<long64>& nominator, const vector<long64>& denominator) {
	nom   = nominator;
	denom = denominator;
}

// add another HilbertSeries to this
HilbertSeries& HilbertSeries::operator+=(const HilbertSeries& other) {
	vector<long64> other_nom = other.nom;
	vector<long64> other_denom = other.denom;

	// adjust denominators
	if (denom.size() < other_denom.size()) 
		denom.resize(other_denom.size());
	else if (denom.size() > other_denom.size())
		other_denom.resize(denom.size());
	size_t d_size = denom.size();
	long64 diff;
	for (size_t i=1; i<d_size; ++i) {
		diff = denom[i] - other_denom[i];
		if (diff > 0) {        // augment other
			other_denom[i] += diff;
			poly_mult_to(other_nom, i, diff);
		} else if (diff < 0) { // augment this
			diff = -diff;
			denom[i] += diff;
			poly_mult_to(nom, i, diff);
		}
	}
	assert (denom == other_denom);

	// now just add the nominators
	poly_add_to(nom,other_nom);

	return (*this);
}


// simplify, see class description
void HilbertSeries::simplify() {

	remove_zeros(denom);
	vector<long64> q, r, poly; //polynomials
	// In cyclPolyExp we collect cyclotomic polynomials in the denominator,
	// During this method the Hilbert series is given by nom/(denom*denom_cyclo)
	// where denom and denom_cyclo are exponent vectors of (1-t^i), i-th cyclotminc poly.
	vector<long64> denom_cyclo = vector<long64>(denom.size());

	for (int i=denom.size()-1; i>0; --i) {
		// check if we can divide the nominator by (1-t^i)
		poly = coeff_vector(i);
		while(denom[i]>0) {
			poly_div(q, r, nom, poly);
			if (r.size() == 0) { // nominator is divisable by poly
				nom = q;
				denom[i]--;
			}
			else {
				break;
			}
		}

		// check if we can divide the nominator by i-th cyclotomic polynomial
		poly = cyclotomicPoly(i);
		while(denom_cyclo[i]>0 || denom[i]>0) {
			poly_div(q, r, nom, poly);
			if (r.size() == 0) { // nominator is divisable by poly
				nom = q;
				if (denom_cyclo[i]>0) {
					denom_cyclo[i]--;
				} else { // decompose (1-t^i)
					denom[i]--;
					// put the remaining factors of (1-t^i) in denom_cyclo
					for(int d=1; d<=i/2; ++d) {
						if (i % d == 0) denom_cyclo[d]++;
					}
				}
			}
			else {
				break;
			}
		}
	} // end for
	// done with canceling
	// now collect the remaining cyclotomic polynomials in (1-t^i) factors
	for (int i=denom.size()-1; i>0; --i) {
		while(denom_cyclo[i]>0) {
			denom[i]++;
			denom_cyclo[i]--;
			for(int d=1; d<=i/2; ++d) {
				if (i % d == 0) {
					if (denom_cyclo[d]>0) {
						denom_cyclo[d]--;
					} else {
						nom = poly_mult(nom, cyclotomicPoly(d));
					}
				}
			}
		}
	}
	remove_zeros(denom);
}


// returns the nominator, repr. as vector of coefficients, the h-vector
const vector<long64>& HilbertSeries::getNominator() const {
	return nom;
}
// returns the denominator, repr. as a vector of the exponents of (1-t^i)^e
const vector<long64>& HilbertSeries::getDenominator() const {
	return denom;
}

ostream& operator<< (ostream& out, const HilbertSeries& HS) {
	out << "(";
	if (HS.nom.size()>0) out << " " << HS.nom[0];
    for (size_t i=1; i<HS.nom.size(); ++i) {
		     if ( HS.nom[i]== 1 ) out << " +t^"<<i;
		else if ( HS.nom[i]==-1 ) out << " -t^"<<i;
		else if ( HS.nom[i] > 0 ) out << " +"<<HS.nom[i]<<"*t^"<<i;
		else if ( HS.nom[i] < 0 ) out << " -"<<-HS.nom[i]<<"*t^"<<i;
	}
	out << " ) / (";
    for (size_t i=1; i<HS.denom.size(); ++i) {
		if ( HS.denom[i]!=0 ) out << " (1-t^"<<i<<")^"<<HS.denom[i];
	}
	out << " )" << std::endl;
	return out;
}



//---------------------------------------------------------------------------
// polynomial operations, for polynomials repr. as vector of coefficients
//---------------------------------------------------------------------------

// returns the coefficient vector of 1-t^i
vector<long64> coeff_vector(size_t i) {
//	cout << "create 1-t^"<<i<<" :   ";
	vector<long64> p(i+1,0);
	p[0] =  1;
	p[i] = -1;
//	cout<<p[0]<<p[i]<<"*t^"<<i<<"   size="<<p.size() <<endl;
	return p;
}

void remove_zeros(vector<long64>& a) {
	size_t i=a.size();
	while ( i>0 && a[i-1]==0 ) --i;

	if (i < a.size()) {
		a.resize(i);
	}
}

// a += b  (also possible to define the += op for vector)
void poly_add_to (vector<long64>& a, const vector<long64>& b) {
	size_t b_size = b.size();
	if (a.size() < b_size) {
		a.resize(b_size);
	}
	for (size_t i=0; i<b_size; ++i) {
		a[i]+=b[i];
	}
	remove_zeros(a);
}
// a -= b  (also possible to define the += op for vector)
void poly_sub_to (vector<long64>& a, const vector<long64>& b) {
	size_t b_size = b.size();
	if (a.size() < b_size) {
		a.resize(b_size);
	}
	for (size_t i=0; i<b_size; ++i) {
		a[i]-=b[i];
	}
	remove_zeros(a);
}

// a * b
vector<long64> poly_mult(const vector<long64>& a, const vector<long64>& b) {
	size_t a_size = a.size();
	size_t b_size = b.size();
	cout << a_size<<"+"<<b_size<<" oder "<<b.size()<<endl;
	vector<long64> p( a_size + b_size - 1 );
	size_t i,j;
	for (i=0; i<a_size; ++i) {
		if (a[i] == 0) continue;
		for (j=0; j<b_size; ++j) {
			if (b[j] == 0) continue;
			p[i+j] += a[i]*b[j];
		}
	}
	return p;
}

// a *= (1-t^i)^e
void poly_mult_to(vector<long64>& a, size_t d, size_t e) {
	assert(d>0);
	int i;
	a.reserve(a.size() + d*e);
	while (e>0) {
		a.resize(a.size() + d);
		for (i=a.size(); i>=d; --i) {
			a[i] -= a[i-d];
		}
		e--;
	}
}

// division with remainder, a = q * b + r, deg(r) < deg(b), needs |leadcoef(b)| = 1
void poly_div(vector<long64>& q, vector<long64>& r, const vector<long64>& a, const vector<long64>&b) {
	assert(b.back()!=0); // no unneeded zeros
	assert(b.back()==1 || b.back()==-1); // then division is always possible
	r = a;
	remove_zeros(r);
	size_t b_size = b.size();
	int degdiff = r.size()-b_size; // degree differenz
	if (r.size() < b_size) {
		q = vector<long64>();
	} else {
		q = vector<long64>(degdiff+1);
	}
	long64 divisor;
	size_t i=0;

//	cout<<"poly_div: r_start="<<r;
//	cout<<"poly_div: div="<<b;
	while (r.size() >= b_size) {
		
		divisor = r.back()/b.back();
		q[degdiff] = divisor;
		// r -= divisor * t^degdiff * b
		for (i=0; i<b_size; ++i) {
			r[i+degdiff] -= divisor * b[i];
		}
		remove_zeros(r);
//		cout<<"r = "<<r;
//		cout<<"q = "<<q;
		degdiff = r.size()-b_size;
	}

	return;
}

vector<long64> cyclotomicPoly(size_t n) {
	// the static variable is initialized only once and then stored
	static vector< vector<long64> > CyclotomicPoly = vector< vector<long64> >();
	int computed = CyclotomicPoly.size();
	if (computed < n) {
		CyclotomicPoly.resize(n);
		vector<long64> poly, q, r;
		for (int i = computed+1; i <= n; ++i) {
			// compute the i-th poly by dividing X^i-1 by the 
			// d-th cycl.poly. with d divides i
			poly = vector<long64>(i+1);
			poly[0] = -1; poly[i] = 1;  // X^i - 1
			for (int d = 1; d < i; ++d) { // <= i/2 should be ok
				if( i % d == 0) {
					poly_div(q, r, poly, CyclotomicPoly[d-1]);
					assert(r.size()==0);
					poly = q;
				}
			}
			CyclotomicPoly[i-1] = poly;
			cout << i << "-th cycl. pol.: " << CyclotomicPoly[i-1];
		}
	}
	return CyclotomicPoly[n-1];
}

} //end namespace libnormaliz
