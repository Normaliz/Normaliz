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
#include "vector_operations.h"
#include "integer.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using std::cout; using std::endl;

// Constructor, creates 0/1
HilbertSeries::HilbertSeries() {
    num   = vector<num_t>(1,0);
    denom = vector<denom_t>(1,0);
}

// Constructor, creates num/denom, see class description for format
HilbertSeries::HilbertSeries(const vector<num_t>& numerator, const vector<denom_t>& denominator) {
    num   = numerator;
    denom = denominator;
}

// add another HilbertSeries to this
HilbertSeries& HilbertSeries::operator+=(const HilbertSeries& other) {
    vector<num_t> other_num = other.num;
    vector<denom_t> other_denom = other.denom;

    // adjust denominators
    if (denom.size() < other_denom.size()) 
        denom.resize(other_denom.size());
    else if (denom.size() > other_denom.size())
        other_denom.resize(denom.size());
    size_t d_size = denom.size();
    denom_t diff;
    for (size_t i=1; i<d_size; ++i) {
        diff = denom[i] - other_denom[i];
        if (diff > 0) {        // augment other
            other_denom[i] += diff;
            poly_mult_to(other_num, i, diff);
        } else if (diff < 0) { // augment this
            diff = -diff;
            denom[i] += diff;
            poly_mult_to(num, i, diff);
        }
    }
    assert (denom == other_denom);

    // now just add the numerators
    poly_add_to(num,other_num);

    return (*this);
}


// simplify, see class description
void HilbertSeries::simplify() {

    remove_zeros(denom);
    vector<num_t> q, r, poly; //polynomials
    // In denom_cyclo we collect cyclotomic polynomials in the denominator.
    // During this method the Hilbert series is given by num/(denom*denom_cyclo)
    // where denom | denom_cyclo are exponent vectors of (1-t^i) | i-th cyclotminc poly.
    vector<denom_t> denom_cyclo = vector<denom_t>(denom.size());

    for (long i=denom.size()-1; i>0; --i) {
        // check if we can divide the numerator by (1-t^i)
        poly = coeff_vector<num_t>(i);
        while(denom[i]>0) {
            poly_div(q, r, num, poly);
            if (r.size() == 0) { // numerator is divisable by poly
                num = q;
                denom[i]--;
            }
            else {
                break;
            }
        }

        // check if we can divide the numerator by i-th cyclotomic polynomial
        poly = cyclotomicPoly<num_t>(i);
        while(denom_cyclo[i]>0 || denom[i]>0) {
            poly_div(q, r, num, poly);
            if (r.size() == 0) { // numerator is divisable by poly
                num = q;
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
    for (int i = denom.size()-1; i > 0; --i) {
        while(denom_cyclo[i]>0) {
            denom[i]++;
            denom_cyclo[i]--;
            for(int d = 1; d <= i/2; ++d) {
                if (i % d == 0) {
                    if (denom_cyclo[d]>0) {
                        denom_cyclo[d]--;
                    } else {
                        num = poly_mult(num, cyclotomicPoly<num_t>(d));
                    }
                }
            }
        }
    }
    remove_zeros(denom);
}


vector< vector<mpz_class> > HilbertSeries::getHilbertQuasiPolynomial() {
    if(quasi_poly.size()==0) {
        computeHilbertQuasiPolynomial();
    }
    return quasi_poly;
}


void HilbertSeries::computeHilbertQuasiPolynomial() {
    //TODO simplify first?
    long periode = 1; //least common multiple of the degrees of t in the denominator
    long dim = 0;
    long i,j;
    for (long d = denom.size()-1; d > 0; --d) {
        dim += denom[d];
        if (denom[d] > 0 && (periode % d) != 0) {
            periode = lcm<long>(periode, d);
        }
    }
    //periode und dim encode the denominator
    //now adjust the numerator
    long num_size = num.size();
    vector<mpz_class> norm_num(num_size);  //normalized numerator
    for (i = 0;  i < num_size;  ++i) {
        norm_num[i] = to_mpz(num[i]);
    }
    for (long d = denom.size()-1; d > 0; --d) {
        vector<mpz_class> factor, r;
        //nothing to do if it already has the correct t-power or exponent is 0
        if (d != periode && denom[d] > 0) {
            //n_num *= (1-t^p / 1-t^d)^denom[d]
            poly_div(factor, r, coeff_vector<mpz_class>(periode), coeff_vector<mpz_class>(d));
            assert(r.size()==0); //assert rest r is 0
            //TODO more efficient method *=
            //TODO Exponentiation by squaring of factor, then *=
            for (i=0; i<denom[d]; ++i) {
                norm_num = poly_mult(norm_num, factor);
            }
        }
    }
    //cut numerator into periode many pieces and apply standard method
    quasi_poly = vector< vector<mpz_class> >(periode);
    long nn_size = norm_num.size();
    for (j=0; j<periode; ++j) {
        quasi_poly[j].reserve(dim);
    }
    for (i=0; i<nn_size; ++i) {
        quasi_poly[i%periode].push_back((norm_num[i]));
    }

    for (j=0; j<periode; ++j) {
        quasi_poly[j] = compute_polynomial(quasi_poly[j], dim);
    }
    
    //substitute t by t/periode:
    //dividing by periode^dim and multipling the coeff with powers of periode
    mpz_class pp=1;
    for (i = dim-2; i >= 0; --i) {
        pp *= periode; //p^i   ok, it is p^(dim-1-i)
        for (j=0; j<periode; ++j) {
            quasi_poly[j][i] *= pp;
        }
    } //at the end pp=p^dim-1
    //the common denominator for all coefficients
    mpz_class common_denom = permutations<mpz_class>(1,dim) * pp;
    //substitute t by t-j
    for (j=0; j<periode; ++j) {
        linear_substitution<mpz_class>(quasi_poly[j], j); // replaces quasi_poly[j]
    }

}

// returns the numerator, repr. as vector of coefficients, the h-vector
const vector<num_t>& HilbertSeries::getNumerator() const {
    return num;
}
// returns the denominator, repr. as a vector of the exponents of (1-t^i)^e
const vector<denom_t>& HilbertSeries::getDenominator() const {
    return denom;
}

ostream& operator<< (ostream& out, const HilbertSeries& HS) {
    out << "(";
    if (HS.num.size()>0) out << " " << HS.num[0];
    for (size_t i=1; i<HS.num.size(); ++i) {
             if ( HS.num[i]== 1 ) out << " +t^"<<i;
        else if ( HS.num[i]==-1 ) out << " -t^"<<i;
        else if ( HS.num[i] > 0 ) out << " +"<<HS.num[i]<<"*t^"<<i;
        else if ( HS.num[i] < 0 ) out << " -"<<-HS.num[i]<<"*t^"<<i;
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
template<typename Integer>
vector<Integer> coeff_vector(size_t i) {
    vector<Integer> p(i+1,0);
    p[0] =  1;
    p[i] = -1;
    return p;
}
template<typename Integer>
void remove_zeros(vector<Integer>& a) {
    size_t i=a.size();
    while ( i>0 && a[i-1]==0 ) --i;

    if (i < a.size()) {
        a.resize(i);
    }
}

// a += b  (also possible to define the += op for vector)
template<typename Integer>
void poly_add_to (vector<Integer>& a, const vector<Integer>& b) {
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
template<typename Integer>
void poly_sub_to (vector<Integer>& a, const vector<Integer>& b) {
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
template<typename Integer>
vector<Integer> poly_mult(const vector<Integer>& a, const vector<Integer>& b) {
    size_t a_size = a.size();
    size_t b_size = b.size();
    vector<Integer> p( a_size + b_size - 1 );
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
template<typename Integer>
void poly_mult_to(vector<Integer>& a, long d, long e) {
    assert(d>0);
    long i;
    a.reserve(a.size() + d*e);
    while (e>0) {
        a.resize(a.size() + d);
        for (i=a.size()-1; i>=d; --i) {
            a[i] -= a[i-d];
        }
        e--;
    }
}

// division with remainder, a = q * b + r, deg(r) < deg(b), needs |leadcoef(b)| = 1
template<typename Integer>
void poly_div(vector<Integer>& q, vector<Integer>& r, const vector<Integer>& a, const vector<Integer>&b) {
    assert(b.back()!=0); // no unneeded zeros
    assert(b.back()==1 || b.back()==-1); // then division is always possible
    r = a;
    remove_zeros(r);
    size_t b_size = b.size();
    int degdiff = r.size()-b_size; // degree differenz
    if (r.size() < b_size) {
        q = vector<Integer>();
    } else {
        q = vector<Integer>(degdiff+1);
    }
    Integer divisor;
    size_t i=0;

//  cout<<"poly_div: r_start="<<r;
//  cout<<"poly_div: div="<<b;
    while (r.size() >= b_size) {
        
        divisor = r.back()/b.back();
        q[degdiff] = divisor;
        // r -= divisor * t^degdiff * b
        for (i=0; i<b_size; ++i) {
            r[i+degdiff] -= divisor * b[i];
        }
        remove_zeros(r);
//      cout<<"r = "<<r;
//      cout<<"q = "<<q;
        degdiff = r.size()-b_size;
    }

    return;
}

template<typename Integer>
vector<Integer> cyclotomicPoly(long n) {
    // the static variable is initialized only once and then stored
    static vector< vector<Integer> > CyclotomicPoly = vector< vector<Integer> >();
    long computed = CyclotomicPoly.size();
    if (computed < n) {
        CyclotomicPoly.resize(n);
        vector<Integer> poly, q, r;
        for (long i = computed+1; i <= n; ++i) {
            // compute the i-th poly by dividing X^i-1 by the 
            // d-th cycl.poly. with d divides i
            poly = vector<Integer>(i+1);
            poly[0] = -1; poly[i] = 1;  // X^i - 1
            for (long d = 1; d < i; ++d) { // <= i/2 should be ok
                if( i % d == 0) {
                    poly_div(q, r, poly, CyclotomicPoly[d-1]);
                    assert(r.size()==0);
                    poly = q;
                }
            }
            CyclotomicPoly[i-1] = poly;
//          cout << i << "-th cycl. pol.: " << CyclotomicPoly[i-1];
        }
    }
    return CyclotomicPoly[n-1];
}



//---------------------------------------------------------------------------
// computing the Hilbert polynomial from h-vector
//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> compute_e_vector(vector<Integer> Q, int dim){
    int i,j;
    vector <Integer> E_Vector(dim,0);
    Q.resize(dim+1);
    for (i = 0; i <dim; i++) {
        for (j = 0; j <dim; j++) {
            E_Vector[i] += Q[j];
        }
        E_Vector[i]/=permutations<Integer>(1,i);
        for (j = 1; j <=dim; j++) {
            Q[j-1]=j*Q[j];
        }
    }
    return E_Vector;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> compute_polynomial(vector<Integer> h_vector, int dim) {
    vector<Integer> Hilbert_Polynomial = vector<Integer>(dim);
    int i,j;
    
    Integer mult_factor;
    vector <Integer> E_Vector=compute_e_vector(h_vector, dim);
    vector <Integer> C(dim,0);
    C[0]=1;
    for (i = 0; i <dim; i++) {
        mult_factor=permutations<Integer>(i,dim);
        if (((dim-1-i)%2)==0) {
            for (j = 0; j <dim; j++) {
                Hilbert_Polynomial[j]+=mult_factor*E_Vector[dim-1-i]*C[j];
            }
        }
        else {
            for (j = 0; j <dim; j++) {
                Hilbert_Polynomial[j]-=mult_factor*E_Vector[dim-1-i]*C[j];
            }
        }
        for (j = dim-1; 0 <j; j--) {
            C[j]=(unsigned long)(i+1)*C[j]+C[j-1];
        }
        C[0]=permutations<Integer>(1,i+1);
    }

    return Hilbert_Polynomial;
}

//template vector<num_t> compute_polynomial(vector<num_t>, int);


//---------------------------------------------------------------------------

// substitutes t by (t-a), overwrites the polynomial!
template<typename Integer>
void linear_substitution(vector<Integer>& poly, const Integer& a) {
    int dim = poly.size();
    // Iterated division by (t+a)
    for (int step=0; step<dim-1; ++step) {
        for (int i = dim-2; i >= step; --i) {
            poly[i] -= a * poly[i+1];
        }
        //the remainders are the coefficients of the transformed polynomial
    }
}

} //end namespace libnormaliz
