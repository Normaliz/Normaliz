/*
 * Copyright (C) 2007-2022  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */


#include <fstream>
#include <sstream>
#include <string>
#include <gmpxx.h>

#ifdef NMZ_COCOA
#include "libnormaliz/nmz_integrate.h"
#endif

#include "libnormaliz/nmz_polynomial.h"
#include "libnormaliz/vector_operations.h"

namespace libnormaliz {

using namespace std;

//-------------------------------------------------------------------
//       OurTerm
//-------------------------------------------------------------------



template<typename Number>
OurTerm<Number>::OurTerm(){

}

template<typename Number>
OurTerm<Number>::OurTerm(const Number& c, const map<key_t, long>& mon, const dynamic_bitset& supp){
    coeff = c;
    monomial = mon;
    support = supp;
}

template<typename Number>
Number OurTerm<Number>::evaluate(const vector<Number>& argument) const{

    Number value = coeff;
    for(auto& M: monomial){
        for(long i = 0; i < M.second; ++i){
         value *= argument[M.first];
         //if(!check_range(value)) -- check done at addition in polynomial
         //    throw ArithmeticException("Overflow in evaluation of polynomial");
        }
    }
    return value;
}

template<typename Number>
void OurTerm<Number>::shift_coordinates(const int& shift){

    OurTerm<Number> transformed;
    transformed.support = dynamic_bitset(support.size() + shift);
    for(auto F: monomial){
        key_t cc = F.first;
        if(shift < 0)
            assert(cc >= -shift);
        cc += shift;
        transformed.support[cc] = 1;
        (transformed.monomial)[cc] = F.second;
    }
    transformed.coeff = coeff;
    *this = transformed;
}

template<typename Number>
void OurTerm<Number>::swap_coordinates(const key_t& first, const key_t& second){

    OurTerm<Number> transformed;
    transformed.support = dynamic_bitset(support.size());
    transformed.coeff = coeff;
    for(auto F: monomial){
        key_t cc = F.first;

        if(cc == first){
            cc = second;
        }
        else
            if(cc == second)
                cc = first;
        transformed.monomial[cc] = F.second;
        transformed.support[cc] = 1;
    }
    *this = transformed;

}

template<typename Number>
void OurTerm<Number>::cyclic_shift_right(const key_t& col){

    v_cyclic_shift_right(support, col);
    vector<long> expo_vec(support.size());
    for(auto& E: monomial)
        expo_vec[E.first] = E.second;
    v_cyclic_shift_right(expo_vec, col);
    monomial.clear();
    for(int i = 0; i< expo_vec.size(); ++i){
        if(expo_vec[i] >0 )
            monomial[i] = expo_vec[i];
    }
}



template <typename Number>
ostream& operator<<(ostream& out, const OurTerm<Number> & T) {
    out << "coeff " << T.coeff << " --- " << T.support << " ---";
    for(auto& F: T.monomial)
        out << F.first << ":" << F.second << "  ";
    out << endl;
    return out;
}

template<typename Number>
void OurTerm<Number>::multiply_by_constant(const Number& factor){
    coeff *= factor;
}

//-------------------------------------------------------------------
//             OurPolynomial
//-------------------------------------------------------------------


template<typename Number>
OurPolynomial<Number>::OurPolynomial(){

}

template<typename Number>
key_t OurPolynomial<Number>::get_highest_indet() const{
        return highest_indet;
}

template<typename Number>
void OurPolynomial<Number>::shift_coordinates(const int& shift){

    support = dynamic_bitset(support.size() + shift);
    for(auto& M: *this){
        M.shift_coordinates(shift);
        support |= M.support;
    }
    highest_indet +=shift;
}

template<typename Number>
void OurPolynomial<Number>::swap_coordinates(const key_t& first, const key_t& second){

    for(auto& M: *this){
        M.swap_coordinates(first, second);
    }

    bool temp = support[first];
    support[first] = support[second];
    support[second] = temp;
    for(size_t i = 0; i < support.size(); ++i){
        if(support[i])
            highest_indet = i;
    }
}

template <typename Number>
ostream& operator<<(ostream& out, const OurPolynomial<Number> & P) {
    out << "terms" << endl;
    for(auto& T: P)
        out << T;
    out << "highest indet " << P.highest_indet << " support " << P.support << endl;
    return out;
}

template<typename Number>
Number OurPolynomial<Number>::evaluate(const vector<Number>& argument) const{

    Number value = 0;
    for(auto& T: *this){
        value += T.evaluate(argument);
        if(!check_range(value))
             throw ArithmeticException("Overflow in evaluation of polynomial");
    }
    return value;
}

template<typename Number>
void OurPolynomial<Number>::cyclic_shift_right(const key_t& col){
    for(auto& T: *this)
        T.cyclic_shift_right(col);

    v_cyclic_shift_right(support, col);
    for(size_t i = 0; i < support.size(); ++i){
        if(support[i])
            highest_indet = i;
    }
}

template<typename Number>
void OurPolynomial<Number>::multiply_by_constant(const Number& factor){
    for(auto& T: *this)
        T.multiply_by_constant(factor);
}


//-------------------------------------------------------------------
//             OurPolynomialSystem
//-------------------------------------------------------------------


template<typename Number>
OurPolynomialSystem<Number>::OurPolynomialSystem(){

}

template<typename Number>
OurPolynomialSystem<Number>::OurPolynomialSystem(const vector<string>& poly_strings, size_t dim, bool verb){

    verbose = verb;
    for(auto& S: poly_strings){
        OurPolynomial<Number> poly(S,dim,verbose);
        this->push_back(poly);
    }
}

template<typename Number>
bool OurPolynomialSystem<Number>::check(const vector<Number>& argument, const bool is_equations, const bool exact_length) const{

    Number test;
    for(auto& P: *this){
        if(P.highest_indet > argument.size() -1)
            continue;
        if(P.highest_indet < argument.size() - 1 && exact_length)
            continue;
        test = P.evaluate(argument);
        if(is_equations && test != 0)
            return false;
        if(!is_equations && test < 0)
            return false;
    }
    return true;
}

template<typename Number>
void OurPolynomialSystem<Number>::shift_coordinates(const int& shift){
    for(auto& P: *this)
        P.shift_coordinates(shift);
}

template<typename Number>
void OurPolynomialSystem<Number>::swap_coordinates(const key_t& first, const key_t& second){
    for(auto& P: *this)
        P.swap_coordinates(first, second);
}

template<typename Number>
void OurPolynomialSystem<Number>::cyclic_shift_right(const key_t& col){
    for(auto& P: *this)
        P.cyclic_shift_right(col);
}

template<typename Number>
void OurPolynomialSystem<Number>::multiply_by_constant(const Number& factor){
    for(auto& P: *this)
        P.multiply_by_constant(factor);
}

template <typename Number>
ostream& operator<<(ostream& out, const OurPolynomialSystem<Number> & S) {
    out << "*****************************" << endl;
    out << "system" << endl;
    for(auto& P: S){
        cout << "************" << endl;
        out << P;
    }
    out << "*****************************" << endl;
    return out;
}



#ifdef NMZ_COCOA


template<typename Number>
OurPolynomial<Number>::OurPolynomial(const string& poly_string, const size_t dim, const bool verbose){
    GlobalManager CoCoAFoundations;

    /*SparsePolyRing RQQ = NewPolyRing_DMPI(RingQQ(), dim + 1, lex);
    string poly_string_new("x[1]^2/2+1/3");
    RingElem FQQ = ReadExpr(RQQ, poly_string_new);
    cout << "FFF " << FQQ << endl;*/

    if(verbose)
        verboseOutput() << poly_string << endl;

    SparsePolyRing RQQ = NewPolyRing_DMPI(RingQQ(), dim + 1, lex);
    RingElem FQQ = ReadExpr(RQQ, poly_string);

    // cout << "DDDD " << FQQ << endl;
    FQQ = ClearDenom(FQQ);

    // cout << "FFFF " << FQQ << endl;

    SparsePolyRing R = NewPolyRing_DMPI(RingZZ(), dim + 1, lex); // in the input shift_coordinates numbered from 1
    RingElem F = makeZZCoeff(FQQ, R);

    // cout << "ZZZZ " << F << endl;

    vector<long> v(NumIndets(R));
    BigInt BI_coeff;
    mpz_class mpz_coeff;
    key_t max_indet = 0;
    support = dynamic_bitset(dim +1);

    INTERRUPT_COMPUTATION_BY_EXCEPTION

    SparsePolyIter mon = BeginIter(F);

    for (; !IsEnded(mon); ++mon) {
        OurTerm<Number> T;

        IsInteger(BI_coeff, coeff(mon)); // in two steps from the coefficient of the term
        mpz_coeff = mpz(BI_coeff);    // to mpz_class
        T.coeff = convertTo<Number>(mpz_coeff); // and one more conversion

        exponents(v, PP(mon));  // this function gives the exponent vector back as v
        T.support = v_support(v);
        for(key_t i = 0; i < v.size(); ++i){
            if(v[i] != 0){
                if(i > max_indet)
                    max_indet = i;
                T.monomial[i] = v[i];
            }
        }
        this->push_back(T);
        support |= T.support;
    }
    highest_indet = max_indet;
}
#else

template<typename Number>
OurPolynomial<Number>::OurPolynomial(const string& poly_string, const size_t dim, const bool verbose){
    assert(false);
}

#endif // NMZ_COCOA

template class OurTerm<long>;
template class OurTerm<long long>;
template class OurTerm<mpz_class>;
#ifdef ENFNORMALIZ
template class OurTerm<renf_elem_class>;
#endif

template class OurPolynomial<long>;
template class OurPolynomial<long long>;
template class OurPolynomial<mpz_class>;
#ifdef ENFNORMALIZ
template class OurPolynomial<renf_elem_class>;
#endif

template class OurPolynomialSystem<long>;
template class OurPolynomialSystem<long long>;
template class OurPolynomialSystem<mpz_class>;
#ifdef ENFNORMALIZ
template class OurPolynomialSystem<renf_elem_class>;
#endif

} // namespace


