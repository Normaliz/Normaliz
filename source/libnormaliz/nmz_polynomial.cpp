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
#include "libnormaliz/list_and_map_operations.h"

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
    mon2vars_expos();
}


template<typename Number>
OurTerm<Number>::OurTerm(const pair<vector<key_t>, Number>& t, size_t dim){
    coeff = t.second;
    monomial = count_in_map<key_t, long>(t.first);
    support = dynamic_bitset(dim);
    for(auto& m: monomial)
        support[m.first] = 1;
    mon2vars_expos();
}

template<typename Number>
void OurTerm<Number>::mon2vars_expos(){
    vars.clear();
    for(auto& M: monomial){
        for(size_t i = 0; i < M.second; ++i)
            vars.push_back(M.first);
    }
}

template<typename Number>
Number OurTerm<Number>::evaluate(const vector<Number>& argument) const{

    Number value = coeff;
    for(size_t i = 0; i < vars.size(); ++i)
        value *= argument[vars[i]];
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
    mon2vars_expos();
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
    mon2vars_expos();
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
    mon2vars_expos();
}

template<typename Number>
void OurTerm<Number>::permute_variables(const vector<key_t>& perm){
    vector<long> expo_vec(support.size());
    map<key_t, long> new_mon;
    for( auto& E: monomial)
        expo_vec[E.first] = E.second;
    // cout << "EEEEEEEEEEEEEEEE " << expo_vec;
    expo_vec = v_permute_coordinates(expo_vec, perm);
    // cout << "FFFFFFFFFFFFFFFF " << expo_vec;
    for(size_t i = 0; i < perm.size(); ++i){
        if(expo_vec[i] != 0)
            new_mon[i] = expo_vec[i];
    }
    monomial = new_mon;
    support = v_permute_coordinates(support, perm);
    mon2vars_expos();
}

template<typename Number>
void OurTerm<Number>::multiply_by_constant(const Number& factor){
    coeff *= factor;
}

/*
template<typename Number>
bool OurTerm<Number>::check_restriction(const dynamic_bitset& set_of_var)  const{
    return support.is_subset_of(set_of_var);
} */

template<typename Number>
bool OurTerm<Number>::is_restrictable_inequ(const dynamic_bitset& set_of_var)  const{
    return support.is_subset_of(set_of_var) || (coeff <= 0);
}

//-------------------------------------------------------------------
//             OurPolynomial
//-------------------------------------------------------------------


template<typename Number>
OurPolynomial<Number>::OurPolynomial(){
    vectorized = false;
}


template<typename Number>
OurPolynomial<Number>::OurPolynomial(const map<vector<key_t>, Number>& poly, size_t dim){

    vectorized = false;
    support = dynamic_bitset(dim);
    for(auto& t: poly){
        pair<vector<key_t>, Number> t_0 = make_pair(t.first, t.second);
        this->push_back(OurTerm<Number>(t_0,dim));
        support |= this->back().support;
    }
    highest_indet = -1;
    for(size_t i = 0; i < support.size(); ++i){
        if(support[i])
            highest_indet = i;
    }
}

template<typename Number>
OurPolynomial<Number>::OurPolynomial(const vector<Number>& linear_form){

    vectorized = false;
    for(size_t i = 0; i < linear_form.size(); ++i){
        if(linear_form[i] == 0)
            continue;
        dynamic_bitset term_supp(linear_form.size());
        term_supp[i] = true;
        map<key_t, long> term_mon;
        term_mon[i] = 1;
        this->push_back(OurTerm<Number>(linear_form[i], term_mon, term_supp));
    }
    this->support = v_support(linear_form);
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
    if(highest_indet >0){
        highest_indet +=shift;
        assert(highest_indet >= 0);
    }
}

template<typename Number>
void OurPolynomial<Number>::swap_coordinates(const key_t& first, const key_t& second){

    for(auto& M: *this){
        M.swap_coordinates(first, second);
    }

    bool temp = support[first];
    support[first] = support[second];
    support[second] = temp;
    highest_indet = -1;
    for(size_t i = 0; i < support.size(); ++i){
        if(support[i])
            highest_indet = i;
    }
}

template<typename Number>
Number OurPolynomial<Number>::evaluate(const vector<Number>& argument) const{

    Number value = 0;
    if(vectorized){
        return evaluate_vectorized(argument);
    }
    for(auto& T: *this){
        value += T.evaluate(argument);
        if(!check_range(value))
             throw ArithmeticException("Overflow in evaluation of polynomial");
    }
    return value;
}

template<typename Number>
Number OurPolynomial<Number>::evaluate_vectorized(const vector<Number>& argument) const{

    Number value = const_term;
    for(size_t i = 0; i <  expo_1_pos.size(); ++i){
        value += argument[expo_1_pos[i]] *  argument[expo_2_pos[i]];
    }
    for(size_t i = 0; i <  expo_1_neg.size(); ++i){
        value -= argument[expo_1_neg[i]] *  argument[expo_2_neg[i]];
    }
    return value;
}


template<typename Number>
OurPolynomial<Number> OurPolynomial<Number>::restrict_to(const dynamic_bitset& variables) const{
    OurPolynomial<Number> Rest;
    for(auto& T: *this){
     if(T.support.is_subset_of(variables))
         Rest.push_back(T);
    }
    return Rest;
}


// splits the polynomial into two parts: terms whose support is contained in support_variables
// and the remaining terms
template<typename Number>
pair<OurPolynomial<Number>, OurPolynomial<Number> > OurPolynomial<Number>::split(const dynamic_bitset& support_variables) const{
    OurPolynomial<Number> Rest;
    OurPolynomial<Number> LeftOver;
    for(auto& T: *this){
     if(T.support.is_subset_of(support_variables))
         Rest.push_back(T);
    else
        LeftOver.push_back(T);
    }
    return make_pair(Rest, LeftOver);
}


template<typename Number>
bool OurPolynomial<Number>::check_linearity(const dynamic_bitset& critical_variables, dynamic_bitset& support_linear) const{
    for(auto& T: *this){
        dynamic_bitset common = T.support & critical_variables;
        if(common.count() == 0)
            return false;
        support_linear |= common;
    }
    return true;
}

template<typename Number>
void OurPolynomial<Number>::vectorize_deg_2(){
    vector<key_t> fact_1_pos, fact_2_pos;
    vector<key_t> fact_1_neg, fact_2_neg;
    // vector<Number> coe;
    Number ct = 0;
    for(auto& T: *this){
        if(T.vars.size() != 2 && T.vars.size() != 0)
            return;
        if(T.vars.size() == 0){
            ct += T.coeff;
            continue;
        }
        if(T.vars.size() == 2){
            if(T.coeff !=1 && T.coeff != -1)
                return;
            if(T.coeff == 1){
                fact_1_pos.push_back(T.vars[0]);
                fact_2_pos.push_back(T.vars[1]);
            }
            if(T.coeff == -1){
                fact_1_neg.push_back(T.vars[0]);
                fact_2_neg.push_back(T.vars[1]);
            }
            // coe.push_back(T.coeff);
        }
    }
    expo_1_pos = fact_1_pos;
    expo_2_pos = fact_2_pos;
    expo_1_neg = fact_1_neg;
    expo_2_neg = fact_2_neg;
    // coeffs = coe;
    const_term = ct;
    vectorized = true;
    (*this).clear();
}



template<typename Number>
Number OurPolynomial<Number>::evaluate_restricted(const vector<Number>& argument, const dynamic_bitset& set_of_var) const{
    Number value = 0;
    for(auto& T: *this){
        if(T.support.is_subset_of(set_of_var))
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
    highest_indet = -1;
    for(size_t i = 0; i < support.size(); ++i){
        if(support[i])
            highest_indet = i;
    }
}


template<typename Number>
void OurPolynomial<Number>::permute_variables(const vector<key_t>& perm){
    for(auto& T: *this)
        T.permute_variables(perm);
    support = v_permute_coordinates(support, perm);
    highest_indet = -1;
    for(size_t i = 0; i < support.size(); ++i)
        if(support[i])
            highest_indet = i;
}

template<typename Number>
void OurPolynomial<Number>::multiply_by_constant(const Number& factor){
    for(auto& T: *this)
        T.multiply_by_constant(factor);
}

template<typename Number>
bool OurPolynomial<Number>::is_restrictable_inequ(const dynamic_bitset& set_of_var)  const{
    size_t nr_negative = 0;
    for(auto& T: *this){
        if(!T.is_restrictable_inequ(set_of_var))
            return false;
        if(T.support.is_subset_of(set_of_var) && T.coeff < 0)
            nr_negative++;
    }
    return nr_negative >= 4;
}

/*
template<typename Number>
bool OurPolynomial<Number>::check_restriction(const dynamic_bitset& set_of_var)  const{
    for(auto& T: *this){
        if(!T.check_restriction(set_of_var))
            return false;
    }
    return true;
}*/

//-------------------------------------------------------------------
//             OurPolynomialCong
//-------------------------------------------------------------------

template<typename Number>
OurPolynomialCong<Number>::OurPolynomialCong(){

}

template<typename Number>
OurPolynomialCong<Number>::OurPolynomialCong(const OurPolynomial<Number>& pol, const Number& mod){
        poly = pol;
        modulus = mod;
}

template<typename Number>
OurPolynomialCong<Number>::OurPolynomialCong(vector<Number> cong){
        modulus = cong.back();
        cong.pop_back();
        poly = OurPolynomial<Number>(cong);
}

template<typename Number>
bool OurPolynomialCong<Number>::check(const vector<Number>& v) const{
    if(poly.evaluate(v) % modulus != 0)
        return false;
    return true;
}

template<>
bool OurPolynomialCong<renf_elem_class>::check(const vector<renf_elem_class>& v) const{
    assert(false);
    return false;
}

//-------------------------------------------------------------------
//             OurPolynomialSystem
//-------------------------------------------------------------------


template<typename Number>
OurPolynomialSystem<Number>::OurPolynomialSystem(){

}

template<typename Number>
OurPolynomialSystem<Number>::OurPolynomialSystem(const set<map<vector<key_t>, Number> >& Polys, size_t dim){
    for(auto& p: Polys)
        this->push_back(OurPolynomial<Number>(p, dim + 1)); // sme conventions for dim as in making from strings by CoCoA below
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
void OurPolynomialSystem<Number>::permute_variables(const vector<key_t>& perm){
    for(auto& P: *this)
        P.permute_variables(perm);
}

template<typename Number>
void OurPolynomialSystem<Number>::multiply_by_constant(const Number& factor){
    for(auto& P: *this)
        P.multiply_by_constant(factor);
}

#ifdef NMZ_COCOA

template<typename Number>
OurPolynomial<Number>::OurPolynomial(const string& poly_string, const size_t dim, const bool verbose){
    GlobalManager CoCoAFoundations;

    /*SparsePolyRing RQQ = NewPolyRing_DMPI(RingQQ(), dim + 1, lex);
    string poly_string_new("x[1]^2/2+1/3");
    RingElem FQQ = ReadExpr(RQQ, poly_string_new);
    cout << "FFF " << FQQ << endl;*/

    vectorized = false;

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
    long max_indet = -1;
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
        for(long i = 0; i < v.size(); ++i){
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

template<typename Number>
RingElem OurTerm<Number>::ToCoCoA(SparsePolyRing R) const{

    mpq_class c;
    c = convertTo<mpq_class>(coeff);
    BigRat ccc = BigRatFromMPQ(c.get_mpq_t());
    RingElem h(R,ccc);
    for(auto& v:vars)
        h *= indet(R,v);
    return h;
}

/*
// We need the special version for long long to avoid a conversion problem
template<>
RingElem OurTerm<long long>::ToCoCoA(SparsePolyRing R) const{

    mpz_class c_mpz = convertTo<mpz_class>(coeff);
    mpq_class c = c_mpz;
    BigRat ccc = BigRatFromMPQ(c.get_mpq_t());
    RingElem h(R,ccc);
    for(auto& v:vars)
        h *= indet(R,v);
    return h;
}

// Another special version ...
template<>
RingElem OurTerm<mpz_class>::ToCoCoA(SparsePolyRing R) const{

    mpq_class c = coeff;
    BigRat ccc = BigRatFromMPQ(c.get_mpq_t());
    RingElem h(R,ccc);
    for(auto& v:vars)
        h *= indet(R,v);
    return h;
}
*/


template<typename Number>
RingElem OurPolynomial<Number>::ToCoCoA(SparsePolyRing R) const{

    RingElem p = zero(R);
    for(auto& T:*this)
        p += T.ToCoCoA(R);
    return p;
}

template<typename Number>
vector<RingElem> OurPolynomialSystem<Number>::ToCoCoA(SparsePolyRing R) const{

    vector<RingElem> CS;
    for(auto& P:*this)
        CS.push_back(P.ToCoCoA(R));
    return CS;
}

bool poly_reduce(RingElem& r, const RingElem&g, const PPMonoidElem& ini_g, const RingElem& h){

    if(!IsDivisible(ini_g, LPP(h)))
        return false;

    PPMonoidElem quot_PP = ini_g/LPP(h);
    RingElem quot_coeff = LC(g)/LC(h);
    RingElem quot = monomial(owner(h), quot_coeff, quot_PP);
    r = g - quot*h;
    return true;
}

bool GB_reduce(const RingElem& f, vector<RingElem>& GB){

    RingElem g = f;
    RingElem r;

    while(g != 0){
        PPMonoidElem ini_g = LPP(g);
        bool reducible = false;
        for(auto& h: GB){
            if(poly_reduce(r,g, ini_g, h)){
                reducible = true;
                g = r;
                break;
            }
        }
        if(!reducible)
            break;
    }
    if(g != 0){
        GB.push_back(g);
        return false;
    }
    else
        return true;
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
OurPolynomialSystem<Number> OurPolynomialSystem<Number>::minimize_equations(const Matrix<Number>& LinEqus) const {

    // cout << "RRRRRRRRRRRRRR  " << LinEqus.rank() << endl;

    size_t EmbDim = LinEqus.nr_of_columns();

    CoCoA::GlobalManager CoCoAFoundations;
    CoCoA::SparsePolyRing R = CoCoA::NewPolyRing_DMPI(CoCoA::RingQQ(), EmbDim , CoCoA::lex);

    /* for(auto& p: HomPol)
        cout << p << endl;*/

    OurPolynomialSystem<Number> LinPolys;
    for(size_t i = 0; i < LinEqus.nr_of_rows(); ++i){ // automatically homogeneous
        LinPolys.push_back(OurPolynomial<Number>(LinEqus[i]));    }
    vector<RingElem> CLin = LinPolys.ToCoCoA(R);

    vector<RingElem> CPol = ToCoCoA(R);
    vector<RingElem> HomPol;
    for(auto& p: CPol){
        if(deg(p) > 2)
            throw BadInputException("Minimization of polynomial is_equations only possible for degree <= 2");
        if(deg(p) == 1)
            HomPol.insert(HomPol.begin(), homogenize(p)); // make sure fegree 1 is inserted into GB before degree 2
        else
            HomPol.push_back(homogenize(p));
    }

    vector<RingElem> GB;
    OurPolynomialSystem<Number> Minis;

    for(auto& lf: CLin){
        GB_reduce(lf,GB);
    }
    // size_t LGGGG = GB.size();
    // cout << "GGGGGGGGGGGGGGGGGGGGGG " << GB.size() << endl;

    for(size_t j = 0; j < this->size(); ++j){
        if(!GB_reduce(HomPol[j],GB))
            Minis.push_back((*this)[j]);
    }

    /* for(auto& p: Minis){
        cout << endl;
        cout << p.ToCoCoA(R) << endl;
    }

    cout << "MMMMMMMM  " << Minis.size()<< " GGGGGGGGGGGG " << GB.size() - LGGGG << endl;*/

    return Minis;
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

template class OurPolynomialCong<long>;
template class OurPolynomialCong<long long>;
template class OurPolynomialCong<mpz_class>;
#ifdef ENFNORMALIZ
template class OurPolynomialCong<renf_elem_class>;
#endif

template class OurPolynomialSystem<long>;
template class OurPolynomialSystem<long long>;
template class OurPolynomialSystem<mpz_class>;
#ifdef ENFNORMALIZ
template class OurPolynomialSystem<renf_elem_class>;
#endif

} // namespace


