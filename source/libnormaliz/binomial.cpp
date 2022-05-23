/*
 * Normaliz
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

#include "libnormaliz/binomial.h"

using std::cout;
using std::endl;
using std::string;
using std::to_string;

unsigned long long winf_ini_coprime = 0;
unsigned long long winf_gm_left = 0;
unsigned long long winf_tail_not_coprime = 0;
unsigned long long winf_s_poly = 0;
unsigned long long winf_red = 0;
unsigned long long winf_red_tail = 0;
unsigned long long winf_red_zero = 0;
unsigned long long winf_red_steps = 0;
unsigned long long winf_gm_steps = 0;
unsigned long long winf_entered_nodes = 0;

void reset_statistics(){
    winf_ini_coprime = 0;
    winf_gm_left = 0;
    winf_tail_not_coprime = 0;
    winf_s_poly = 0;
    winf_red = 0;
    winf_red_tail = 0;
    winf_red_zero = 0;
    winf_red_steps = 0;
    winf_gm_steps = 0;
    winf_entered_nodes = 0;
}


struct timeval OUR_TIME_begin, OUR_TIME_end;

void OURStartTime() {
    gettimeofday(&OUR_TIME_begin, 0);
}

void OURMeasureTime(bool verbose, const std::string& step) {
    gettimeofday(&OUR_TIME_end, 0);
    long seconds = OUR_TIME_end.tv_sec - OUR_TIME_begin.tv_sec;
    long microseconds = OUR_TIME_end.tv_usec - OUR_TIME_begin.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    if (verbose)
        std::cout << step << ": " << elapsed << " sec" << std::endl;
    OUR_TIME_begin = OUR_TIME_end;
}

bool revlex(const exponent_vec& lhs, const exponent_vec& rhs) {
    assert(lhs.size() == rhs.size());
    for (size_t i = 1; i <= lhs.size(); ++i) {
        if (lhs[lhs.size()-i] > rhs[lhs.size()-i])
            return true;
        if (lhs[lhs.size()-i] < rhs[lhs.size()-i])
            return false;
    }
    return false; // equality
}

bool revlex_nonstrict(const exponent_vec& lhs, const exponent_vec& rhs) {
    assert(lhs.size() == rhs.size());
    for (size_t i = 1; i <= lhs.size(); ++i) {
        if (lhs[lhs.size()-i] > rhs[lhs.size()-i])
            return true;
        if (lhs[lhs.size()-i] < rhs[lhs.size()-i])
            return false;
    }
    return true; // equality, and we are doing nonstrict comparison
}

monomial_order::monomial_order(const bool t,
                               const exponent_vec& g) :
exponent_vec(g),
type(t) {}

monomial_order::monomial_order(const std::string& type_string,
                               const exponent_vec& g) :
    exponent_vec(g) {
    set_type(type_string);
}

void monomial_order::set_type(const std::string& type_string) {
    if ("deglex" == type_string) {
        type = false;
    } else if ("degrevlex" == type_string) {
        type = true;
    } else {
        std::cout << "Error: Monomial order \""
                  << type_string
                  << "\" unknown; possible values: \"deglex\", \"degrevlex\"."
                  << std::endl;
        exit(1);
    }
}

void monomial_order::set_weight(const exponent_vec& g) {
    exponent_vec::operator =(g);
}

bool monomial_order::get_type() const {
    return type;
}

std::string monomial_order::get_type_string() const {
    return (type ? "degrevlex" : "deglex");
}

exponent_vec monomial_order::get_weight() const {
    return (*this);
}

bool monomial_order::compare(const exponent_vec& lhs,
                             const exponent_vec& rhs) const {
    assert(size() == lhs.size());
    assert(size() == rhs.size());
    exponent_t wdeg_lhs(libnormaliz::v_scalar_product(*this, lhs));
    exponent_t wdeg_rhs(libnormaliz::v_scalar_product(*this, rhs));

    if (wdeg_lhs != wdeg_rhs)
        return (wdeg_lhs < wdeg_rhs);
    return (type ? revlex(lhs, rhs) : (lhs < rhs));
}

bool monomial_order::compare_nonstrict(const exponent_vec& lhs,
                                       const exponent_vec& rhs) const {
    assert(size() == lhs.size());
    assert(size() == rhs.size());
    exponent_t wdeg_lhs(libnormaliz::v_scalar_product(*this, lhs));
    exponent_t wdeg_rhs(libnormaliz::v_scalar_product(*this, rhs));
    if (wdeg_lhs != wdeg_rhs)
        return (wdeg_lhs < wdeg_rhs);
    return (type ? revlex_nonstrict(lhs, rhs) : (lhs <= rhs));
}


bool exp_vec_compare_componentwise(const exponent_vec& lhs,
                                   const exponent_vec& rhs) {
    assert(lhs.size() == rhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
        if (lhs[i] > rhs[i])
            return false;
    return true;
}


void binomial::set_mo_degrees(const monomial_order& mo) {
    mo_degree_pos = libnormaliz::v_scalar_product(mo, get_exponent_pos());
    mo_degree_neg = libnormaliz::v_scalar_product(mo, get_exponent_neg());
}

/*void binomial::compute_exponent_pos() const {
    for (size_t i = 0; i < size(); ++i)
        exponent_pos[i] = ((*this)[i] > 0 ? (*this)[i] : 0);
}
*/

exponent_vec binomial::get_exponent_pos() const {
    exponent_vec exponent_pos(size());
    for (size_t i = 0; i < size(); ++i)
        exponent_pos[i] = ((*this)[i] > 0 ? (*this)[i] : 0);
    return exponent_pos;
}

exponent_vec binomial::get_exponent_neg() const {
    exponent_vec neg_vec(size());
    for (size_t i = 0; i < size(); ++i)
        neg_vec[i] = ((*this)[i] < 0 ? -(*this)[i] : 0);
    return neg_vec;
}

void binomial::clear() {
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] = 0;
    }
    mo_degree_pos = 0;
    mo_degree_neg = 0;
}

// void binomial::compute_total_degrees() const {
//     total_degree_pos = std::accumulate(begin(), end(), 0,
//                                [](const exponent_t& e1, const exponent_t& e2)
//                                { return (0 < e2 ? e1 + e2 : e1); });
//     total_degrees_computed = true;
// }

// exponent_t binomial::get_total_degree_pos() const {
//     if (!total_degrees_computed)
//         compute_total_degrees();
//     return total_degree_pos;
// }

// exponent_t binomial::get_total_degree_neg() const {
//     exponent_vec neg_vec = get_exponent_neg();
//     return std::accumulate(neg_vec.begin(), neg_vec.end(), 0);
// }

bool binomial::operator ==(const exponent_vec& rhs) const {
    // for (size_t i = 0; i < size(); ++i)
    //     if ((*this)[i] != rhs[i])
    //         return false;
    // return true;
    return (static_cast<exponent_vec>(*this) == rhs);
}

binomial binomial::operator -(const binomial& rhs) const {
    assert(size() == rhs.size());
    binomial w(size());
    for (size_t i = 0; i < size(); ++i)
        w[i] = (*this)[i] - rhs[i];
    return w;
}

binomial binomial::operator *(const exponent_t rhs) const {
    binomial w(size());
    for (size_t i = 0; i < size(); ++i)
        w[i] = rhs * (*this)[i];
    return w;
}

void binomial::operator -=(const binomial& rhs) {
    assert(size() == rhs.size());
    for (size_t i = 0; i < size(); ++i)
        (*this)[i] -= rhs[i];
    mo_degree_pos = -1;
    mo_degree_neg = -1;
}

void binomial::operator *=(const exponent_t rhs) {
    for (size_t i = 0; i < size(); ++i)
        (*this)[i] *= rhs;
    mo_degree_pos = -1;
    mo_degree_neg = -1;
}

// Compare *this with the binomial rhs as follows:
// *this < rhs iff one of the following holds:
// (a) totdeg(exponent_pos) < totdeg(rhs.exponent_pos)
// (b) totdeg(exponent_pos) = totdeg(rhs.exponent_pos)
//     and    exponent_pos  < rhs.exponent_pos (lexicographically)
// (c) totdeg(exponent_pos) = totdeg(rhs.exponent_pos)
//     and    exponent_pos  = rhs.exponent_pos (lexicographically)
//     and totdeg(exponent_neg) < totdeg(rhs.exponent_neg)
// (d) totdeg(exponent_pos) = totdeg(rhs.exponent_pos)
//     and    exponent_pos  = rhs.exponent_pos (lexicographically)
//     and totdeg(exponent_neg) = totdeg(rhs.exponent_neg)
//     and        exponent_neg  < rhs.exponent_neg (lexicographically)
// bool binomial::operator <(const binomial& rhs) const {
//     assert(size() == rhs.size());
//     // if total degree of positive monomial decides,
//     // avoid lexicographic comparison
//     if (get_total_degree_pos() != rhs.get_total_degree_pos())
//         return (get_total_degree_pos() < rhs.get_total_degree_pos());
//     // Now we are in above case (b), (c) or (d).
//     // compare positive monomials (lexicographic comparison):
//     if (get_exponent_pos() != rhs.get_exponent_pos())
//         return (get_exponent_pos() < rhs.get_exponent_pos());
//     // if total degree of negative monomial decides,
//     // avoid lexicographic comparison:
//     if (get_total_degree_neg() != rhs.get_total_degree_neg())
//         return (get_total_degree_neg() < rhs.get_total_degree_neg());
//     // compare negative monomials:
//     return (get_exponent_neg() < rhs.get_exponent_neg());
// }

// bool binomial::operator <=(const binomial& rhs) const {
//     assert(size() == rhs.size());
//     return !(*this < rhs); // total order
// }
// Reduce nonnegative (!) exponent_vec to_reduce
// by binomial_list [begin, end[
// Returns true iff to_reduce is changed
// criterion_true is set to true iff "criterion tail" applies
// It is the responsibility of the caller to handle "criterion_true" !

// bool binomial::compare_mo(const monomial_order& mo,
//                           const binomial& rhs) const {
//     assert(size() == mo.size());
//     assert(size() == rhs.size());
//     if (get_exponent_pos() != rhs.get_exponent_pos())
//         return mo.compare(get_exponent_pos(), rhs.get_exponent_pos());
//     return mo.compare(get_exponent_neg(), rhs.get_exponent_neg());
// }

// bool binomial::compare_mo_nonstrict(const monomial_order& mo,
//                                     const binomial& rhs) const {
//     return !rhs.compare_mo(mo, *this); // total order
// }

// For monomials, | denotes divisibility.
// Compare exponent vectors by product order, i.e. component wise
// true iff v1[i] <= v2[i] for all i
bool binomial::operator |(const exponent_vec& rhs) const {
    assert(size() == rhs.size());
    // assert(rhs.nonnegative());
    assert(std::all_of(rhs.begin(), rhs.end(),
                       [](const exponent_t& e) { return (0 <= e); }));
    for (size_t i = 0; i < size(); ++i)
        if ((*this)[i] > rhs[i])
            return false;
    return true;
}

binomial binomial::lcm(const exponent_vec& rhs) const {
    assert(size() == rhs.size());
    binomial w(size());
    for (size_t i = 0; i < size(); ++i)
        w[i] = std::max((*this)[i], rhs[i]);
    return w;
}

bool binomial::zero() const {
    return std::all_of(begin(), end(),
                       [](const exponent_t& e) { return (0 == e); });
}

// bool binomial::nonnegative() const {
//     return std::all_of(begin(), end(),
//                        [](const exponent_t& e) { return (0 <= e); });
// }

bool binomial::normal(const monomial_order& mo) const {
    return (mo.compare(get_exponent_neg(), get_exponent_pos()));
}

void binomial::invert() {
    *this *= -1;
    std::swap(mo_degree_pos, mo_degree_neg);
    // total_degrees_computed = false;
}

void binomial::normalize(const monomial_order& mo) {
    if (!normal(mo))
        invert();
    // set_mo_degrees(mo);
    set_mo_degrees(mo);
}

//////////////////////////////////////////////////////////////////////////////
// S-pair criteria:

// Test if positive monomials have an indeterminate in common.
// For normalized input vectors this is one of the S-vector criterions
bool binomial::positive_coprime(const binomial& rhs) const {
    for (auto& i: pos_support_key)
        if(0 < rhs[i])
            return false;
    winf_ini_coprime++;
    return true;
}

bool binomial::criterion_tail(const binomial& rhs) const {
    for (auto& i: neg_support_key)
        if ( 0 > rhs[i]) {
            winf_tail_not_coprime++;
            return true;
        }

    return false;
}

//////////////////////////////////////////////////////////////////////////////

string binomial::to_polystring() const {
    // if (0 == size())
    //     return "";

    string ps_pos;
    string ps_neg;
    bool found_pos = false;
    bool found_neg = false;
    for (size_t i = 0; i < size(); ++i) {
        if (0 < (*this)[i]) {
            if (found_pos)
                ps_pos += "*";
            else
                found_pos = true;
            ps_pos += "x" + to_string(i + 1) + "^" + to_string((*this)[i]);
        } else if (0 > (*this)[i]) {
            if (found_neg)
                ps_neg += "*";
            else
                found_neg = true;
            ps_neg += "x" + to_string(i + 1) + "^" + to_string(-(*this)[i]);
        }
    }

    if (!(found_pos || found_neg))
        return "0";
    if (!found_pos)
        ps_pos = "1";
    if (!found_neg)
        ps_neg = "1";

    return (ps_pos + " - " + ps_neg);
}

void binomial::pretty_print(std::ostream& out) const {
    static_cast<matrix_t>(*this).pretty_print(out);
}

void binomial::set_support_keys(const dynamic_bitset& sat_support){
    neg_support_key.clear();
    pos_support_key.clear();
    for(int i = 0; i < size(); ++i){
        if((*this)[i] < 0 && sat_support[i])
            neg_support_key.push_back(i);
        if((*this)[i] > 0)
            pos_support_key.push_back(i);
    }
}

