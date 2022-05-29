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

#ifndef LATTICE_BINOMIAL_H
#define LATTICE_BINOMIAL_H

#include <sys/time.h>
#include <string>
#include "libnormaliz/matrix.h"

typedef libnormaliz::dynamic_bitset dynamic_bitset;
//ypedef std::bitset<64> dynamic_bitset;

// #include <bitset>

// entries of exponent vectors:
typedef long long exponent_t;

// exponent vectors:
typedef std::vector<exponent_t> exponent_vec;

// matrices:
typedef libnormaliz::Matrix<exponent_t> matrix_t;

extern unsigned long long winf_ini_coprime;
extern unsigned long long winf_tail_not_coprime;
extern unsigned long long winf_gm_left;
extern unsigned long long winf_s_poly;
extern unsigned long long winf_red;
extern unsigned long long winf_red_tail;
extern unsigned long long winf_red_zero;
extern unsigned long long winf_red_steps;
extern unsigned long long winf_gm_steps;
extern unsigned long long winf_entered_nodes;

void reset_statistics();

bool revlex(const exponent_vec& lhs, const exponent_vec& rhs);
bool revlex_nonstrict(const exponent_vec& lhs, const exponent_vec& rhs);

class monomial_order : public exponent_vec {
public:
    monomial_order() = default;
    monomial_order(const bool t,
                   const exponent_vec& g);
    monomial_order(const std::string& type_string,
                   const exponent_vec& g);

    void set_type(const std::string& type_string);
    void set_weight(const exponent_vec& g);

    bool get_type() const; // false: deglex; true: degrevlex
    std::string get_type_string() const;
    exponent_vec get_weight() const;

    bool compare(const exponent_vec& lhs,
                 const exponent_vec& rhs) const;
    bool compare_nonstrict(const exponent_vec& lhs,
                           const exponent_vec& rhs) const;

private:
    bool type{false}; // false: deglex; true: degrevlex
};

bool exp_vec_compare_componentwise(const exponent_vec& lhs,
                                   const exponent_vec& rhs);

class binomial : public exponent_vec {
    // inherit all of exponent_vec's constructors:
    // using exponent_vec::exponent_vec;

public:
    // Constructors:
    // Construct from number of indeterminates:
    explicit binomial(const size_t num_ind) :
    exponent_vec(num_ind){}

    binomial() :
    exponent_vec(0){}

    explicit binomial(const size_t num_ind, const monomial_order& mo) :
    exponent_vec(num_ind){
        set_mo_degrees(mo);
    }

    // Construct from exponent vector:
    explicit binomial(const exponent_vec& e) :
    exponent_vec(e){}


    void set_mo_degrees(const monomial_order& mo);

    exponent_vec get_exponent_pos() const;
    exponent_vec get_exponent_neg() const;

    // exponent_t get_total_degree_pos() const;
    // exponent_t get_total_degree_neg() const;

    exponent_t get_mo_degree_pos() const {
        return mo_degree_pos;
    }
    exponent_t get_mo_degree_neg() const {
        return mo_degree_neg;
    }

    void clear();

    // Operators:
    binomial operator -(const binomial& rhs) const;
    binomial operator *(const exponent_t rhs) const; // scalar multiplication

    void operator -=(const binomial& rhs);
    void operator *=(const exponent_t rhs); // scalar multiplication

    bool operator ==(const exponent_vec& rhs) const;
    bool operator |(const exponent_vec& rhs) const;

    // General member functions:
    binomial lcm(const exponent_vec& rhs) const;

    bool zero() const;
    // bool nonnegative() const;

    bool normal(const monomial_order& mo) const;
    void invert();
    void normalize(const monomial_order& mo);

    bool positive_coprime(const binomial& rhs) const;
    bool criterion_tail(const binomial& rhs) const;
    std::string to_polystring() const;
    void pretty_print(std::ostream& out) const;

    vector<int> neg_support_key;
    vector<int> pos_support_key;
    void set_support_keys(const dynamic_bitset& sat_support);

    void compute_exponent_pos() const;

private:
    exponent_t mo_degree_pos{-1};
    exponent_t mo_degree_neg{-1};
}; // class binomial

// used for sorting a binomial_list w.r.t. "weighted deglex":
class binomial_compare_wdeglex_class {
public:
    bool operator ()(const binomial& lhs, const binomial& rhs) const {
        assert(lhs.size() == rhs.size());
        assert(-1 != lhs.get_mo_degree_pos());
        assert(-1 != lhs.get_mo_degree_neg());
        assert(-1 != rhs.get_mo_degree_pos());
        assert(-1 != rhs.get_mo_degree_neg());
        if (lhs.get_mo_degree_pos() != rhs.get_mo_degree_pos())
            return (lhs.get_mo_degree_pos() < rhs.get_mo_degree_pos()); // modeg
        if (lhs.get_exponent_pos() != rhs.get_exponent_pos())
            return (lhs.get_exponent_pos() < rhs.get_exponent_pos()); // lex

        if (lhs.get_mo_degree_neg() != rhs.get_mo_degree_neg())
            return (lhs.get_mo_degree_neg() < rhs.get_mo_degree_neg()); // modeg
        return (lhs.get_exponent_neg() < rhs.get_exponent_neg()); // lex
    }
};

// used for sorting a binomial_list w.r.t. "weighted degrevlex":
class binomial_compare_wdegrevlex_class {
public:
    bool operator ()(const binomial& lhs, const binomial& rhs) const {
        assert(lhs.size() == rhs.size());
        assert(-1 != lhs.get_mo_degree_pos());
        assert(-1 != lhs.get_mo_degree_neg());
        assert(-1 != rhs.get_mo_degree_pos());
        assert(-1 != rhs.get_mo_degree_neg());
        if (lhs.get_mo_degree_pos() != rhs.get_mo_degree_pos())
            return (lhs.get_mo_degree_pos() < rhs.get_mo_degree_pos()); // modeg
        if (lhs.get_exponent_pos() != rhs.get_exponent_pos())
            return revlex(lhs.get_exponent_pos(), rhs.get_exponent_pos());

        if (lhs.get_mo_degree_neg() != rhs.get_mo_degree_neg())
            return (lhs.get_mo_degree_neg() < rhs.get_mo_degree_neg()); // modeg
        return revlex(lhs.get_exponent_neg(), rhs.get_exponent_neg());
    }
};

#endif // include guard
