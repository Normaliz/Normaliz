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

#ifndef LATTICE_BINOMIAL_CONTAINERS_H
#define LATTICE_BINOMIAL_CONTAINERS_H

#include <set>
#include "libnormaliz/binomial.h"

using std::vector;
using std::list;
using std::pair;
using std::cout;
using std::endl;
using std::set;

namespace libnormaliz{

class binomial_tree_node {
public:
    // constructors:
    binomial_tree_node();
    explicit binomial_tree_node(const binomial& b);
    // copy constructor:
    explicit binomial_tree_node(const binomial_tree_node& rhs);
    // destructor:
    ~binomial_tree_node();

    void pretty_print(std::ostream& out);

    // void set_data(const int val);
    // void insert(const binomial& b);

    // data fields:
    // binomial data{0};
    binomial node_binomial;
    bool has_binomial;
    vector<binomial> minimization_binomials;

    // children:
    // std::list<std::pair<size_t, binomial_tree_node*>> children;
    vector<pair < pair<size_t, exponent_t>, binomial_tree_node*> > children{};

    // vector<key_t> support_key;

    bool reduce(exponent_vec& to_reduce, bool auto_reduce);
    bool reduce_by_list(exponent_vec& to_reduce, bool auto_reduce);
    bool collect_neighbors(const exponent_vec& mon_start, const exponent_vec&  mon_goal,
                           const set<exponent_vec>& old_neighbors, set<exponent_vec>& new_neighbors);

}; // class binomial_tree_node

class binomial_tree {
public:
    // constructor:
    binomial_tree();
    binomial_tree(const monomial_order& mo,const dynamic_bitset& sat_supp);
    // destructor:
    ~binomial_tree();
    // copy constructor (to avoid "double frees"):
    binomial_tree(const binomial_tree& rhs);
    // copy assignment operator (using "copy and swap" idiom):
    binomial_tree& operator =(const binomial_tree& rhs);

    // assignment operator?
    // void set_data(const int val);

    void insert(const binomial& b);

    void clear();

    bool is_trivial() const;

    void pretty_print(std::ostream& out) const;

    bool reduce(binomial& to_reduce, bool& tail_criterion);
    bool collect_neighbors(const exponent_vec& mon_start, const exponent_vec& mon_goal,
                           const set<exponent_vec>& old_neighbors, set<exponent_vec>& new_neighbors);

    // root node:
    binomial_tree_node* root{};

    monomial_order mon_ord;
    dynamic_bitset sat_support;
    bool auto_reduce;
    bool minimization_tree;
    void set_minimization_tree();

private:
    // for copy-and-swap idiom:
    void swap(binomial_tree& rhs);
}; // class binomial_tree


class binomial_list : public std::list<binomial> {
public:
    binomial_list() = default;
    explicit binomial_list(const matrix_t& binomial_matrix);

    size_t get_number_indets() const;

    void mo_sort();

    void normalize();
    void insert_back(const binomial& b);
    void customize(binomial& b);


    void auto_reduce(binomial_tree& red_ree, const bool = false);
    template<typename Iterator>
    void intermediate_auto_reduce(binomial_tree& red_tree, Iterator& new_binom);
    void buchberger(const exponent_vec& weight_vec,
                               const bool degrevlex_mode,
                               const dynamic_bitset& sat_supp);
    // void buchberger(const monomial_order& mo);

    matrix_t to_matrix() const;
    void pretty_print(std::ostream& out,
                      const bool with_row_nr = true) const;
    std::string to_polystring() const;

    void set_degree_bound(const long deg_bound);

    template<typename Iterator>
    bool make_and_reduce_s_poly(binomial& s_poly, const Iterator match,
                                const Iterator new_binom,
                                binomial_tree& red_tree);

    // combinatorial minimization of Markov
    binomial_list  graph_minimize(bool& success); // const vector<long long>& grading);
    // minimization by degree controlled Buchberger
    binomial_list  bb_and_minimize(const vector<long long>& weight); // const vector<long long>& grading);
    vector<mpz_class> compute_HilbertSeries(const vector<long long>& given_grading);

    // binomial_list bb_and_minimize(const vector<long long>& grading, bool starting_from_GB, binomial_list& G);

    mutable monomial_order mon_ord;
    mutable dynamic_bitset sat_support;
    vector<long long> grading;
    long degree_bound;
    bool degree_bound_set;
    // bool only_monomials= false;
    // Premise: *b < *c, i.e. *b = min{*b, *c}
    template<typename Iterator>
    bool criterion_gm_left(const Iterator& b,
                           const Iterator& c) const;
    template<typename Iterator>
    bool criterion_gm_middle(const Iterator& b,
                            const Iterator& c) const;
    template<typename Iterator>
    bool criterion_gm_right(const Iterator& b,
                        const Iterator& c) const;

    void start_bb(binomial_tree& red_tree);
    void sort_by_nonzero_weight_and_normalize();

    void set_grading(const vector<long long>& grad);

    bool verbose;
    void set_verbose(bool verb);
};

class monomial_list : public std::list<exponent_vec> {
public:

    monomial_list() = default;
    monomial_list(const binomial_list& BL);
    mutable dynamic_bitset appearing_at_least_twice;

    void minimize_generating_monomials();

    monomial_list colon_by_monmial(const int& indet, const int& power);
    monomial_list add_monmial(const int& indet, const int& power) const;
    int find_pivot(int& indet) const;
    bool check_complete_intersection() const;
    vector<mpz_class> compute_HilbertSeries_inner(int level, const vector<long long>& grading);

};

class binomial_list_by_degrees : public set<pair<size_t, binomial> > {

public:

    binomial_list_by_degrees(const binomial_list& BL);
    binomial_list_by_degrees(const vector<long long>& grad);
    vector<long long> grading;
    void bin_insert(const binomial& b);
    void s_poly_insert(binomial_list& BL, binomial_tree& red_tree);

private:

    dynamic_bitset sat_support;
    monomial_order mon_ord;

};

} // namespace
#endif // include guard
