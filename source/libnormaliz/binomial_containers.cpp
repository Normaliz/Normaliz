/*
 * Normaliz
 * Copyright (C) 2007-2025  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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

#include <set>
#include <map>
#include <iostream>

#include "libnormaliz/binomial_containers.h"
#include "libnormaliz/list_and_map_operations.h"

namespace libnormaliz{

using std::cout;
using std::endl;
using std::pair;
using std::list;
using std::string;
using std::set;
using std::map;

bool do_test = false;

// -----------------------------------------------------
// binomial tree with winf_entered_nodes
// -----------------------------------------------------

// constructors:
binomial_tree_node::binomial_tree_node() {
    has_binomial = false;
}

/*
binomial_tree_node::binomial_tree_node(const binomial& b) {
    data.push_back(b);
    sorted = true;
}
*/


// copy constructor:
binomial_tree_node::binomial_tree_node(const binomial_tree_node& rhs) :
node_binomial(rhs.node_binomial) {
    for (auto child : rhs.children) {
        if (nullptr == child.second) {
            pair < pair<size_t, exponent_t>, binomial_tree_node*> n(child.first, nullptr);
            children.push_back(n);
        } else { // copy recursively
            binomial_tree_node* copy = new binomial_tree_node(*child.second);
            pair < pair<size_t, exponent_t>, binomial_tree_node*> n(child.first, copy);
            children.push_back(n);
        }
    }
}


// destructor:
binomial_tree_node::~binomial_tree_node() {
    for (auto child : children)
        delete child.second;
}

/*
bool binomial_tree_node::reduce_by_list(exponent_vec& to_reduce, const monomial_order mon_ord, bool auto_reduce){

        for(auto& B: data){
            winf_red_steps++;
            bool reduces = true;
            for(auto& i:support_key){
                if(to_reduce[i] < B[i]){
                    reduces = false;
                    break;
                }
            }

            if(auto_reduce){
                if(to_reduce == B.get_exponent_pos())
                    continue;
            }

            if(reduces){
                for(size_t j = 0; j < to_reduce.size(); ++j)
                    to_reduce[j] -= B[j];
                return true;
            }
        }
        return false;
}
*/

bool binomial_tree_node::reduce(exponent_vec& to_reduce, bool auto_reduce){

    winf_entered_nodes++;

    if(has_binomial){
        if(auto_reduce){
            if(to_reduce == node_binomial.get_exponent_pos()){
                return false;
            }
        }
        // cout << "TTTTTT " << to_reduce;
        // cout << "RRRRRR " << node_binomial;
        for(size_t i = 0; i < to_reduce.size(); ++i)
            to_reduce[i] -= node_binomial[i];
        winf_red_steps++;
        return true;
    }

    for(auto& C: children){
        if(to_reduce[C.first.first] >=  C.first.second && C.second->reduce(to_reduce, auto_reduce)){
            return true;
        }
    }
    return false;
}

bool binomial_tree_node::collect_neighbors(const exponent_vec& mon_start, const exponent_vec& mon_goal, const set<exponent_vec>& old_neighbors, set<exponent_vec>& new_neighbors){

    exponent_vec candidate;

    if(has_binomial){
        for(auto& min_bin: minimization_binomials){
            // cout << "In Schleife " << minimization_binomials.size() << endl;
            // if(minimization_binomials.size() > 1)
            //    cout << "In Schleife " << minimization_binomials.size() << endl;
            candidate = mon_start;
            for(size_t i = 0; i < candidate.size(); ++i){
                candidate[i] -= min_bin[i];
                assert(candidate[i] >= 0);
            }
            if(candidate == mon_goal)
                return true;
            if(old_neighbors.find(candidate) == old_neighbors.end())
                new_neighbors.insert(candidate);
        }
    }
    for(auto& C: children){
        if(mon_start[C.first.first] >=  C.first.second &&
                C.second->collect_neighbors(mon_start, mon_goal, old_neighbors, new_neighbors)){
            return true;
        }
    }
    return false;


}

void binomial_tree_node::pretty_print(std::ostream& out) {
    out << "begin node" << endl;
    node_binomial.pretty_print(cout);
    // out << "(";
    for (auto child : children) {
        if (nullptr == child.second)
            out << "nullptr";
        else {
            out << "| " << child.first.first  << " " << child.first.second << endl;
            child.second->pretty_print(out);
        }
    }
    // out << ")";
    out << "end node" << endl;
}


//////// class binomial_tree:

// constructor:
binomial_tree::binomial_tree() {
    // cout << "binomial_tree() called" << endl;
    root = new binomial_tree_node;
    root->has_binomial = false;
    minimization_tree= false;
}

binomial_tree::binomial_tree(const monomial_order& mo,const dynamic_bitset& sat_supp) {
    // cout << "binomial_tree() called" << endl;
    root = new binomial_tree_node;
    mon_ord = mo;
    sat_support = sat_supp;
    auto_reduce = false;
    minimization_tree= false;
}

void binomial_tree::set_minimization_tree(){
    minimization_tree = true;
}


// destructor:
binomial_tree::~binomial_tree() {
    // cout << "~binomial_tree() called" << endl;
    delete root;
}

// copy constructor (to avoid "double frees"):
binomial_tree::binomial_tree(const binomial_tree& rhs) {
    root = new binomial_tree_node(*(rhs.root));
}

// copy assignment operator (using "copy and swap" idiom):
binomial_tree& binomial_tree::operator =(const binomial_tree& rhs) {
    // cout << "operator =() called" << endl;
    binomial_tree copy = rhs;
    copy.swap(*this);
    return (*this);
}

// void binomial_tree::set_data(const int val) {
//     root->set_data(val);
// }

void binomial_tree::insert(const binomial& b) {
    binomial_tree_node* cur_node = root;
    for (size_t i = 0; i < b.size(); ++i) {
        if (0 < b[i]) { // only positive part matters
            size_t j = 0;
            while (   cur_node->children.size() > j
                   && (cur_node->children[j].first.first != i|| cur_node->children[j].first.second != b[i]) ) {
                ++j;
            }
            if (cur_node->children.size() > j) {
                // Child with first == (i, b[i]) exists and children[j] is that child.
                // (We have reached edge with label (i, b[i]).)
                cur_node = cur_node->children[j].second; // traverse
            } else { // There is no child with first == i yet. We create one.
                binomial_tree_node* next = new binomial_tree_node;
                cur_node->children.push_back(std::make_pair(std::make_pair(i,b[i]), next));
                cur_node = next;
                cur_node->has_binomial = false;
            }
        }
    }
    // now cur_node points to correct node (possibly newly created)
    cur_node -> has_binomial = true;
    if(!minimization_tree){
        cur_node-> node_binomial = b;
    }
    else{
        cur_node-> minimization_binomials.push_back(b);
    }
}

bool binomial_tree::reduce(binomial& to_reduce, bool& tail_criterion){

/*    if(to_reduce == test_vec){
        do_test = true;
        test_pos = to_reduce.get_exponent_pos();
        cout << "Aktiviere Test" << endl;
        to_reduce.pretty_print(cout);
        cout << "$$$$$$$$$$ " << endl;
    }
    else
        do_test = false; */

    // pos and neg are reduced separately against the tree
    exponent_vec pos(to_reduce.get_exponent_pos());
    exponent_vec neg(to_reduce.get_exponent_neg());

    exponent_vec pos_ori;
    if(auto_reduce)
        pos_ori = pos;

    tail_criterion = false;
    bool pos_changed = false;
    while(true){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        bool changed = root->reduce(pos, auto_reduce);
        if(changed)
            pos_changed = true;
        for (size_t i = 0; i < to_reduce.size(); ++i) {
            if (sat_support[i] && pos[i] != 0 && neg[i] != 0){
                tail_criterion = true;
                break;
            }
        }
        if(tail_criterion || !changed)
            break;
    }
    if(tail_criterion)
        return true;

    bool neg_changed = false;
    while(true){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        bool changed = root->reduce(neg, false);
        if(changed)
            neg_changed = true;
        for (size_t i = 0; i < to_reduce.size(); ++i) {
            if (sat_support[i] && pos[i] != 0 && neg[i] != 0){
                tail_criterion = true;
                break;
            }
        }
        if(tail_criterion || !changed)
            break;
    }
    if(tail_criterion)
        return true;

    if(neg_changed || pos_changed){
        for(size_t i = 0; i < to_reduce.size(); ++i)
            to_reduce[i] = pos[i] - neg[i];
        to_reduce.normalize(mon_ord);
    }

    return neg_changed || pos_changed;
}

bool binomial_tree::collect_neighbors(const exponent_vec& mon_start, const exponent_vec& mon_goal,
                        const set<exponent_vec>& old_neighbors, set<exponent_vec>& new_neighbors){
        return root-> collect_neighbors(mon_start, mon_goal, old_neighbors, new_neighbors);
}

void binomial_tree::clear() {
    root->node_binomial.clear();
    root->children.clear();
}

bool binomial_tree::is_trivial() const { // trivial == only root
    return std::all_of(root->children.begin(), root->children.end(),
                       [](const pair<pair<size_t, exponent_t>, binomial_tree_node*>& child) {
                           return (nullptr == child.second);
                       });
}

void binomial_tree::pretty_print(std::ostream& out) const {
    if (nullptr == root)
        out << "()";
    else
        root->pretty_print(out);
}

void binomial_tree::swap(binomial_tree& rhs) {
    // cout << "swap() called" << endl;
    std::swap(root, rhs.root);
}



// -----------------------------------------------------
//monomioal list
// -----------------------------------------------------

size_t nr_branches = 0;
size_t max_level;

void monomial_list::minimize_generating_monomials(){

    /*cout << "In MIN " << endl;
    pretty_print(cout);
    cout << "===== " << endl;*/

    if(size() <= 1)
        return;

    sort(); // ivisors precede potentail multiples
    for(auto M = begin(); M!= end(); ++M){
        for(auto N = std::next(M); N != end();){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            bool M_div_N = true;
            for(size_t k = 0; k < M->size(); ++k){
                if( (*M)[k] > (*N)[k]){
                    M_div_N = false;
                    break;
                }
            }
            if(M_div_N)
                N = erase(N);
            else
                ++N;
        }
    }
    /*cout << "Nach MIN " << endl;
    pretty_print(cout);
    cout << "===== " << endl;*/
}

monomial_list monomial_list::add_monmial(const int& indet, const int& power) const{

    monomial_list new_gen_set;
    for(auto& M: *this){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if(M[indet] < power)
            new_gen_set.push_back(M);
    }

    exponent_vec add_gen(front().size());
    add_gen[indet] = power;
    new_gen_set.push_back(binomial(add_gen));
    new_gen_set.appearing_at_least_twice = appearing_at_least_twice;
    return new_gen_set;
}

bool mon_divides(const vector<long long>& M1, const vector<long long> M2){
    for(size_t i = 0; i< M1.size(); ++i){
        if(M1[i] > M2[i])
            return false;
    }
    return true;

}

monomial_list monomial_list::colon_by_monmial(const int& indet, const int& power){

    /* monomial_list test_gen_set = *this;
    for(auto& M: test_gen_set){
        if(M[indet] > power)
            M[indet] -= power;
        else
            M[indet] = 0;
    }
    test_gen_set.appearing_vars = appearing_vars;
    test_gen_set.minimize_generating_monomials();
    // return test_gen_set;
    */



    map<int, list<list<exponent_vec>::iterator> > by_degrees;
    for(auto it = begin(); it != end(); ++it){
        by_degrees[(*it)[indet]].push_back(it);

    }
    int previous = -1;
    vector<int> degrees_prsent;
    for(auto& BD: by_degrees){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if(BD.first != previous){
            degrees_prsent.push_back(BD.first);
            previous = BD.first;
        }
    }

    for(int j = 0; j < degrees_prsent.size(); ++j){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        int m = degrees_prsent[j];
        for(auto& M: by_degrees[m]){
            if((*M)[indet] > power)
                (*M)[indet] -= power;
            else
                (*M)[indet] = 0;
        }
    }


    for(int j = 0; j < degrees_prsent.size(); ++j){
        int l  = degrees_prsent[j];
        if(l > power)
            break;
        for(auto& M1: by_degrees[l]){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            for(int k = 0; k < j; k++){
                int m = degrees_prsent[k];
                for(auto M2 = by_degrees[m].begin();
                                    M2 != by_degrees[m].end(); ){
                    if( mon_divides(*M1,*(*M2))){
                     M2 =  by_degrees[m].erase(M2);
                    }
                    else
                        M2++;
                }
            }
        }
    }

    monomial_list new_gen_set;
    new_gen_set.appearing_at_least_twice = appearing_at_least_twice;

    for(int j = 0; j < degrees_prsent.size(); ++j){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        int m = degrees_prsent[j];
        for(auto& P: by_degrees[m])
            new_gen_set.splice(new_gen_set.end(), *this, P);
    }

    new_gen_set.sort();

    /* if(test_gen_set.size() != new_gen_set.size()){
        cout << "indet " << indet << " power " << power << endl;
        cout << "TTTT " << test_gen_set.size() << " NNNN " << new_gen_set.size() << endl;
        cout << "degrees present " << degrees_prsent;
        for(auto& M: test_gen_set)
            cout << M;
        cout << "-------------" << endl;
        for(auto& M: new_gen_set)
                cout << M;
        cout << "-------------" << endl;

        for(int j = 0; j < degrees_prsent.size(); ++j){
            int m = degrees_prsent[j];
            cout << "mmmmm " << m << "ssssss " << by_degrees[m].size() << endl;
            for(auto& P: by_degrees[m])
                cout<< *P;;
            cout << "-------------" << endl;
        }


        assert(false);
    }*/


    return new_gen_set;
}

int  monomial_list::find_pivot(int& indet) const{

    if(empty()){
        return -1;
    }

    size_t N = front().size();
    int max_nr_hits = 0;
    int max_hits_indet;
    int max_appear_power = 0;
    int min_appear_power = 0;
    for(size_t k = 0; k < N; ++k){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if(!appearing_at_least_twice[k])
            continue;
        int min_power = 0;
        int max_power = 0;
        int number_hits = 0;
        for(auto& M: *this){
            if(M[k] == 0)
                continue;
            number_hits++;
            if(M[k] < min_power || min_power == 0)
                min_power = M[k];
            if(M[k] > max_power)
                max_power = M[k];
        }
        if(number_hits <= 1)
            appearing_at_least_twice[k] = false;
        if(number_hits > max_nr_hits){
            max_hits_indet = k;
            max_nr_hits = number_hits;
            min_appear_power = min_power;
            max_appear_power = max_power;

        }
    }
    if(max_nr_hits <= 1)
        return -1;

    indet = max_hits_indet;
    return (max_appear_power + min_appear_power)/2;
}

monomial_list::monomial_list(const binomial_list& BL){
    for(auto& B: BL){
        push_back(B.get_exponent_pos());
    }
    if(!BL.empty())
        appearing_at_least_twice.resize(BL.get_number_indets());
    appearing_at_least_twice.flip();
}

int level_bound_for_omp = 0;

vector<mpz_class> monomial_list::compute_HilbertSeries_inner(int level, const vector<long long>& grading){

    // cout << "LEVEL " << level << endl;

    if(level > max_level)
        max_level = level;
    nr_branches++;

    int indet;
    int d = find_pivot(indet);

    INTERRUPT_COMPUTATION_BY_EXCEPTION

    // if(check_complete_intersection()){
    if(d < 0){
        mpz_class One = 1;
        vector<mpz_class> numerator(1,One);
        for(auto& M: *this){
            long long deg = v_scalar_product(grading, M);
            mpz_class deg_mpz = convertTo<mpz_class>(deg);
            vector<mpz_class> shifted(numerator.size()+deg);
            for(size_t i = 0; i< numerator.size(); ++i)
                shifted[i+ deg] = numerator[i];
            numerator.resize(shifted.size());
            for(size_t i = 0; i< numerator.size(); ++i)
                numerator[i] = numerator[i] - shifted[i];
        }
        /* cout << "INTER " << " Level " << level << endl;
        pretty_print(cout);
        cout << "Inter num " << numerator << endl;*/
        return numerator;
    }

    //cout << "PIVOT " << indet << " -- " << d << endl;

    monomial_list sum = add_monmial(indet, d);

    // IMPORTANT: the elements of colon get spliced in from *this.
    // After this operation *this is eddentially destroyed.
    monomial_list colon = colon_by_monmial(indet, d);
    clear();
    // cout << "Vor Sum " << endl;
    // sum.pretty_print(cout);

    vector<mpz_class> numerator_sum;
    vector<mpz_class> numerator_colon;

#pragma omp parallel sections if(level <= level_bound_for_omp)
    {
#pragma omp section
    numerator_sum = sum.compute_HilbertSeries_inner(level +1,grading);
    /* cout << "Nach sum" << indet << " Power " << d  << " Level " << level << endl;
    sum.pretty_print(cout);
    cout << "SUM num " << numerator_sum;*/
#pragma omp section
    numerator_colon = colon.compute_HilbertSeries_inner(level +1,grading);
    }

    long long deg = d*grading[indet];
    /* cout << "COLON " << indet << " Power " << d << "deg " << deg  << " Level " << level << endl;
    colon.pretty_print(cout);*/
    vector<mpz_class> shifted(numerator_colon.size()+deg);
    for(size_t i = 0; i< numerator_colon.size(); ++i)
        shifted[i+ deg] = numerator_colon[i];
    // cout << "shifted " << shifted;

    int max_size = std::max(shifted.size(), numerator_sum.size());
    shifted.resize(max_size);
    numerator_sum.resize(max_size);
    for(size_t i = 0 ; i < numerator_sum.size(); ++i)
        numerator_sum[i] += shifted[i];

    /* cout << "TOTAL " << endl;
    pretty_print(cout);*/

    // cout << " total  num " << numerator_sum;

    return numerator_sum;
}

// -----------------------------------------------------
//binomioal list
// -----------------------------------------------------

binomial_list::binomial_list(const matrix_t& binomial_matrix) {
    degree_bound = -1;
    degree_bound_set = false;
    for (size_t i = 0; i < binomial_matrix.nr_of_rows(); ++i) {
        binomial bi(binomial_matrix[i]);
        push_back(bi);
    }
}

void binomial_list::set_degree_bound(const long deg_bound){
    assert(grading.size() > 0);
    degree_bound = deg_bound;
}

void binomial_list::set_grading(const vector<long long>& grad){
    grading = grad;
}

void binomial_list::set_verbose(bool verb){
    verbose = verb;
}

size_t binomial_list::get_number_indets() const {
    return (empty() ? 0 : front().size());
}

vector<mpz_class> binomial_list::compute_HilbertSeries(const vector<long long>& grad){

    grading = grad;

    // cout << grading.size() << " -- " <<grading;
    monomial_list the_monomials(*this);
    /* cout << "START" << endl;
    the_monomials.pretty_print(cout);
    cout << "$$$$$$$$$$$$$$$" << endl;*/

    int mt = omp_get_max_threads(); // must limit the nested parallelization
    while(mt > 0){                  // otherwise we risk a crash
        level_bound_for_omp++;
        mt /= 2;
    }
    level_bound_for_omp++;
    // cout << "LLLLLLLLLLL " << level_bound_for_omp << endl;

    omp_set_nested(1);
    vector<mpz_class> Num = the_monomials.compute_HilbertSeries_inner(0,grading);
    omp_set_nested(0);
    // cout << "max level " << max_level << " branches " << nr_branches << endl;
    return Num;

}

void binomial_list::mo_sort() {
    // sort(binomial_compare_class()); // to reduce against changed elements
    if (mon_ord.get_type())
        sort(binomial_compare_wdegrevlex_class());
    else
        sort(binomial_compare_wdeglex_class());
}

void binomial_list::normalize() {
    for (auto b = begin(); end() != b; ++b)
        b->normalize(mon_ord);
}

void binomial_list::customize(binomial& b) {
    b.normalize(mon_ord);
    b.set_support_keys(sat_support);
}

void binomial_list::insert_back(const binomial& b) {
    push_back(b);
    customize(back());
}

// not in use at present
template<typename Iterator>
void binomial_list::intermediate_auto_reduce(binomial_tree& red_tree, Iterator& new_binom) {
    red_tree.auto_reduce = true;

    auto b = begin();
    while (end() != b) {
        binomial b_ori(*b);
        bool tail_criterion = false;
        bool changed = red_tree.reduce(*b, tail_criterion);
        if (!changed && !tail_criterion) { // *b irreducible
            ++b;
            continue;     // nothing changed
        }
        bool zero = tail_criterion || b->zero();
        if(!zero){
            insert_back(*b);
        }
        if(b == new_binom){
            new_binom++;
        }
        b = erase(b);
    }

    red_tree.auto_reduce = false;
}

void binomial_list::sort_by_nonzero_weight_and_normalize(){
    bool weight_temporarily_added = false;
    size_t nr_vars = get_number_indets();
    exponent_vec zero_test = vector<exponent_t>(nr_vars);
    if(mon_ord == zero_test){
        weight_temporarily_added = true;
        exponent_vec total_degree = vector<exponent_t>(nr_vars,1);
        mon_ord.set_weight(total_degree);
        normalize();
        mo_sort();
    }
    if(weight_temporarily_added){
        mon_ord.set_weight(zero_test);
        normalize();
    }
    else{
        normalize();
        mo_sort();
    }
}


void binomial_list::auto_reduce(binomial_tree& red_tree, const bool initial) {
    red_tree.auto_reduce = true;
    // list<binomial> new_bins;
    bool changed;
    do {
        // mo_sort(mo); // done in buchberger
        // unique();
        changed = false;
        auto b = begin();
        while (end() != b) {

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            binomial b_ori(*b);
            bool tail_criterion = false;
            changed = red_tree.reduce(*b, tail_criterion);
            if (!changed && !tail_criterion) { // *b irreducible
                ++b;
                continue;     // nothing changed
            }
            if (tail_criterion || b->zero()) {  // reduction to zero
                b = erase(b); // leaves b pointing at next element
                continue;     // no need to set changed
            }
            changed = (b_ori != *b);
            if(!initial){ // we do not change the b_original input
                   // *b changed (and not erased)
                    b->set_support_keys(sat_support);
            }
            else
                *b = b_ori;
            ++b;
        }
    } while (changed && !initial);
    red_tree.auto_reduce = false;

    /*for(auto& B: new_bins){
        cout << "insert " << endl;
        red_tree.insert(B);
    }
    cout << "Done " << endl;*/

    //mo_sort();
    sort_by_nonzero_weight_and_normalize();
    unique();
}

// Premise: *b < *c, i.e. *b = min{*b, *c}
template<typename Iterator>
bool binomial_list::criterion_gm_left(const Iterator& b,
                        const Iterator& c) const {

    binomial lcm = c->lcm(b->get_exponent_pos());
    // cout << "LCM "; lcm.pretty_print(cout);
    for (auto it = begin(); it != b; ++it){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        winf_gm_steps++;
        bool divides = true;
        for(auto& i: it->pos_support_key){
            if((*it)[i] > lcm[i]){
                divides = false;
                break;
            }
        }
        if(divides){
            // cout << "DIV "; it->pretty_print(cout);
            winf_gm_left++;
            return true;
        }
    }
    return false;
}


template<typename Iterator>
bool binomial_list::make_and_reduce_s_poly(binomial& s_poly, const Iterator match,
                                const Iterator new_binom,
                                binomial_tree& red_tree){

    INTERRUPT_COMPUTATION_BY_EXCEPTION

    winf_s_poly++;

    /* if(match->criterion_tail(*new_binom)){
        cout << "TAIL " << endl;
    }

    if(match->positive_coprime(*new_binom)){
        cout << "COPRIME" << endl;
    }

    if(criterion_gm_left(match, new_binom)){
        cout << "GEB MÖL" << endl;
    } */

    if ( (match->criterion_tail(*new_binom)) // non-coprime tails
                || (match->positive_coprime(*new_binom)) // coprime heads
                || (criterion_gm_left(match, new_binom))  )         // GM "left"
        return true;

    s_poly = *match - *new_binom;
    if(degree_bound_set && pos_degree(s_poly, grading) > degree_bound)
        return true;

    winf_red++;
    s_poly.normalize(mon_ord);

    bool tail_criterion = false;
    red_tree.reduce(s_poly, tail_criterion);

    if(tail_criterion)
        winf_red_tail++;
    if(s_poly.zero())
        winf_red_zero ++;
    if (!tail_criterion &&  !s_poly.zero()){
        return false;
    }
    return true;
}

void binomial_list::start_bb(binomial_tree& red_tree){



    sort_by_nonzero_weight_and_normalize();
    for(auto& B: *this){
        B.set_support_keys(sat_support);
        red_tree.insert(B);
    }
    /* pretty_print(cout);
    cout << "--------------" << endl;
    red_tree.pretty_print(cout);
    cout << "--------------" << endl;
    exit(0);*/
    auto_reduce(red_tree, true); // true == initialb_ori

        if(verbose)
            verboseOutput() << "After initial auto-reduction " << size() << endl;
}


void binomial_list::buchberger(const exponent_vec& weight_vec,
                               const bool degrevlex_mode,
                               const dynamic_bitset& sat_supp) {

    mon_ord = monomial_order(degrevlex_mode, weight_vec);
    sat_support  = sat_supp;

    if(degree_bound >= 0){
        degree_bound_set = true;
        assert(grading.size() > 0);

        for(auto b = begin(); b != end(); ){
            if(pos_degree(*b, grading) > degree_bound)
                b = erase(b);
            else
                ++b;
        }
    }

    /* size_t Bind = 0;
    bool too_many = false;
    for(auto&B: *this){
        size_t Cind = 0;
        for(auto& C: *this){
            if(Cind == Bind)
                continue;
            bool does_not_divide = false;
            for(size_t i = 0; i< B.size(); ++i){
                if(C[i] <= 0)
                    continue;
                if(C[i] > B[i]){
                    does_not_divide = true;
                    break;
                }
            }
            if(!does_not_divide){
                too_many = true;
                std::cout << endl;
                B.pretty_print(std::cout);
                std::cout << "divisible " << Bind << " -- " << Cind;
                if(Bind < Cind && B == C)
                    std::cout << "equal";
                std::cout << std::endl;
                break;
            }
            Cind++;
        }
        Bind++;
    }*/

    StartTime();

    size_t inserted = 0;
    binomial_tree red_tree(mon_ord, sat_supp); // our reduction tree
    start_bb(red_tree); // inserts the input into red_tree and auto-reduces it

    size_t nr_vars = sat_supp.size();

    binomial s_poly(nr_vars);

    auto new_binom = begin(); // to be matched with the preceding binomials
    while(true){

        // cout << " new_binom "; new_binom -> pretty_print(cout);
        for(auto match = begin(); match != new_binom; ++match){
        // cout << " match "; match -> pretty_print(cout);
            bool is_zero = make_and_reduce_s_poly(s_poly,match, new_binom, red_tree);
            if(!is_zero){
                // s_poly.set_support_keys(sat_supp);
                // cout << "s_poly "; s_poly.pretty_print(cout);
                red_tree.insert(s_poly);
                insert_back(s_poly);
                inserted++;
            }
        }
        ++new_binom;
        // takes much time and does little
        /* if(inserted >= 1000){
            inserted = 0;
            intermediate_auto_reduce(red_tree, new_binom);
        }*/
        if(new_binom == end())
            break;
    }

        if(verbose)
            verboseOutput() << "Before final auto-reduction " << size() << endl;

    auto_reduce(red_tree);
    mo_sort();

    /*  binomial_list lead_mons = extract_pos_monomials();
    lead_mons.minimize_generating_monomials();
    assert(lead_mons.size() == size());*/

    /* too_many =false;
    Bind = 0;
    for(auto&B: *this){
        size_t Cind = 0;
        for(auto& C: *this){
            if(Cind == Bind)
                continue;
            bool does_not_divide = false;
            for(size_t i = 0; i< B.size(); ++i){
                if(C[i] <= 0)
                    continue;
                if(C[i] > B[i]){
                    does_not_divide = true;
                    break;
                }
            }
            if(!does_not_divide){
                too_many = true;
                std::cout << endl;
                B.pretty_print(std::cout);
                C.pretty_print(std::cout);
                std::cout << "divisible " << Bind << " -- " << Cind;
                if(Bind < Cind && B == C)
                    std::cout << "equal";
                std::cout << std::endl;
                break;
            }
            Cind++;
        }
        Bind++;
    }

    if(too_many){
        pretty_print(std::cout);
        red_tree.pretty_print(cout);
        exit(0);

    }*/

    MeasureTime(verbose, "Buchberger");
}

matrix_t binomial_list::to_matrix() const {
    matrix_t bmat(0, get_number_indets());
    for (auto b : *this)
        bmat.append(b);
    return bmat;
}

void binomial_list::pretty_print(std::ostream& out,
                                 const bool with_row_nr) const {
    to_matrix().pretty_print(out, with_row_nr);
}

string binomial_list::to_polystring() const {
    string ps;
    for (auto b = begin(); end() != b; ++b) {
        ps += b->to_polystring();
        if (end() != std::next(b))
            ps += ",\n";
    }
    return ps;
}

void s_poly_insert(binomial_list& G, binomial_list_by_degrees& B){

    if(G.size() <=1)
        return;

    binomial s_poly(G.get_number_indets());

    auto last = G.end();
    last --;
    binomial last_bin =G.back();
    last_bin.set_support_keys(G.sat_support);

    for(auto match = G.begin(); match != last; ++match){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        winf_s_poly++;

        if ( (match->criterion_tail(last_bin)) // non-coprime tails
            || (match->positive_coprime(last_bin)) // coprime heads
            || (G.criterion_gm_left(match, last))  )         // GM "left"
            continue;

        s_poly = last_bin - *match;
        if(G.degree_bound_set && pos_degree(s_poly, G.grading) > G.degree_bound)
            continue;
        s_poly.normalize(G.mon_ord);
        size_t deg = libnormaliz::v_scalar_product(B.grading, s_poly.get_exponent_pos());
        s_poly.set_support_keys(G.sat_support);
        B.insert(make_pair(deg, s_poly));
    }
}

// The two minimization routines work without degree bound presently.
// The reson is that the grading in PandL is not handled correctly.:
// its coordinates are not permuted correctly.
// Therefore we always compute the full minimal Markov basis.
// The grading used in the minimization is the lifted weight.

// minimizatiom by trying to connect lead and tail in tghe graph whose edges are defined by
// previous minimal geneartors
binomial_list binomial_list::graph_minimize(bool& success){

    assert(grading.size() > 0);
    vector<long long> weight = grading;
    StartTime();
    success = true;

    if(size() <= 1)
        return *this;

    // settings for *this
    sat_support = dynamic_bitset(weight.size());
    sat_support.flip(); // we have a positive weight and can take revlex
    mon_ord = monomial_order(true, weight);

    binomial_tree Min_red_tree(mon_ord, sat_support);
    Min_red_tree.set_minimization_tree();

    binomial_list Vmin; // minimal Markov
    binomial_list Vstopped; //returned in case of failure

    binomial_list_by_degrees  W(*this); // ordered version of *this
    set<exponent_vec> G_set;

    size_t min_degree = W.begin()->first;

    // We startb from a reduced Gröbner basis
    // Binomials of minimal degree belong to a minimal system of generators
    // For any binomial b we need b and -b in the reduction tree.
    binomial b1;
    for(auto& b: W){
        if(b.first > min_degree)
            break;
        b1 = b.second;
        Vmin.push_back(b1);
        b1.set_support_keys(sat_support);
        Min_red_tree.insert(b1);
        for(size_t i = 0; i < b1.size(); ++i)
            b1[i] = -b1[i];
        b1.set_support_keys(sat_support);
        Min_red_tree.insert(b1);
    }

    size_t nnn = 0;
    for(auto& b: W){

        INTERRUPT_COMPUTATION_BY_EXCEPTION
        nnn++;
        // cout << "********************************** " << nnn << endl;
        bool is_minimal = true;
        if(b.first == min_degree)
            continue;
        b1 = b.second;
        b1.set_support_keys(sat_support);
        exponent_vec bpos = b1.get_exponent_pos();
        exponent_vec bneg = b1.get_exponent_neg();
        set<exponent_vec> old_neighbors;
        set<exponent_vec> new_neighbors;
        // old_neighbors.insert(bpos);
        new_neighbors.insert(bpos);
        while(!new_neighbors.empty()){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            exponent_vec next = *new_neighbors.begin();
            old_neighbors.insert(next);
            if(old_neighbors.size() > 2000){
                success = false;
                MeasureTime(verbose, "graph_minimize stopped");
                return Vstopped;
            }
            /* if( old_neighbors.size() % 100000 == 0)
                cout << "OOOOO " << old_neighbors.size() << endl;*/
            new_neighbors.erase(next);
            if(Min_red_tree.collect_neighbors(next, bneg, old_neighbors, new_neighbors)){
                is_minimal =false;
                break;
            }
            if(new_neighbors.size() > 2000){
                MeasureTime(verbose, "graph_minimize stopped");
                success = false;
                return Vstopped;
            }
        }
        if(is_minimal){
                Vmin.push_back(b1);
            b1.set_support_keys(sat_support);
            Min_red_tree.insert(b1);
            for(size_t i = 0; i < b1.size(); ++i)
                b1[i] = -b1[i];
            b1.set_support_keys(sat_support);
            Min_red_tree.insert(b1);
        }
    }

    MeasureTime(verbose, "graph_minimize");

    return Vmin;
}


// realizes the algorithm in Kreuzer-Robbiano CCA II, Theorem 4.6.7
// notation as used there
// not used anymore
 binomial_list binomial_list::bb_and_minimize(const vector<long long>& weight){
 //binomial_list binomial_list::bb_and_minimize(const vector<long long>& weight, bool starting_from_GB, binomial_list& G){

    StartTime();

    // assert(G.empty());

    if(size() <= 1)
        return *this;

    bool starting_from_GB = true;

    // settings for *this
    sat_support = dynamic_bitset(weight.size());
    sat_support.flip(); // we have a positive weight and can take revlex
    mon_ord = monomial_order(true, weight);
    normalize();

    binomial_list G; // the GB built by degrees
    G.mon_ord = mon_ord;
    G.sat_support = sat_support;
    G.grading = grading;
    G.degree_bound = degree_bound;
    G.degree_bound_set = false;
    // The following is not used at present
    //see comment preceding graph_minimize
    if(degree_bound >= 0){
        G.degree_bound_set = true;
        // cout << "RRRRRRRR " << grading;
    }

    binomial_tree G_red_tree(mon_ord, sat_support);

    binomial_list Vmin; // minimal Markov

    binomial_list_by_degrees  W(*this); // ordered version of *this
    binomial_list_by_degrees B(weight);
    set<exponent_vec> G_set;

    long long min_degree;

    while(!W.empty()){

        if(B.empty()){
            min_degree = W.begin()->first;
        }
        if(!B.empty()){
            min_degree = std::min(W.begin()->first,B.begin()->first);
        }

        binomial b;

        while(!B.empty() && min_degree == B.begin()->first){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            b = B.begin()->second;
            B.erase(B.begin());
            b.set_support_keys(sat_support);

            bool tail_criterion = false;
            G_red_tree.reduce(b, tail_criterion);
            if(tail_criterion)
                winf_red_tail++;
            if(!tail_criterion && b.zero())
                winf_red_zero ++;
            if (tail_criterion || b.zero())
                continue;
            G.insert_back(b);
            G_red_tree.insert(b);
            G_set.insert(b.get_exponent_pos());
            s_poly_insert(G, B);
        }

        while(!W.empty() && min_degree == W.begin()->first){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            b = W.begin()->second;
            W.erase(W.begin());

            if(starting_from_GB && G_set.find(b.get_exponent_pos())!= G_set.end())
                continue;
            G.insert_back(b);
            G_red_tree.insert(b);
            Vmin.push_back(b);
            s_poly_insert(G, B);
        }
    }

    MeasureTime(verbose, "bb_and_minimize");

    return Vmin;
}

// -----------------------------------------------------
// binomial list by degrees
// -----------------------------------------------------

binomial_list_by_degrees::binomial_list_by_degrees(const binomial_list& BL){

    grading = BL.mon_ord.get_weight();  // we want to keep the order in BL as much as possible
    vector<long long> bounding_grad = BL.grading;
    long long degree_bound = BL.degree_bound;

    bool no_degree_bound = true;
    if(degree_bound > -1)
        no_degree_bound = false;
    if(!BL.empty())
        assert(grading.size() == BL.front().size());
    for(auto& b: BL){
        if(no_degree_bound || pos_degree(b, bounding_grad) <= degree_bound)
            bin_insert(b);
    }
}


binomial_list_by_degrees::binomial_list_by_degrees(const vector<long long>& grad){

    grading = grad;
}

void binomial_list_by_degrees::bin_insert(const binomial& b){
    int deg = libnormaliz::v_scalar_product(grading,b.get_exponent_pos());
    insert(make_pair(deg,b));
}

} // namespace


