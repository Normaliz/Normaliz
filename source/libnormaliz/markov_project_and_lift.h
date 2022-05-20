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

#ifndef PROJECT_ASND_LIFT_H
#define PROJECT_ASND_LIFT_H

#include "libnormaliz/matrix.h"
#include "libnormaliz/HilbertSeries.h"
#include "libnormaliz/binomial_containers.h"

using namespace libnormaliz;

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::transform;

class groebner_project {
public:
    groebner_project() {}
    explicit groebner_project(const matrix_t& binomial_matrix,
                              const monomial_order& mo);
    explicit groebner_project(const matrix_t& binomial_matrix,
                              const exponent_vec& weight_vec,
                              const bool degrevlex_mode);
    explicit groebner_project(const matrix_t& binomial_matrix,
                              const monomial_order& mo,
                              const dynamic_bitset& sat_supp);
    explicit groebner_project(const matrix_t& binomial_matrix,
                              const exponent_vec& weight_vec,
                              const bool degrevlex_mode,
                              const dynamic_bitset& sat_supp);

    size_t get_number_indets() const;
    // lattice get_lattice() const;
    binomial_list get_binomials() const;
    monomial_order get_monomial_order() const;

    // void set_lattice(const lattice& l);
    // void set_binomials(const binomial_list& blst);
    void set_monomial_order(const monomial_order& mo);

    binomial_list get_groebner_basis() const;
    binomial_list get_minimal_markov() const;

    void write_gb() const;
    void pretty_print(std::ostream& out,
                      const bool with_row_nr = false) const;

private:
    void print_usage() const;
    std::vector<std::string> command_line{};
    std::string input_filename{};
    std::string output_filename{};

    // lattice lat{};
    binomial_list binomials{};
    monomial_order mon_ord{};
    dynamic_bitset saturation_support{};

    // mutable bool lattice_binomials_sync{false};

    mutable binomial_list gb{};
    void compute_gb() const;
    // void compute_minimal_markov() const;
    mutable bool gb_computed{false};
    mutable bool min_computed{false};
};

using namespace libnormaliz;
typedef long long Integer;

class MarkovProjectAndLift {
public:

    MarkovProjectAndLift(Matrix<Integer>& LatticeIdeal, const vector<Integer>& given_grading, const bool verb);
    Matrix<Integer> getMarkovBasis(const bool minimize_markov);

private:
    vector<Integer> grading;    
    bool is_positively_graded;
    bool verbose;

    size_t rank;
    size_t nr_vars;

    Matrix<Integer> LattiiceIdealInput; // the original argument from the constructor (witout any reordering)
    Matrix<Integer> LatticeBasis;  // as computed by row echelon reduce
    Matrix<Integer> LatticeBasisReordered; // built successively vy ading lifted columns
    Matrix<Integer> LatticeBasisTranspose;
    Matrix<Integer> LatticeBasisReorderedTranspose;
    vector<libnormaliz::key_t> StartPerm; // key vector for first permutation by L_1 norm;
    vector<libnormaliz::key_t> ColumnKey; // key vector for projection/lifting order of columns
    // First columns reordered according to StartPerm, then reordered columns added to matrix as recorded in Columnkey
    libnormaliz::dynamic_bitset TestedUnbounded; // register coordinates tested for unbounded lift
    libnormaliz::dynamic_bitset Lifted; // register lifted coordinates

    Matrix<Integer> CurrentMarkov;
    Matrix<Integer> MinimalMarkov;
    vector<Integer> CurrentWeight; // used for Buchberger
    vector<Integer> LiftedWeight; // same as CurrentWeight, but with last coordinate
    vector<Integer> PreComputedFinalGrading; // to check positively graded a priori
    libnormaliz::dynamic_bitset CurrentSatturationSupport; // for tail criterion
    bool CurrentOrder; // true for revles, false for kex

    // Note: this is a rational map with denominator denom
    // in all entries
    Matrix<Integer> TransformToTop;
    Integer denom;

    void start_column_key();
    void columns_to_old_order();
    void Compute_lift_map();
    void make_normal_form();
    void find_projection();
    void add_new_coordinate_to_Markov();
    bool compute_current_weight();
    bool lift_next_not_yet_lifted(bool allow_revlex);
    bool find_and_lift_next_unbounded();
    void lift_not_yet_lifted(bool allow_revlex);
    void lift_single_unbounded(const vector<Integer>& new_vector);
    vector<Integer> find_new_element_for_unbounded();
    void update_bookkeeping(const size_t& coord_to_lift);
    void lift_unbounded();
};

class LatticeIdeal {
    
public:

    LatticeIdeal(Matrix<Integer>& Input, const vector<Integer>& given_grading, const bool verb);
    ConeProperties compute(ConeProperties ToCompute);
    
    Matrix<Integer> getMarkovBasis();
    Matrix<Integer> getGroebnerBasis();
    HilbertSeries getHilbertSeries();


private:
    
    ConeProperties is_Computed;
    
    bool isComputed(ConeProperty::Enum prop) const;
    // returns true, when ALL properties in CheckComputed are computed
    void setComputed(ConeProperty::Enum prop);
    void setComputed(ConeProperty::Enum prop, bool value);
    
    HilbertSeries HilbSer;
    Matrix<Integer> Markov;
    Matrix<Integer> MinimalMarkov; // minimal Markov basis
    Matrix<Integer> Groebner;  

    vector<Integer> grading;
    bool is_positively_graded;
    bool verbose;

    size_t nr_vars;

};


#endif // include guard
