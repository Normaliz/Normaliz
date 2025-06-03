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

#include <fstream>

#include "libnormaliz/lattice_ideal.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/list_and_map_operations.h"
#include "libnormaliz/project_and_lift.h"

namespace  libnormaliz{

typedef long long Integer;
typedef mpz_class BigInt;

Integer find_nonzero_degree(const Matrix<Integer>& M,
                                       const vector<Integer>& grading, const long min_degree){
    bool first = true;
    Integer degree_found = -1;
    for(size_t i = 0; i < M.nr_of_rows(); ++i){
        Integer PD = pos_degree(M[i], grading);
        if(PD < min_degree)
            continue;
        if(first || PD < degree_found){
            first = false;
            degree_found = PD;
        }
    }
    return degree_found;
}

void sort_by_pos_degree(Matrix<Integer>& M, const vector<Integer>& grading){
    list<pair<long long, size_t> > to_be_sorted;
    for(size_t i = 0; i < M.nr_of_rows(); ++i){
        // cout << i << "  " << pos_degree(M[i], grading) << endl;
        to_be_sorted.push_back(make_pair(pos_degree(M[i], grading),i));
    }
    to_be_sorted.sort();
    vector<key_t> perm;
    for(auto& s:to_be_sorted)
        perm.push_back(s.second);
    // cout << perm;
    M.order_rows_by_perm(perm);
}


// degree_bound = - 2: find first degree with elements of degree >= min_degree
// not used at present
Matrix<Integer> select_by_degree(const Matrix<Integer>& M,
                                       const vector<Integer>& grading, long degree_bound, const long min_degree){
    if(degree_bound == -2){
        degree_bound = find_nonzero_degree(M, grading, min_degree);
    }
    vector<key_t> satisfies_degree_bound;
    for(size_t i = 0; i < M.nr_of_rows(); ++i){
        Integer PD = pos_degree(M[i], grading);
        if(degree_bound != -1 && PD > degree_bound)
            continue;
        if(PD  < min_degree)
            continue;
        satisfies_degree_bound.push_back(i);
    }
    return M.submatrix(satisfies_degree_bound);
}


//---------------------------------------------------------------
//         MarkovProjectAndLift
//---------------------------------------------------------------

MarkovProjectAndLift::MarkovProjectAndLift(Matrix<Integer>& LatticeIdeal, const bool verb){

    verbose = verb;
    // degree_bound = -1;

    // cout << "Given lattice ideal in Laurent polynomial ring" << endl;
    // LatticeIdeal.pretty_print(cout);
    LattiiceIdealInput = LatticeIdeal;

    Matrix<Integer> LItranspose = LatticeIdeal.transpose();
    Matrix<Integer> Weights(0,LItranspose.nr_of_columns());
    Weights.append(vector<Integer> (LItranspose.nr_of_columns(),1));
    vector<bool> absolute(1,1);

    StartPerm = LItranspose.perm_by_weights(Weights, absolute);
    LItranspose.order_rows_by_perm(StartPerm);
    if(verbose){
        verboseOutput() << "---------------------------------------------------" << endl;
        verboseOutput() << "Starting project-and-lift for Markov basis" << endl << endl;
        verboseOutput() << "Columns reordered "<< StartPerm << endl;
    }
    LatticeIdeal = LItranspose.transpose();

    // LatticeIdeal.pretty_print(cout);

    LatticeBasis = LatticeIdeal;
    nr_vars = LatticeBasis.nr_of_columns();
    rank = LatticeBasis.row_echelon_reduce();
    // cout << "Row echelon form of lattice basis" << endl;
    //LatticeBasis.debug_print();
    // cout << "rank " << rank << endl;
    LatticeBasis.resize(rank);
    start_column_key();
    make_normal_form();
    Compute_lift_map();
    /* cout << "Transformation to top" << endl;
    TransformToTop.pretty_print(cout);
    cout << "denom " << denom << endl; */
}

void MarkovProjectAndLift::start_column_key(){

    Lifted.resize(nr_vars);
    TestedUnbounded.resize(nr_vars);

    for(size_t i = 0; i < rank; ++i){
        for(size_t j = 0; j < nr_vars; ++j){
            if(LatticeBasis[i][j] != 0){
                ColumnKey.push_back(j);
                Lifted[j] = true;
                TestedUnbounded[j] = true;
                break;
            }
        }
    }

    if(verbose){
        verboseOutput() << "Projection to new coordinates" << endl;
    verboseOutput() << ColumnKey;;
    }
}

// makes kind of a Hermite normal form
void MarkovProjectAndLift::make_normal_form(){

    for(size_t i=1; i < rank; ++i){
        for(size_t j = 0; j < i; ++j){
            if(LatticeBasis[j][ColumnKey[i]] <= 0)
                continue;
            Integer fact = LatticeBasis[j][ColumnKey[i]]/LatticeBasis[i][ColumnKey[i]];
            if(LatticeBasis[j][ColumnKey[i]] % LatticeBasis[i][ColumnKey[i]] != 0)
                fact++;
            for(size_t k = i; k < nr_vars; ++k)
                LatticeBasis[j][k] -= fact* LatticeBasis[i][k];
        }
    }
    // LatticeBasis.debug_print();
    // cout << "Lattice basis after transformation " << endl;
    // LatticeBasis.pretty_print(cout);

    LatticeBasisTranspose = LatticeBasis.transpose();
}

void MarkovProjectAndLift::columns_to_old_order(){

    Matrix<Integer> Copy = CurrentMarkov;
    for(size_t i = 0; i < Copy.nr_of_rows(); ++i){
        for(size_t j = 0; j < nr_vars; ++j)
            CurrentMarkov[i] [ StartPerm[ ColumnKey[j]] ] = Copy[i][j];
    }

    if(MinimalMarkov.nr_of_rows() == 0){
        MinimalMarkov.resize(0, nr_vars); // should have correc t format
        return;
    }

    Copy = MinimalMarkov;
    for(size_t i = 0; i < Copy.nr_of_rows(); ++i){
        for(size_t j = 0; j < nr_vars; ++j)
            MinimalMarkov[i] [ StartPerm[ ColumnKey[j]] ] = Copy[i][j];
    }
}

void MarkovProjectAndLift::Compute_lift_map(){

    Matrix<Integer> rxr(rank,rank);
    for(size_t i = 0; i< rank; ++i)
        for(size_t j = 0; j < rank; ++j)
        rxr[i][j] = LatticeBasis[i][ColumnKey[j]];

    TransformToTop = rxr.solve(LatticeBasis, denom);

    /* Matrix<Integer> Check = rxr.multiplication(TransformToTop);
    Check.scalar_division(denom);
    cout << "========================7" << endl;
    Check.pretty_print(cout);*/

}

void MarkovProjectAndLift::add_new_coordinate_to_Markov(){

    size_t new_coord = LatticeBasisReorderedTranspose.nr_of_rows() - 1;

    vector<Integer> new_column(CurrentMarkov.nr_of_rows());
    for(size_t i = 0; i < new_column.size(); ++i){
        Integer new_entry = 0;
        for(size_t j = 0; j < rank; ++j){
            new_entry += CurrentMarkov[i][j]*TransformToTop[j][ColumnKey[new_coord]];
        }
        new_column[i] = new_entry/denom;
    }
    CurrentMarkov.insert_column(new_coord, new_column);
}

bool MarkovProjectAndLift::compute_current_weight(){

    /* cout << "Truncated lattice basis for coordinate " << new_coord << endl;
    TruncBasis.pretty_print(cout); */

    size_t new_coord = LatticeBasisReordered.nr_of_columns()-1;

    // bool save_global_verbose = libnormaliz::verbose;
    // libnormaliz::verbose = false;
    Matrix<BigInt> LBR_Big;
    convert(LBR_Big, LatticeBasisReordered);
    suppressNextConstructorVerbose();
    Cone<BigInt> WeightCone(Type::equations, LBR_Big); // intersects with positive orthant
    WeightCone.setVerbose(false);
    // libnormaliz::verbose = save_global_verbose;
    Matrix<BigInt> ER_big = WeightCone.getExtremeRaysMatrix();
    Matrix<Integer> ExtRays;
    convert(ExtRays, ER_big);
    vector<Integer> GradingOnCurrentQuotient(new_coord+1,0);
    CurrentWeight = vector<Integer>(new_coord+1,0);
    for(size_t i = 0; i < ExtRays.nr_of_rows(); ++i){
        CurrentWeight = v_add(CurrentWeight, ExtRays[i]);
        if(ExtRays[i].back() == 0)
            GradingOnCurrentQuotient = v_add(GradingOnCurrentQuotient,ExtRays[i]);
    }
    v_make_prime(CurrentWeight);
        bool good_for_bounded = (CurrentWeight.back() > 0);
    LiftedWeight = CurrentWeight;
    CurrentWeight.resize(new_coord); // shotened by lat coordinate

    GradingOnCurrentQuotient.resize(new_coord);

    // cout << "Grading on current quotient for saturation support" << endl;
    // cout << GradingOnCurrentQuotient;

    // cout << "Current weight for monomial order" << endl;
    // cout << CurrentWeight;

    CurrentSatturationSupport.resize(new_coord);
    for(size_t i = 0; i < new_coord; ++i){
        if(GradingOnCurrentQuotient[i] > 0)
            CurrentSatturationSupport[i] = true;
        else
            CurrentSatturationSupport[i] = false; // to be on the safe side
    }

    // must remove the last entry for Buchberger
    return good_for_bounded;
}

void MarkovProjectAndLift::update_bookkeeping(const size_t& coord_to_lift){

    Lifted[coord_to_lift] = true;
    ColumnKey.push_back(coord_to_lift);
    LatticeBasisReordered.append_column(LatticeBasisTranspose[coord_to_lift]);
    LatticeBasisReorderedTranspose.append(LatticeBasisTranspose[coord_to_lift]);
}

bool MarkovProjectAndLift::lift_next_not_yet_lifted(bool allow_revlex){

    dynamic_bitset NotLifted = ~Lifted;
    if(!NotLifted.any())
        return false;

    update_bookkeeping(NotLifted.find_first());

    bool good_for_bounded = compute_current_weight();
    if(!good_for_bounded){ // can happen in "straight" mode
        lift_single_unbounded(vector<Integer>());
        return true;
    }
    if(verbose)
        verboseOutput() << "Lift step " << ColumnKey.size() - 1 << " bounded to sorted coordinate "
                << ColumnKey.back() <<  ", original coordinate "
                << StartPerm[ ColumnKey.back() ] << endl;

    bool full_support_of_weight = allow_revlex;
    if(allow_revlex){
        for(size_t i = 0; i< CurrentWeight.size(); ++ i){
            if(CurrentWeight[i] == 0){
                full_support_of_weight = false;
                break;
            }
        }
    }
    CurrentOrder = full_support_of_weight;
    binomial_list grp(CurrentMarkov);
    grp.set_verbose(verbose);
    // cout << CurrentWeight; // *****************
    grp.buchberger(CurrentWeight, full_support_of_weight, CurrentSatturationSupport);
    CurrentMarkov = grp.to_matrix();
    if(verbose)
        verboseOutput() << "Size of current Markov after Buchberger " << CurrentMarkov.nr_of_rows() << endl;
    add_new_coordinate_to_Markov();

    if(verbose)
        verboseOutput() << "Dim reached " << CurrentMarkov.nr_of_columns() << endl;
    if(verbose)
        verboseOutput() << "---------------------------------------------------" << endl;


    if(CurrentMarkov.nr_of_columns() < nr_vars)
        return true;

    for(size_t i = 0; i < LiftedWeight.size(); ++i){
        if(LiftedWeight[i] == 0)
            return true;
    }


    // now full dimension and positively graded
    // can minimize

    if(verbose)
        verboseOutput() << "Computing minimal Markov basis" << endl;
    binomial_list gr(CurrentMarkov);
    gr.set_verbose(verbose);
    bool graph_success;
    // if(degree_bound >= 0){
        // gr.set_grading(grading);
        // gr.set_degree_bound(degree_bound); // only by selection now
    // }
    // else{
        gr.set_grading(LiftedWeight);
    // }
    binomial_list min_markov = gr.graph_minimize(graph_success);
    if(!graph_success){
        min_markov = gr.bb_and_minimize(LiftedWeight);
    }
    MinimalMarkov = min_markov.to_matrix();

    if(verbose)
        verboseOutput() << "Size of minimal Markov basis " << MinimalMarkov.nr_of_rows() << endl;

    if(verbose)
        verboseOutput() << "---------------------------------------------------" << endl;

    return true;
}

bool MarkovProjectAndLift::find_and_lift_next_unbounded(){

    dynamic_bitset NotLifted = ~TestedUnbounded;
    if(!NotLifted.any())
        return false;

    size_t first_coord_to_test = NotLifted.find_first();

    //  Find  cone of coefficient vectors that give nonnegative linear combination of alltice basis
    // within the already done columns
    // bool save_global_verbose = libnormaliz::verbose;
    // libnormaliz::verbose = false;
    Matrix<BigInt> LBRT_Big;
    convert(LBRT_Big, LatticeBasisReorderedTranspose);
    suppressNextConstructorVerbose();
    Cone<BigInt> CheckBounded(Type::inequalities, LBRT_Big);  // TODO Use Normaliz dynamic -- we must add
    //  one inequality only to last inequalities
    CheckBounded.setVerbose(false);
    //libnormaliz::verbose = save_global_verbose;
    Matrix<BigInt> ER_big = CheckBounded.getExtremeRaysMatrix();
    Matrix<Integer> ExtRays;
    convert(ExtRays, ER_big);

    // Now find next column that admits positive value under one of the extreme rays
    size_t good_ext_ray = ExtRays.nr_of_rows();
    size_t new_column;

    for(size_t i = first_coord_to_test; i < nr_vars; ++i){ // we need one extreme ray that gives positive value on next column

        if(Lifted[i])
            continue;
        TestedUnbounded[i] = true;

        if(verbose)
            verboseOutput() << "checking coordinate " << i << endl;
        for(size_t k = 0; k < ExtRays.nr_of_rows(); ++k){
            if(v_scalar_product(ExtRays[k], LatticeBasisTranspose[i]) > 0){
                good_ext_ray = k;
                break;
            }
        }
        if(good_ext_ray < ExtRays.nr_of_rows()){
            new_column = i;
            break;
        }
    }

    if(good_ext_ray == ExtRays.nr_of_rows()){ // no unbounded lift possible
        return false;
    }

    update_bookkeeping(new_column);

    if(verbose)
        verboseOutput() << "Lift step " << ColumnKey.size() - 1 <<  " un-bounded to sorted coordinate "
                << ColumnKey.back()  <<  ", original coordinate " << StartPerm[ ColumnKey.back() ] << endl;

    vector<Integer> new_vector = LatticeBasisReorderedTranspose.MxV(ExtRays[good_ext_ray]);
    lift_single_unbounded(new_vector);

    return true;
}

// we need a second method for finding the new Markov element for unbounded lifting
// if this comes out of the blue (and not from find_and_lift_next_unbounded)
vector<Integer> MarkovProjectAndLift::find_new_element_for_unbounded(){

    Matrix<BigInt> UnitMat(LatticeBasisReordered.nr_of_columns());
    Matrix<BigInt> LBR_Big;
    convert(LBR_Big, LatticeBasisReordered);
    suppressNextConstructorVerbose();
    Cone<BigInt> WeightCone(Type::cone, LBR_Big, Type::inequalities, UnitMat);
    WeightCone.setVerbose(false);
    Matrix<BigInt> ER_big = WeightCone.getExtremeRaysMatrix();
    Matrix<Integer> ExtRays;
    convert(ExtRays, ER_big);
    assert(ExtRays.nr_of_rows()> 0);
    size_t good_ext_ray = ExtRays.nr_of_rows();
    for(size_t i=0; i < ExtRays.nr_of_rows(); ++i){
        if(ExtRays[i].back() > 0){
            good_ext_ray = i;
            break;
        }
    }
    assert(good_ext_ray < ExtRays.nr_of_rows());

    return ExtRays[good_ext_ray];
}

void MarkovProjectAndLift::lift_single_unbounded(const vector<Integer>& new_vector){

    // add new coordinate to CurrentMarkov
    add_new_coordinate_to_Markov();

    vector<Integer> vector_to_add;

    if(new_vector.size() > 0){
        vector_to_add = new_vector;
    }
    else{
        vector_to_add = find_new_element_for_unbounded();
    }

    // extend CurrentMarkov by new vector
    CurrentMarkov.append(vector_to_add);

    if(verbose)
        verboseOutput()<< "Size of current Markov after unbounded lift " << CurrentMarkov.nr_of_rows() << endl;
    if(verbose)
        verboseOutput() << "---------------------------------------------------" << endl;
}

void MarkovProjectAndLift::lift_unbounded(){

    if(verbose)
        verboseOutput() << "searching unbounded coordinates" << endl;

    while(find_and_lift_next_unbounded()){
    }
}

void MarkovProjectAndLift::lift_not_yet_lifted(bool allow_revlex){

    while(lift_next_not_yet_lifted(allow_revlex)){
    }
}

void MarkovProjectAndLift::find_projection(){

    bool diagonal_is_ones = true;

    for(size_t i = 0; i< rank; ++i){  // TODO doe we need this condition?
        if(LatticeBasis[i][ColumnKey[i]] != 1){
            diagonal_is_ones = false;
            break;
        }
    }

    vector<int> ExtensionKey;

    if(diagonal_is_ones){
        for(size_t j = 0; j < nr_vars; ++j){
            bool column_non_positive = true;
            for(size_t i = 0; i < rank; ++i){
                if(LatticeBasis[i][j] > 0){
                    column_non_positive = false;
                    break;
                }
            }
            if(column_non_positive){
                ColumnKey.push_back(j);
                Lifted[j] = true;
                TestedUnbounded[j] = true;
                ExtensionKey.push_back(j);
            }
        }
    }

    if(ExtensionKey.size() > 0)
        if(verbose)
            verboseOutput() << "Extending projection to new coordinates " << ExtensionKey;

    CurrentMarkov = LatticeBasisTranspose.submatrix(ColumnKey).transpose();
    LatticeBasisReordered = CurrentMarkov;
    LatticeBasisReorderedTranspose = LatticeBasisReordered.transpose();

    // cout << "Markov basis at projrection " << endl;
    // CurrentMarkov.pretty_print(cout);
}



void MarkovProjectAndLift::compute(Matrix<long long>& Mark, Matrix<long long>& MinMark){

    find_projection();
    lift_unbounded(); // straight no longer used
    lift_not_yet_lifted(true); // revlex allowed
    columns_to_old_order();

    // cout << "Pand L " << CurrentMarkov.nr_of_rows() << " -- " << MinimalMarkov.nr_of_rows() << endl;

    swap(CurrentMarkov, Mark);
    swap(MinimalMarkov, MinMark);
}

/*
void MarkovProjectAndLift::set_degree_bound(const long deg_bound) {
    assert(grading.size() > 0);
    degree_bound = deg_bound;
}


void MarkovProjectAndLift::set_grading(const vector<long long>& grad) {
    grading = grad;
}
*/

//---------------------------------------------------------------
// lattice ideal
//---------------------------------------------------------------

LatticeIdeal::LatticeIdeal(const Matrix<long long>& Input, const vector<Integer>& given_grading, const bool verb){
    verbose = verb;
    Grading = given_grading;
    OurInput = Input;
    is_positively_graded = false;
    nr_vars = Input.nr_of_columns();
    degree_bound= -1;
    min_degree = -1;
}

bool LatticeIdeal::isComputed(ConeProperty::Enum prop) const {
    return is_Computed.test(prop);
}

void LatticeIdeal::set_degree_bound(const long deg_bound) {
    assert(Grading.size() > 0); // make sonly sense with grading
    degree_bound = deg_bound;
    setComputed(ConeProperty::MarkovBasis, false);
    setComputed(ConeProperty::GroebnerBasis, false);
}

void LatticeIdeal::set_gb_weight(const vector<long long>& given_weight){
    gb_weight = given_weight;
}

void LatticeIdeal::set_min_degree(const long deg) {
    assert(Grading.size() > 0); // make sonly sense with grading
    min_degree = deg;
    setComputed(ConeProperty::MarkovBasis, false);
    setComputed(ConeProperty::GroebnerBasis, false);
}

void LatticeIdeal::setComputed(ConeProperty::Enum prop) {
    is_Computed.set(prop);
}

void LatticeIdeal::setComputed(ConeProperty::Enum prop, bool value) {
    is_Computed.set(prop, value);
}

Matrix<Integer>  LatticeIdeal::getMarkovBasis(){
    if(!isComputed(ConeProperty::MarkovBasis))
        compute(ConeProperty::MarkovBasis);
    vector<long long> OurDegrees;
    for(size_t i = 0; i< MinimalMarkov.nr_of_rows(); ++i){
        OurDegrees.push_back(pos_degree(MinimalMarkov[i], Grading));
    }
    if(verbose){
        map<long long, size_t> DegMap = count_in_map<long long, size_t>(OurDegrees);
        verboseOutput() << "Grading " << Grading;
        verboseOutput() << "Degrees in minimal Markovbasis " << DegMap;
        verboseOutput() << "---------------------------------------------------" << endl;
    }

    if(MinimalMarkov.nr_of_rows() >0 ){
        if(degree_bound >= 0 || min_degree >= 0){
            sort_by_pos_degree(MinimalMarkov, Grading);
            if(verbose){
                verboseOutput() << "Selecting max deg " << degree_bound << " min deg " << min_degree << endl;
            }
            return select_by_degree(MinimalMarkov, Grading, degree_bound, min_degree);
        }
        else
            return MinimalMarkov;
    }
    else
        return Markov;
}

Matrix<Integer>  LatticeIdeal::getGroebnerBasis(){
    if(!isComputed(ConeProperty::GroebnerBasis))
        compute(ConeProperty::GroebnerBasis);
    if(degree_bound >= 0  || min_degree >= 0){
        sort_by_pos_degree(Groebner, Grading);
        return select_by_degree(Groebner, Grading, degree_bound, min_degree);
    }
    else
        return Groebner;
}

HilbertSeries  LatticeIdeal::getHilbertSeries(){
    if(!isComputed(ConeProperty::HilbertSeries))
        compute(ConeProperty::HilbertSeries);
    return HilbSer;
}

void LatticeIdeal::computeMarkov(){

    // cout << "IIIIIIIIIIIIIIII " << OurInput.nr_of_rows() << endl;
    // OurInput.debug_print();

    MarkovProjectAndLift PandL(OurInput, verbose);
    // if(Grading.size() > 0 && degree_bound != -1){
        // PandL.set_grading(Grading);
        // cout << "BBBBBBBBBBBBB " << degree_bound << endl;
        //PandL.set_degree_bound(degree_bound);
    // }
    PandL.compute(Markov, MinimalMarkov);
    if(MinimalMarkov.nr_of_rows() > 0){
        is_positively_graded = true;
    }
    // cout << "Mark " << Markov.nr_of_rows() << " MinMark " << MinimalMarkov.nr_of_rows() << endl;
    // Markov.pretty_print(cout);
}

void LatticeIdeal::computeGroebner(ConeProperties ToCompute){

    // cout << "GRÖBNER " << ToCompute << endl;

    string FinalGB = "RevLex";
    vector<Integer> all_one(Markov.nr_of_columns(),1);
    if(gb_weight.size() > 0){
        all_one = gb_weight;
        FinalGB = "weighted " + FinalGB;
    }
    bool use_rev_lex = true;

    if(ToCompute.test(ConeProperty::Lex)){
        FinalGB = "Lex";
        use_rev_lex = false;
        all_one = vector<Integer> (nr_vars,0);
        if(gb_weight.size() > 0){
            all_one = gb_weight;
          FinalGB = "weighted " + FinalGB;
        }
    }
    if(ToCompute.test(ConeProperty::DegLex)){
        use_rev_lex = false;
        FinalGB = "Deglex";
    }
    if(verbose)
        verboseOutput()<< "Final Gröbner basis " << FinalGB << endl;

    dynamic_bitset CurrentSatturationSupport(nr_vars);
    if(is_positively_graded)
        CurrentSatturationSupport.flip();

    // cout << CurrentSatturationSupport.size() << "   " << CurrentSatturationSupport  << endl;
    reset_statistics();

    binomial_list grp(Markov);
    grp.set_verbose(verbose);
    if(degree_bound != -1){ // so far no effect
        assert(Grading.size() > 0);
        grp.set_grading(Grading);
        grp.set_degree_bound(degree_bound);
    }
    grp.buchberger(all_one, use_rev_lex, CurrentSatturationSupport);

    Groebner = grp.to_matrix();

    // Groebner = select_by_degree(Groebner, Grading, degree_bound, min_degree);
    if(verbose)
        verboseOutput() << "Gröbner basis elements " << Groebner.nr_of_rows() << endl;
    // cout << "GGGGGG " << Groebner.nr_of_rows() << endl;
    if(verbose)
        verboseOutput()<<"---------------------------------------------------" << endl;
}

void LatticeIdeal::computeHilbertSeries(const ConeProperties& ToCompute){

    assert(degree_bound == -1);
    assert(Grading.size() > 0);

    StartTime();
    // cout << "Final quotient psoitively graded" << endl;
    binomial_list bl_HilbertSeries(Markov);
    vector<mpz_class> numerator = bl_HilbertSeries.compute_HilbertSeries(Grading);
    vector<long> Grading_long;
    convert(Grading_long, Grading);
    long save_nr_coeff_quasipol = HilbSer.get_nr_coeff_quasipol();
    long save_expansion_degree = HilbSer.get_expansion_degree();
    HilbSer = HilbertSeries(numerator, Grading_long);
    HilbSer.set_expansion_degree(save_expansion_degree);
    HilbSer.set_nr_coeff_quasipol(save_nr_coeff_quasipol);
    if(ToCompute.test(ConeProperty::OnlyCyclotomicHilbSer)){
        HilbSer.set_only_cyclotomic(true);
        HilbSer.forbid_quasipol(true);
    }
    if(ToCompute.test(ConeProperty::NoQuasiPolynomial)){
        HilbSer.forbid_quasipol(true);
    }
    HilbSer.simplify();
    /* if(verbose){
        verboseOutput() << "Hilbert series numerator  " << HilbSer.getNum();
        verboseOutput() << "Hilbert series denominator " <<  HilbSer.getDenom();
    }*/
    MeasureTime(verbose, "Hilbert series");

    if(verbose)
        verboseOutput() << "---------------------------------------------------" << endl;
}

ConeProperties LatticeIdeal::compute(ConeProperties ToCompute){

    ToCompute.reset(is_Computed);
    if(!ToCompute.any())
        return ToCompute;

    if(ToCompute.test(ConeProperty::HilbertSeries))
        ToCompute.set(ConeProperty::MarkovBasis);

    if(ToCompute.test(ConeProperty::GroebnerBasis))
        ToCompute.set(ConeProperty::MarkovBasis);

    ToCompute.reset(is_Computed);
    if(!ToCompute.any())
        return ToCompute;

    if(ToCompute.test(ConeProperty::MarkovBasis)){
        computeMarkov();
        setComputed(ConeProperty::MarkovBasis);
        ToCompute.reset(is_Computed);
    }

    if(ToCompute.test(ConeProperty::GroebnerBasis)){
        computeGroebner(ToCompute);
        setComputed(ConeProperty::GroebnerBasis);
        ToCompute.reset(is_Computed);
    }

    if(ToCompute.test(ConeProperty::HilbertSeries)){
        computeHilbertSeries(ToCompute);
        setComputed(ConeProperty::HilbertSeries);
        ToCompute.reset(is_Computed);
    }

    return ToCompute;
}

//----------------------------------------------------------------
// HilbertBasisMonoid
//----------------------------------------------------------------

HilbertBasisMonoid::HilbertBasisMonoid(const Matrix<long long>& Gens, const Matrix<long long>& Supps){

    vector< pair< vector <long long>, vector<long long> > > GensWithValues;

    dim = Gens.nr_of_columns();
    nr_supps = Supps.nr_of_rows();
    nr_gens = Gens.nr_of_rows();
    GensWithValues.resize(nr_gens);


    // We order the generators by the lexicographic order of their values under supps
    for(size_t i = 0; i< nr_gens; ++i){
        GensWithValues[i].second = Gens[i];
        GensWithValues[i].first.resize(nr_supps+1);
        for(size_t j = 0; j < nr_supps; ++j){
            GensWithValues[i].first[j] = v_scalar_product(Supps[j], Gens[i]);
        }
        GensWithValues[i].first[nr_supps] = i; // we register the index w.r.t. Gens
    }

     /* for(size_t i = 0; i< nr_gens; ++i){
        cout << GensWithValues[i].first;
        cout << "               " << GensWithValues[i].second;
    }
    cout << "---------------" << endl; */

    sort(GensWithValues.begin(), GensWithValues.end());
    // we register the permutation
    for(size_t u = 0; u < nr_gens; ++u){
        ExternalKey.push_back(GensWithValues[u].first.back());

    }

    Gens_ordered.resize(0,dim);
    GensVal_ordered.resize(0, nr_supps);
    vector<long long> transfer;
    for(size_t u = 0; u < GensWithValues.size(); ++u){
        Gens_ordered.append(GensWithValues[u].second);
        transfer = GensWithValues[u].first;
        // remove the last component used for registeing the order
        transfer.resize(nr_supps);
        GensVal_ordered.append(transfer);
    }

    HilbertBasis.resize(0, dim);
    Representations.resize(0,nr_gens);
    internal_max_deg_ind.resize(nr_gens);
}

// compute Hilbert bbasis and representations via equation method
// not used at present
void HilbertBasisMonoid::computeHB_Equ(){

    // we skip zero vectors and put the first nonzero into Hilbert basis
    size_t u = 0;
    for(size_t i = 0; i < nr_gens; ++i){
        if(Gens_ordered[i] != vector<long long>(dim,0)){
            HilbertBasis.append(Gens_ordered[i]);
            InternalHilbBasKey.push_back(i);
            HilbertBasisKey.push_back(ExternalKey[u]);
            u++;
            break;
        }
        u++;
    }

    for(; u < nr_gens; ++u){

        // we assemble the inequalities for project_and_lift
        Matrix<long long> Help = Gens_ordered.submatrix(InternalHilbBasKey);
        Help = Help.transpose();
        Help.insert_column(0,0);
        Matrix<long long> Inequs = Help;
        Help.scalar_multiplication(-1); // equations split into uwo inequalities
        Inequs.append(Help);
        Inequs.append(Matrix<long long>(Inequs.nr_of_columns()) ); // nonnegativity

        // Now we compute the representations
        // we must insert element and -element into column 0
        // since we have split the nequations into ineqiualities

        for(size_t j = 0; j< dim; ++j){
            Inequs[j][0] = -Gens_ordered[u][j];
            Inequs[j+ dim][0] = Gens_ordered[u][j];
        }

        vector<dynamic_bitset> dummy_Ind;
        size_t dummy_rank = 0;

        ProjectAndLift<long long, long long> PL(Inequs, dummy_Ind, dummy_rank);
        PL.set_primitive();
        PL.set_LLL(false);
        PL.set_verbose(false);
        PL.compute(false,false,false); // single point, no float, not only counting
        vector<long long> sol;
        PL.put_single_point_into(sol);
        if(sol.size() == 0){
            HilbertBasis.append(Gens_ordered[u]);
            InternalHilbBasKey.push_back(u);
            HilbertBasisKey.push_back(ExternalKey[u]);
        }
        else{
            vector<long long> rel(nr_gens);
            for(size_t j=1; j < sol.size(); ++j){
                rel[HilbertBasisKey[j-1]] = - sol[j];
            }
            rel[ExternalKey[u]] = 1;
            Representations.append(rel);
        }
    }
}

// compute Hilbert bbasis and representations by the subtrction  method
// with backtracking

void HilbertBasisMonoid::computeHB_Sub(){

    pair<bool, vector<long long> > answer;
    vector<long long> rep(nr_gens);
    for(size_t u = 0; u < nr_gens; ++u){
        if(Gens_ordered[u] == vector<long long>(dim,0))
            continue;
        //if(!internal_max_deg_ind[u])
        answer = subtract_recursively(GensVal_ordered[u],0, rep,0);
        // if(!internal_max_deg_ind[u] && !answer.first){ // an element of the Hilbert basis
        if(!answer.first){
            InternalHilbBasKey.push_back(u);
            HilbertBasisKey.push_back(ExternalKey[u]);
            HilbertBasis.append(Gens_ordered[u]);
        }
        else{ // reducibloe
            vector<long long> rep_ext(nr_gens);
            for(size_t j = 0; j < nr_gens; ++j)
                rep_ext[ExternalKey[j]] = answer.second[j];
            rep_ext[ExternalKey[u]] = 1;
            Representations.append(rep_ext);
        }
    }
}

pair<bool, vector<long long> > HilbertBasisMonoid::subtract_recursively(vector<long long> val, size_t start, vector<long long> rep, int level){

    if(val == vector<long long>(nr_supps))
        return make_pair(true,rep);
    for(size_t uu = start; uu < InternalHilbBasKey.size(); ++uu){
        key_t i = InternalHilbBasKey[uu];
        bool subtractible = true;
        for(size_t j = 0; j < nr_supps; ++j){
            if(val[j] - GensVal_ordered[i][j] < 0){
                    subtractible = false;;
                    break;
            }
            if(!subtractible)
                continue;
        }
        if(subtractible){
            vector<long long> new_val = val;
            vector<long long> new_rep = rep;
            for(size_t j = 0; j < nr_supps; ++j){
                new_val[j] -= GensVal_ordered[i][j];
            }
            new_rep[i]--;
            pair<bool, vector<long long> > answer = subtract_recursively(new_val,uu,new_rep, level + 1);
                        if(answer.first){
                return answer;
            }
        }
    }
    return make_pair(false, rep);
}

void HilbertBasisMonoid::set_max_deg_ind(const dynamic_bitset& mdi){
        max_deg_ind = mdi;
}


void HilbertBasisMonoid::compute_HilbertBasis(){
    if(max_deg_ind.size() > 0){
        assert(max_deg_ind.size() == nr_gens);
        for(size_t i = 0; i < max_deg_ind.size(); ++i)
            internal_max_deg_ind[i] = max_deg_ind[ExternalKey[i]];
    }

    computeHB_Sub();
    // computeHB_Equ();
}

void HilbertBasisMonoid::put_HilbertBasis_into(Matrix<long long>& HB){
    swap(HB, HilbertBasis);
}

void HilbertBasisMonoid::put_Representations_into(Matrix<long long>& Rep){
    swap(Rep, Representations);
}

void HilbertBasisMonoid::put_HilbertBasisKey_into(vector<key_t>& Ind){
        sort(HilbertBasisKey.begin(), HilbertBasisKey.end());
        swap(Ind, HilbertBasisKey);
}

} // namespace
