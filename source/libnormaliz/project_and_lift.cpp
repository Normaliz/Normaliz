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

#include "libnormaliz/project_and_lift.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/input.h"

namespace libnormaliz{
using std::vector;
using std::string;
using std::list;
using std::pair;
//---------------------------------------------------------------------------

// classless functions

void write_control_file(const size_t split_level, const size_t nr_vectors){

    if(verbose)
        verboseOutput() << "split_level " << split_level << endl;

    SplitData def_split(global_project, split_level, nr_vectors);
    def_split.write_data();
}

template <typename Integer>
void write_local_solutions(const size_t level, const vector<vector<Integer> >& Sols){
    string file_name = global_project;
    file_name += "." + to_string(level) + ".sls";
    ofstream sols_out(file_name);
    sols_out << Sols.size() << endl;
    sols_out << Sols[0].size() << endl;
    sols_out << Sols;
    if(verbose)
        verboseOutput() << Sols.size()<< " local solutions stored on level " << level << endl;
}

template <typename Integer>
Matrix<Integer> reconstruct_equations(const Matrix<Integer>& Inequalities){

    Matrix<Integer> Equations(0,Inequalities.nr_of_columns());
    if(using_renf<Integer>())
        return Equations;
    if(Inequalities.nr_of_rows() == 0)
        return Equations;

    vector<Integer> test(Inequalities.nr_of_columns());
    set<vector<Integer> > Ineq;
    for(size_t i = 0; i < Inequalities.nr_of_rows(); ++i){
        Ineq.insert(Inequalities[i]);
    }

    Integer MinusOne = -1;
    for(size_t i = 0; i < Inequalities.nr_of_rows(); ++i){
        test = Inequalities[i];
        v_scalar_multiplication(test, MinusOne);
        if(Ineq.find(test) != Ineq.end()){
            Equations.append(Inequalities[i]);
            Ineq.erase(test); // we don't want it twice
            Ineq.erase(Inequalities[i]);
        }
    }

    Equations.remove_duplicate_and_zero_rows();
    return Equations;
}

template<typename Integer>
void coarsen_this_cong(const vector<Integer>& cong, const size_t k, set<vector<Integer> >& CongSet){

    for(size_t i = k; i < cong.size() - 1; ++i){
        if(cong[i] == 0)
            continue;
        Integer new_mod = libnormaliz::gcd(cong[i], cong.back());
        if(new_mod == 1)
            break;
        vector<Integer> coarser_cong(cong.size());
        for(size_t j = 0; j < cong.size() - 1; ++j)
            coarser_cong[j] = cong[j] % new_mod;
        coarser_cong.back() = new_mod;
        CongSet.insert(coarser_cong);
        coarsen_this_cong(coarser_cong, i + 1, CongSet);
    }
}

template<>
void coarsen_this_cong(const vector<renf_elem_class>& cong, const size_t k, set<vector<renf_elem_class> >& CongSet){
    assert(false);
}

template<>
void coarsen_this_cong(const vector<nmz_float>& cong, const size_t k, set<vector<nmz_float> >& CongSet){
    assert(false);
}

//--------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::add_congruences_from_equations() {

    if(using_renf<IntegerPL>() || using_float<IntegerPL>()){
            return;
    }

    set<vector<IntegerRet> > CongSet;
    for(size_t i = 0; i < Congs.nr_of_rows(); ++i)
        CongSet.insert(Congs[i]);
    for(size_t i = 0; i < Congs.nr_of_rows(); ++i){
        coarsen_this_cong(Congs[i], 0, CongSet);
    }
    Matrix<IntegerPL> RecosntructedEquations = reconstruct_equations(AllSupps[EmbDim]);
    for(size_t i = 0; i< RecosntructedEquations.nr_of_rows(); ++i){
        vector<IntegerRet> cong_candidate;
        convert(cong_candidate, RecosntructedEquations[i]);
        cong_candidate.resize(RecosntructedEquations.nr_of_columns() + 1);
        coarsen_this_cong(cong_candidate, 0, CongSet);
    }

    Congs.resize(0);
    for(auto& c: CongSet)
        Congs.append(c);

    // Congs.debug_print('$');
}

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::restrict_congruences() {

    if(using_renf<IntegerPL>() && using_float<IntegerPL>()){
            return;
    }

    // Congs.debug_print();
    for(size_t i = 1; i < AllCongs.size(); ++i){
        AllCongs[i] = Matrix<IntegerRet>(0, i+1);
        for(size_t j = 0; j <Congs.nr_of_rows(); ++j){
            if(Congs[j][i-1] == 0)
                continue;
            bool restrictable = true;
            for(size_t k = i; k < EmbDim; ++k){
                if(Congs[j][k] != 0){
                    restrictable = false;
                    break;
                }
            }
            if(!restrictable)
                continue;
            vector<IntegerRet> new_cong = Congs[j];
            new_cong.resize(i+1);
            new_cong.back() = Congs[j].back();
            AllCongs[i].append(new_cong);
        }
    }
}

// computes c1*v1-c2*v2
template <typename Integer>
vector<Integer> FM_comb(Integer c1, const vector<Integer>& v1, Integer c2, const vector<Integer>& v2, bool& is_zero) {
    size_t dim = v1.size();
    vector<Integer> new_supp(dim);
    is_zero = false;
    size_t k = 0;
    for (; k < dim; ++k) {
        new_supp[k] = c1 * v1[k] - c2 * v2[k];
        if (!check_range(new_supp[k]))
            break;
    }
    Integer g = 0;
    if (k == dim)
        g = v_make_prime(new_supp);
    else {  // redo in GMP if necessary
#pragma omp atomic
        GMP_hyp++;
        vector<mpz_class> mpz_neg(dim), mpz_pos(dim), mpz_sum(dim);
        convert(mpz_neg, v1);
        convert(mpz_pos, v2);
        for (k = 0; k < dim; k++)
            mpz_sum[k] = convertTo<mpz_class>(c1) * mpz_neg[k] - convertTo<mpz_class>(c2) * mpz_pos[k];
        mpz_class GG = v_make_prime(mpz_sum);
        convert(new_supp, mpz_sum);
        convert(g, GG);
    }
    if (g == 0)
        is_zero = true;

    return new_supp;
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::check_and_prepare_sparse() {


    size_t nr_all_supps = AllSupps[EmbDim].nr_of_rows();

    Indicator.resize(nr_all_supps); // indicaor of nonzero coordinates in inequality
    upper_bounds.resize(nr_all_supps); // indicator of inequalities giving upper boounds
    max_sparse.resize(nr_all_supps); // indicator of inequalities used in covering by "sparse" inequalities

    // Congs.debug_print();

    for(size_t i = 0; i< nr_all_supps; ++i){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        bool is_upper_bound = true;
        Indicator[i].resize(AllSupps[EmbDim][i].size());
        Indicator[i][0] = 1; // zeroeth coordinate always in the support (even if zero)
        for(size_t j = 1; j < AllSupps[EmbDim][i].size(); ++j){
            if(AllSupps[EmbDim][i][j] != 0){
                Indicator[i][j] = 1;
                if(AllSupps[EmbDim][i][j] > 0)
                    is_upper_bound = false;
            }
        }
        if(is_upper_bound)
            upper_bounds[i] = 1;
    }
    dynamic_bitset union_upper_bounds(EmbDim); // to check whether sparse upper bound inequ cover
    dynamic_bitset sparse_bounds(nr_all_supps);

    for(size_t i = 0; i< nr_all_supps; ++i){
        if(upper_bounds[i] && Indicator[i].count() < EmbDim){ // our criterion for sparseness
            union_upper_bounds |= Indicator[i];
            sparse_bounds[i] = 1;
        }
    }

    sparse = (union_upper_bounds.count() == EmbDim);

    if(!sparse){
        if(verbose)
            verboseOutput() << "System not sparse or only 1 patch" << endl;
        return;
    }

    if(verbose)
        verboseOutput() << "Preparing data for patching algorithm " << endl;

    // now we want to find the sparse upper bounds with maximal support
    vector<dynamic_bitset> help(nr_all_supps); // TODO find them, not only preparation
    for(size_t i = 0; i < nr_all_supps; ++i){
        help[i].resize(EmbDim);
        if(sparse_bounds[i])
            help[i] = Indicator[i];
    }
    dynamic_bitset max_sparse = sparse_bounds; // TODO here they must be found. So far: all

    dynamic_bitset covered(EmbDim);  // registers covered coordinates
    covered[0] = 1; // the 0-th coordinate is covered by all local PL
    AllLocalPL.resize(EmbDim);
    AllIntersections_key.resize(EmbDim);
    AllNew_coords_key.resize(EmbDim);
    AllCovered.resize(EmbDim);
    AllCoveredKey.resize(EmbDim);
    AllCoveredKeyInverse.resize(EmbDim);
    AllPatches.resize(EmbDim);
    active_coords.resize(EmbDim);

    AllCongsRestricted.resize(EmbDim);
    AllPolyEqus.resize(EmbDim);
    AllPolyEqusThread.resize(EmbDim);
    for(auto& T: AllPolyEqusThread)
        T.resize(omp_get_max_threads());

    AllPolyInequs.resize(EmbDim);
    AllPolyInequsThread.resize(EmbDim);
    for(auto& T: AllPolyInequsThread)
        T.resize(omp_get_max_threads());
    AllAutoms.resize(EmbDim);

    NrRemainingLP.resize(EmbDim,0);
    NrDoneLP.resize(EmbDim,0);
    AllLocalSolutions_by_intersection_and_cong.resize(EmbDim);
    AllLocalSolutions.resize(EmbDim);
    AllShortLocalSolutions.resize(EmbDim);

    DefiningSupps.resize(EmbDim);

    if(fusion.check_simplicity){
        if(fusion.candidate_given)
            check_simplicity_cand = true;
        else
            check_simplicity_all = true;
    }

    if(check_simplicity_all)
        fusion.all_critical_coords_keys.resize(EmbDim);

    poly_equs_minimized.resize(EmbDim);
    poly_inequs_minimized.resize(EmbDim);
    // poly_congs_minimized.resize(EmbDim);
    automs_minimized.resize(EmbDim);

    // First we compute the patches
    for(size_t coord = 1; coord < EmbDim; coord++){
        if(covered[coord] == 1)
            continue;

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        active_coords[coord] = 1;

        key_t next_supp = 0; //  the next inequality we will use // = 0 to makemgcc happy

        // try to find an inequality with coord in its support and minimal extension
        size_t min_card_extension = 0;
        bool found_extension = false;
        for(size_t i = 0; i < nr_all_supps; i++){
            if(used_supps[i] || !max_sparse[i])
                continue;
            if(Indicator[i][coord]){
                size_t card_extension = 0;
                for(size_t k = 0; k < Indicator[i].size(); ++k){
                    if(Indicator[i][k] && !covered[k])
                        card_extension++;
                }
                if(!found_extension || card_extension < min_card_extension){
                    next_supp = i;
                    min_card_extension = card_extension;
                }
                found_extension = true;
            }
        }
        assert(found_extension);

        // The coordinates of the patch found
        AllPatches[coord] = Indicator[next_supp];
        DefiningSupps[coord] = AllSupps[EmbDim][next_supp];
        covered |= AllPatches[coord];

    }

    // Wec want to restrict global congruences
    for(size_t k = 0; k < Congs.nr_of_rows(); ++k){
        dynamic_bitset cong_support(EmbDim);  // modulus not indicated
        for(size_t i = 0; i < EmbDim; ++i){
            if(Congs[k][i] != 0)
                cong_support[i] = true;
        }
        CongIndicator.push_back(cong_support);
    }

    // we compute basic data for simplicity test
    if(check_simplicity_cand){
        critical_coord_simplicity = -1;
    }

    compute_covers();
    // exit(0);

    covered = dynamic_bitset(covered.size()); // reset covered to empty set
    covered[0] = true;

    // now the main work: find "local" solutions and patach them
    // we extend the set of covered coordinates by at least
    // the first noncovered

    // Now wer go over the insertion order

    //  Indicator[next_supp] ersetzen durch AllPatches

    vector<key_t> old_coord_key;

    for(auto& coord: InsertionOrderPatches){

        // now we want to build "local" systems
        // First we add was is is contained in AllPatches[coord]
        vector<key_t> relevant_supps_now;
        for(size_t i = 0; i < nr_all_supps; ++i){
            if(Indicator[i].is_subset_of(AllPatches[coord])){
                relevant_supps_now.push_back(i);
                used_supps[i] = 1;
            }
        }

        // next we add what is implied by other upper bounds
        for(size_t i = 0; i < nr_all_supps; ++i){
            if(!used_supps[i] && upper_bounds[i]){
                relevant_supps_now.push_back(i);
            }
        }

        // now the intersections and new_coords
        dynamic_bitset intersection_coods = covered & AllPatches[coord];
        dynamic_bitset new_coords(AllPatches[coord]);
        for(size_t i = 0; i< EmbDim; ++i)
            if(new_coords[i] && covered[i])
                new_coords[i] = 0;

        vector<key_t> intersection_key = bitset_to_key(intersection_coods); // w.r.t. full coordinates
        vector<key_t> new_coords_key = bitset_to_key(new_coords); // w.r.t. to full coordinates
        AllIntersections_key[coord] = intersection_key;
        AllNew_coords_key[coord] = new_coords_key;

        if(talkative){
                verboseOutput() << "level " << LevelPatches[coord] << endl;
                verboseOutput() << "new coords " << new_coords_key;
        }
        // for the "local" project-and-lift we need their suport hyperplanes
        vector<key_t> LocalKey = bitset_to_key(AllPatches[coord]);
        Matrix<IntegerPL> LocalSuppsRaw; // could be avoided, only for convenience
        LocalSuppsRaw = AllSupps[EmbDim].submatrix(relevant_supps_now);
        /* vector<IntegerPL> test_v(EmbDim);
        test_v[0] = -1;
        for(size_t kkn = 0; kkn < LocalSuppsRaw.nr_of_rows(); ++kkn)
            if(test_v == LocalSuppsRaw[kkn])
                assert(false);*/
        // LocalSuppsRaw.debug_print('+');
        // convert(LocalSuppsRaw, AllSupps[EmbDim].submatrix(relevant_supps_now));
        // Matrix<IntegerPL> Localsupps = LocalSuppsRaw.transpose().submatrix(LocalKey).transpose(); // select columns

        // in For the "local" project-and-lift we must put the intersection coordinates first
        // then the new coordinates
        size_t nr_coordinates = intersection_key.size() + new_coords_key.size();
        vector<key_t> OrderedCoordinates(nr_coordinates);
        for(size_t i = 0; i< nr_coordinates; ++i){
            if(i < intersection_key.size())
                OrderedCoordinates[i] = intersection_key[i];
            else
                OrderedCoordinates[i] = new_coords_key[i- intersection_key.size()];
        }

        // for the "local" project-and-lift we must correspondingly reorder the support hyperplanes.
        Matrix<IntegerPL> LocalSuppsReordered(LocalSuppsRaw.nr_of_rows(), nr_coordinates);
        for(size_t i = 0; i < LocalSuppsRaw.nr_of_rows(); ++i){
            for(size_t j = 0; j < nr_coordinates; ++j)
                LocalSuppsReordered[i][j]= LocalSuppsRaw[i][OrderedCoordinates[j]];
        }

        // Now we selct the congruences that apply locally

        dynamic_bitset relevant_congs_now(Congs.nr_of_rows());
        for(size_t k = 0; k < relevant_congs_now.size(); ++k){
            if(CongIndicator[k].is_subset_of((AllPatches[coord])))
                relevant_congs_now[k] = true;
        }

        Matrix<IntegerRet> LocalCongsRaw;
        LocalCongsRaw = Congs.submatrix(bitset_to_key(relevant_congs_now));

        Matrix<IntegerRet> LocalCongsReordered(LocalCongsRaw.nr_of_rows(), nr_coordinates + 1); // +1 for modulus
        for(size_t i = 0; i < LocalCongsRaw.nr_of_rows(); ++i){
            for(size_t j = 0; j < nr_coordinates; ++j)
                LocalCongsReordered[i][j]= LocalCongsRaw[i][OrderedCoordinates[j]];
            LocalCongsReordered[i][nr_coordinates] = LocalCongsRaw[i].back(); // modulus transferred
        }
        // Congs done


        // Now we can set up the local project-and-lift
        vector<dynamic_bitset> DummyInd;
        LocalSuppsReordered.remove_duplicate_and_zero_rows();
        ProjectAndLift<IntegerPL, IntegerRet> PL(LocalSuppsReordered, DummyInd, 0); // 0 is dummy
        PL.set_LLL(false);
        PL.set_primitive();
        PL.set_verbose(false);
        PL.set_short_int(use_short_int);
        PL.set_congruences(LocalCongsReordered);
        PL.compute_projections_primitive(LocalKey.size());
        PL.add_congruences_from_equations();
        PL.restrict_congruences();
        if(false){ // verbose){
            if(PL.Congs.nr_of_rows() >0 ){
                verboseOutput() << "coord " << coord << endl;
                PL.Congs.debug_print();
            }
        }
        AllLocalPL[coord] = PL;

        dynamic_bitset new_covered = covered | AllPatches[coord];

        // first we check whether the simplicity test can be done at this point
        if(check_simplicity_cand){
            if(critical_coord_simplicity == -1){
                if(fusion.coords_to_check_ind[0].is_subset_of(new_covered) ){
                    critical_coord_simplicity = coord;
                }
            }
        }
        if(check_simplicity_all){
            for(size_t i = 0; i < fusion.coords_to_check_ind.size(); ++i){
                if(fusion.coords_to_check_ind[i].is_subset_of(new_covered) && !fusion.coords_to_check_ind[i].is_subset_of(covered))
                    fusion.all_critical_coords_keys[coord].push_back(fusion.coords_to_check_key[i]);
            }
        }

        // now we find the congruences that can be restricted to new_covered, but not to
        // covered or AllPatches[coord]

        for(size_t i = 0; i < Congs.nr_of_rows(); ++i){
            if(!CongIndicator[i].is_subset_of(new_covered))
                continue;
            if(CongIndicator[i].is_subset_of(covered))
                continue;
            if(CongIndicator[i].is_subset_of(AllPatches[coord]))
                continue;
            AllCongsRestricted[coord].push_back(OurPolynomialCong<IntegerRet>(Congs[i]));
        }

         // now the extra constraints, in perticular the polynomial ones
         // First the linear inequalities
         // will be added to thepolynomial ones TODO simplify
         Matrix<IntegerRet> ExtraInequalities(0, EmbDim);
         for(size_t i = 0; i <nr_all_supps ; ++i){
            if(used_supps[i])
                continue;
            if(Indicator[i].is_subset_of(covered))
                continue;
            if(Indicator[i].is_subset_of(AllPatches[coord]))
                continue;
            if((Indicator[i] & new_covered).count() == 0)
                continue;
             if(!using_renf<IntegerPL>()){
                bool can_be_restricted = true;
                for(size_t j = 0; j< EmbDim; ++j){
                    if(new_covered[j])
                        continue;
                    if(AllSupps[EmbDim][i][j] > 0){
                        can_be_restricted = false;
                        break;
                    }
                }
                if(can_be_restricted){
                    vector<IntegerRet> inequ;
                    convert(inequ, AllSupps[EmbDim][i]);
                    ExtraInequalities.append(inequ);
                }
             }
        }

        // Collect relevant polynomial constraints
        // The keys are used for terminal output
        vector<key_t> PolyEqusKey, PolyInequsKey, RestrictablePolyInequsKey;

        /*
         * The following is the first step towards exploiting polynomial equations as
         * congruiences.
         * Not yet implemented further since it is unclear whether it makes sense.
         *  TODO don't delete

        size_t linear_in_new_covered = 0;
        for(size_t i = 0; i < PolyEquations.size(); ++i){
            if(PolyEquations[i].support.is_subset_of(new_covered)) // not yet usable
                continue;
            OurPolynomial<IntegerRet> new_part = PolyEquations[i].split(new_covered).second;
            dynamic_bitset support_linear(EmbDim);
            bool is_linear_in_new_part = new_part.check_linearity(new_covered, support_linear);
            if(is_linear_in_new_part && support_linear.count() <= 3){
                linear_in_new_covered++;
                cout << support_linear.count() << " ";
            }
        }
        cout << " Linears " << linear_in_new_covered << endl;
        */

        // first the equations
        // bool first_rest = true;
        for(size_t i = 0; i < PolyEquations.size(); ++i){
            if(!PolyEquations[i].support.is_subset_of(new_covered)) // not yet usable
                continue;
            if(PolyEquations[i].support.is_subset_of(covered)) // already used
                continue;
            OurPolynomial<IntegerRet> Restrict = PolyEquations[i].restrict_to(covered);
            /* if(first_rest){
                cout << "Rest " << i << " --- " <<  Restrict.size() << " of " << PolyEquations[i].size() << endl;
                first_rest = false;
            }*/
            AllPolyEqus[coord].push_back(PolyEquations[i].split(covered));
            AllPolyEqus[coord].back().first.vectorize_deg_2();
            AllPolyEqus[coord].back().second.vectorize_deg_2();
            PolyEqusKey.push_back(i);
        }
        for(auto& T: AllPolyEqusThread[coord]){ // vcopy for each thread
            T = AllPolyEqus[coord];
        }

        // now the full inequalities
        for(size_t i = 0; i < PolyInequalities.size(); ++i){

            if(PolyInequalities[i].support.is_subset_of(covered)) // already used
                continue;
            if(!(PolyInequalities[i]).support.is_subset_of(new_covered))
                continue;
            AllPolyInequs[coord].push_back(PolyInequalities[i]);
            PolyInequsKey.push_back(i);
        }

        // now the inequalities which can be restricted
        // treated separately to avoid double evaluations of equations
        for(size_t i = 0; i < RestrictablePolyInequs.size(); ++i){
            if(RestrictablePolyInequs[i].support.is_subset_of(new_covered))
                continue;
            if(!(RestrictablePolyInequs[i]).is_restrictable_inequ(new_covered))
                continue;
            AllPolyInequs[coord].push_back(RestrictablePolyInequs[i]);
            RestrictablePolyInequsKey.push_back(i);
        }

        /* Matrix<IntegerRet> ExtraEquations = reconstruct_equations(ExtraInequalities);

        ExtraEquations.debug_print('E');*/

        for(size_t i = 0; i < ExtraInequalities.nr_of_rows(); ++i){
            AllPolyInequs[coord].push_back(OurPolynomial<IntegerRet>(ExtraInequalities[i]));
        }

        for(auto& T: AllPolyInequsThread[coord]){ // vcopy for each thread
            T = AllPolyInequs[coord];
        }

        AllAutoms[coord] =identity_key(fusion.Automorphisms.size()); // initialized with all automorphisms in the given order

        if(talkative){
                verboseOutput() << "index coord " << coord << " nr covered coordinates " << new_covered.count() << endl;
            if(coord == critical_coord_simplicity)
                verboseOutput() << "simplicity check at this coordinate" << endl;
            if(verbose && PolyEqusKey.size() > 0)
                verboseOutput() << "nr poly equations " << PolyEqusKey.size() << endl;
            if(verbose && PolyInequsKey.size() > 0)
                verboseOutput() << " nr poly inequalities " << PolyInequsKey.size() << endl;
            if(verbose && RestrictablePolyInequsKey.size() > 0)
                verboseOutput() <<  "nr restrictable poly inequalities " << RestrictablePolyInequsKey.size() << endl;
            if(verbose && AllCongsRestricted[coord].size() > 0){
                verboseOutput() <<  "nr congruences " << AllCongsRestricted[coord].size() << endl;
            }
            if(verbose)
                verboseOutput() << "---------------------------------------------------------------" << endl;
        }

        if(old_coord_key.size() == 0){
            vector<key_t> addition = bitset_to_key(new_covered);
            AllCoveredKey[coord].insert(AllCoveredKey[coord].end(), addition.begin(),addition.end());
        }
        else{
            AllCoveredKey[coord] = old_coord_key;
            for(key_t i = 0; i< new_covered.size(); ++i){
                if(!new_covered[i] || covered[i])
                    continue;
                AllCoveredKey[coord].push_back(i);
            }
        }
        old_coord_key = AllCoveredKey[coord];

        AllCoveredKeyInverse[coord].resize(new_covered.size());
        for(size_t i = 0; i < AllCoveredKey[coord].size(); ++i)
            AllCoveredKeyInverse[coord][AllCoveredKey[coord][i]] = i;

        AllCovered[coord] = new_covered;
        covered = new_covered;

        if(fusion.total_FPdim == 0 || no_heuristic_minimization) // no discarding of equations if not fusion
            poly_equs_minimized[coord] = true;

    } // coord

    max_nr_new_latt_points_total = libnormaliz::max_nr_new_latt_points_total;
    max_nr_new_latt_points_patching = libnormaliz::max_nr_new_latt_points_patching;
    nr_extensions_for_elimination_equs = libnormaliz::nr_extensions_for_elimination_equs;
    nr_extensions_for_elimination_inequs = libnormaliz:: nr_extensions_for_elimination_inequs;
    nr_extensions_for_elimination_automs = libnormaliz::nr_extensions_for_elimination_automs;

    // to reduce the danger of disaster by premature discarding of poly equations
    if(fusion.total_FPdim >= 1000)
        nr_extensions_for_elimination_equs *= 10;
    if(fusion.total_FPdim >= 10000)
        nr_extensions_for_elimination_equs *= 10;

    FreeVectThread.resize(omp_get_max_threads());
}

//---------------------------------------------------------------------------

// First we cover the supports of polynomial equations by patches
// using as few patches as possioble
// The results are stored in covering_patches
template <typename IntegerPL, typename IntegerRet>
vector<pair<size_t, vector<key_t> > > ProjectAndLift<IntegerPL,IntegerRet>::
              cover_supports(const vector<dynamic_bitset>& supports) {

    vector<pair<size_t, vector<key_t> > > covering_equations;
    size_t dim = EmbDim;

    for(size_t i = 0; i< supports.size(); ++i){
        dynamic_bitset poly_supp = supports[i];
        dynamic_bitset already_covered(dim);
        size_t nr_patches_needed = 0;
        vector<key_t> patches_used;
        while(true){
            bool first = true;
            size_t max_nr_new_covered = 0; // = 0 to make gcc happy
            size_t max_covering = 0; // to make gcc happy
            for(size_t k = 0; k < EmbDim; ++k){
                if(AllPatches[k].size() == 0)
                    continue;
                size_t nr_new_covered = ((already_covered | AllPatches[k]) & poly_supp).count();
                if(first || nr_new_covered > max_nr_new_covered){
                    first= false;
                    max_nr_new_covered = nr_new_covered;
                    max_covering = k;
                }
            }
            already_covered = (already_covered |AllPatches[max_covering]) & poly_supp;
            patches_used.push_back(max_covering);
            nr_patches_needed++;
            if(already_covered == poly_supp){
                sort(patches_used.begin(), patches_used.end());
                covering_equations.push_back(make_pair(patches_used.size(),
                                                       patches_used));
                break;
            }
        }
    }

    return covering_equations;
}

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::find_order_congruences() {

    dynamic_bitset used_covering_congruences(Congs.nr_of_rows());
    dynamic_bitset covered_coords(EmbDim); // = AllPatches[InsertionOrderPatches[0]];
    dynamic_bitset used_patches(EmbDim);
    // used_patches[InsertionOrderPatches[0]] = true;
    while(covered_coords.count() < EmbDim && used_covering_congruences.count() < Congs.nr_of_rows()){
        bool first = true;
        bool first_cong = true;
        size_t max_at = 0;
        size_t max_nr_congs = 0;
        bool have_congs = false;
        double min_weight = 0;
        double cong_min_weight = 0;
        dynamic_bitset max_covered(EmbDim);
        for(size_t i = 0; i < AllPatches.size(); ++i){  // i is the index of the patch tested for optimality
            if(AllPatches[i].size() == 0 || used_patches[i]) // used or not existent
                continue;

            dynamic_bitset test_covered = covered_coords | AllPatches[i];

            double test_weight = 0;
            for(size_t w = 0; w < EmbDim; ++w){
                if(covered_coords[w] || !test_covered[w])
                    continue;
                test_weight += WeightOfCoord[i][w];
            }
            size_t nr_congs = 0;
            for(size_t j = 0; j < Congs.nr_of_rows(); ++j){
                if(CongIndicator[j].is_subset_of(covered_coords))
                    continue;
                if(CongIndicator[j].is_subset_of(AllPatches[i]))
                    continue;
                if(!CongIndicator[j].is_subset_of(test_covered))
                    continue;
                if(!use_coord_weights){ // in this case we optimize the number of "global" congruences
                    nr_congs++;
                    continue;
                }
                // with use_coord_weights we only count congruences that are not contained
                // in lower weght patches
                bool optimal = true; // for the tested congruence
                for(size_t k = 0; k < EmbDim; ++k){  // k is index of alternative patch
                    if(AllPatches[k].size() == 0 || used_patches[k])
                        continue;
                    if(!CongIndicator[j].is_subset_of(AllPatches[k])) // we count only congruences
                        continue;                                    // that are not gotten with smaller weight
                    double alternative_weight = 0; // weight for alternative patch
                    for(size_t w = 0; w < EmbDim; ++w){
                        if(covered_coords[w] || !AllPatches[k][w])
                            continue;
                        alternative_weight += WeightOfCoord[w][k];
                    }
                    if(alternative_weight < test_weight){
                        optimal = false;
                        break;
                    }
                }
                if(optimal)
                    nr_congs++;
            }
            if(!use_coord_weights &&  (first || (nr_congs > max_nr_congs))){
                first = false;
                max_at = i;
                max_nr_congs = nr_congs;
                max_covered = test_covered;
            }
            if(use_coord_weights){
                if(nr_congs >0)
                    have_congs = true;
                if(first){
                    first = false;
                    min_weight = test_weight;
                    max_at = i;
                }
                if(have_congs && first_cong){
                    first_cong = false;
                    cong_min_weight = test_weight;
                    max_at = i;
                }
                if(have_congs){
                    if(test_weight < cong_min_weight){
                        max_at = i;
                        cong_min_weight = test_weight;
                    }
                }
                if(!have_congs && first_cong){ // still waiting for a congruence
                    if(test_weight < min_weight){
                        min_weight = test_weight;
                        max_at = i;
                    }
                }
            }
        }
        InsertionOrderPatches.push_back(max_at);
        used_patches[max_at] = true;
        covered_coords |= AllPatches[max_at];
        for(size_t j = 0; j < Congs.nr_of_rows(); ++j){  // register used congs
            if(CongIndicator[j].is_subset_of(max_covered))
                used_covering_congruences[j] = true;
        }
    }

    finalize_order(used_patches);
}

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::find_order_linear() {


    dynamic_bitset covered_coords(EmbDim);
    dynamic_bitset used_patches(EmbDim);
    while(covered_coords.count() < EmbDim){
        bool first = true;
        size_t min_at = 0;
        double min_weight = 0;
        dynamic_bitset min_covered(EmbDim);
        for(size_t i = 0; i < AllPatches.size(); ++i){
            if(AllPatches[i].size() == 0 || used_patches[i])
                continue;

            dynamic_bitset test_covered = covered_coords | AllPatches[i];
            double test_weight = 0;
            for(size_t j = 0;j < test_covered.size(); ++j){
                if(!covered_coords[j] && test_covered[j])
                    test_weight += WeightOfCoord[i][j];
            }
            if(first || test_weight < min_weight){
                first = false;
                min_at = i;
                min_covered = test_covered;
                min_weight = test_weight;
            }
        }
        InsertionOrderPatches.push_back(min_at);
        used_patches[min_at] = true;
        covered_coords |= AllPatches[min_at];
    }

    finalize_order(used_patches);
}

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::finalize_order(const dynamic_bitset& used_patches) {

    // There could be patches not used.
    for(size_t j = 0; j < EmbDim; ++j){
        if(!used_patches[j] && AllPatches[j].size() > 0)
            InsertionOrderPatches.push_back(j);
    }

    if(verbose){
        verboseOutput() << "Insertion order linear patches " << endl;
        verboseOutput() << InsertionOrderPatches << endl;
    }

    for(size_t k = 0; k < InsertionOrderPatches.size(); ++k)
        LevelPatches[InsertionOrderPatches[k]] = k;

    ExpectedNrRounds.resize(InsertionOrderPatches.size());
    TimeToLevel.resize(InsertionOrderPatches.size() +1);
    NrNodes.resize(InsertionOrderPatches.size() +1, 1);

}

template <typename IntegerPL, typename IntegerRet>
bool ProjectAndLift<IntegerPL,IntegerRet>::order_patches_user_defined() {

    string name = global_project + ".order.patches";
    const char* file_name = name.c_str();
    ifstream in_order;
    in_order.open(file_name, ifstream::in);
    if(in_order.is_open()){
        string test;
        in_order >> test;
        if(test != "nr_patches")
            throw BadInputException("<project>.order.patches does not start with nr_patches");
        long nr_patch;
        in_order >> nr_patch;
        dynamic_bitset used_patches(EmbDim);
        for(size_t i = 0; i < nr_patch; ++i){
            size_t j;
            in_order >> j;
            if(j >= EmbDim || AllPatches[j].empty() )
                throw BadInputException("File defining insertion order corrupt");
            if(used_patches[j])
                throw BadInputException("<project>.order.patches contains " + to_string(j) + " more than once");
            used_patches[j] = true;
            InsertionOrderPatches.push_back(j);
        }
        in_order.close();
        finalize_order(used_patches);
        return true;
    }

    return false;
}

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_covers() {

    // Note the indices of thze patches are the coordinates
    // with whom they are associated

    LevelPatches.resize(EmbDim);


    if(talkative){
        for(size_t i = 0; i < EmbDim; ++i)
            verboseOutput() << i << " ----- " << bitset_to_key(AllPatches[i]);
    }


    if(order_patches_user_defined()){
        if(verbose)
            verboseOutput() << "Insertion order user defined" << endl;
        return;
    }

    WeightOfCoord.resize(EmbDim,EmbDim);
    for(size_t i = 1; i< EmbDim; ++i){
        if(AllPatches[i].size() == 0)
            continue;
        WeightOfCoord[i][0] = 0;
        for(size_t j = 1; j < EmbDim; ++j){
            if(!AllPatches[i][j])
                continue;
            double num = convertTo<nmz_float>(DefiningSupps[i][0]);
            double den = convertTo<nmz_float>(-DefiningSupps[i][j]);
            WeightOfCoord[i][j] = num/den;
        }
    }

    vector<double> TotalWeights(EmbDim);
    double min_weight = 0, max_weight= 0;
    bool first = true;
    for(size_t i = 1; i< EmbDim; ++i){
        if(AllPatches[i].size() == 0)
            continue;
        for(size_t j = 1; j < EmbDim; ++j){
            if(!AllPatches[i][j])
                continue;
            TotalWeights[i] += WeightOfCoord[i][j];
        }
        if(first && TotalWeights[i] > 0){
            first = false;
            min_weight = TotalWeights[i];
            max_weight = TotalWeights[i];
        }
        if(TotalWeights[i] > 0){
            if(TotalWeights[i] < min_weight)
                min_weight = TotalWeights[i];
            if(TotalWeights[i] > max_weight)
                max_weight = TotalWeights[i];
        }
    }

    if(max_weight / min_weight > 4 && !no_weights && !use_coord_weights){
        use_coord_weights = true;
        if(verbose)
            verboseOutput() << "Weights activated" << endl;
    }

    if(!use_coord_weights){
        for(size_t i = 1; i< EmbDim; ++i){
            if(AllPatches[i].size() == 0)
                continue;
            for(size_t j = 1; j < EmbDim; ++j){
                if(!AllPatches[i][j])
                    continue;
                WeightOfCoord[i][j] = 1.0;
            }
        }
    }

    assert(!linear_order_patches || !cong_order_patches);

    if(!linear_order_patches && !cong_order_patches){
        if(PolyEquations.empty())
            linear_order_patches = true;
    }

    if(linear_order_patches){
        find_order_linear();
        return;
    }

    if(cong_order_patches){
        find_order_congruences();
        return;
    }

    // now we coompute the order base on polynomial equations

    vector<dynamic_bitset> our_supports;
    for(size_t i = 0; i < PolyEquations.size();++i)
            our_supports.push_back(PolyEquations[i].support);

    vector<pair<size_t, vector<key_t> > > covering_equations =
            cover_supports(our_supports);

    // Now we must find an order of the supports of the polynomial equations
    // in which we insert their patches.
    // The idea is to have as miuch overlap in the covered coordinates as
    // possible.

    sort(covering_equations.begin(), covering_equations.end());

    /* for(auto& c: covering_equations)
        cout << c.second; */

    size_t dim = EmbDim;

    vector<key_t> InsertionOrderEquations;
    dynamic_bitset used_covering_equations(covering_equations.size());
    dynamic_bitset covered_coords(dim);
    while(covered_coords.count() < dim && used_covering_equations.count() < covering_equations.size()){
        dynamic_bitset test_covered_coords(dim);
        dynamic_bitset min_covered_coords(dim);
        bool first = true;
        size_t min_at = 0; //value to make g++ happy
        double min_added_weight = 0; // =0 to make gcc happy

        for( size_t i = 0; i < covering_equations.size(); ++i){
            if(used_covering_equations[i])
                continue;
            test_covered_coords = covered_coords;
            for(auto& c: covering_equations[i].second) // does equation add anything ?
                test_covered_coords |= AllPatches[c];
            if(test_covered_coords == covered_coords){
                used_covering_equations[i] = true; // no new coordinates
                continue;
            }
            vector<double> coord_added_weight(EmbDim);
            dynamic_bitset weight_found(EmbDim);
            for(size_t j = 0; j < EmbDim; ++j){
                if(covered_coords[j] || !test_covered_coords[j])
                    continue;
                for(auto& c: covering_equations[i].second){
                    if(!AllPatches[c][j])
                        continue;
                    if(!weight_found[j])
                        coord_added_weight[j] = WeightOfCoord[c][j];
                    else{
                        coord_added_weight[j] = std::min(WeightOfCoord[c][j],coord_added_weight[j]);
                    }
                }
            }
            double new_added_weight = 0;
            for(auto& w: coord_added_weight)
                new_added_weight += w;

            if(first || new_added_weight < min_added_weight){
                first = false;
                min_added_weight = new_added_weight;
                min_covered_coords = test_covered_coords;
                min_at = i;
            }
        }
        used_covering_equations[min_at] = true;
        InsertionOrderEquations.push_back(min_at);
        for(auto& c:covering_equations[min_at].second){
            covered_coords |= AllPatches[c];
        }
    }

    dynamic_bitset inserted_patches(EmbDim);
    for(auto& i: InsertionOrderEquations){
        for(auto& c: covering_equations[i].second){
            if(!inserted_patches[c])
                InsertionOrderPatches.push_back(c);
            inserted_patches[c] = true;
        }
    }

    finalize_order(inserted_patches);
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_local_solutions_for_saving() {

    vector<IntegerRet> start(EmbDim);
    start[0] = GD;

    for(long this_patch = 0; this_patch <= level_local_solutions; ++this_patch){

        size_t coord = InsertionOrderPatches[this_patch];
        vector<IntegerRet> start(1, GD);
        start_list.push_back(start);
        AllLocalPL[coord].lift_points_to_this_dim(start_list);

        if(use_short_int){
            vector<vector<short> > ShortLocalSolutionsNow;
            AllLocalPL[coord].put_short_deg1Points_into(ShortLocalSolutionsNow);
            write_local_solutions(this_patch, ShortLocalSolutionsNow);
        }
        else{
            vector<vector<IntegerRet> > LocalSolutionsNow;
            AllLocalPL[coord].put_deg1Points_into(LocalSolutionsNow);
            write_local_solutions(this_patch, LocalSolutionsNow);
        }
    }
    return;
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_latt_points_by_patching() {

    if(save_local_solutions){
        compute_local_solutions_for_saving();
        throw NoComputationException("No output with SavedLocalSolutions");
    }

    vector<IntegerRet> start(EmbDim);
    start[0] = GD;
    list<vector<IntegerRet> > start_list;
    start_list.push_back(start);
    extend_points_to_next_coord(start_list, 0);
    NrLP[EmbDim] = TotalNrLP;
    if(verbose){
        verboseOutput() << endl << "=======================================" << endl;
        verboseOutput() << "Final number of lattice points "  << NrLP[EmbDim] << endl;
    }

    if(!only_single_point && !distributed_computation){
        for(auto& n: NrRemainingLP){
                assert(n == 0);
        }
    }

    // we stop other splits if this one was successful
    if(is_split_patching && only_single_point && NrLP[EmbDim] > 0){
        string name = global_project + ".spst";
        ofstream stop_file(name);
        stop_file << " ";
        stop_file.close();
    }
}

//---------------------------------------------------------------------------


// evaluate congruence partially and return negative of result
template <typename IntegerRet>
IntegerRet eval_cong_partially(const OurPolynomialCong<IntegerRet>& cong,
                                  const vector<IntegerRet>& local_solution_new,
                                  const dynamic_bitset& restriction, const bool take_neg){

    IntegerRet res = cong.poly.evaluate_restricted(local_solution_new, restriction);
    res %= cong.modulus;
    // norm res 0 <= res < modulus
    if(res < 0)
        res += cong.modulus;
    // res = -res mod modulus
    if(take_neg && res != 0)
        res = cong.modulus - res;
    return res;
}

//---------------------------------------------------------------------------

vector<key_t> global_intersection_key;

// lexicographic comparison with overlap first
template <typename Integer>
bool intersect_compare(const vector<Integer>& v, const vector<Integer>& w){

    if(v_select_coordinates(v, global_intersection_key) < v_select_coordinates(w, global_intersection_key))
        return true;
    if(v_select_coordinates(v, global_intersection_key) == v_select_coordinates(w, global_intersection_key))
        return (v < w);
    return false;
}


template <typename Integer>
void sort_lattice_points_by_overlap(list<vector<Integer> >& LatticePoints,const vector<key_t>&  intersection_key){
    global_intersection_key = intersection_key;
    LatticePoints.sort(intersect_compare<Integer>);
}
//---------------------------------------------------------------------------

template <typename Integer>
void select_and_split(list<vector<Integer> >& LatticePoints, const key_t& this_patch, const long& split_modulus, const long& split_residue, const size_t& done_indices, const vector<key_t>&  intersection_key){


    if(verbose){
        verboseOutput() << "==========================" << endl;
        verboseOutput() << LatticePoints.size() << " lattice points before splitting and selection" << endl;
        verboseOutput() << "Spilt level " << this_patch << " modulus " << split_modulus << " residue " << split_residue << endl;
    }
    sort_lattice_points_by_overlap(LatticePoints,intersection_key);
    list<vector<Integer> > Selection;


    if(done_indices > 0){
        list<vector<Integer> > PreSelection;
        size_t k = 0;
        for(auto& v: LatticePoints){
            if(k >= done_indices)
                PreSelection.push_back(v);
            k++;
        }
        size_t given_lattice_points = LatticePoints.size();
        swap(LatticePoints, PreSelection);
        if(verbose)
            verboseOutput() << done_indices << " already done lattice points of " << given_lattice_points << " discarded, "
            << LatticePoints.size() << " remaining" << endl;
        if(done_indices > given_lattice_points){
            verboseOutput() << "ALARM" << endl;
            assert(false);
        }
        PreSelection.clear();
    }

    size_t nr_left = LatticePoints.size();
    size_t nr_per_split = nr_left / split_modulus;
    size_t to_add_one = nr_left % split_modulus;

    size_t start, last;
    if(split_residue < to_add_one){
        start = (nr_per_split + 1) * split_residue;
        last = start + nr_per_split + 1;
    }
    else{
        size_t new_start = to_add_one * (nr_per_split +1);
        start = new_start + (split_residue - to_add_one) * nr_per_split;
        last = start + nr_per_split;
    }

    if(split_residue + 1 == split_modulus)
        assert(nr_left == last);

    size_t i = 0;
    for(auto& p: LatticePoints){
        if(i >= start && i < last)
            Selection.push_back(p);
       i++;
    }

    if(verbose)
        verboseOutput() << Selection.size() << " lattice points after splitting" << endl;

    swap(LatticePoints, Selection);
}

//--------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::prepare_split(list<vector<IntegerRet> >& LatticePoints, const key_t& this_patch){

    const size_t coord = InsertionOrderPatches[this_patch];
    auto& intersection_key = AllIntersections_key[coord];

    for(size_t i = 0; i < our_split.nr_split_levels; ++i){
        if(this_patch == our_split.this_split_levels[i]){
            long split_modulus = our_split.split_moduli[i];
            long split_residue = our_split.this_split_residues[i];
            size_t done_indices = 0;
            size_t total_indices = 0;
            if(i >= 1){
                done_indices = our_split.this_split_done_indices[i - 1];
                total_indices = our_split.this_split_total_indices[i - 1];
                assert(LatticePoints.size() == total_indices);
            }
            select_and_split(LatticePoints, this_patch, split_modulus,split_residue, done_indices, intersection_key);
        }
    }
}

//--------------------------------------------------------------------------

template <typename Integer>
bool import_local_solutions(vector<vector<Integer> >& SavedLocalSolutions, const key_t& this_patch){

    bool read_local_solutions = false;

    string file_name = global_project + "." + to_string(this_patch) + ".sls";
    ifstream sls(file_name);
    if(sls.is_open()){
        read_local_solutions = true;
        size_t nr_rows, nr_cols;
        sls >> nr_rows;
        sls >> nr_cols;
        SavedLocalSolutions.resize(nr_rows);
        for(size_t i = 0; i < nr_rows; ++i){
            if(i % 1000000 == 0 && verbose)
                verboseOutput() << i << " local solutions read on level " << this_patch << endl;
            SavedLocalSolutions[i].resize(nr_cols);
            for(size_t j = 0; j < nr_cols; ++j)
                sls >> SavedLocalSolutions[i][j];
        }
        if(sls.fail())
            throw BadInputException("Corrupt file " + file_name);
        if(verbose)
            verboseOutput() <<  SavedLocalSolutions.size() << " local solutions read on level " << this_patch << endl;
    }
    return read_local_solutions;
}

//--------------------------------------------------------------------------
// Recyccling of lattice points in patching algorithm. Unclear whetehr it has much effect
// on time or RAM usage. See 3.10.2 for straight version without recyccling.
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::store_new_vector(const vector<IntegerRet>& new_vect, const int tn){

    bool vector_available = true;

    if(FreeVectThread[tn].empty()){
        if(FreeVect.empty()){
            vector_available = false;
        }
        else{
    #pragma omp critical(FREEVECT)
            {
            if(FreeVect.empty()) // must test again because another thread may have emtied FreeVect
                vector_available = false;
            else{  // vector_available is still true
                // take 1000 vectors from FreeVect or what you can get
                auto F = FreeVect.begin();
                size_t q;
                for (q = 0; q < 1000; ++q, ++F) {
                    if (F == FreeVect.end())
                        break;
                }
                if (q < 1000)
                    FreeVectThread[tn].splice(FreeVectThread[tn].begin(), FreeVect);
                else
                    FreeVectThread[tn].splice(FreeVectThread[tn].begin(), FreeVect, FreeVect.begin(), F);
            } // FreeVect empty in critical
            } // critical
        } // FreeVect empty outer
    } // FreeVectThread empty

    if(vector_available){
        Deg1Thread[tn].splice( Deg1Thread[tn].begin(),FreeVectThread[tn], FreeVectThread[tn].begin());
        Deg1Thread[tn].front() = new_vect;
        if(new_vect[0] == 0)
            assert(false);
    }
    else
        Deg1Thread[tn].push_front(new_vect);
}

//--------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_local_solutions(const key_t this_patch,
                            list<vector<IntegerRet> >& start_list){

    const size_t coord = InsertionOrderPatches[this_patch]; // the coord marking this patch

    // for easier access and simpler code
    auto& intersection_key = AllIntersections_key[coord];
    auto& new_coords_key = AllNew_coords_key[coord];

    auto& LocalPL = AllLocalPL[coord];
    auto&& LocalSolutions_by_intersection_and_cong =
            AllLocalSolutions_by_intersection_and_cong[coord];

    auto& LocalSolutions = AllLocalSolutions[coord];
    auto& ShortLocalSolutions = AllShortLocalSolutions[coord];

    auto& CongsRestricted = AllCongsRestricted[coord];


    // Now we extend the "new" intersection coordinates by the local system
    vector<vector<IntegerRet> > LocalSolutionsNow;
    vector<vector<short> > ShortLocalSolutionsNow;
    if(stored_local_solutions && !ImportedLocalSolutions[coord]){
        if(use_short_int){
            if(verbose)
                verboseOutput() << "Importing short integers" << endl;
            stored_local_solutions = import_local_solutions(ShortLocalSolutionsNow, this_patch);
        }
        else
            stored_local_solutions = import_local_solutions(LocalSolutionsNow, this_patch);
        if(stored_local_solutions)
            ImportedLocalSolutions[coord] = true;
    }
    if(!ImportedLocalSolutions[coord]){
        LocalPL. set_startList(start_list);
        LocalPL.lift_points_to_this_dim(start_list); // computes the extensions
        if(use_short_int)
            LocalPL.put_short_deg1Points_into(ShortLocalSolutionsNow);
        else
            LocalPL.put_deg1Points_into(LocalSolutionsNow);
    }
    size_t nr_old_solutions = LocalSolutions.size();
    if(use_short_int)
        nr_old_solutions = ShortLocalSolutions.size();
    size_t nr_new_solutions = LocalSolutionsNow.size();
    if(use_short_int)
        nr_new_solutions = ShortLocalSolutionsNow.size();
    if(nr_old_solutions == 0){
        swap(LocalSolutions,LocalSolutionsNow);
        swap(ShortLocalSolutions, ShortLocalSolutionsNow);
    }
    else{
        if(use_short_int){
            ShortLocalSolutions.resize(nr_old_solutions + nr_new_solutions);
            for(size_t i = 0; i < nr_new_solutions; ++i)
                swap(ShortLocalSolutions[nr_old_solutions +i], ShortLocalSolutionsNow[i]);
        }
        else{
            LocalSolutions.resize(nr_old_solutions + nr_new_solutions);
            for(size_t i = 0; i < nr_new_solutions; ++i)
                swap(LocalSolutions[nr_old_solutions +i], LocalSolutionsNow[i]);
        }
    }
    LocalSolutionsNow.clear();
    ShortLocalSolutionsNow.clear();

    if(talkative){
            // verbose_0 = "Local solutions total " + to_string(LocalSolutions.size()) + " new " < to_string(nr_new_solutions);
            verboseOutput() << "--" << endl;
            verboseOutput() << this_patch << " / " << InsertionOrderPatches[this_patch] << " Local solutions total "
                        << nr_old_solutions + nr_new_solutions << " new " << nr_new_solutions << endl;
            // LocalSolutions.debug_print('+');
    }

    // Next the newly computed extensions are registered
    // based on overlap,
    // and for each overlap by the partial values of the congruences.
    // This allows only patchings that satisfy the congr5uences
    vector<IntegerRet> overlap(intersection_key.size());
    size_t nr_intersect = intersection_key.size();
    size_t nr_cong = CongsRestricted.size();
    vector<IntegerRet> partial_cong_values(nr_cong);
    dynamic_bitset new_coords_ind = key_to_bitset(new_coords_key, EmbDim);

    for(size_t i = nr_old_solutions; i < nr_old_solutions + nr_new_solutions; i++){ // take only the new solutions
        for(size_t j = 0; j < intersection_key.size(); ++j){
            if(use_short_int)
                overlap[j] = ShortLocalSolutions[i][j];
            else
                overlap[j] = LocalSolutions[i][j];
        }
        // insert "new" coordinates of local solution into place in full vector
        // On overlap the congruences will be evaluated below
        // For patching we need that the eval on the new coordinates = - eval on the
        // other coordiantes
        vector<IntegerRet> local_solution_new(EmbDim);
        if(use_short_int){
            for(size_t k = 0; k < new_coords_key.size(); ++k)
                local_solution_new[new_coords_key[k]] = ShortLocalSolutions[i][nr_intersect + k];
        }
        else{
            for(size_t k = 0; k < new_coords_key.size(); ++k)
                local_solution_new[new_coords_key[k]] = LocalSolutions[i][nr_intersect + k];
        }
        for(size_t k = 0; k < nr_cong; ++k){
            partial_cong_values[k] =
                eval_cong_partially(CongsRestricted[k],local_solution_new, new_coords_ind, true);
        }

        if(LocalSolutions_by_intersection_and_cong.find(overlap) != LocalSolutions_by_intersection_and_cong.end()){
            if(LocalSolutions_by_intersection_and_cong[overlap].find(partial_cong_values)
                        != LocalSolutions_by_intersection_and_cong[overlap].end()) {
                LocalSolutions_by_intersection_and_cong[overlap][partial_cong_values].push_back(i);
            }
            else{
                LocalSolutions_by_intersection_and_cong[overlap][partial_cong_values] ={};
                LocalSolutions_by_intersection_and_cong[overlap][partial_cong_values].push_back(i);
            }
        }
        else{
            LocalSolutions_by_intersection_and_cong[overlap] = {};
            LocalSolutions_by_intersection_and_cong[overlap][partial_cong_values] ={};
            LocalSolutions_by_intersection_and_cong[overlap][partial_cong_values].push_back(i);
        }
    }  // for i
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::extend_points_to_next_coord(list<vector<IntegerRet> >& LatticePoints, const key_t this_patch) {

    if(is_split_patching && only_single_point){
        double time_spent = MeasureTime(stop_ckeck_begin);
        if(time_spent > 1){
            StartTime(stop_ckeck_begin);
            string name = global_project + ".spst";
            ifstream stop(name);
            if(stop.is_open()){
                ofstream out(lat_file_name);
                string lat_type = "fusion_rings";
                if(!fusion_rings_computation)
                    lat_type = "lattice_points";
                out << lat_type << endl;
                out << "0" << endl;
                out << "0 " << endl;
                out.close();
                throw NoComputationException("Single lattice point already found");
            }
        }
    }

    const size_t coord = InsertionOrderPatches[this_patch]; // the coord marking this patch

    // for easier access and simpler code
    auto& intersection_key = AllIntersections_key[coord];
    auto& new_coords_key = AllNew_coords_key[coord];

    auto&& PolyEqusThread = AllPolyEqusThread[coord];
    auto& PolyInequsThread = AllPolyInequsThread[coord];
    auto& Automs = AllAutoms[coord];

    auto&& CoveredKey = AllCoveredKey[coord];
    auto&& CoveredKeyInverse = AllCoveredKeyInverse[coord];
    auto&  Automorphisms = fusion.Automorphisms;

    // auto& LocalPL = AllLocalPL[coord];
    auto&& LocalSolutions_by_intersection_and_cong =
            AllLocalSolutions_by_intersection_and_cong[coord];

    auto& LocalSolutions = AllLocalSolutions[coord];
    auto& ShortLocalSolutions = AllShortLocalSolutions[coord];
    dynamic_bitset covered = AllCovered[coord];

    auto& CongsRestricted = AllCongsRestricted[coord];

    /*
     * size_t max_nr_per_thread = 4000 / omp_get_max_threads(); //max_nr_new_latt_points_total/ omp_get_max_threads();
    if(distributed_computation)
        max_nr_per_thread = 10000; // to avoid very few lattice points at the lowest split level
    */

    size_t max_nr_latt_points_processed = 4000;
    if(distributed_computation)
        max_nr_latt_points_processed = 10000;

    key_t max_split_level = 0;
    if(is_split_patching){
        max_split_level = our_split.this_split_levels.back();
        if(this_patch <= max_split_level){
            prepare_split(LatticePoints, this_patch);
            // max_nr_per_thread /= 10; // we want to force a quick min_return
            max_nr_latt_points_processed /= 10;
        }
    }

    StartTime();

    size_t min_fall_back = 0;
    bool min_found = false;
    size_t max_fall_back = 0;
    for(size_t i = 0; i < this_patch; ++i){
        if(NrRemainingLP[i] > 0){
            max_fall_back = i;
            if(!min_found){
                min_found = true;
                min_fall_back = i;
            }
        }
    }

    if(min_fall_back == 0)
        sort_lattice_points_by_overlap(LatticePoints, intersection_key);

    INTERRUPT_COMPUTATION_BY_EXCEPTION

    // We extract the "intersection coordinates" from the LatticePoints
    // and extend them to solutions of the local system LocalPL
    // but only those whose extension has not yet been computed

    set< vector<IntegerRet> > LatticePoints_restricted_to_intersection;
    for(auto& P: LatticePoints){
        LatticePoints_restricted_to_intersection.insert(
                                v_select_coordinates(P,intersection_key));
    }
    list<vector<IntegerRet> > start_list;
    start_list.insert(start_list.begin(),LatticePoints_restricted_to_intersection.begin(),
                      LatticePoints_restricted_to_intersection.end());

    // Now we have the "intersection coordinates" of all LatticePoints.
    // We must remove those from start_list whose extension has already been computed
    // and initialize the others.

    // Matrix<IntegerRet>(start_list).debug_print();

    for(auto IC = start_list.begin(); IC != start_list.end(); ){
        if(LocalSolutions_by_intersection_and_cong.find(*IC) == LocalSolutions_by_intersection_and_cong.end()){
            LocalSolutions_by_intersection_and_cong[*IC] = {};
            IC++;
        }
        else{
            IC =  start_list.erase(IC);
        }
    }

    dynamic_bitset full_coords_ind(EmbDim);
    full_coords_ind.flip();

    compute_local_solutions(this_patch, start_list);

    bool last_coord = (coord == InsertionOrderPatches.back());

    size_t nr_to_match = LatticePoints.size();
    size_t nr_points_matched = 0;

    // In the following the oputer paralleliozed loop is over the lattice points
    // given to the routine, and the inner is over the extensions along the overlao.
    // One could think about exchanging the order of the loops to get better
    // parallelization. But it is unclear whether this can be achieved.

    struct timeval time_begin;
    StartTime(time_begin);

    size_t nr_rounds = 0;

    bool apply_automorphisms = true;
    if(LatticePoints.size() < 100 && Automorphisms.size() > 1000  && !last_coord) {
        apply_automorphisms = false;
        automs_minimized[coord] = true;
    }

    while (true) {

        /* if(talkative){
            if(nr_rounds > 0 && verbose){
                verboseOutput() << "----" << endl;
            }
        }*/

        nr_rounds++;
        if(GlobalTimeBound > 0 &&  TimeSinceStart() > GlobalTimeBound){
            throw TimeBoundException("while patching");
        }

        Check_Stop();

        struct timeval step_time_begin;
        StartTime(step_time_begin);

        if(talkative){
            verboseOutput() << "----" << endl;
            verboseOutput() <<  LevelPatches[coord] << " / " << coord << " left " << nr_to_match - nr_points_matched;
                if(min_fall_back > 0)
                    verboseOutput() << " min " << min_fall_back << " max " << max_fall_back;
            verboseOutput()    << endl;
        }

        bool skip_remaining;
        std::exception_ptr tmp_exception;

        skip_remaining = false;
        int omp_start_level = omp_get_level();


        size_t nr_new_latt_points = 0;  // counts number of new lattice points produced in the loop
        size_t nr_extensions = 0;  //statistics for this run of the while loop
        size_t nr_caught_by_restricted = 0;  //restricted inequalities
        size_t nr_caught_by_simplicity = 0;  // caught by simplicity check
        size_t nr_caught_by_equations = 0;  //statistics for this run of the while loop
        size_t nr_caught_by_automs = 0;  //statistics for this run of the while loop
        size_t nr_points_done_in_this_round = 0;

        vector<vector<size_t> > poly_equs_stat; // TODO make macro
        vector<size_t> poly_equs_stat_total;
        if(!poly_equs_minimized[coord]){
            poly_equs_stat.resize(omp_get_max_threads());  // counts the number of "successful". i.e. != 0,
            for(auto& p:poly_equs_stat)                    // evaluations of a mpolynomial erquation
                p.resize(PolyEqusThread[0].size());
            poly_equs_stat_total.resize(poly_equs_stat[0].size());
        }

        vector<vector<size_t> > poly_inequs_stat;
        vector<size_t> poly_inequs_stat_total;
        if(!poly_inequs_minimized[coord]){
            poly_inequs_stat.resize(omp_get_max_threads());  // counts the number of "successful". i.e. < 0,
            for(auto& p:poly_inequs_stat)                    // evaluations of a restricted mpolynomial inequalities
                p.resize(PolyInequsThread[0].size());
            poly_inequs_stat_total.resize(poly_inequs_stat[0].size());
        }

        vector<vector<size_t> > automs_stat;
        vector<size_t> automs_stat_total;
        if(!automs_minimized[coord]){
            automs_stat.resize(omp_get_max_threads());  // counts the number of "successful". i.e. < 0,
            for(auto& p:automs_stat)                    // evaluations of a restricted mpolynomial inequalities
                p.resize(Automs.size());
            automs_stat_total.resize(automs_stat[0].size());
        }

        size_t nr_intersect = intersection_key.size();
        dynamic_bitset full_support(EmbDim);
        full_support.flip();

        size_t nr_latt_points_processed = 0;

#pragma omp parallel
        {

        vector<IntegerRet> overlap(nr_intersect);
        vector<IntegerRet> old_cong(CongsRestricted.size());

        vector<IntegerRet> NewLattPoint(EmbDim);

        vector<IntegerRet> restricted(CoveredKey.size());
        vector<IntegerRet> restricted_conjugate(CoveredKey.size());

        int tn;
        if (omp_get_level() == omp_start_level)
            tn = 0;
        else
            tn = omp_get_ancestor_thread_num(omp_start_level + 1);

        // size_t nr_points_in_thread = 0;

        size_t ppos = 0;
        auto P = LatticePoints.begin();

        list<pair<key_t, IntegerRet> > order_poly_equs; // suddessful constraints to the front !!
        for(key_t k = 0; k < PolyEqusThread[tn].size(); ++k) // the nsecond component will contain values
            order_poly_equs.push_back(make_pair(k,0));       // of the restrictions of the polynomials
        list<key_t> order_poly_inequs;
        for(key_t k = 0; k < PolyInequsThread[tn].size(); ++k)
            order_poly_inequs.push_back(k);
        list<key_t> order_automs;
        for(key_t k = 0; k < Automs.size(); ++k)
            order_automs.push_back(k);

#pragma omp for schedule(dynamic)
        for (size_t ppp = 0; ppp < nr_to_match; ++ppp) {

            if (skip_remaining)
                continue;

            for (; ppp > ppos; ++ppos, ++P)
                ;
            for (; ppp < ppos; --ppos, --P)
                ;

            if ((*P)[0] == 0)  // point done
                continue;

#pragma omp atomic
            nr_points_matched++;

#pragma omp atomic
            nr_points_done_in_this_round++;

        NewLattPoint = *P;

        overlap = v_select_coordinates(NewLattPoint, intersection_key);
        for(size_t k = 0; k < CongsRestricted.size(); ++k)
            old_cong[k]  = eval_cong_partially(CongsRestricted[k],NewLattPoint, full_coords_ind, false);

        if(LocalSolutions_by_intersection_and_cong[overlap].find(old_cong)
                    == LocalSolutions_by_intersection_and_cong[overlap].end()){
            (*P)[0] = 0;
            continue;
        }

        try{  // now the given lattice points are extended along their overlaps and congruence values

            if(ppp % 10000 == 0 && GlobalTimeBound > 0 &&  TimeSinceStart() > GlobalTimeBound){
                throw TimeBoundException("extending");
            }

            // the - is important because the sum of the two parts should give 0 <==> equation satisfied
            for(auto pp = order_poly_equs.begin(); pp!= order_poly_equs.end(); ++pp){
                pp-> second = - PolyEqusThread[tn][pp->first].first.evaluate(NewLattPoint);
            }

            // now the extensions of the overlap
            for(auto& i: LocalSolutions_by_intersection_and_cong[overlap][old_cong]){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

                if(use_short_int){
                    for(size_t j = 0; j < new_coords_key.size(); ++j)
                        NewLattPoint[new_coords_key[j]] = ShortLocalSolutions[i][j + intersection_key.size()];
                }
                else{
                    for(size_t j = 0; j < new_coords_key.size(); ++j)
                        NewLattPoint[new_coords_key[j]] = LocalSolutions[i][j + intersection_key.size()];

                }

                bool can_be_inserted = true;
#pragma omp atomic
                nr_extensions++;

                if(check_simplicity_cand && critical_coord_simplicity == coord){
                    if(!fusion.simplicity_check(fusion.coords_to_check_key[0], NewLattPoint)){
                        can_be_inserted = false;
#pragma omp atomic
                        nr_caught_by_simplicity++;
                    }
                }
                if(can_be_inserted && check_simplicity_all){
                    if(!fusion.simplicity_check(fusion.all_critical_coords_keys[coord], NewLattPoint)){
                        can_be_inserted = false;
#pragma omp atomic
                        nr_caught_by_simplicity++;
                    }
                }
                if(can_be_inserted){
                    for(auto pp = order_poly_equs.begin(); pp!= order_poly_equs.end(); ++pp){
                        if(pp->second != PolyEqusThread[tn][pp->first].second.evaluate(NewLattPoint)){
                            can_be_inserted = false;
#pragma omp atomic
                            nr_caught_by_equations++;
                            if(!poly_equs_minimized[coord])
                                poly_equs_stat[tn][pp->first]++;
                           /* else
                                order_poly_equs.splice(order_poly_equs.begin(), order_poly_equs, pp);*/
                            break;
                        }
                    }
                }
                if(can_be_inserted){
                    for(auto pp = order_poly_inequs.begin(); pp!= order_poly_inequs.end(); ++pp){
                        if(PolyInequsThread[tn][*pp].evaluate(NewLattPoint) < 0){
                            can_be_inserted = false;
#pragma omp atomic
                            nr_caught_by_restricted++;
                            if(!poly_inequs_minimized[coord])
                                poly_inequs_stat[tn][*pp]++;
                            else
                                order_poly_inequs.splice(order_poly_inequs.begin(), order_poly_inequs, pp);
                            break;
                        }
                    }
                }
               if(can_be_inserted && fusion.use_automorphisms && apply_automorphisms){
                    for(size_t i = 0; i < CoveredKey.size(); ++i)
                        restricted[i] = NewLattPoint[CoveredKey[i]];
                    for(auto pp = order_automs.begin(); pp!= order_automs.end(); ++pp){
                        for(auto& r: restricted_conjugate)
                            r = 0;
                        for(size_t i = 0; i < CoveredKey.size(); ++i){
                            key_t image = Automorphisms[Automs[*pp]][CoveredKey[i]];
                            if(covered[image]){
                                key_t inverse_image = CoveredKeyInverse[image];
                                restricted_conjugate[inverse_image] = restricted[i];
                            }
                        }
                        if(restricted < restricted_conjugate){
                            can_be_inserted = false;
#pragma omp atomic
                           nr_caught_by_automs++;
                           if(!automs_minimized[coord])
                                automs_stat[tn][*pp]++;
                            else
                                order_automs.splice(order_automs.begin(), order_automs, pp);
                            break;
                        }
                    }
                }
                if(can_be_inserted){
                    // nr_points_in_thread++;
#pragma omp atomic
                    nr_new_latt_points++;
                    if(last_coord)
                        finalize_latt_point(NewLattPoint, tn);
                    else{
                        store_new_vector(NewLattPoint,tn);
                        // Deg1Thread[tn].push_back(NewLattPoint); // straight version
                    }
                }

            } // for i (inner for loop)

            (*P)[0] = 0;  // mark point as done

#pragma omp atomic
            nr_latt_points_processed++;
            if(this_patch >= max_split_level &&  nr_new_latt_points > max_nr_new_latt_points_patching && !last_coord) {  // thread is full and we are allowed to break
                skip_remaining = true;

#pragma omp flush(skip_remaining)
            }

        } catch (const std::exception&) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
        } // for ppp (outer for loop)
        } // parallel

        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);

        list<vector<IntegerRet> > NewLatticePoints;  // lattice points computed in this round

        if(!poly_equs_minimized[coord]){
            for(size_t i = 0; i< poly_equs_stat.size(); ++i){
                    for(size_t j = 0; j < poly_equs_stat[0].size(); ++j)
                        poly_equs_stat_total[j] += poly_equs_stat[i][j];
            }
        }

        if(!poly_inequs_minimized[coord]){
            for(size_t i = 0; i< poly_equs_stat.size(); ++i){
                    for(size_t j = 0; j < poly_inequs_stat[0].size(); ++j)
                        poly_inequs_stat_total[j] += poly_inequs_stat[i][j];
            }
        }

        assert(NewLatticePoints.empty());

        for (size_t i = 0; i < Deg1Thread.size(); ++i)
            NewLatticePoints.splice(NewLatticePoints.end(), Deg1Thread[i]);

        if(!automs_minimized[coord]){
            for(size_t i = 0; i< automs_stat.size(); ++i){
                    for(size_t j = 0; j < automs_stat[0].size(); ++j)
                        automs_stat_total[j] += automs_stat[i][j];
            }
        }

        if(last_coord)
            collect_results(NewLatticePoints); // clears NewLatticePoints

        NrDoneLP[coord] += nr_extensions;

        // *******************************************
        // Heuristic minimization of equations etc.
        // *******************************************
        if(!poly_equs_minimized[coord] &&  nr_extensions > nr_extensions_for_elimination_equs){
            poly_equs_minimized[coord] = true;
             vector < pair<OurPolynomial<IntegerRet>, OurPolynomial<IntegerRet> > > EffectivePolys;
            for(size_t i = 0; i < PolyEqusThread[0].size(); ++i){
                if(poly_equs_stat_total[i] > 0){
                    EffectivePolys.push_back(PolyEqusThread[0][i]);
                }
            }

            if(talkative && PolyEqusThread[0].size() > 0){
                verboseOutput() << LevelPatches[coord] << " / " << coord << " active polynomial equations " << EffectivePolys.size() << " of " <<  PolyEqusThread[0].size() << endl;
            }

            for(size_t thr = 0; thr < PolyEqusThread.size(); ++thr)
                    PolyEqusThread[thr] = EffectivePolys;
        }

        if(!automs_minimized[coord] &&  nr_extensions > nr_extensions_for_elimination_automs &&!last_coord){
            automs_minimized[coord] = true;
            vector<key_t> EffectiveAutoms;
            for(size_t i = 0; i < Automs.size(); ++i){
                if(automs_stat_total[i] > 0)
                    EffectiveAutoms.push_back(Automs[i]);
            }

            if(talkative && Automs.size() > 0){
                verboseOutput() << LevelPatches[coord] << " / " << coord << " active automorphisms "
                  << EffectiveAutoms.size() << " of " << Automorphisms.size() << endl;
            }
             Automs = EffectiveAutoms;
        }

        //Discard ineffective restrictable plolynjomial inequalities
        if(!poly_inequs_minimized[coord] &&  nr_extensions > nr_extensions_for_elimination_inequs){
            poly_inequs_minimized[coord] = true;
            OurPolynomialSystem<IntegerRet> EffectivePolyInequs;
            for(size_t i = 0; i < PolyInequsThread[0].size(); ++i){
                if(poly_inequs_stat_total[i] > 0)
                    EffectivePolyInequs.push_back(PolyInequsThread[0][i]);
            }

            if(talkative && PolyInequsThread[0].size() > 0){
                verboseOutput() << LevelPatches[coord] << " / " << coord << " active restricted polynomial inequalities "
                  << EffectivePolyInequs.size() << " of " <<  PolyInequsThread[0].size() << endl;
            }

             for(size_t thr = 0; thr < PolyEqusThread.size(); ++thr)
                PolyInequsThread[thr] = EffectivePolyInequs;
        }

       NrRemainingLP[this_patch] = nr_to_match - nr_points_matched;

       if(distributed_computation && NrRemainingLP[this_patch] > 0){
            write_control_file(this_patch, LatticePoints.size());
            throw NoComputationException("No output with DistribitedComp for patching");
       }

       if(min_fall_back == 0 && NrRemainingLP[this_patch] == 0){ // no return to this level
            LocalSolutions_by_intersection_and_cong.clear();  // save memory
            LocalSolutions.clear();
            ShortLocalSolutions.clear();
            poly_equs_stat.clear();
            poly_equs_stat_total.clear();
            poly_inequs_stat.clear();
            poly_inequs_stat_total.clear();
            automs_stat.clear();
            automs_stat_total.clear();
       }

        double expected_number_of_rounds = NrRemainingLP[this_patch];
        if(nr_points_done_in_this_round > 0){
            expected_number_of_rounds/= nr_points_done_in_this_round;
        }

        if(talkative){
            verboseOutput() << "done " << nr_points_done_in_this_round;
            verboseOutput() << " ext " << nr_extensions;
            if(PolyEqusThread[0].size() > 0)
                verboseOutput() << " equ " << nr_caught_by_equations;
            if(PolyInequsThread[0].size() > 0)
                verboseOutput() << " ine " << nr_caught_by_restricted;
            if(nr_caught_by_simplicity > 0)
                verboseOutput() << " smp " << nr_caught_by_simplicity;
            if(nr_caught_by_automs > 0)
                verboseOutput() << " aut " << nr_caught_by_automs;
            if(nr_points_done_in_this_round > 0)
                verboseOutput() << " exp rnd " << expected_number_of_rounds;
            ExpectedNrRounds[this_patch] = expected_number_of_rounds;
            verboseOutput() << endl;
        }
        // we write a file that preserves what has been done. First we tecord the lowest patch
        // to which we must return
        if(NrRemainingLP[this_patch] > 0 && min_fall_back == 0 && min_return_patch == 0){  // all previous patches were fully done
            min_return_patch = this_patch;
            if(is_split_patching){
                ofstream prel_data;
                prel_data.open(lat_file_name, ofstream::app);
                if(fusion.use_automorphisms){
                    if(check_simplicity_all)
                        prel_data << "simple_";
                    prel_data << "fusion_rings" << endl;
                }
                else{
                    prel_data << "lattice_points" << endl;
                }
                prel_data << endl << "min_return " << min_return_patch << endl << endl;

                prel_data << "total_indices " << LatticePoints.size() << endl;

                prel_data << "done_indices" << endl;
                prel_data << "0" << endl;
                prel_data << "found_solutions" << endl;
                prel_data << "0" << endl;
                prel_data << EmbDim << endl;
                prel_data.close();
            }
        }

        double time_to_ascent = 0;

        if(!last_coord && NewLatticePoints.size() > 0){ // we must go up
            if(verbose && !talkative){
                verboseOutput() << "+" << flush;
                verb_length++;
                if(verb_length == 50){
                    verboseOutput() << endl;
                    verb_length = 0;
                }
            }

            if(NrNodes[this_patch + 1] == 1 || this_patch == min_return_patch)
                NrNodes[this_patch +1] = NrNodes[this_patch] * (expected_number_of_rounds +1);

            for(size_t i = this_patch + 2; i < NrNodes.size(); ++i)
                NrNodes[i] = 1;

            // ***********************************   ascent to next patch
            time_to_ascent = MeasureTime(time_begin);
            extend_points_to_next_coord(NewLatticePoints, this_patch + 1);
            // ****************************************** down from next patch

            if(verbose && !talkative){
                verboseOutput() << "-" << flush;
                verb_length++;
                if(verb_length == 50){
                    verboseOutput() << endl;
                    verb_length = 0;
                 }
            }
        }

        FreeVect.splice(FreeVect.end(), NewLatticePoints);
        // NewLatticePoints.clear(); // straight version

        if(single_point_found)
            break;

        double expected_time = 0; // to make gcc happy
        double total_expected_time = 0;

        bool time_measured = false;


        if(nr_points_done_in_this_round > 0 && NrRemainingLP[this_patch] > 0 && nr_rounds == 1){
            time_measured = true;
            double time_spent = MeasureTime(time_begin);
            expected_time = time_spent*expected_number_of_rounds;

            // -------------------------------

            TimeToLevel[this_patch + 1] = NrNodes[this_patch] * time_to_ascent + TimeToLevel[this_patch];
            total_expected_time = TimeToLevel[this_patch] + NrNodes[this_patch + 1] * time_spent + TimeSinceStart();
            // -------------------------
        }

        if(time_measured){
            if(((talkative || nr_time_printed <= 10) ||  this_patch == min_return_patch)
                            && verbose && (nr_rounds ==1 || this_patch == min_return_patch) ){
                nr_time_printed++;
                if(verbose){
                    verboseOutput() << "---------" << endl;
                    verboseOutput() << "expected future time on level  " << LevelPatches[coord] << "  " << expected_time;
                    verboseOutput() << " / total " << total_expected_time << endl;
                }
            }
        }

        if(GlobalPredictionTimeBound > 0 && total_expected_time > GlobalPredictionTimeBound){
            errorOutput() << "expected time exceeds bound of " << GlobalPredictionTimeBound << " sec" << endl;
            throw TimeBoundException("patching");
        }

        if(nr_points_matched == nr_to_match)
            break;

        // Here we record what has been done on the lowest level to which we reurn for further work
        // this information is compeletely written whenever qwe return here.
        if(is_split_patching && this_patch == min_return_patch){

            ofstream prel_data;
            prel_data.open(lat_file_name, ofstream::app);

            size_t counter = 0;
            for(auto& p: LatticePoints){  // we only register the first "block"; there may be sporadic further done indices
                if(p[0] != 0){
                    break;
                }
                counter++;
            }
            prel_data << "done_indices " << counter << endl;
            if(verbose){
                verboseOutput() << endl << "done_indices " << counter << " of " << LatticePoints.size() << endl;
            }
            auto check = Deg1Points;
            check.sort();
            check.unique();
            Matrix<IntegerRet> FoundSolutions(check);
            if(FoundSolutions.nr_of_rows() > 0){
                FoundSolutions.cyclic_shift_left(FoundSolutions.nr_of_columns()-1); // to get final coordinates
                /* if(fusion.use_automorphisms)
                    fusion.do_iso_classes(FoundSolutions); // to use the natural order of coordiantes for lex */
            }

            prel_data << "found_solutions" << endl;
            prel_data << FoundSolutions.nr_of_rows() << endl;
            prel_data << EmbDim << endl;
            FoundSolutions.pretty_print(prel_data);
            prel_data.close();
        }

    }  // while not done

    return;
}
//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
vector<size_t> ProjectAndLift<IntegerPL, IntegerRet>::order_supps(const Matrix<IntegerPL>& Supps) {
    assert(Supps.nr_of_rows() > 0);
    size_t dim = Supps.nr_of_columns();

    vector<pair<nmz_float, size_t> > NewPos, NewNeg, NewNeutr;  // to record the order of the support haperplanes
    for (size_t i = 0; i < Supps.nr_of_rows(); ++i) {
        if (Supps[i][dim - 1] == 0) {
            NewNeutr.push_back(make_pair(0.0, i));
            continue;
        }
        nmz_float num, den;
        convert(num, Supps[i][0]);
        convert(den, Supps[i][dim - 1]);
        nmz_float quot = num / den;
        if (Supps[i][dim - 1] > 0)
            NewPos.push_back(make_pair(Iabs(quot), i));
        else
            NewNeg.push_back(make_pair(Iabs(quot), i));
    }
    sort(NewPos.begin(), NewPos.end());
    sort(NewNeg.begin(), NewNeg.end());
    NewPos.insert(NewPos.end(), NewNeutr.begin(), NewNeutr.end());

    size_t min_length = NewNeg.size();
    if (NewPos.size() < min_length)
        min_length = NewPos.size();

    vector<size_t> Order;

    for (size_t i = 0; i < min_length; ++i) {
        Order.push_back(NewPos[i].second);
        Order.push_back(NewNeg[i].second);
    }
    for (size_t i = min_length; i < NewPos.size(); ++i)
        Order.push_back(NewPos[i].second);
    for (size_t i = min_length; i < NewNeg.size(); ++i)
        Order.push_back(NewNeg[i].second);

    assert(Order.size() == Supps.nr_of_rows());

    return Order;
}

//---------------------------------------------------------------------------

// not used
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::reorder_coordinates(){

    bool has_poly_equs = (PolyEquations.size() > 0);
    bool has_poly_inequs = (PolyInequalities.size() > 0);
    if(!has_poly_equs && !has_poly_inequs)
        return;

    size_t dim = AllSupps[EmbDim][0].size();

    dynamic_bitset covered(dim);
    covered[0] = true; // coordinate 0 is fixed
    vector<key_t> NewOrder(1,0);
    dynamic_bitset new_covered(dim);

    while(true){
        bool first = true;
        dynamic_bitset test_covered(dim);
        size_t nr_new_covered = 0;
        for(size_t i = 0; i < PolyEquations.size(); ++i){
            if(PolyEquations[i].support.is_subset_of(covered))
                continue;
            test_covered = covered | PolyEquations[i].support;
            if(first || test_covered.count() < nr_new_covered){
                new_covered = test_covered;
                nr_new_covered = new_covered.count();
            }
        }
        if(nr_new_covered == 0)
            break;
        for(size_t j = 0; j < dim; ++j){
            if(covered[j] || !new_covered[j])
                continue;
            NewOrder.push_back(j);
        }
        covered = new_covered;
    }
    for(size_t i = 0; i < dim; ++i){
        if(!new_covered[i])
            NewOrder.push_back(i);
    }
    AllSupps[EmbDim].permute_columns(NewOrder);
    PolyEquations.permute_variables(NewOrder);
    PolyInequalities.permute_variables(NewOrder);
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_projections_primitive(size_t dim){

    size_t dim1 = dim - 1;

    if (dim == 1)
        return;

    Matrix<IntegerPL> SuppsProj(0,dim1);

    for(size_t i = 0; i< AllSupps[EmbDim].nr_of_rows(); ++i){
        if(AllSupps[EmbDim][i][0] < 0){
            bool unsolvable = true;
            for(size_t j = 1; j < AllSupps[EmbDim][i].size(); ++j){
                if(AllSupps[EmbDim][i][j] > 0){
                    unsolvable = false;
                    break;
                }
            }
            if(unsolvable){
                system_unsolvable = true;
                return;
            }
        }
    }


    for(size_t i = 0; i< AllSupps[EmbDim].nr_of_rows(); ++i){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        bool can_be_restricted = true;
        for(size_t j= dim1; j <= EmbDim-1; ++j){
            if(AllSupps[EmbDim][i][j] >0){
                can_be_restricted = false;
                break;
            }
        }
        if(can_be_restricted){
            vector<IntegerPL> Restriction = AllSupps[EmbDim][i];
            Restriction.resize(dim1);
            SuppsProj.append(Restriction);
        }
    }

    SuppsProj.remove_duplicate_and_zero_rows();

    if (verbose)
        verboseOutput() << "embdim " << dim << " inequalities " << SuppsProj.nr_of_rows() << endl;

    AllOrders[dim1] = order_supps(SuppsProj);
    swap(AllSupps[dim1], SuppsProj);
    compute_projections_primitive(dim1);
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_projections(size_t dim,
                                                                size_t down_to,
                                                                vector<dynamic_bitset>& Ind,
                                                                vector<dynamic_bitset>& Pair,
                                                                vector<dynamic_bitset>& ParaInPair,
                                                                size_t rank,
                                                                bool only_projections) {
    INTERRUPT_COMPUTATION_BY_EXCEPTION

    const Matrix<IntegerPL>& Supps = AllSupps[dim];

    size_t dim1 = dim - 1;

    if (verbose)
        verboseOutput() << "embdim " << dim << " inequalities " << Supps.nr_of_rows() << endl;

    if (dim == down_to)
        return;

    // We now augment the given cone by the last basis vector and its negative
    // Afterwards we project modulo the subspace spanned by them

    vector<key_t> Neg, Pos;               // for the Fourier-Motzkin elimination of inequalities
    Matrix<IntegerPL> SuppsProj(0, dim);  // for the support hyperplanes of the projection
    Matrix<IntegerPL> EqusProj(0, dim);   // for the equations (both later minimized)

    // First we make incidence vectors with the given generators
    vector<dynamic_bitset> NewInd;         // for the incidence vectors of the new hyperplanes
    vector<dynamic_bitset> NewPair;        // for the incidence vectors of the new hyperplanes
    vector<dynamic_bitset> NewParaInPair;  // for the incidence vectors of the new hyperplanes

    dynamic_bitset TRUE;
    if (!is_parallelotope && !primitive) {
        TRUE.resize(Ind[0].size());
        TRUE.set();
    }

    vector<bool> IsEquation(Supps.nr_of_rows());

    bool rank_goes_up = false;  // if we add the last unit vector
    size_t PosEquAt = 0;        // we memorize the positions of pos/neg equations if rank goes up
    size_t NegEquAt = 0;

    for (size_t i = 0; i < Supps.nr_of_rows(); ++i) {
        if (!is_parallelotope && Ind[i] == TRUE)
            IsEquation[i] = true;

        if (Supps[i][dim1] == 0) {  // already independent of last coordinate
            no_crunch = false;
            if (IsEquation[i])
                EqusProj.append(Supps[i]);  // is equation
            else {
                SuppsProj.append(Supps[i]);  // neutral support hyperplane
                if (!is_parallelotope)
                    NewInd.push_back(Ind[i]);
                else {
                    NewPair.push_back(Pair[i]);
                    NewParaInPair.push_back(ParaInPair[i]);
                }
            }
            continue;
        }
        if (IsEquation[i])
            rank_goes_up = true;
        if (Supps[i][dim1] > 0) {
            if (IsEquation[i])
                PosEquAt = i;
            Pos.push_back(static_cast<key_t>(i));
            continue;
        }
        Neg.push_back(static_cast<key_t>(i));
        if (IsEquation[i])
            NegEquAt = i;
    }

    // now the elimination, patching Pos and Neg

    bool skip_remaining;
    std::exception_ptr tmp_exception;

    if (rank_goes_up) {
        assert(!is_parallelotope);

        for (size_t p : Pos) {  // match pos and neg equations
            if (!IsEquation[p])
                continue;
            IntegerPL PosVal = Supps[p][dim1];
            for (size_t n : Neg) {
                if (!IsEquation[n])
                    continue;
                IntegerPL NegVal = Supps[n][dim1];
                bool is_zero;
                vector<IntegerPL> new_equ = FM_comb(PosVal, Supps[n], NegVal, Supps[p], is_zero);
                if (is_zero)
                    continue;
                EqusProj.append(new_equ);
            }
        }

        for (size_t p : Pos) {  // match pos inequalities with a negative equation
            if (IsEquation[p])
                continue;
            IntegerPL PosVal = Supps[p][dim1];
            IntegerPL NegVal = Supps[NegEquAt][dim1];
            vector<IntegerPL> new_supp(dim);
            bool is_zero;
            new_supp = FM_comb(PosVal, Supps[NegEquAt], NegVal, Supps[p], is_zero);
            if (is_zero)  // cannot happen, but included for analogy
                continue;
            SuppsProj.append(new_supp);
            NewInd.push_back(Ind[p]);
        }

        for (size_t n : Neg) {  // match neg inequalities with a positive equation
            if (IsEquation[n])
                continue;
            IntegerPL PosVal = Supps[PosEquAt][dim1];
            IntegerPL NegVal = Supps[n][dim1];
            vector<IntegerPL> new_supp(dim);
            bool is_zero;
            new_supp = FM_comb(PosVal, Supps[n], NegVal, Supps[PosEquAt], is_zero);

            if (is_zero)  // cannot happen, but included for analogy
                continue;
            SuppsProj.append(new_supp);
            NewInd.push_back(Ind[n]);
        }
    }

    if (!rank_goes_up && !is_parallelotope) {  // must match pos and neg hyperplanes

        skip_remaining = false;

        size_t min_nr_vertices = rank - 2;

#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < Pos.size(); ++i) {
            if (skip_remaining)
                continue;

            try {
                size_t p = Pos[i];
                IntegerPL PosVal = Supps[p][dim1];
                vector<key_t> PosKey;
                for (size_t k = 0; k < Ind[i].size(); ++k)
                    if (Ind[p][k])
                        PosKey.push_back(static_cast<key_t>(k));

                for (size_t n : Neg) {
                    INTERRUPT_COMPUTATION_BY_EXCEPTION

                    // // to give a facet of the extended cone
                    // match incidence vectors
                    dynamic_bitset incidence(TRUE.size());
                    size_t nr_match = 0;
                    vector<key_t> CommonKey;
                    for (unsigned int k : PosKey)
                        if (Ind[n][k]) {
                            incidence[k] = true;
                            CommonKey.push_back(k);
                            nr_match++;
                        }
                    if (rank >= 2 && nr_match < min_nr_vertices)  // cannot make subfacet of augmented cone
                        continue;

                    bool IsSubfacet = true;
                    for (size_t k = 0; k < Supps.nr_of_rows(); ++k) {
                        if (k == p || k == n || IsEquation[k])
                            continue;
                        bool contained = true;
                        for (unsigned int j : CommonKey) {
                            if (!Ind[k][j]) {
                                contained = false;
                                break;
                            }
                        }
                        if (contained) {
                            IsSubfacet = false;
                            break;
                        }
                    }
                    if (!IsSubfacet)
                        continue;
                    //}

                    IntegerPL NegVal = Supps[n][dim1];
                    vector<IntegerPL> new_supp(dim);
                    bool is_zero;
                    new_supp = FM_comb(PosVal, Supps[n], NegVal, Supps[p], is_zero);
                    if (is_zero)  // linear combination is 0
                        continue;

                    if (nr_match == TRUE.size()) {  // gives an equation
#pragma omp critical(NEWEQ)
                        EqusProj.append(new_supp);
                        continue;
                    }
#pragma omp critical(NEWSUPP)
                    {
                        SuppsProj.append(new_supp);
                        NewInd.push_back(incidence);
                    }
                }

            } catch (const std::exception&) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
        }

    }  // !rank_goes_up && !is_parallelotope

    if (!(tmp_exception == 0))
        std::rethrow_exception(tmp_exception);

    if (!rank_goes_up && is_parallelotope) {  // must match pos and neg hyperplanes

        size_t codim = dim1 - 1;  // the minimal codim a face of the original cone must have
                                  // in order to project to a subfacet of the current one
        size_t original_dim = Pair[0].size() + 1;
        size_t max_number_containing_factes = original_dim - codim;

        skip_remaining = false;

        size_t nr_pos = Pos.size();
        size_t nr_neg = Neg.size();

#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < nr_pos; ++i) {
            if (skip_remaining)
                continue;

            try {
                INTERRUPT_COMPUTATION_BY_EXCEPTION

                size_t p = Pos[i];
                IntegerPL PosVal = Supps[p][dim1];

                for (size_t j = 0; j < nr_neg; ++j) {
                    size_t n = Neg[j];
                    dynamic_bitset IntersectionPair(Pair[p].size());
                    size_t nr_hyp_intersection = 0;
                    bool in_parallel_hyperplanes = false;
                    bool codim_too_small = false;

                    for (size_t k = 0; k < Pair[p].size(); ++k) {  // run over all pairs
                        if (Pair[p][k] || Pair[n][k]) {
                            nr_hyp_intersection++;
                            IntersectionPair[k] = true;
                            if (nr_hyp_intersection > max_number_containing_factes) {
                                codim_too_small = true;
                                break;
                            }
                        }
                        if (Pair[p][k] && Pair[n][k]) {
                            if (ParaInPair[p][k] != ParaInPair[n][k]) {
                                in_parallel_hyperplanes = true;
                                break;
                            }
                        }
                    }
                    if (in_parallel_hyperplanes || codim_too_small)
                        continue;

                    dynamic_bitset IntersectionParaInPair(Pair[p].size());
                    for (size_t k = 0; k < ParaInPair[p].size(); ++k) {
                        if (Pair[p][k])
                            IntersectionParaInPair[k] = ParaInPair[p][k];
                        else if (Pair[n][k])
                            IntersectionParaInPair[k] = ParaInPair[n][k];
                    }

                    // we must nevertheless use the comparison test
                    bool IsSubfacet = true;
                    if (!no_crunch) {
                        for (size_t k = 0; k < Supps.nr_of_rows(); ++k) {
                            if (k == p || k == n || IsEquation[k])
                                continue;
                            bool contained = true;

                            for (size_t u = 0; u < IntersectionPair.size(); ++u) {
                                if (Pair[k][u] && !IntersectionPair[u]) {  // hyperplane k contains facet of Supp
                                    contained = false;                     // not our intersection
                                    continue;
                                }
                                if (Pair[k][u] && IntersectionPair[u]) {
                                    if (ParaInPair[k][u] != IntersectionParaInPair[u]) {  // they are contained in parallel
                                        contained = false;                                // original facets
                                        continue;
                                    }
                                }
                            }

                            if (contained) {
                                IsSubfacet = false;
                                break;
                            }
                        }
                    }
                    if (!IsSubfacet)
                        continue;

                    IntegerPL NegVal = Supps[n][dim1];
                    bool dummy;
                    vector<IntegerPL> new_supp = FM_comb(PosVal, Supps[n], NegVal, Supps[p], dummy);
#pragma omp critical(NEWSUPP)
                    {
                        SuppsProj.append(new_supp);
                        NewPair.push_back(IntersectionPair);
                        NewParaInPair.push_back(IntersectionParaInPair);
                    }
                }

            } catch (const std::exception&) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
        }

        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);

    }  // !rank_goes_up && is_parallelotope

    Ind.clear();  // no longer needed

    EqusProj.resize_columns(dim1);   // cut off the trailing 0
    SuppsProj.resize_columns(dim1);  // project hyperplanes

    // Equations have not yet been appended to support hypwerplanes
    EqusProj.row_echelon();  // reduce equations
    SuppsProj.append(EqusProj);  // append them as pairs of inequalities
    EqusProj.scalar_multiplication(-1);
    SuppsProj.append(EqusProj);
    AllNrEqus[dim1] = EqusProj.nr_of_rows();
    // We must add indicator vectors for the equations
    for (size_t i = 0; i < 2 * EqusProj.nr_of_rows(); ++i)
        NewInd.push_back(TRUE);

    if (dim1 > 1 && !only_projections)
        AllOrders[dim1] = order_supps(SuppsProj);
    swap(AllSupps[dim1], SuppsProj);

    size_t new_rank = dim1 - EqusProj.nr_of_rows();

    compute_projections(dim - 1, down_to, NewInd, NewPair, NewParaInPair, new_rank);
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
bool ProjectAndLift<IntegerPL, IntegerRet>::fiber_interval(IntegerRet& MinInterval,
                                                           IntegerRet& MaxInterval,
                                                           const vector<IntegerRet>& base_point) {

    size_t dim = base_point.size() + 1;
    Matrix<IntegerPL>& Supps = AllSupps[dim];
    vector<size_t>& Order = AllOrders[dim];
    assert(Order.size() == Supps.nr_of_rows());

    bool FirstMin = true, FirstMax = true;
    vector<IntegerPL> LiftedGen;
    convert(LiftedGen, base_point);
    size_t check_supps = Supps.nr_of_rows();
    if (check_supps > 1000 && dim < EmbDim && !no_relax)
        check_supps = 1000;
    for (size_t j = 0; j < check_supps; ++j) {

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        IntegerPL Den = Supps[Order[j]].back();
        if (Den == 0)
            continue;
        IntegerPL Num = -v_scalar_product_vectors_unequal_lungth(LiftedGen, Supps[Order[j]]);
        IntegerRet Bound = 0;
        if (Den > 0) {  // we must produce a lower bound of the interval
            Bound = ceil_quot<IntegerRet, IntegerPL>(Num, Den);
            if (FirstMin || Bound > MinInterval) {
                MinInterval = Bound;
                FirstMin = false;
            }
        }
        if (Den < 0) {  // we must produce an upper bound of the interval
            Bound = floor_quot<IntegerRet, IntegerPL>(Num, Den);
            if (FirstMax || Bound < MaxInterval) {
                MaxInterval = Bound;
                FirstMax = false;
            }
        }
        if (!FirstMax && !FirstMin && MaxInterval < MinInterval)
            return false;  // interval empty
    }
    return true;  // interval nonempty
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::finalize_latt_point(vector<IntegerRet> NewPoint, const int tn) {

    if(only_single_point && single_point_found)
        return;

    vector<IntegerPL> NewPointPL;
    if(sparse){ // we must make sure that all inequalities are applied to our lattice point
        convert(NewPointPL, NewPoint);
        for(size_t i = 0; i < AllSupps[EmbDim].nr_of_rows(); ++i){
            if(used_supps[i])
                continue;
            if(v_scalar_product(NewPointPL, AllSupps[EmbDim][i]) < 0){
                return;
            }
        }
        // and we check all polynomial equations because we have suppressed the "noneffective" ones
        if(!PolyEquations.check(NewPoint, true, false)) // true = equations, fasle = any length
            return;
        if(!PolyInequalities.check(NewPoint, false, false)) // false = inequlities, fasle = any length
            return;
    }

    if(fusion.use_automorphisms)
        NewPoint = fusion.normal_form_of(NewPoint);

    if(only_single_point || !first_solution_printed){
#pragma omp critical(FINALSOL)
        {
        if(!first_solution_printed){
            if(verbose)
                verboseOutput() << endl << "Final solution 1 (preliminary format)-----  "  << NewPoint;
            verb_length = 0;
            }
            SingleDeg1Point = NewPoint;
        }
        first_solution_printed = true;
        if(only_single_point){
            TotalNrLP = 1;
            single_point_found = true;
        }
    }

    if(only_single_point && single_point_found)
        return;

#pragma omp atomic
    TotalNrLP++;

    if (!count_only)
        Deg1Thread[tn].push_back(NewPoint);

    if (Grading.size() > 0) {
        long deg = convertToLong(v_scalar_product(Grading, NewPoint));
        if (deg >= 0) {
            if (deg >= (long)h_vec_pos_thread[tn].size())
                h_vec_pos_thread[tn].resize(deg + 1);
            h_vec_pos_thread[tn][deg]++;
        }
        else {
            deg *= -1;
            if (deg >= (long)h_vec_neg_thread[tn].size())
                h_vec_neg_thread[tn].resize(deg + 1);
            h_vec_neg_thread[tn][deg]++;
        }
    }
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::splice_into_short_deg1_points(list<vector<IntegerRet> >& Deg1PointsComputed) {

    if(Deg1PointsComputed.empty())
        return;

    size_t dim = Deg1PointsComputed.front().size();

    vector<short> short_sol(dim);
    long long bridge;

    while(!Deg1PointsComputed.empty()){
        for(size_t i = 0; i < dim; ++i){
            bridge = convertTo<long long>(Deg1PointsComputed.front()[i]);
            if(bridge > 32767 || bridge < -32768)
                throw NoComputationException("Range short int not sufficient");
            short_sol[i] = bridge;
        }
        ShortDeg1Points.push_back(short_sol);
        Deg1PointsComputed.pop_front();
    }
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::collect_results(list<vector<IntegerRet> >& Deg1PointsComputed) {

    if(use_short_int){
        splice_into_short_deg1_points(Deg1PointsComputed);
    }
    else{
        Deg1Points.splice(Deg1Points.end(), Deg1PointsComputed);
    }

    for (size_t i = 0; i < Deg1Thread.size(); ++i) {
        if (h_vec_pos_thread[i].size() > h_vec_pos.size())
            h_vec_pos.resize(h_vec_pos_thread[i].size());
        for (size_t j = 0; j < h_vec_pos_thread[i].size(); ++j)
            h_vec_pos[j] += h_vec_pos_thread[i][j];
        h_vec_pos_thread[i].clear();
    }

    for (size_t i = 0; i < Deg1Thread.size(); ++i) {
        if (h_vec_neg_thread[i].size() > h_vec_neg.size())
            h_vec_neg.resize(h_vec_neg_thread[i].size());
        for (size_t j = 0; j < h_vec_neg_thread[i].size(); ++j)
            h_vec_neg[j] += h_vec_neg_thread[i][j];
        h_vec_neg_thread[i].clear();
    }
}

size_t our_counter = 0;

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::lift_points_to_this_dim(list<vector<IntegerRet> >& Deg1Proj) {
    if (Deg1Proj.empty())
        return;

    size_t dim1 = Deg1Proj.front().size(); // length iof vectors so far
    size_t dim = dim1 + 1; // old length // extended length

    // AllCongs[dim].debug_print();

    if(dim > EmbDim){ // can happen in patching algorithm
        used_supps.reset(); // to make finalize_latt_point work
        sparse = true; // ditto
        for(auto& P: Deg1Proj)
            finalize_latt_point(P, 0);
        if(use_short_int){
            splice_into_short_deg1_points(Deg1Thread[0]);
        }
        else{
            Deg1Points.splice(Deg1Points.begin(), Deg1Thread[0]);
        }
        return;
    }


    list<vector<IntegerRet> > Deg1Lifted;  // to this dimension if < EmbDim

    max_nr_new_latt_points_total = libnormaliz::max_nr_new_latt_points_total;

    size_t max_nr_per_thread = max_nr_new_latt_points_total / omp_get_max_threads();

    size_t nr_to_lift = Deg1Proj.size();
    NrLP[dim1] += nr_to_lift;
    size_t already_lifted = 0;

    bool not_done = true;
    bool has_poly_equs = (PolyEquations.size() > 0);
    bool has_poly_inequs = (PolyInequalities.size() > 0);

    while (not_done) {


        if(GlobalTimeBound > 0 &&  TimeSinceStart() > GlobalTimeBound){
            throw TimeBoundException("project-and-lift");
        }


        not_done = false;
        bool message_printed = false;

        bool skip_remaining;
        std::exception_ptr tmp_exception;

        skip_remaining = false;
        int omp_start_level = omp_get_level();

        // size_t nr_cong_killed = 0;


#pragma omp parallel
        {
            int tn;
            if (omp_get_level() == omp_start_level)
                tn = 0;
            else
                tn = omp_get_ancestor_thread_num(omp_start_level + 1);

            size_t nr_points_in_thread = 0;

            size_t ppos = 0;
            auto p = Deg1Proj.begin();
#pragma omp for schedule(dynamic)
            for (size_t i = 0; i < nr_to_lift; ++i) {
                if (skip_remaining)
                    continue;

                for (; i > ppos; ++ppos, ++p)
                    ;
                for (; i < ppos; --ppos, --p)
                    ;

                if ((*p)[0] == 0)  // point done
                    continue;

                if (!not_done && verbose) {
#pragma omp critical
                    {
                        if (!message_printed)
                            verboseOutput() << "Lifting to dimension " << dim << endl;
                        message_printed = true;
                    }
                }

                not_done = true;

#pragma omp atomic
                already_lifted ++;

                try {
                    IntegerRet MinInterval = 0, MaxInterval = 0;  // the fiber over *p is an interval -- 0 to make gcc happy
                    fiber_interval(MinInterval, MaxInterval, *p);;
                    IntegerRet add_nr_Int = 0;
                    if (MaxInterval >= MinInterval)
                        add_nr_Int = 1 + MaxInterval - MinInterval;
                    long long add_nr = convertToLongLong(add_nr_Int);
                    if (dim == EmbDim && count_only && add_nr >= 1 && !primitive
                           && Congs.nr_of_rows() == 0 && Grading.size() == 0 && PolyEquations.size() == 0
                       && PolyInequalities.size() == 0){
#pragma omp atomic
                        TotalNrLP += add_nr;
                    }
                    else {  // lift ppoint
                        for (IntegerRet k = MinInterval; k <= MaxInterval; ++k) {
                            INTERRUPT_COMPUTATION_BY_EXCEPTION

                            vector<IntegerRet> NewPoint(dim);
                            for (size_t j = 0; j < dim1; ++j)
                                NewPoint[j] = (*p)[j];
                            NewPoint[dim1] = k;

                            if(has_poly_equs && !PolyEquations.check(NewPoint, true, true)) // true = equations, true  = exact length
                                continue;
                            if(has_poly_inequs && !PolyInequalities.check(NewPoint, false, true)){ // false = inequalities, true  = exact length
                                continue;
                            }

                            if (!AllCongs[dim].check_congruences(NewPoint)){
                                continue;
                            }

                            if (dim == EmbDim) {
                                finalize_latt_point(NewPoint, tn);
                            }
                            else{
                                Deg1Thread[tn].push_back(NewPoint);
                            }
                        }
                    }

                    (*p)[0] = 0;  // mark point as done
                    if (dim < EmbDim)
                        nr_points_in_thread += add_nr;
                    if (nr_points_in_thread > max_nr_per_thread) {  // thread is full
                        skip_remaining = true;
#pragma omp flush(skip_remaining)
                    }

                } catch (const std::exception&) {
                    tmp_exception = std::current_exception();
                    skip_remaining = true;
#pragma omp flush(skip_remaining)
                }

            }  // lifting

        }  // pararllel

        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);

        for (size_t i = 0; i < Deg1Thread.size(); ++i)
            Deg1Lifted.splice(Deg1Lifted.begin(), Deg1Thread[i]);

        if (dim == EmbDim){
            collect_results(Deg1Lifted);
        }

        if(already_lifted == nr_to_lift){
            if(dim1 <= 1){
                if(DoneWithDim.size() > 1)
                    DoneWithDim[1] = true;
                DoneWithDim[0] = true;
            }

            if(dim1 >=1 && DoneWithDim[dim1-1]){
                if(verbose && !DoneWithDim[dim1])
                    verboseOutput() << "Done with dim " << dim1 << " LatticePoints " << NrLP[dim1] << endl;
                DoneWithDim[dim1] = true;
            }
        }
        lift_points_to_this_dim(Deg1Lifted);
        Deg1Lifted.clear();

    }  // not_done

    if(verbose && dim == EmbDim){
        verboseOutput() << "Complete lattice points so far " << TotalNrLP << endl;
    }

    return;
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::lift_point_recursively(vector<IntegerRet>& final_latt_point,
                                                                   const vector<IntegerRet>& latt_point_proj) {
    size_t dim1 = latt_point_proj.size();
    size_t dim = dim1 + 1;
    size_t final_dim = AllSupps.size() - 1;

    IntegerRet MinInterval = 0, MaxInterval = 0;  // the fiber over Deg1Proj[i] is an interval -- 0 to make gcc happy
    fiber_interval(MinInterval, MaxInterval, latt_point_proj);
    for (IntegerRet k = MinInterval; k <= MaxInterval; ++k) {
        INTERRUPT_COMPUTATION_BY_EXCEPTION

        vector<IntegerRet> NewPoint(dim);
        for (size_t j = 0; j < dim1; ++j)
            NewPoint[j] = latt_point_proj[j];
        NewPoint[dim1] = k;

        if (!AllCongs[dim].check_congruences(NewPoint))
            continue;

        if (dim == final_dim && NewPoint != excluded_point) {
            final_latt_point = NewPoint;
            break;
        }
        if (dim < final_dim) {
            lift_point_recursively(final_latt_point, NewPoint);
            if (final_latt_point.size() > 0)
                break;
        }
    }
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::find_single_point() {
    size_t dim = AllSupps.size() - 1;
    assert(dim >= 2);

    vector<IntegerRet> start(1, GD);
    vector<IntegerRet> final_latt_point;
    lift_point_recursively(final_latt_point, start);
    if (final_latt_point.size() > 0) {
        SingleDeg1Point = final_latt_point;
        if (verbose)
            verboseOutput() << "Found point" << endl;
    }
    else {
        if (verbose)
            verboseOutput() << "No point found" << endl;
    }
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_latt_points() {

    size_t dim = AllSupps.size() - 1;
    assert(dim >= 2);

    if(start_list.empty()){
        vector<IntegerRet> start(1, GD);
        start_list.push_back(start);
    }
    lift_points_to_this_dim(start_list);
    NrLP[EmbDim] = TotalNrLP;
    if(verbose){
        verboseOutput() << endl << "=======================================" << endl;
        verboseOutput() << "Final number of lattice points "  << NrLP[EmbDim] << endl;
    }
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_latt_points_float() {

    ProjectAndLift<nmz_float, IntegerRet> FloatLift(*this);
    FloatLift.compute_latt_points();
    Deg1Points.swap(FloatLift.Deg1Points);
    TotalNrLP = FloatLift.TotalNrLP;
    h_vec_pos = FloatLift.h_vec_pos;
    h_vec_neg = FloatLift.h_vec_neg;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::initialize(const Matrix<IntegerPL>& Supps, size_t rank) {
    EmbDim = Supps.nr_of_columns();  // our embedding dimension
    AllSupps.resize(EmbDim + 1);
    AllCongs.resize(EmbDim + 1);
    AllOrders.resize(EmbDim + 1);
    AllNrEqus.resize(EmbDim + 1);
    AllSupps[EmbDim] = Supps;
    Congs.resize(0, EmbDim+1);
    ImportedLocalSolutions.resize(EmbDim + 1);
    AllSupps[EmbDim].remove_duplicate_and_zero_rows();
    AllOrders[EmbDim] = order_supps(AllSupps[EmbDim]);
    DoneWithDim.resize(EmbDim+1);
    used_supps.resize(AllSupps[EmbDim].nr_of_rows());
    StartRank = rank;
    GD = 1;  // the default choice
    verbose = true;
    is_parallelotope = false;
    no_crunch = true;
    use_LLL = false;
    no_relax = false;
    primitive = false;
    sparse = false;
    patching_allowed = true;
    count_only = false;
    system_unsolvable = false;
    use_coord_weights = false;
    no_weights = false;
    fusion_rings_computation = false;
    single_point_found = false;
    first_solution_printed = false;
    only_single_point = false; // if in case don't come via compute
    linear_order_patches = false;
    cong_order_patches = false;
    distributed_computation = false;
    check_simplicity_all = false;
    check_simplicity_cand = false;
    stored_local_solutions = false;
    use_short_int = false;
    no_heuristic_minimization = false;
    TotalNrLP = 0;
    min_return_patch = 0;
    NrLP.resize(EmbDim + 1);
    nr_time_printed = 0;

    Congs = Matrix<IntegerRet>(0, EmbDim + 1);

    Deg1Thread.resize(omp_get_max_threads());
    h_vec_pos_thread.resize(omp_get_max_threads());
    h_vec_neg_thread.resize(omp_get_max_threads());

    LLL_Coordinates = Sublattice_Representation<IntegerRet>(EmbDim);  // identity
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
ProjectAndLift<IntegerPL, IntegerRet>::ProjectAndLift() {
}
//---------------------------------------------------------------------------
// General constructor
template <typename IntegerPL, typename IntegerRet>
ProjectAndLift<IntegerPL, IntegerRet>::ProjectAndLift(const Matrix<IntegerPL>& Supps,
                                                      const vector<dynamic_bitset>& Ind,
                                                      size_t rank) {
    initialize(Supps, rank);
    StartInd = Ind;
}

//---------------------------------------------------------------------------
// Constructor for parallelotopes
template <typename IntegerPL, typename IntegerRet>
ProjectAndLift<IntegerPL, IntegerRet>::ProjectAndLift(const Matrix<IntegerPL>& Supps,
                                                      const vector<dynamic_bitset>& Pair,
                                                      const vector<dynamic_bitset>& ParaInPair,
                                                      size_t rank) {
    initialize(Supps, rank);
    is_parallelotope = true;
    StartPair = Pair;
    StartParaInPair = ParaInPair;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_congruences(const Matrix<IntegerRet>& congruences) {
    Congs = congruences;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_PolyEquations(const OurPolynomialSystem<IntegerRet>& PolyEqus, const bool minimize) {
    PolyEquations = PolyEqus;
    OurPolynomialSystem<IntegerRet> DerivedPolyInequs = PolyEquations;
    RestrictablePolyInequs.insert(RestrictablePolyInequs.begin(), DerivedPolyInequs.begin(), DerivedPolyInequs.end());
    IntegerRet MinusOne = -1;
    DerivedPolyInequs.multiply_by_constant(MinusOne);
    RestrictablePolyInequs.insert(RestrictablePolyInequs.begin(), DerivedPolyInequs.begin(), DerivedPolyInequs.end());
    Matrix<IntegerPL> LinEqusPL = reconstruct_equations(AllSupps[EmbDim]);
    Matrix<IntegerRet> LinEqus;
    convert(LinEqus, LinEqusPL);
    if(minimize){
        if(verbose){
            verboseOutput() << "Minimizing polynomial equations (may take long time)"<< endl;
            verboseOutput() << "System has " << PolyEquations.size() << " equations" << endl;
        }
#ifdef NMZ_COCOA
        PolyEquations = PolyEquations.minimize_equations(LinEqus);
#else
        assert(false);
#endif
        if(verbose){
            verboseOutput() << "Minimal system has " << PolyEquations.size() << " equations" << endl;
        }
    }
}
//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_PolyInequalities(const OurPolynomialSystem<IntegerRet>& PolyInequs) {
    PolyInequalities = PolyInequs;
    RestrictablePolyInequs.insert(RestrictablePolyInequs.begin(), PolyInequs.begin(), PolyInequs.end());
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_startList(const list<vector<IntegerRet> >& start_from) {
    start_list = start_from;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_verbose(bool on_off) {
    verbose = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_short_int(bool on_off) {
    use_short_int = on_off;
}
//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_no_heuristic_minimization(bool on_off) {
    no_heuristic_minimization = on_off;
}
//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_LLL(bool on_off) {
    use_LLL = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_no_relax(bool on_off) {
    no_relax = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_coord_weights(bool on_off) {
    use_coord_weights = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_no_weights(bool on_off) {
    no_weights = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_linear_order_patches(bool on_off) {
    linear_order_patches = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_cong_order_patches(bool on_off) {
    cong_order_patches = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_primitive() {
    primitive = true;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_distributed_computation(const bool on_off) {
    distributed_computation = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_patching_allowed(bool on_off) {
    patching_allowed = on_off;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_grading_denom(const IntegerRet GradingDenom) {
    GD = GradingDenom;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_grading(const vector<IntegerRet>& grad) {
    Grading = grad;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_excluded_point(const vector<IntegerRet>& excl_point) {
    excluded_point = excl_point;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_vertices(const Matrix<IntegerPL>& Verts) {
    Vertices = Verts;
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::read_split_data() {

    distributed_computation = false;
    our_split.read_data(global_project);
    our_split.set_this_split(split_index_option);
    split_refinement = our_split.this_refinement; // needed in cone for output of lat file
    if(verbose){
        verboseOutput() << "split levels " << our_split.this_split_levels;
        verboseOutput() << "split moduli " << our_split.split_moduli;
        verboseOutput() << "split residues " << our_split.this_split_residues;
        verboseOutput() << "done indices " << our_split.this_split_done_indices;
        verboseOutput() << "refinement " << our_split.this_refinement << endl;
        if(split_refinement >0)
            verboseOutput() << "split residues " << our_split.this_split_min_returns;
    }
    // starting the new lat file
    lat_file_name = global_project + "." + to_string(split_refinement) + "." + to_string(split_index_rounds) + ".lat";
    if(verbose)
        verboseOutput() << "Writing " << lat_file_name << endl;
    ofstream prel_data(lat_file_name);
    prel_data << "preliminary_stage" << endl;
    prel_data.close();
}
//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute(bool all_points, bool lifting_float, bool do_only_count) {
    // Project-and-lift for lattice points in a polytope.
    // The first coordinate is homogenizing. Its value for polytope points ism set by GD so that
    // a grading denominator 1=1 can be accommodated.
    // We need only the support hyperplanes Supps and the facet-vertex incidence matrix Ind.
    // Its rows correspond to facets.

#ifdef NMZ_EXTENDED_TESTS
    if (!using_GMP<IntegerRet>() && !using_renf<IntegerRet>() && test_arith_overflow_proj_and_lift)
        throw ArithmeticException(0);
#endif

    if(is_split_patching){
        read_split_data();
    }

    if(is_split_patching || distributed_computation){
        stored_local_solutions = true; // will of course be tested
    }

    if(fusion.nr_coordinates > 0 && fusion.nr_coordinates != EmbDim -1){
        throw BadInputException("Wrong number of coordinates in fusion data. Mismatch of duality or commutativity.");
    }

    assert(all_points || !lifting_float);  // only all points allowed with float

    assert(all_points || !do_only_count);  // counting maks only sense for all points

    //used in patching
    only_single_point = !all_points;

    if (use_LLL) {
        LLL_coordinates_without_1st_col(LLL_Coordinates, AllSupps[EmbDim], Vertices, verbose);
        // Note: LLL_Coordinates is of type IntegerRet.
        Matrix<IntegerPL>
            Aconv;  // we cannot use to_sublattice_dual directly (not even with convert) since the integer types may not match
        convert(Aconv, LLL_Coordinates.getEmbeddingMatrix());
        AllSupps[EmbDim] = AllSupps[EmbDim].multiplication(Aconv.transpose());

        if (Congs.nr_of_rows() > 0) {  // must also transform congruences
            vector<IntegerRet> Moduli(Congs.nr_of_rows());
            for (size_t i = 0; i < Congs.nr_of_rows(); ++i)
                Moduli[i] = Congs[i][Congs.nr_of_columns() - 1];
            Matrix<IntegerRet> WithoutModuli(0, Congs.nr_of_columns() - 1);
            for (size_t i = 0; i < Congs.nr_of_rows(); ++i) {
                vector<IntegerRet> trans = Congs[i];
                trans.resize(trans.size() - 1);
                WithoutModuli.append(trans);
            }
            Congs = LLL_Coordinates.to_sublattice_dual(WithoutModuli);
            Congs.insert_column(Congs.nr_of_columns(), Moduli);
        }
        if (Grading.size() > 0)
            Grading = LLL_Coordinates.to_sublattice_dual_no_div(Grading);
    }

    add_congruences_from_equations();
    restrict_congruences();

    count_only = do_only_count;  // count_only belongs to *this

    if(primitive && patching_allowed){
        if(verbose)
            verboseOutput() << "Checking if patching possible" << endl;
        check_and_prepare_sparse();
    }

    if(!sparse){
        if (verbose)
            verboseOutput() << "Projection";
        if(primitive){
            if(verbose)
                verboseOutput() << " with relaxation for positive system " << endl;
            compute_projections_primitive(EmbDim);
        }
        else{
            if(verbose)
                verboseOutput() << "for general system" << endl;
            compute_projections(EmbDim, 1, StartInd, StartPair, StartParaInPair, StartRank);
        }
        //restrict_congruences();
    }

    if(system_unsolvable)
        return;

    if(all_points){
        if(sparse){ // all & sparse
            if(verbose)
                verboseOutput() << "Patching for all points" << endl;
            compute_latt_points_by_patching();
        }
        else{ // al and mot sparse
            if (verbose)
                verboseOutput() << "Lifting" << endl;
            if (!lifting_float || (lifting_float && using_float<IntegerPL>())) {
                compute_latt_points();
            }
            else {
                compute_latt_points_float();  // with intermediate conversion to float
            }
        }
    }
    else{ // single
        if(sparse){ // single and sparse
            if(verbose)
                verboseOutput() << "Patching for a single point" << endl;
            compute_latt_points_by_patching();
        }
        else{ // single and not sparse
            if (verbose)
                verboseOutput() << "Try finding a lattice point" << endl;
            find_single_point();
        }
    }
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_only_projection(size_t down_to) {
    assert(down_to >= 1);
    compute_projections(EmbDim, down_to, StartInd, StartPair, StartParaInPair, StartRank, true);
}

//---------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::put_short_deg1Points_into(vector<vector< short> >& LattPoints) {

    assert(!use_LLL);

    while (!ShortDeg1Points.empty()) {
        LattPoints.push_back(ShortDeg1Points.front());
        ShortDeg1Points.pop_front();
    }
}

//---------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::put_deg1Points_into(vector<vector< IntegerRet> >& LattPoints) {

    while (!Deg1Points.empty()) {
        if (use_LLL) {
            LattPoints.push_back(LLL_Coordinates.from_sublattice(Deg1Points.front()));
        }
        else
            LattPoints.push_back(Deg1Points.front());
        Deg1Points.pop_front();
    }
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::put_eg1Points_into(Matrix<IntegerRet>& LattPoints) {

    if(Deg1Points.empty() && !ShortDeg1Points.empty()){
        vector<IntegerRet> bridge(ShortDeg1Points.front().size());
        for(auto& p: ShortDeg1Points){
            for(size_t i= 0; i < bridge.size(); ++i)
                bridge[i] = p[i];
            Deg1Points.push_back(bridge);
        }
    }

    while (!Deg1Points.empty()) {
        if (use_LLL) {
            LattPoints.append(LLL_Coordinates.from_sublattice(Deg1Points.front()));
        }
        else
            LattPoints.append(Deg1Points.front());
        Deg1Points.pop_front();
    }
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::put_single_point_into(vector<IntegerRet>& LattPoint) {
    if (use_LLL && SingleDeg1Point.size() > 0)
        LattPoint = LLL_Coordinates.from_sublattice(SingleDeg1Point);
    else
        LattPoint = SingleDeg1Point;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
size_t ProjectAndLift<IntegerPL, IntegerRet>::getNumberLatticePoints() const {
    return TotalNrLP;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::get_h_vectors(vector<num_t>& pos, vector<num_t>& neg) const {
    pos = h_vec_pos;
    neg = h_vec_neg;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::setFusion(const  FusionBasic& FC){
    fusion = FusionComp<IntegerRet>(FC);
    /*
    if(fusion.fusion_type.size() == 0){
        FusionBasic basic;
        basic.data_from_string(global_project,false); // We want the final data
        fusion = FusionComp<IntegerRet>(basic);
    }
    if(fusion.fusion_type.size() == 0)
        throw BadInputException("Fusion rings asked for, but fusion data not available");
    */
}


//---------------------------------------------------------------------------
// For projection of cones
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::putSuppsAndEqus(Matrix<IntegerPL>& SuppsRet,
                                                            Matrix<IntegerPL>& EqusRet,
                                                            size_t in_dim) {
    assert(in_dim < EmbDim);
    assert(in_dim > 0);

    EqusRet.resize(0, in_dim);  // to make it well-defined
    size_t equs_start_in_row = AllSupps[in_dim].nr_of_rows() - 2 * AllNrEqus[in_dim];
    for (size_t i = equs_start_in_row; i < AllSupps[in_dim].nr_of_rows(); i += 2)  // equations come in +- pairs
        EqusRet.append(AllSupps[in_dim][i]);
    AllSupps[in_dim].swap(SuppsRet);
    SuppsRet.resize(equs_start_in_row);  // we must delete the superfluous rows because the transformation
                                         // to vector<vector> could else fail.
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::setOptions(const ConeProperties& ToCompute, const bool primitive, const bool our_verbose){

    if(is_split_patching)
        StartTime(stop_ckeck_begin);

    if(ToCompute.test(ConeProperty::FusionRings) || ToCompute.test(ConeProperty::SimpleFusionRings)){
        fusion_rings_computation = true;
        fusion.set_options(ToCompute, our_verbose);
    }

    if(ToCompute.test(ConeProperty::ShortInt))
        use_short_int = true;

    if(ToCompute.test(ConeProperty::NoHeuristicMinimization))
        no_heuristic_minimization = true;

    if(primitive){
        set_primitive();
        set_LLL(false);
        set_patching_allowed(!ToCompute.test(ConeProperty::NoPatching));
        set_cong_order_patches(ToCompute.test(ConeProperty::CongOrderPatches));
        set_linear_order_patches(ToCompute.test(ConeProperty::LinearOrderPatches));
        set_coord_weights(ToCompute.test(ConeProperty::UseWeightsPatching));
        set_no_weights(ToCompute.test(ConeProperty::NoWeights));
        if(!is_split_patching) // must exclude Split and DistributedComp simultaneously
            set_distributed_computation(ToCompute.test(ConeProperty::DistributedComp));
    }
    set_verbose(our_verbose);
    set_no_relax(ToCompute.test(ConeProperty::NoRelax));
    if(using_renf<IntegerPL>())
        set_LLL(false);
    else{
        if(!primitive)
            set_LLL(!ToCompute.test(ConeProperty::NoLLL));
    }
}

//---------------------------------------------------------------------------

template class ProjectAndLift<mpz_class, mpz_class>;
template class ProjectAndLift<long, long long>;
template class ProjectAndLift<mpz_class, long long>;
template class ProjectAndLift<long long, long long>;
template class ProjectAndLift<nmz_float, mpz_class>;
template class ProjectAndLift<nmz_float, long long>;
// template class ProjectAndLift<nmz_float, nmz_float>;
#ifndef NMZ_MIC_OFFLOAD  // offload with long is not supported
template class ProjectAndLift<long, long>;
template class ProjectAndLift<nmz_float, long>;
#endif

#ifdef ENFNORMALIZ
template class ProjectAndLift<renf_elem_class, mpz_class>;
#endif

//---------------------------------------------------------------------------

// interface cone/project-and-lift

#ifdef ENFNORMALIZ
// special version to avoid problems with machine integer etc.
template <>
void project_and_lift(Cone<renf_elem_class>&  C, const ConeProperties& ToCompute,
                                             Matrix<renf_elem_class>& Deg1,
                                             const Matrix<renf_elem_class>& Gens,
                                             const Matrix<renf_elem_class>& Supps,
                                             const Matrix<renf_elem_class>& Congs,
                                             const vector<renf_elem_class>& GradingOnPolytope,
                                             const bool primitive,
                                             const OurPolynomialSystem<renf_elem_class>& PolyEqus,
                                             const OurPolynomialSystem<renf_elem_class>& PolyInequs) {   // no primitive vgersion yet for renf
    bool count_only = ToCompute.test(ConeProperty::NumberLatticePoints);
    bool all_points = !( ToCompute.test(ConeProperty::SingleLatticePoint) ||ToCompute.test(ConeProperty::SingleFusionRing) );

    vector<dynamic_bitset> Ind;
    if(!primitive){
        Ind = vector<dynamic_bitset>(Supps.nr_of_rows(), dynamic_bitset(Gens.nr_of_rows()));
        for (size_t i = 0; i < Supps.nr_of_rows(); ++i)
            for (size_t j = 0; j < Gens.nr_of_rows(); ++j)
                if (v_scalar_product(Supps[i], Gens[j]) == 0)
                    Ind[i][j] = true;
    }

    size_t rank = C.getRankRaw();

    Matrix<renf_elem_class> Verts;
    if(!primitive){
        if (C.isComputed(ConeProperty::Generators)) {
            vector<key_t> choice = identity_key(Gens.nr_of_rows());  // Gens.max_rank_submatrix_lex();
            if (choice.size() >= C.getEmbeddingDim())
                Verts = Gens.submatrix(choice);
        }
    }

    // Matrix<mpz_class> Raw(0, Gens.nr_of_columns());

    vector<renf_elem_class> Dummy;
    ProjectAndLift<renf_elem_class, mpz_class> PL;
    PL = ProjectAndLift<renf_elem_class, mpz_class>(Supps, Ind, rank);

    PL.setFusion(C.getFusionBasicCone());
    PL.setOptions(ToCompute, primitive, C.getVerbose());

    PL.set_grading_denom(1);
    PL.set_vertices(Verts);
    OurPolynomialSystem<mpz_class> PolyEqus_mpz;
    convert(PolyEqus_mpz, PolyEqus);
    PL.set_PolyEquations(PolyEqus_mpz, ToCompute.test(ConeProperty::MinimizePolyEquations));
    OurPolynomialSystem<mpz_class> PolyInequs_mpz;
    convert(PolyInequs_mpz, PolyInequs);
    PL.set_PolyInequalities(PolyInequs_mpz);
    PL.compute(all_points, false, count_only);

    Matrix<mpz_class> Deg1_mpz(0,Supps.nr_of_columns());

    if(all_points){
            PL.put_eg1Points_into(Deg1_mpz);
            C.setNumberLatticePoints(PL.getNumberLatticePoints());
        }
    else{
        vector<mpz_class> SLP;
        PL.put_single_point_into(SLP);
        if(SLP.size() > 0){
            Deg1_mpz.append(SLP);
        }
    }

    convert(Deg1, Deg1_mpz);

}

#endif

//---------------------------------------------------------------------------

template <typename Integer>
void project_and_lift(Cone<Integer>&  C, const ConeProperties& ToCompute,
                                     Matrix<Integer>& Deg1,
                                     const Matrix<Integer>& Gens,
                                     const Matrix<Integer>& Supps,
                                     const Matrix<Integer>& Congs,
                                     const vector<Integer>& GradingOnPolytope,
                                     const bool primitive,
                                     const OurPolynomialSystem<Integer>& PolyEqus,
                                     const OurPolynomialSystem<Integer>& PolyInequs ) {
    bool float_projection = ToCompute.test(ConeProperty::ProjectionFloat);
    bool count_only = ToCompute.test(ConeProperty::NumberLatticePoints);
    bool all_points = !( ToCompute.test(ConeProperty::SingleLatticePoint) ||ToCompute.test(ConeProperty::SingleFusionRing) );
    vector<dynamic_bitset> Ind;

    if (!primitive && !C.isParallelotope()) {
        Ind = vector<dynamic_bitset>(Supps.nr_of_rows(), dynamic_bitset(Gens.nr_of_rows()));
        for (size_t i = 0; i < Supps.nr_of_rows(); ++i)
            for (size_t j = 0; j < Gens.nr_of_rows(); ++j)
                if (v_scalar_product(Supps[i], Gens[j]) == 0)
                    Ind[i][j] = true;
    }

    size_t rank = C.getRankRaw();

    Matrix<Integer> Verts;
        if(!primitive){
        if (C.isComputed(ConeProperty::Generators)) {
            vector<key_t> choice = identity_key(Gens.nr_of_rows());  // Gens.max_rank_submatrix_lex();
            if (choice.size() >= C.getEmbeddingDim())
                Verts = Gens.submatrix(choice);
        }
    }

    vector<num_t> h_vec_pos, h_vec_neg;

    if (float_projection) {  // conversion to float inside project-and-lift
        // vector<Integer> Dummy;
        ProjectAndLift<Integer, MachineInteger> PL;
        if (!C.isParallelotope())
            PL = ProjectAndLift<Integer, MachineInteger>(Supps, Ind, rank);
        else
            PL = ProjectAndLift<Integer, MachineInteger>(Supps, C.getPair(), C.getParaInPair(), rank);
        Matrix<MachineInteger> CongsMI;
        convert(CongsMI, Congs);
        PL.set_congruences(CongsMI);
        PL.set_grading_denom(convertTo<MachineInteger>(C.getGradingDenomRaw()));
        vector<MachineInteger> GOPMI;
        convert(GOPMI, GradingOnPolytope);
        PL.set_grading(GOPMI);
        PL.set_verbose(verbose);
        PL.set_LLL(!ToCompute.test(ConeProperty::NoLLL));
        PL.set_no_relax(ToCompute.test(ConeProperty::NoRelax));
        PL.set_vertices(Verts);

        PL.compute(true, true, count_only);  // the first true for all_points, the second for float
        Matrix<MachineInteger> Deg1MI(0, Deg1.nr_of_columns());
        PL.put_eg1Points_into(Deg1MI);
        convert(Deg1, Deg1MI);
        C.setNumberLatticePoints(PL.getNumberLatticePoints());
        PL.get_h_vectors(h_vec_pos, h_vec_neg);
    }
    else {
        if (C.getChangeIntegerType()) {
            Matrix<MachineInteger> Deg1MI(0, Deg1.nr_of_columns());
            // Matrix<MachineInteger> GensMI;
            Matrix<MachineInteger> SuppsMI;
            try {
                // convert(GensMI,Gens);
                convert(SuppsMI, Supps);
                MachineInteger GDMI = convertTo<MachineInteger>(C.getGradingDenomRaw());
                ProjectAndLift<MachineInteger, MachineInteger> PL;
                if (!C.isParallelotope() || primitive)
                    PL = ProjectAndLift<MachineInteger, MachineInteger>(SuppsMI, Ind, rank);
                else
                    PL = ProjectAndLift<MachineInteger, MachineInteger>(SuppsMI, C.getPair(), C.getParaInPair(), rank);
                Matrix<MachineInteger> CongsMI;
                convert(CongsMI, Congs);
                PL.set_congruences(CongsMI);

                PL.setFusion(C.getFusionBasicCone());
                PL.setOptions(ToCompute, primitive, C.getVerbose());

                PL.set_grading_denom(GDMI);
                vector<MachineInteger> GOPMI;
                convert(GOPMI, GradingOnPolytope);
                PL.set_grading(GOPMI);
                Matrix<MachineInteger> VertsMI;
                convert(VertsMI, Verts);
                PL.set_vertices(VertsMI);

                OurPolynomialSystem<MachineInteger> PolyEqus_MI;
                OurPolynomialSystem<MachineInteger> PolyInequs_MI;
                convert(PolyEqus_MI, PolyEqus);
                convert(PolyInequs_MI, PolyInequs);
                PL.set_PolyEquations(PolyEqus_MI, ToCompute.test(ConeProperty::MinimizePolyEquations));
                PL.set_PolyInequalities(PolyInequs_MI);
                if(PolyInequs.size() > 0 || PolyEqus.size() > 0)
                    PL.set_LLL(false);

                PL.compute(all_points, false, count_only);

                if(all_points){
                    PL.put_eg1Points_into(Deg1MI);
                    C.setNumberLatticePoints(PL.getNumberLatticePoints());
                }
                else{
                    vector<MachineInteger> SLP_MI;
                    PL.put_single_point_into(SLP_MI);
                    if(SLP_MI.size() > 0){
                        Deg1MI.append(SLP_MI);
                    }
                }

                PL.get_h_vectors(h_vec_pos, h_vec_neg);
            } catch (const ArithmeticException& e) {
                if (verbose) {
                    verboseOutput() << e.what() << endl;
                    verboseOutput() << "Restarting with a bigger type." << endl;
                }
                C.setChangeIntegerType(false);
            }
            if (C.getChangeIntegerType()) {
                convert(Deg1, Deg1MI);
            }
        }

        if (!C.getChangeIntegerType()) {
            ProjectAndLift<Integer, Integer> PL;
            if (!C.isParallelotope() || primitive)
                PL = ProjectAndLift<Integer, Integer>(Supps, Ind, rank);
            else
                PL = ProjectAndLift<Integer, Integer>(Supps, C.getPair(), C.getParaInPair(), rank);
            PL.set_congruences(Congs);

            PL.setFusion(C.getFusionBasicCone());
            PL.setOptions(ToCompute, primitive, C.getVerbose());

            PL.set_grading_denom(C.getGradingDenomRaw());
            PL.set_grading(GradingOnPolytope);
            PL.set_vertices(Verts);
            PL.set_PolyEquations(PolyEqus, ToCompute.test(ConeProperty::MinimizePolyEquations));
            PL.set_PolyInequalities(PolyInequs);
            if(PolyInequs.size() > 0 || PolyEqus.size() > 0)
                PL.set_LLL(false);

            PL.compute(all_points, false, count_only);

            if(all_points){
                    PL.put_eg1Points_into(Deg1);
                    C.setNumberLatticePoints(PL.getNumberLatticePoints());
                }
            else{
                vector<Integer> SLP;
                PL.put_single_point_into(SLP);
                if(SLP.size() > 0){
                    Deg1.append(SLP);
                }
            }

            PL.get_h_vectors(h_vec_pos, h_vec_neg);
        }
    }

    if (ToCompute.test(ConeProperty::HilbertSeries) && C.isComputed(ConeProperty::Grading)) {
        C.make_Hilbert_series_from_pos_and_neg(h_vec_pos, h_vec_neg);
    }

    if(C.getVerbose())
        verboseOutput() << "Project-and-lift complete" << endl
                        << "------------------------------------------------------------" << endl;
}

//---------------------------------------------------------------------------

template void project_and_lift<long long>(Cone<long long>&  C, const ConeProperties& ToCompute,
                                     Matrix<long long>& Deg1,
                                     const Matrix<long long>& Gens,
                                     const Matrix<long long>& Supps,
                                     const Matrix<long long>& Congs,
                                     const vector<long long>& GradingOnPolytope,
                                     const bool primitive,
                                     const OurPolynomialSystem<long long>& PolyEqus,
                                     const OurPolynomialSystem<long long>& PolyInequs );

template void project_and_lift<long>(Cone<long>&  C, const ConeProperties& ToCompute,
                                     Matrix<long>& Deg1,
                                     const Matrix<long>& Gens,
                                     const Matrix<long>& Supps,
                                     const Matrix<long>& Congs,
                                     const vector<long>& GradingOnPolytope,
                                     const bool primitive,
                                     const OurPolynomialSystem<long>& PolyEqus,
                                     const OurPolynomialSystem<long>& PolyInequs );

template void project_and_lift<mpz_class>(Cone<mpz_class>&  C, const ConeProperties& ToCompute,
                                     Matrix<mpz_class>& Deg1,
                                     const Matrix<mpz_class>& Gens,
                                     const Matrix<mpz_class>& Supps,
                                     const Matrix<mpz_class>& Congs,
                                     const vector<mpz_class>& GradingOnPolytope,
                                     const bool primitive,
                                     const OurPolynomialSystem<mpz_class>& PolyEqus,
                                     const OurPolynomialSystem<mpz_class>& PolyInequs );


#ifdef ENFNORMALIZ

template void project_and_lift<renf_elem_class>(Cone<renf_elem_class>&  C, const ConeProperties& ToCompute,
                                     Matrix<renf_elem_class>& Deg1,
                                     const Matrix<renf_elem_class>& Gens,
                                     const Matrix<renf_elem_class>& Supps,
                                     const Matrix<renf_elem_class>& Congs,
                                     const vector<renf_elem_class>& GradingOnPolytope,
                                     const bool primitive,
                                     const OurPolynomialSystem<renf_elem_class>& PolyEqus,
                                     const OurPolynomialSystem<renf_elem_class>& PolyInequs );

#endif



}  // end namespace libnormaliz
