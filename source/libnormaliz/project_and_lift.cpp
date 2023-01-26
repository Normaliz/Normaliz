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

#include "libnormaliz/project_and_lift.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/sublattice_representation.h"
#include "libnormaliz/cone.h"

namespace libnormaliz {
using std::vector;
using std::string;
using std::list;
using std::pair;

//---------------------------------------------------------------------------
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

    // reorder_coordinates();

    size_t nr_all_supps = AllSupps[EmbDim].nr_of_rows();

    Indicator.resize(nr_all_supps); // indicaor of nonzero coordinates in inequality
    upper_bounds.resize(nr_all_supps); // indicator of inequalities giving upper boounds
    max_sparse.resize(nr_all_supps); // indicator of inequalities used in covering by "sparse" inequalities

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
        if(upper_bounds[i] /* && Indicator[i].count() <= EmbDim/5*/){ // our criterion for sparseness
            union_upper_bounds |= Indicator[i];
            sparse_bounds[i] = 1;
        }
    }

    sparse = (union_upper_bounds.count() == EmbDim);

    if(!sparse){
        if(verbose)
            verboseOutput() << "System not sparse" << endl;
        return;
    }

    if(verbose)
        verboseOutput() << "Preparing data for patching algorithm " << endl;

    // now we want to find the sparse upper bounds with maximal support
    vector<dynamic_bitset> help(nr_all_supps);
    for(size_t i = 0; i < nr_all_supps; ++i){
        help[i].resize(EmbDim);
        if(sparse_bounds[i])
            help[i] = Indicator[i];
    }
    dynamic_bitset max_sparse = sparse_bounds;

    dynamic_bitset covered(EmbDim);  // registers covered coordinates
    covered[0] = 1; // the 0-th coordinate is covered by all local PL
    AllLocalPL.resize(EmbDim);
    AllIntersections_key.resize(EmbDim);
    AllNew_coords_key.resize(EmbDim);
    AllCovered.resize(EmbDim);
    AllPatches.resize(EmbDim);
    active_coords.resize(EmbDim);

    AllPolyEqus.resize(EmbDim);
    AllPolyEqusThread.resize(EmbDim);
    for(auto& T: AllPolyEqusThread)
        T.resize(omp_get_max_threads());

    AllPolyInequs.resize(EmbDim);
    AllPolyInequsThread.resize(EmbDim);
    for(auto& T: AllPolyInequsThread)
        T.resize(omp_get_max_threads());

    NrRemainingLP.resize(EmbDim,0);
    NrDoneLP.resize(EmbDim,0);
    AllLocalSolutions_by_intersecion.resize(EmbDim);
    AllLocalSolutions.resize(EmbDim);

    poly_equs_minimized.resize(EmbDim);

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
        covered |= AllPatches[coord];

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
            if(!used_supps[i] && upper_bounds[i])
                relevant_supps_now.push_back(i);
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

#ifdef NMZ_DEVELOP
        if(verbose){
            verboseOutput() << "level " << LevelPatches[coord] << endl;
            verboseOutput() << "new coords " << new_coords_key;
        }
#endif
        // for the "local" project-and-lift we need their suport hyperplanes
        vector<key_t> LocalKey = bitset_to_key(AllPatches[coord]);
        Matrix<IntegerPL> LocalSuppsRaw;
        LocalSuppsRaw = AllSupps[EmbDim].submatrix(relevant_supps_now);
        // convert(LocalSuppsRaw, AllSupps[EmbDim].submatrix(relevant_supps_now));
        Matrix<IntegerPL> Localsupps = LocalSuppsRaw.transpose().submatrix(LocalKey).transpose(); // select columns

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
        Matrix<IntegerPL> LocalSuppsReordered(Localsupps.nr_of_rows(), nr_coordinates);
        for(size_t i = 0; i < Localsupps.nr_of_rows(); ++i){
            for(size_t j = 0; j < nr_coordinates; ++j)
                LocalSuppsReordered[i][j]= LocalSuppsRaw[i][OrderedCoordinates[j]];
        }

        // Now we can set up the local project-and-lift
        vector<dynamic_bitset> DummyInd;
        ProjectAndLift<IntegerPL, IntegerRet> PL(LocalSuppsReordered, DummyInd, 0); // 0 is dummy
        PL.set_LLL(false);
        PL.set_primitive();
        PL.set_verbose(false);
        PL.compute_projections_primitive(LocalKey.size());
        AllLocalPL[coord] = PL;

        dynamic_bitset new_covered = covered | AllPatches[coord];

        // Collect relevant polynomial constraints
        vector<key_t> PolyEqusKey, PolyInequsKey;

        // first the equations
        for(size_t i = 0; i < PolyEquations.size(); ++i){
            if(!PolyEquations[i].support.is_subset_of(new_covered)) // not yet usable
                continue;
            if(PolyEquations[i].support.is_subset_of(covered)) // already used
                continue;
            AllPolyEqus[coord].push_back(PolyEquations[i]);
            PolyEqusKey.push_back(i);
        }
        for(auto& T: AllPolyEqusThread[coord]){ // vcopy for each thread
            T = AllPolyEqus[coord];
        }


        // next the inequalities for which all coordinates are covered
        for(size_t i = 0; i < PolyInequalities.size(); ++i){
            if(!(PolyInequalities[i]).support.is_subset_of(new_covered))
                continue;
            if(PolyInequalities[i].support.is_subset_of(covered))
                continue;
            AllPolyInequs[coord].push_back(PolyInequalities[i]);
            PolyInequsKey.push_back(i);
        }
        for(auto& T: AllPolyInequsThread[coord]){ // vcopy for each thread
            T = AllPolyInequs[coord];
        }

#ifdef NMZ_DEVELOP
        if(verbose)
            verboseOutput() << endl << "index coord " << coord << " nr covered coordinates " << new_covered.count() << " coordinates " << bitset_to_key(new_covered);
        if(verbose && AllPolyEqus.size() > 0)
            verboseOutput() << endl << "poly equations " << PolyEqusKey;
        if(verbose && PolyInequsKey.size() > 0)
            verboseOutput() << endl << coord << " poly inequalities " << PolyInequsKey;
        if(verbose)
            verboseOutput() << "---------------------------------------------------------------" << endl;
#endif

        covered = new_covered;

        AllCovered[coord] = covered;

    } // coord


}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_covers() {

    // Note the indices of thze patches are the coordinates
    // with whom they are associated

    LevelPatches.resize(EmbDim);

    if(PolyEquations.empty()){
        for(size_t i = 0; i < EmbDim; ++i){
            if(AllPatches[i].size() > 0)
                InsertionOrderPatches.push_back(i);
        }

        for(size_t k = 0; k < InsertionOrderPatches.size(); ++k)
            LevelPatches[InsertionOrderPatches[k]] = k;
        return;
    }

#ifdef NMZ_DEVELOP
    if(verbose){
        for(size_t i = 0; i < EmbDim; ++i){
            verboseOutput() << i << " ----- " << bitset_to_key(AllPatches[i]);
        }
    }
#endif

    // First we cover the supports of polynom,ial equations by patches
    // using as few patches as possioble
    // The results are stored in covering_patches
    vector<pair<size_t, vector<key_t> > > covering_equations;

    size_t dim = PolyEquations[0].support.size();

    for(size_t i = 0; i< PolyEquations.size(); ++i){
        dynamic_bitset poly_supp = PolyEquations[i].support;
        dynamic_bitset already_covered(dim);
        size_t nr_patches_needed = 0;
        vector<key_t> patches_used;
        while(true){
            bool first = true;
            size_t max_nr_new_covered = 0; // = 0 to make gcc happy
            size_t max_covering;
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
                covering_equations.push_back(make_pair(patches_used.size(), patches_used));
                break;
            }
        }
    }

    // Now we must find an order of the supports of the polynomial equations
    // in which we insert their patches.
    // The idea is to have as miuch overlap in the covered coordinates as
    // possible.

    sort(covering_equations.begin(), covering_equations.end());

#ifdef NMZ_DEVELOP
    if(verbose){
        for(auto& N: covering_equations){
            verboseOutput() << N.second;
        }
        verboseOutput() << endl;
    }
#endif

    vector<key_t> InsertionOrderEquations;
    dynamic_bitset used_covering_equations(covering_equations.size());
    dynamic_bitset covered_coords(dim);
    while(covered_coords.count() < dim && used_covering_equations.count() < covering_equations.size()){
        dynamic_bitset test_covered_coords(dim);
        dynamic_bitset min_covered_coords(dim);
        bool first = true;
        size_t min_at;
        size_t nr_min_covered_coords = 0; // =0 to make gcc happy
        for( size_t i = 0; i < covering_equations.size(); ++i){
            if(used_covering_equations[i])
                continue;
            test_covered_coords = covered_coords;
            for(auto& c: covering_equations[i].second) // finding the coordinates covered by this patch
                test_covered_coords |= AllPatches[c];
            if(test_covered_coords == covered_coords){
                used_covering_equations[i] = true; // no new coordinates, not used explicitly , but not needed
                continue;
            }
            if(first || test_covered_coords.count() < nr_min_covered_coords){
                first = false;
                nr_min_covered_coords = test_covered_coords.count();
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

    // There could be patches not involved in polynomial equations
    for(size_t j = 0; j < EmbDim; ++j){
        if(!inserted_patches[j] && AllPatches[j].size() > 0)
            InsertionOrderPatches.push_back(j);
    }

    if(verbose){
        verboseOutput() << "Insertion order linear patches " << endl;
        verboseOutput() << InsertionOrderPatches << endl;
    }

    for(size_t k = 0; k < InsertionOrderPatches.size(); ++k)
        LevelPatches[InsertionOrderPatches[k]] = k;

    /*
    set<key_t> HI;
    for(auto& P: PolyEquations)
        HI.insert(P.highest_indet);
    for(auto& H: HI)
        cout << H << " ";
    cout << endl;*/
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_latt_points_by_patching() {

    vector<IntegerRet> start(EmbDim);
    start[0] = GD;
    list<vector<IntegerRet> > start_list;
    start_list.push_back(start);
    extend_points_to_next_coord(start_list, 0);
    NrLP[EmbDim] = TotalNrLP;
    if(verbose)
        verboseOutput() << "Final number of lattice points "  << NrLP[EmbDim] << endl;

    for(auto& n: NrRemainingLP){
            assert(n == 0);
    }
}

const size_t max_nr_new_latt_points_total = 1000000;
const size_t nr_new_latt_points_for_elimination_equs = 10000;

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::extend_points_to_next_coord(list<vector<IntegerRet> >& LatticePoints, const key_t this_patch) {

    size_t max_nr_per_thread =  max_nr_new_latt_points_total/ omp_get_max_threads();

    size_t coord = InsertionOrderPatches[this_patch]; // the coord marking the next patch

    StartTime();

    /* size_t min_fall_back = 0;
    bool min_found = false;
    size_t max_fall_back = 0;
    for(size_t i = 0; i < coord; ++i){
        if(NrRemiaaingLP[i] > 0){
            max_fall_back = i;
            if(!min_found){
                min_found = true;
                min_fall_back = i;
            }
        }
    }*/

    /* if(verbose){
        verboseOutput() << "coord " << coord;
        if(min_found){
            verboseOutput() << " back min " << min_fall_back << " max " << max_fall_back;
            verboseOutput() << endl;
            // verboseOutput() << NrRemainingLP << endl;
        }
        else{
            verboseOutput() << endl;
        }
    }*/

    INTERRUPT_COMPUTATION_BY_EXCEPTION

    auto& intersection_key = AllIntersections_key[coord];
    auto& new_coords_key = AllNew_coords_key[coord];

    vector<OurPolynomialSystem<IntegerRet>>& PolyEqusThread = AllPolyEqusThread[coord];
    vector<OurPolynomialSystem<IntegerRet>>& PolyInequsThread = AllPolyInequsThread[coord];

    ProjectAndLift<IntegerPL, IntegerRet>& LocalPL = AllLocalPL[coord];
    map<vector<IntegerRet>, vector<key_t> >& LocalSolutions_by_intersecion = AllLocalSolutions_by_intersecion[coord];
    Matrix<IntegerRet>& LocalSolutions = AllLocalSolutions[coord];
    dynamic_bitset covered = AllCovered[coord];

    vector<key_t> CoveredKey = bitset_to_key(covered);
    // cout << "Covered " << CoveredKey;

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

    for(auto IC = start_list.begin(); IC != start_list.end(); ){
        if(LocalSolutions_by_intersecion.find(*IC) == LocalSolutions_by_intersecion.end()){
            LocalSolutions_by_intersecion[*IC] = {};
            IC++;
        }
        else{
            IC =  start_list.erase(IC);
        }
    }

    // Now we extend the "new" intersection coordinates" by the local system
    LocalPL. set_startList(start_list);
    LocalPL.lift_points_to_this_dim(start_list); // computes the extensions
    Matrix<IntegerRet> LocalSolutionsNow(0, intersection_key.size() + new_coords_key.size());
    LocalPL.put_eg1Points_into(LocalSolutionsNow);
    size_t nr_old_solutions = LocalSolutions.nr_of_rows();
    if(nr_old_solutions == 0)
        LocalSolutions = LocalSolutionsNow;
    else
        LocalSolutions.append(LocalSolutionsNow);

    // Next the newly computed extensions are registered
    vector<IntegerRet> overlap(intersection_key.size());
    for(size_t i = nr_old_solutions; i < LocalSolutions.nr_of_rows(); i++){
        for(size_t j = 0; j < intersection_key.size(); ++j)
            overlap[j] = LocalSolutions[i][j];
        LocalSolutions_by_intersecion[overlap].push_back(i);
    }

    bool last_coord = (coord == InsertionOrderPatches.back());

    size_t nr_to_match = LatticePoints.size();
    size_t nr_points_matched = 0;

    // In the following the oputer paralleliozed loop is over the lattice points
    // given to the routine, and the inner is over the extensions along the overlao.
    // One could think about exchanging the order of the loops to get better
    // parallelization. But it is unclear whether this can be achieved.

    while (true) {
        if(verbose)
            verboseOutput() <<  LevelPatches[coord] << " / " << coord << " left " << nr_to_match - nr_points_matched << endl;

        bool skip_remaining;
        std::exception_ptr tmp_exception;

        skip_remaining = false;
        int omp_start_level = omp_get_level();

        size_t nr_latt_points_total = 0;  //statistics for this run of the while loop

        vector<vector<size_t> > poly_stat;
        vector<size_t> poly_stat_total;
        if(!poly_equs_minimized[coord]){
            poly_stat.resize(omp_get_max_threads());  // counts the number of "successful". iu.e. != 0,
            for(auto& p:poly_stat)                    // evaluations of a mpolynomial erquationj
                p.resize(PolyEqusThread[0].size());
            poly_stat_total.resize(poly_stat[0].size());
        }

#pragma omp parallel
        {

        vector<IntegerRet> overlap(intersection_key.size());
        vector<IntegerRet> NewLattPoint(EmbDim);

        int tn;
        if (omp_get_level() == omp_start_level)
            tn = 0;
        else
            tn = omp_get_ancestor_thread_num(omp_start_level + 1);

        size_t nr_points_in_thread = 0;

        size_t ppos = 0;
        auto P = LatticePoints.begin();

        list<key_t> order_poly_equs; // suddessful constraints to the front !!
        for(key_t k = 0; k < PolyEqusThread[tn].size(); ++k)
            order_poly_equs.push_back(k);
        list<key_t> order_poly_inequs;
        for(key_t k = 0; k < PolyInequsThread[tn].size(); ++k)
            order_poly_inequs.push_back(k);

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

            /*if(ppp % 2 != 0){
                (*P)[0] = 0;
                continue;
            }*/

        try{  // now the given lattice points are extended along their overlaps
            overlap = v_select_coordinates(*P, intersection_key);

            // now the extensions of the overlap
            for(auto& i: LocalSolutions_by_intersecion[overlap]){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

                NewLattPoint = *P;
                for(size_t j = 0; j < new_coords_key.size(); ++j)
                    NewLattPoint[new_coords_key[j]] = LocalSolutions[i][j + intersection_key.size()];

                bool can_be_inserted = true;

#pragma omp atomic
                nr_latt_points_total++;
                if(can_be_inserted){
                    for(auto pp = order_poly_equs.begin(); pp!= order_poly_equs.end(); ++pp){
                        if(PolyEqusThread[tn][*pp].evaluate(NewLattPoint) != 0){
                            can_be_inserted = false;
                            if(!poly_equs_minimized[coord])
                                poly_stat[tn][*pp]++;
                            else
                                order_poly_equs.splice(order_poly_equs.begin(), order_poly_equs, pp);
                            break;
                        }
                    }
                }
                if(can_be_inserted){
                    for(auto pp = order_poly_inequs.begin(); pp!= order_poly_inequs.end(); ++pp){
                        if(PolyInequsThread[tn][*pp].evaluate(NewLattPoint)< 0){
                            can_be_inserted = false;
                            order_poly_inequs.splice(order_poly_inequs.begin(), order_poly_inequs, pp);
                            break;
                        }
                    }
                }
                if(can_be_inserted){
                    nr_points_in_thread++;
                    if(last_coord)
                        finalize_latt_point(NewLattPoint, tn);
                    else
                        Deg1Thread[tn].push_back(NewLattPoint);
                }
            } // for i (inner for loop)

            (*P)[0] = 0;  // mark point as done
            if (nr_points_in_thread > max_nr_per_thread && !last_coord) {  // thread is full
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
            for(size_t i = 0; i< poly_stat.size(); ++i){
                    for(size_t j = 0; j < poly_stat[0].size(); ++j)
                        poly_stat_total[j] += poly_stat[i][j];
            }
            // cout << "PPPP " << poly_stat_total;
        }

        for (size_t i = 0; i < Deg1Thread.size(); ++i)
            NewLatticePoints.splice(NewLatticePoints.end(), Deg1Thread[i]);

        // size_t nr_new_latt_points = NewLatticePoints.size();

        if(last_coord)
            collect_results(NewLatticePoints); // clears NewLatticePoints

        NrDoneLP[coord] += nr_latt_points_total;

        // Discard ineffective polynomial equations
        if(!poly_equs_minimized[coord] &&  nr_latt_points_total > nr_new_latt_points_for_elimination_equs){
            poly_equs_minimized[coord] = true;
            OurPolynomialSystem<IntegerRet> EffectivePolys;
            for(size_t i = 0; i < PolyEqusThread[0].size(); ++i){
                if(poly_stat_total[i] > 0)
                    EffectivePolys.push_back(PolyEqusThread[0][i]);
            }

            for(size_t thr = 0; thr < PolyEqusThread.size(); ++thr)
                PolyEqusThread[thr] = EffectivePolys;
        }

        /*if(verbose){
            // verboseOutput() << " --- ext " << nr_latt_points_total << endl;
            // verboseOutput() << " cst " << nr_new_latt_points << endl;
            // verboseOutput() << "rtd " << nr_caught_by_restricted << endl;
        }*/

        MeasureTime(false, "Elapsed ");

        /* for(auto& P: NewLatticePoints)
            cout << P;

        if(verbose)
            verboseOutput() << "----------" << endl;*/

        NrRemainingLP[coord] = nr_to_match - nr_points_matched;

        if(!last_coord && NewLatticePoints.size() > 0)
            extend_points_to_next_coord(NewLatticePoints, this_patch + 1);
        NewLatticePoints.clear();
        if(nr_points_matched == nr_to_match)
            break;

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
        dynamic_bitset test_coevered(dim);
        size_t nr_new_covered = 0;
        for(size_t i = 0; i < PolyEquations.size(); ++i){
            if(PolyEquations[i].support.is_subset_of(covered))
                continue;
            test_coevered = covered | PolyEquations[i].support;
            if(first || test_coevered.count() < nr_new_covered){
                new_covered = test_coevered;
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
    // cout << "NNNNNNNNNNN " << NewOrder;

    // NewOrder = identity_key(dim);

    // convert(TestPL, TestV);
    /* cout << "Equations Old " << endl;
    for(size_t i = 0; i < AllSupps[EmbDim].nr_of_rows(); ++i){
        cout << i << "    " << v_scalar_product(AllSupps[EmbDim][i], TestPL) << endl;
    }*/

    // cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

    AllSupps[EmbDim].permute_columns(NewOrder);

    // cout << PolyEquations[0] << endl;

    // convert(TestRet, TestV);



    PolyEquations.permute_variables(NewOrder);
    // cout << PolyEquations[0] << endl;;

    PolyInequalities.permute_variables(NewOrder);

    // TestRet = v_permute_coordinates(TestRet, NewOrder);
    //for(size_t i = 0; i < PolyEquations.size(); ++i)
    //     cout << "New " << PolyEquations[i].evaluate(TestRet) << endl;

    // convert(TestPL, TestRet);
     /* cout << "Equations New" << endl;
    for(size_t i = 0; i < AllSupps[EmbDim].nr_of_rows(); ++i){
        cout << i << "    " << v_scalar_product(AllSupps[EmbDim][i], TestPL) << endl;
    }*/

    // xit(0);

    /* PolyEquations.clear();
    PolyInequalities.clear();*/
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

        // cout << "SSSSS " << Supps[Order[j]];

        IntegerPL Den = Supps[Order[j]].back();
        if (Den == 0)
            continue;
        IntegerPL Num = -v_scalar_product_vectors_unequal_lungth(LiftedGen, Supps[Order[j]]);
        // cout << "Num " << Num << endl;
        IntegerRet Bound = 0;
        // frac=(Num % Den !=0);
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
void ProjectAndLift<IntegerPL, IntegerRet>::finalize_latt_point(const vector<IntegerRet>& NewPoint, const int tn) {

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
    }

    if (!Congs.check_congruences(NewPoint))
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

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::collect_results(list<vector<IntegerRet> >& Deg1PointsComputed) {

    Deg1Points.splice(Deg1Points.end(), Deg1PointsComputed);

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

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::lift_points_to_this_dim(list<vector<IntegerRet> >& Deg1Proj) {
    if (Deg1Proj.empty())
        return;
    size_t dim1 = Deg1Proj.front().size(); // length iof vectors so far
    size_t dim = dim1 + 1; // old length // extended length

    if(dim > EmbDim){ // can happen in patching algorithm
        used_supps.reset(); // to make finalize_latt_point work
        sparse = true; // ditto
        for(auto& P: Deg1Proj)
            finalize_latt_point(P, 0);
        Deg1Points.splice(Deg1Points.begin(), Deg1Thread[0]);
        return;
    }


    list<vector<IntegerRet> > Deg1Lifted;  // to this dimension if < EmbDim

    size_t max_nr_per_thread = 1000000 / omp_get_max_threads();

    size_t nr_to_lift = Deg1Proj.size();
    NrLP[dim1] += nr_to_lift;
    size_t already_lifted = 0;

    bool not_done = true;
    bool has_poly_equs = (PolyEquations.size() > 0);
    bool has_poly_inequs = (PolyInequalities.size() > 0);

    while (not_done) {
        not_done = false;
        bool message_printed = false;

        bool skip_remaining;
        std::exception_ptr tmp_exception;

        skip_remaining = false;
        int omp_start_level = omp_get_level();


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
                    fiber_interval(MinInterval, MaxInterval, *p);
                    // cout << "Min " << MinInterval << " Max " << MaxInterval << endl;
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
    if(verbose)
        verboseOutput() << "Final number of lattice points "  << NrLP[EmbDim] << endl;
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
    AllOrders.resize(EmbDim + 1);
    AllNrEqus.resize(EmbDim + 1);
    AllSupps[EmbDim] = Supps;
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
    TotalNrLP = 0;
    NrLP.resize(EmbDim + 1);

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
void ProjectAndLift<IntegerPL, IntegerRet>::set_PolyEquations(const OurPolynomialSystem<IntegerRet>& PolyEqus) {
    PolyEquations = PolyEqus;
}
//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::set_PolyInequalities(const OurPolynomialSystem<IntegerRet>& PolyInequs) {
    PolyInequalities = PolyInequs;
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
void ProjectAndLift<IntegerPL, IntegerRet>::set_primitive() {
    primitive = true;
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

    assert(all_points || !lifting_float);  // only all points allowed with float

    assert(all_points || !do_only_count);  // counting maks only sense for all points

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

    count_only = do_only_count;  // count_only belongs to *this

    if(primitive && patching_allowed && all_points){
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
    }

    if(system_unsolvable)
        return;

    if (all_points) {
        if(sparse){
            if(verbose)
                verboseOutput() << "Patching" << endl;
            compute_latt_points_by_patching();
        }
        else{
            if (verbose)
                verboseOutput() << "Lifting" << endl;
            if (!lifting_float || (lifting_float && using_float<IntegerPL>())) {
                compute_latt_points();
            }
            else {
                compute_latt_points_float();  // with intermediate conversion to float
            }
            /* if(verbose)
                verboseOutput() << "Number of lattice points " << TotalNrLP << endl;*/
        }
    }
    else {
        if (verbose)
            verboseOutput() << "Try finding a lattice point" << endl;
        find_single_point();
    }
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_only_projection(size_t down_to) {
    assert(down_to >= 1);
    compute_projections(EmbDim, down_to, StartInd, StartPair, StartParaInPair, StartRank, true);
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::put_eg1Points_into(Matrix<IntegerRet>& LattPoints) {

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
    if (use_LLL)
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

}  // end namespace libnormaliz
