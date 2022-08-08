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
    dynamic_bitset union_upper_bounds(EmbDim); // to check whether sparse upper bound inequ covere
    dynamic_bitset sparse_bounds(nr_all_supps);

    for(size_t i = 0; i< nr_all_supps; ++i){
        if(upper_bounds[i] && Indicator[i].count() <= EmbDim/5){ // our criterion for sparseness
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
        verboseOutput() << "Preparing data for matching algorithm " << endl;
    
    // now we want to find the sparse upper bounds with maximal support
    vector<dynamic_bitset> help(nr_all_supps);
    for(size_t i = 0; i < nr_all_supps; ++i){
        help[i].resize(EmbDim);
        if(sparse_bounds[i])
            help[i] = Indicator[i];
    }
    dynamic_bitset max_sparse(nr_all_supps);
    max_sparse.flip();
    maximal_subsets(help, max_sparse);

    // now the main work: find "local" solutions and patach them
    // we extend the set of covered coordinates by at least
    // the first noncovered
    
    dynamic_bitset covered(EmbDim);  // registers covered coordinates 
    AllLocalSolutions.resize(EmbDim);
    AllIntersections_key.resize(EmbDim);
    AllNew_coords_key.resize(EmbDim);
    active_coords.resize(EmbDim);
    AllExtraInequalities.resize(EmbDim);
    AllPolyEqusKey.resize(EmbDim);
    AllPolyInequsKey.resize(EmbDim);
    
    dynamic_bitset used_supps(nr_all_supps);

    // main loop for preparation of coord dependent data
    for(size_t coord = 1; coord < EmbDim; coord++){
        if(covered[coord] == 1)
            continue;
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        active_coords[coord] = 1;
        
        key_t next_supp = 0; //  the next inequality we will use // = 0 to makemgcc happy
 
        // try to find an inequality with coord in its support and maximal intersection
        size_t max_card_intersect = 0;  // with the set of covered coordinates // = 0 to make gcc happy
        bool found_extension = false;
        for(size_t i = 0; i < nr_all_supps; i++){
            if(used_supps[i] || !max_sparse[i])
                continue;
            if(Indicator[i][coord]){
                size_t card_intersect = (covered & Indicator[i]).count();
                if(!found_extension || card_intersect > max_card_intersect){
                    next_supp = i;
                    max_card_intersect = card_intersect;
                }
                found_extension = true;
            }
        }
        assert(found_extension);

        // now we want to compute the "local" solutions for inequalities whose support
        // is contained in Indicator[next_supp]
        vector<key_t> relevant_supps_now;  
        for(size_t i = 0; i < nr_all_supps; ++i){
            if(Indicator[i].is_subset_of(Indicator[next_supp])){
                relevant_supps_now.push_back(i);
                used_supps[i] = 1;  // also used_supps[next_supp] set here
            }
        }

        vector<key_t> LocalKey = bitset_to_key(Indicator[next_supp]);
        Matrix<IntegerRet> LocalSuppsRaw;
        convert(LocalSuppsRaw, AllSupps[EmbDim].submatrix(relevant_supps_now));
        Matrix<IntegerRet> Localsupps = LocalSuppsRaw.transpose().submatrix(LocalKey).transpose();
        vector<IntegerRet> LocalDehom(Localsupps.nr_of_columns());
        LocalDehom[0] = 1;
        Matrix<IntegerRet> DehomMat(LocalDehom);
        // Localsupps.pretty_print(cout);
        Cone<IntegerRet> LocalCone(Type::inequalities,Localsupps,Type::dehomogenization, DehomMat);
        LocalCone.setVerbose(false);
        Matrix<IntegerRet> LocalSolutions = LocalCone.getLatticePointsMatrix();
        // LocalSolutions.pretty_print(cout);
        // cout << bitset_to_key(Indicator[next_supp]);
        // extending LocalSolutions to full coordinates
        Matrix<IntegerRet> LocalSolutionsGlobal(LocalSolutions.nr_of_rows(), EmbDim);
        for(size_t i = 0; i < LocalSolutions.nr_of_rows(); ++i){
            for(size_t j = 0; j < LocalKey.size(); ++j)
                LocalSolutionsGlobal[i][LocalKey[j]] = LocalSolutions[i][j];
        }
        AllLocalSolutions[coord] = LocalSolutionsGlobal;
        
        // now the intersections and new_coords
        dynamic_bitset intersection_coods = covered & Indicator[next_supp];
        dynamic_bitset new_coords(Indicator[next_supp]);
        for(size_t i = 0; i< EmbDim; ++i)
            if(new_coords[i] && covered[i])
                new_coords[i] = 0;

        vector<key_t> intersection_key = bitset_to_key(intersection_coods); // w.r.t. full coordinates
        vector<key_t> new_coords_key = bitset_to_key(new_coords); // w.r.t. to full coordinates
        // cout << "intersection " << intersection_key;
        // cout << "New coords   " << new_coords_key;
        AllIntersections_key[coord] = intersection_key;
        AllNew_coords_key[coord] = new_coords_key;
        
        dynamic_bitset new_covered = covered | Indicator[next_supp];
        Matrix<IntegerRet> ExtraInequalities(0, EmbDim);
        for(size_t i = 0; i <nr_all_supps ; ++i){
            if(used_supps[i])
                continue;
            if(Indicator[i].is_subset_of(covered))
                continue;
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
        AllExtraInequalities[coord] = ExtraInequalities;    
        
        vector<key_t> PolyEqusKey, PolyInequsKey;
        for(size_t i = 0; i < PolyEquations.size(); ++i){
            if(!PolyEquations[i].support.is_subset_of(new_covered))
                continue;
            if(PolyEquations[i].support.is_subset_of(covered))
                continue;
            PolyEqusKey.push_back(i);
        }
        AllPolyEqusKey[coord] = PolyEqusKey;
        
        for(size_t i = 0; i < PolyInequalities.size(); ++i){
            if(!(PolyInequalities[i]).support.is_subset_of(new_covered))
                continue;
            if(PolyInequalities[i].support.is_subset_of(covered))
                continue;
            PolyInequsKey.push_back(i);
        }
        AllPolyInequsKey[coord] = PolyInequsKey;
        
        covered = new_covered;
 
        if(verbose)
            verboseOutput() << "nr_covered coordinates " << covered.count() << " coordinates " << bitset_to_key(covered);
        
    } // coord
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::compute_latt_points_by_matching() {
    
    vector<IntegerRet> start(EmbDim);
    start[0] = GD;
    list<vector<IntegerRet> > start_list;
    start_list.push_back(start);
    extend_points_to_next_coord(start_list, 1);
    // cout << "TTTT " << TotalNrLP << endl;
    NrLP[EmbDim] = TotalNrLP;
    /* if (verbose) {
        for (size_t i = 2; i < NrLP.size(); ++i)
            verboseOutput() << "embdim " << i << " LatticePoints " << NrLP[i] << endl;
    }*/
    if(verbose)
        verboseOutput() << "Final number of lattice points "  << NrLP[EmbDim] << endl;
}

//---------------------------------------------------------------------------

template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL,IntegerRet>::extend_points_to_next_coord(list<vector<IntegerRet> >& LatticePoints, const key_t start_coord) {

    size_t last_active_coord = 0;
    for(size_t i = 0; i < active_coords.size(); ++i){
        if(active_coords[i])
            last_active_coord = i;
    }
    
    size_t max_nr_per_thread = 1000000 / omp_get_max_threads();
    
    size_t coord; // the coord to which we want to extend (at least)
    for(coord = start_coord; coord < EmbDim; coord++)
        if(active_coords[coord])
            break;
        
    StartTime();
    size_t nr_latt_points_linear = 0;
    
    if(verbose)
        verboseOutput() << "coordinate " << coord << endl;
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    vector<key_t>& intersection_key = AllIntersections_key[coord];
    vector<key_t>& new_coords_key = AllNew_coords_key[coord];
    vector<key_t>& PolyEqusKey = AllPolyEqusKey[coord];
    vector<key_t>& PolyInequsKey = AllPolyInequsKey[coord];
    Matrix<IntegerRet>& ExtraInequalities = AllExtraInequalities[coord];
    Matrix<IntegerRet>& LocalSolutionsGlobal = AllLocalSolutions[coord];
    
    map<vector<IntegerRet>, vector<key_t> >LocalSolutions_by_intersecion;
    vector<IntegerRet> overlap(intersection_key.size());
    for(size_t i = 0; i < LocalSolutionsGlobal.nr_of_rows(); i++){
        overlap = v_select_coordinates(LocalSolutionsGlobal[i], intersection_key);
        LocalSolutions_by_intersecion[overlap].push_back(i);
    }
    
    bool last_coord = (coord == last_active_coord);
    
    if(PolyEqusKey.size() > 0 && verbose)
        verboseOutput() << "Pplynomial equations " << PolyEqusKey.size() << endl;
    if(PolyInequsKey.size() > 0 && verbose)
        verboseOutput() << "Polynomial inequalities " << PolyEqusKey.size() << endl;
    
    size_t nr_to_match = LatticePoints.size();
    size_t nr_points_matched = 0;
        
    while (true) {
    // cout << "Durchgang dim " << dim << endl;

        // bool message_printed = false;
        
        if(verbose)
            verboseOutput() << "lattice points left " << nr_to_match - nr_points_matched << endl;

        bool skip_remaining;
        std::exception_ptr tmp_exception;

        skip_remaining = false;
        int omp_start_level = omp_get_level();
        
        
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

        try{
            overlap = v_select_coordinates(*P, intersection_key);
            if(LocalSolutions_by_intersecion.find(overlap) == LocalSolutions_by_intersecion.end()){
                (*P)[0] = 0;
                continue; // vector cannot be extended
            }
            // now the extensions
            for(auto& i: LocalSolutions_by_intersecion[overlap]){
                
                INTERRUPT_COMPUTATION_BY_EXCEPTION
                
                NewLattPoint = *P;
                for(auto& j: new_coords_key)
                    NewLattPoint[j] = LocalSolutionsGlobal[i][j];
                bool can_be_inserted = true;
                for(size_t k = 0; k < ExtraInequalities.nr_of_rows(); ++k){
                    if(v_scalar_product(NewLattPoint, ExtraInequalities[k]) < 0){
                        can_be_inserted = false;
                        break;
                    }
                }
                if(can_be_inserted){
#pragma omp atomic
                    nr_latt_points_linear++;
                    for(auto& k: PolyEqusKey){
                        if((PolyEquations[k]).evaluate(NewLattPoint) != 0){
                            can_be_inserted = false;
                            // poly_equs_hit[k] = true;
                        }
                    }
                }
                if(can_be_inserted){
                    for(auto& k: PolyInequsKey){
                        if((PolyInequalities[k]).evaluate(NewLattPoint)< 0){
                            can_be_inserted = false;
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
        
        for (size_t i = 0; i < Deg1Thread.size(); ++i)
            NewLatticePoints.splice(NewLatticePoints.end(), Deg1Thread[i]);
        
        size_t nr_new_latt_points = NewLatticePoints.size();
        
        if(last_coord)
            collect_results(NewLatticePoints); // clears NewLatticePoints

        /* cout << "PolyEqus effective ";
        for(auto& k: PolyEqusKey){
            if(poly_equs_hit[k])
                cout << k << " "; 
        }
        cout << endl;*/
        
        if(verbose){
            verboseOutput() << "latt points satisfying linear constraints " << nr_latt_points_linear << endl;
            verboseOutput() << "latt points satisfying   all  constraints " << nr_new_latt_points << endl;
        }
        
        MeasureTime(verbose, "Elapsed ");
        
        if(verbose)
            verboseOutput() << "------------------" << endl;
        
        // cout << "LC LC LC " << last_active_coord << " " << coord << " " << last_coord << " not done " << not_done << endl;;
        
        if(!last_coord && NewLatticePoints.size() > 0)
            extend_points_to_next_coord(NewLatticePoints, coord + 1);
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

    /* for(size_t i=0;i<Order.size();++i)
        cout << Supps[Order[i]][dim-1] << " ";
    cout << endl;*/

    return Order;
}

//---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_projections_primitive(size_t dim){

    size_t dim1 = dim - 1;

    if (dim == 1)
        return;

    Matrix<IntegerPL> SuppsProj(0,dim1);

    // cout << "dim1 dim1 " << dim1 << endl;
    // AllSupps[EmbDim].debug_print('$');

    // InEqusByDim[dim1].resize(0,dim1);
    // cout << "AAAAAAA " << AllSupps[EmbDim].nr_of_rows() << endl;
    for(size_t i = 0; i< AllSupps[EmbDim].nr_of_rows(); ++i){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
            
        bool can_be_restricted = true;
        // cout << "ddd " << dim1  << " eeee " << EmbDim << " sss " << AllSupps[EmbDim][i].size() << endl;
        // cout << AllSupps[EmbDim][i];
        for(size_t j= dim1; j <= EmbDim-1; ++j){
            if(AllSupps[EmbDim][i][j] >0){
                // cout << "pos at " << j << " val " << AllSupps[EmbDim][i][j] << endl;
                can_be_restricted = false;
                break;
            }
        }
        if(can_be_restricted){
            vector<IntegerPL> Restriction = AllSupps[EmbDim][i];
            Restriction.resize(dim1);
            // cout << "rest " << Restriction;
            SuppsProj.append(Restriction);
        }
    }

    SuppsProj.remove_duplicate_and_zero_rows();

    if (verbose)
        verboseOutput() << "embdim " << dim << " inequalities " << SuppsProj.nr_of_rows() << endl;

    // SuppsProj.debug_print();

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

    // Supps.debug_print();

    if (dim == down_to)
        return;

    // Supps.pretty_print(cout);
    // cout << Ind;

    // cout << "SSS" << Ind.size() << " " << Ind;

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

    // cout << "Nach Pos/Neg " << EqusProj.nr_of_rows() << " " << Pos.size() << " " << Neg.size() << endl;

    // now the elimination, matching Pos and Neg

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
                // cout << Supps[p];
                // cout << Supps[n];
                vector<IntegerPL> new_equ = FM_comb(PosVal, Supps[n], NegVal, Supps[p], is_zero);
                // cout << "zero " << is_zero << endl;
                // cout << "=====================" << endl;
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
            /* cout << Supps[NegEquAt];
            cout << Supps[p];
            cout << new_supp;
            cout << "zero " << is_zero << endl;
            cout << "+++++++++++++++++++++" << endl; */
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
            /* cout << Supps[PosEquAt];
            cout << Supps[n];
            cout << new_supp;
            cout << "zero " << is_zero << endl;
            cout << "=====================" << endl;*/

            if (is_zero)  // cannot happen, but included for analogy
                continue;
            SuppsProj.append(new_supp);
            NewInd.push_back(Ind[n]);
        }
    }

    // cout << "Nach RGU " << EqusProj.nr_of_rows() << " " << SuppsProj.nr_of_rows() << endl;

    if (!rank_goes_up && !is_parallelotope) {  // must match pos and neg hyperplanes

        // cout << "Pos " << Pos.size() << " Neg " << Neg.size() << " Supps " << SuppsProj.nr_of_rows() << endl;

        skip_remaining = false;

        size_t min_nr_vertices = rank - 2;
        /*if(rank>=3){
            min_nr_vertices=1;
            for(long i=0;i<(long) rank -3;++i)
                min_nr_vertices*=2;

        }*/

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
        // if(nr_pos>10000)
        //    nr_pos=10000;
        size_t nr_neg = Neg.size();
        // if(nr_neg>10000)
        //    nr_neg=10000;

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

    // cout << "Nach FM " << EqusProj.nr_of_rows() << " " << SuppsProj.nr_of_rows() << endl;

    Ind.clear();  // no longer needed

    EqusProj.resize_columns(dim1);   // cut off the trailing 0
    SuppsProj.resize_columns(dim1);  // project hyperplanes

    // Equations have not yet been appended to support hypwerplanes
    EqusProj.row_echelon();  // reduce equations
    // cout << "Nach eche " << EqusProj.nr_of_rows() << endl;
    /* for(size_t i=0;i<EqusProj.nr_of_rows(); ++i)
        cout << EqusProj[i]; */
    SuppsProj.append(EqusProj);  // append them as pairs of inequalities
    EqusProj.scalar_multiplication(-1);
    SuppsProj.append(EqusProj);
    AllNrEqus[dim1] = EqusProj.nr_of_rows();
    // We must add indictor vectors for the equations
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
    // cout << LiftedGen;
    size_t check_supps = Supps.nr_of_rows();
    if (check_supps > 1000 && dim < EmbDim && !no_relax)
        check_supps = 1000;
    for (size_t j = 0; j < check_supps; ++j) {
        INTERRUPT_COMPUTATION_BY_EXCEPTION

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
    size_t dim = Deg1Proj.front().size() + 1;
    size_t dim1 = dim - 1;

    /*if (dim == EmbDim) {
        InEqus.debug_print('+');
    }*/

    list<vector<IntegerRet> > Deg1Lifted;  // to this dimension if < EmbDim   
    
    size_t max_nr_per_thread = 1000000 / omp_get_max_threads();

    size_t nr_to_lift = Deg1Proj.size();
    NrLP[dim1] += nr_to_lift;
    size_t already_lifted = 0;

    bool not_done = true;
    bool has_poly_equs = (PolyEquations.size() > 0);
    bool has_poly_inequs = (PolyInequalities.size() > 0);

    while (not_done) {
        // cout << "Durchgang dim " << dim << endl;

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

                            /*if(primitive){ // in this case we must check equations and true inequalities
                                if(InEqusByDim[EmbDim].nr_of_rows() > 0){
                                    if(!v_non_negative(InEqusByDim[dim].MxV(NewPoint)))
                                        continue;
                                }
                            }*/

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

        // cout << nr_to_lift << " " << already_lifted << endl;
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

    // cout << "RET FROM " << dim1 << endl;

    if(verbose && dim == EmbDim){
        verboseOutput() << "Complete lattice points so far " << TotalNrLP << endl;
    }

    return;

    /* Deg1.pretty_print(cout);
    cout << "*******************" << endl; */
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

    vector<IntegerRet> start(1, GD);
    list<vector<IntegerRet> > start_list;
    start_list.push_back(start);
    lift_points_to_this_dim(start_list);
    // cout << "TTTT " << TotalNrLP << endl;
    NrLP[EmbDim] = TotalNrLP;
    /* if (verbose) {
        for (size_t i = 2; i < NrLP.size(); ++i)
            verboseOutput() << "embdim " << i << " LatticePoints " << NrLP[i] << endl;
    }*/
    if(verbose)
        verboseOutput() << "Final number of lattice points "  << NrLP[EmbDim] << endl;
}

///---------------------------------------------------------------------------
template <typename IntegerPL, typename IntegerRet>
void ProjectAndLift<IntegerPL, IntegerRet>::compute_latt_points_float() {
    
    cout << "FFFFFFF " << Deg1Thread.size() << endl;
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
    StartRank = rank;
    GD = 1;  // the default choice
    verbose = true;
    is_parallelotope = false;
    no_crunch = true;
    use_LLL = false;
    no_relax = false;
    primitive = false;
    sparse = false;
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
    // a grading denominator 1=1 can be accomodated.
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
        // Aconv.transpose().pretty_print(cout);
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
    
    if(primitive){
        if(verbose)
            verboseOutput() << "Checking if matching possible" << endl; 
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

    /* for(size_t i = 0; i <= EmbDim; ++i){
        cout << "iiiiii " << i << endl;
        AllSupps[i].debug_print();
    }*/
    if (all_points) {
        if(sparse){
            if(verbose)
                verboseOutput() << "Matching" << endl;
            compute_latt_points_by_matching();            
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

    // cout << " POS " << h_vec_pos;
    // cout << " NEG " << h_vec_neg;
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
            /* cout << "ori  " << Deg1Points.front();
            cout << "tra  " << LLL_Coordinates.from_sublattice(Deg1Points.front());
            cout << "tra1 " << LLL_Coordinates.A.VxM(Deg1Points.front());*/
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
    // SuppsRet.resize_colums(equs_start_in_row,in_dim);
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
