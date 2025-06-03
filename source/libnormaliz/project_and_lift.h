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

#ifndef LIBNORMALIZ_PROJECT_AND_LIFT_H_
#define LIBNORMALIZ_PROJECT_AND_LIFT_H_

#include <vector>
#include <list>

#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/sublattice_representation.h"
#include "libnormaliz/HilbertSeries.h"
#include "libnormaliz/dynamic_bitset.h"
#include "libnormaliz/nmz_polynomial.h"
#include "libnormaliz/fusion.h"
#include "libnormaliz/collect_lat.h"

namespace libnormaliz {
using std::vector;

const size_t max_nr_new_latt_points_total = 100000; // before we go one level up in the recursion
const size_t max_nr_new_latt_points_patching = 4000; // the same for patching
const size_t nr_extensions_for_elimination_equs = 10000;
const size_t nr_extensions_for_elimination_inequs = 10000;
const size_t nr_extensions_for_elimination_automs = 10000;

// the project-and-lift algorithm for lattice points in a polytope

template <typename Integer>
class Cone;

template <typename IntegerPL, typename IntegerRet>
class ProjectAndLift {
    template <typename, typename>
    friend class ProjectAndLift;

    size_t max_nr_new_latt_points_total;
    size_t max_nr_new_latt_points_patching;
    size_t nr_extensions_for_elimination_equs;
    size_t nr_extensions_for_elimination_inequs;
    size_t nr_extensions_for_elimination_automs;

    list<vector<IntegerRet> > start_list; // list of lattice points to start from
                                          // and to be lifted to full dimension

    vector<Matrix<IntegerPL> > AllSupps; // thwe support hyperplanes in all dimensions
    vector<Matrix<IntegerRet> > AllCongs; // congruences derived from Congs and equations by coarsening
    vector<vector<size_t> > AllOrders;
    vector<size_t> AllNrEqus;  // the numbers of equations --- well defined
                               // in dimensions < start dimension !!!!
    Matrix<IntegerRet> Congs;  // congruences used to sieve out points
    Matrix<IntegerPL> Vertices;  // only used for LLL coordinates

    Sublattice_Representation<IntegerRet> LLL_Coordinates;

    vector<dynamic_bitset> StartInd;
    vector<dynamic_bitset> StartPair;
    vector<dynamic_bitset> StartParaInPair;

    size_t StartRank;

    list<vector<IntegerRet> > FreeVect;
    vector<list<vector<IntegerRet> > > FreeVectThread;
    void store_new_vector(const vector<IntegerRet>& new_vect, const int tn);

    vector<list<vector<IntegerRet> > > Deg1Thread;
    vector<vector<num_t> > h_vec_pos_thread;
    vector<vector<num_t> > h_vec_neg_thread;

    list<vector<IntegerRet> > Deg1Points; // all Deg1Points found
    list<vector<short> > ShortDeg1Points; // all Deg1Points found
    vector<IntegerRet> SingleDeg1Point;
    vector<IntegerRet> excluded_point;
    IntegerRet GD;

    OurPolynomialSystem<IntegerRet> PolyEquations;
    OurPolynomialSystem<IntegerRet> PolyInequalities;
    OurPolynomialSystem<IntegerRet> RestrictablePolyInequs;

    vector<IntegerRet> Grading;
    size_t TotalNrLP;
    vector<size_t> NrLP;  // number of lattice points by dimension

    dynamic_bitset DoneWithDim;

    vector<num_t> h_vec_pos;
    vector<num_t> h_vec_neg;

    size_t EmbDim;
    bool verbose;

    bool is_parallelotope;
    bool no_crunch;  // indicates that the projection vector is never parallel to a facet of
                     // the parallelotope (in all dimensions)
    bool use_LLL;
    bool no_relax;
    bool count_only;
    bool use_coord_weights;
    bool no_weights;
    bool linear_order_patches;
    bool cong_order_patches;
    bool first_solution_printed;
    bool fusion_rings_computation;
    bool stored_local_solutions;
    bool use_short_int;
    bool no_heuristic_minimization;

    bool system_unsolvable;

    long min_return_patch;  // saet to the lowest patch that is not completely done before going up

    bool primitive; // true = using positive_bounded (a priori x >= 0 and upper bounds)
    bool sparse; // true = using the patching method
    bool patching_allowed; // if true blocks patching
    bool only_single_point;
    bool single_point_found;
    bool distributed_computation;

    long critical_coord_simplicity;
    bool check_simplicity_cand;
    bool check_simplicity_all;

    size_t nr_time_printed;

    struct timeval stop_ckeck_begin;

    // data for patching method
    vector<dynamic_bitset> Indicator; // indicaor of nonzero coordinates in inequality
    dynamic_bitset upper_bounds; // indicator of inequalities giving upper boounds
    dynamic_bitset max_sparse; // indicator of inequalities used in covering by "sparse" inequalities

    // data for patching depending on coordinates
    vector<key_t> InsertionOrderPatches;
    vector<key_t> LevelPatches; // index of coord in InsertionOrderPatches
    vector<vector<vector<IntegerRet> > > AllLocalSolutions; // "local" solutions that will be patched
    vector<vector<vector<short> > > AllShortLocalSolutions; // the same in short int format

    vector< map< vector<IntegerRet>, map< vector<IntegerRet>, vector<key_t>  > > >  AllLocalSolutions_by_intersection_and_cong;

    vector<vector<key_t> > AllIntersections_key;
    vector<vector<key_t> > AllNew_coords_key;
    vector<dynamic_bitset > AllCovered;
    vector<vector<key_t > > AllCoveredKey;  // key in the order of insertion of patches
    vector<vector<key_t > > AllCoveredKeyInverse;  // ionverse key
    vector<dynamic_bitset > AllPatches; // patches associated with the coordinates
    vector<ProjectAndLift<IntegerPL, IntegerRet> > AllLocalPL;
    dynamic_bitset active_coords;

    vector<vector<OurPolynomialCong<IntegerRet> > > AllCongsRestricted; // the congruences that can be restrictd to the unions of patches built successively
    vector< vector < pair<OurPolynomial<IntegerRet>, OurPolynomial<IntegerRet> > > > AllPolyEqus; // indexed by coord, poly equs (split into two components) applied with this coord
    vector<OurPolynomialSystem<IntegerRet> > AllPolyInequs; // ditto for inequalities
    vector<vector<key_t> > AllAutoms; // ditto for automorphisms where the key is the index

    vector<vector<vector < pair<OurPolynomial<IntegerRet>, OurPolynomial<IntegerRet> > > > > AllPolyEqusThread; // a copy for each thread
    vector<vector<OurPolynomialSystem<IntegerRet> > > AllPolyInequsThread; // ditto for inequalities

    dynamic_bitset used_supps; // registers which inequalities are used in the patching process
    dynamic_bitset poly_equs_minimized; // redisters whether the polynomial equations at a coord have been minimized
    dynamic_bitset poly_inequs_minimized; //  ditto for the inequalities from equations
    dynamic_bitset automs_minimized; //  ditto for automorphisms
    vector<dynamic_bitset> CongIndicator; // stores supports of congruences
    Matrix<double> WeightOfCoord;
    Matrix<IntegerPL> DefiningSupps;

    // Matrix<IntegerRet> SavedLocalSolutions;

    vector<double> ExpectedNrRounds;
    vector<double> TimeToLevel;
    vector<size_t> NrNodes;

    SplitData our_split;

    vector<size_t> NrRemainingLP;
    vector<size_t> NrDoneLP;

    vector<size_t> order_supps(const Matrix<IntegerPL>& Supps);
    bool fiber_interval(IntegerRet& MinInterval, IntegerRet& MaxInterval, const vector<IntegerRet>& base_point);

    void lift_point_recursively(vector<IntegerRet>& final_latt_point, const vector<IntegerRet>& latt_point_proj); // single point
    void lift_points_to_this_dim(list<vector<IntegerRet> >& Deg1Proj);

    void compute_latt_points_by_patching();
    void extend_points_to_next_coord(list<vector<IntegerRet> >& LatticePoints, const key_t this_patch);

    void find_single_point();
    void compute_latt_points();
    void compute_latt_points_float();
    void finalize_latt_point(vector<IntegerRet> NewPoint, const int tn);
    void collect_results(list<vector<IntegerRet> >& Deg1PointsComputed);
    void splice_into_short_deg1_points(list<vector<IntegerRet> >& Deg1PointsComputed);

    void compute_projections(size_t dim,
                             size_t down_to,
                             vector<dynamic_bitset>& Ind,
                             vector<dynamic_bitset>& Pair,
                             vector<dynamic_bitset>& ParaInPair,
                             size_t rank,
                             bool only_projections = false);

    void compute_projections_primitive(size_t dim);

    void initialize(const Matrix<IntegerPL>& Supps, size_t rank);

    void make_PolyEquations();
    bool check_PolyEquations(const vector<IntegerRet>& point, const size_t dim) const;
    void check_and_prepare_sparse();

    void reorder_coordinates(); // to use polynomial constraints as early as possible
    vector<pair<size_t, vector<key_t> > > cover_supports(const vector<dynamic_bitset>& supports);
    void compute_covers();

    void add_congruences_from_equations();
    void restrict_congruences();
    void find_order_congruences();
    void find_order_linear();
    bool order_patches_user_defined();
    void finalize_order(const dynamic_bitset& used_patches);

    void prepare_split(list<vector<IntegerRet> >& LatticePoints, const key_t& this_patch);
    void read_split_data();


    void compute_local_solutions(const key_t this_patch, list<vector<IntegerRet> >& start_list);
    void compute_local_solutions_for_saving();
    dynamic_bitset ImportedLocalSolutions; // registers the levels with imported local solutions

    // void make_LLL_coordinates();

   public:
    ProjectAndLift();
    ProjectAndLift(const Matrix<IntegerPL>& Supps, const vector<dynamic_bitset>& Ind, size_t rank);
    ProjectAndLift(const Matrix<IntegerPL>& Supps,
                   const vector<dynamic_bitset>& Pair,
                   const vector<dynamic_bitset>& ParaInPair,
                   size_t rank);
    template <typename IntegerPLOri, typename IntegerRetOri>
    ProjectAndLift(const ProjectAndLift<IntegerPLOri, IntegerRetOri>& Original);

    void set_excluded_point(const vector<IntegerRet>& excl_point);
    void set_grading_denom(const IntegerRet GradingDenom);
    void set_verbose(bool on_off);
    void set_short_int(bool on_off);
    void set_LLL(bool on_off);
    void set_no_relax(bool on_off);
    void set_primitive();
    void set_coord_weights(bool on_off);
    void set_no_weights(bool on_off);
    void set_patching_allowed(bool on_off);
    void set_vertices(const Matrix<IntegerPL>& Verts);
    void set_congruences(const Matrix<IntegerRet>& congruences);
    void set_grading(const vector<IntegerRet>& grad);
    void set_PolyEquations(const OurPolynomialSystem<IntegerRet>& PolyEqs, const bool minimize = true);
    void set_PolyInequalities(const OurPolynomialSystem<IntegerRet>& PolyInequs);
    void set_startList(const list<vector<IntegerRet> >& start_from);
    void set_linear_order_patches(const bool on_off);
    void set_cong_order_patches(const bool on_off);
    void set_distributed_computation(const bool on_off);
    void set_no_heuristic_minimization(bool on_off);

    void compute(bool do_all_points = true, bool lifting_float = false, bool count_only = false);
    void compute_only_projection(size_t down_to);

    void putSuppsAndEqus(Matrix<IntegerPL>& SuppsRet, Matrix<IntegerPL>& EqusRet, size_t in_dim);
    void put_eg1Points_into(Matrix<IntegerRet>& LattPoints);
    void put_deg1Points_into(vector<vector< IntegerRet> >& LattPoints); // we have no matrices for sjort
    void put_short_deg1Points_into(vector<vector<short> >& ShortLattPoints); // therefore vector<vector >
    void put_single_point_into(vector<IntegerRet>& LattPoint);
    void get_h_vectors(vector<num_t>& pos, vector<num_t>& neg) const;

    size_t getNumberLatticePoints() const;

    FusionComp<IntegerRet> fusion;
    void read_subring_data();

    void setOptions(const ConeProperties& ToCompute, const bool primitive, const bool our_verbose);
    void setFusion(const FusionBasic& FC);
};



// constructor by conversion

template <typename IntegerPL, typename IntegerRet>
template <typename IntegerPLOri, typename IntegerRetOri>
ProjectAndLift<IntegerPL, IntegerRet>::ProjectAndLift(const ProjectAndLift<IntegerPLOri, IntegerRetOri>& Original) {
    // The constructed PL is only good for lifting!!
    // Don't apply initialize to it

    EmbDim = Original.EmbDim;
    AllOrders = Original.AllOrders;
    verbose = Original.verbose;
    no_relax = Original.no_relax;
    only_single_point = false;
    use_short_int = false;
    convert(GD, Original.GD);
    AllSupps.resize(EmbDim + 1);
    AllCongs.resize(EmbDim + 1);
    for (size_t i = 0; i < AllSupps.size(); ++i)
        convert(AllSupps[i], Original.AllSupps[i]);
    convert(Congs, Original.Congs);
    for (size_t i = 0; i < Original.AllCongs.size(); ++i)
        convert(AllCongs[i], Original.AllCongs[i]);
    TotalNrLP = 0;
    Grading = Original.Grading;
    count_only = Original.count_only;
    NrLP.resize(EmbDim + 1);
    DoneWithDim.resize(EmbDim + 1);
    used_supps.resize(AllSupps[EmbDim].nr_of_rows());

    Deg1Thread.resize(omp_get_max_threads());
    h_vec_pos_thread.resize(omp_get_max_threads());
    h_vec_neg_thread.resize(omp_get_max_threads());

}

// computes c1*v1-c2*v2
template <typename Integer>
vector<Integer> FM_comb(Integer c1, const vector<Integer>& v1, Integer c2, const vector<Integer>& v2, bool& is_zero);

// setup function

template <typename Integer>
void project_and_lift(Cone<Integer>&  C, const ConeProperties& ToCompute,
                                     Matrix<Integer>& Deg1,
                                     const Matrix<Integer>& Gens,
                                     const Matrix<Integer>& Supps,
                                     const Matrix<Integer>& Congs,
                                     const vector<Integer>& GradingOnPolytope,
                                     const bool primitive,
                                     const OurPolynomialSystem<Integer>& PolyEqus,
                                     const OurPolynomialSystem<Integer>& PolyInequs );


}  // end namespace libnormaliz

#endif /* PROJECT_AND_LIFT_H_ */
