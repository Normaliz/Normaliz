/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#ifndef CONE_H_
#define CONE_H_

#include <vector>
#include <map>
#include <utility> //for pair
//#include <boost/dynamic_bitset.hpp>

#include <libQnormaliz/libQnormaliz.h>
#include <libQnormaliz/Qcone_property.h>
#include <libQnormaliz/Qsublattice_representation.h>
#include <libQnormaliz/Qmatrix.h>
// #include <libQnormaliz/QHilbertSeries.h>

namespace libQnormaliz {
using std::vector;
using std::map;
using std::pair;

template<typename Number> class Full_Cone;
//template<typename Number> class Matrix;

// type for simplex, short in contrast to class Simplex
template<typename Number> struct SHORTSIMPLEX {
    vector<key_t> key;                // full key of simplex
    Number height;                   // height of last vertex over opposite facet
    Number vol;                      // volume if computed, 0 else
    vector<bool> Excluded;           // for disjoint decomposition of cone
                                      // true in position i indictate sthat the facet 
                                      // opposite of generator i must be excluded
};


template<typename Number>
class Cone {

//---------------------------------------------------------------------------
//                               public methods
//---------------------------------------------------------------------------
public:

//---------------------------------------------------------------------------
//                    Constructors, they preprocess the input
//---------------------------------------------------------------------------

    /* give up to 3 matrices as input
     * the types must be pairwise different
     */
    Cone(InputType type, const vector< vector<Number> >& input_data);

    Cone(InputType type1, const vector< vector<Number> >& input_data1,
         InputType type2, const vector< vector<Number> >& input_data2);

    Cone(InputType type1, const vector< vector<Number> >& input_data1,
         InputType type2, const vector< vector<Number> >& input_data2,
         InputType type3, const vector< vector<Number> >& input_data3);

    /* give multiple input */
    Cone(const map< InputType , vector< vector<Number> > >& multi_input_data);
    
    // Now with Matrix
    Cone(InputType type, const Matrix<Number>& input_data);

    Cone(InputType type1, const Matrix<Number>& input_data1,
         InputType type2, const Matrix<Number>& input_data2);

    Cone(InputType type1, const Matrix<Number>& input_data1,
         InputType type2, const Matrix<Number>& input_data2,
         InputType type3, const Matrix<Number>& input_data3);

    /* give multiple input */
    Cone(const map< InputType , Matrix<Number> >& multi_input_data);
//---------------------------------------------------------------------------
//                                Destructor
//---------------------------------------------------------------------------

    ~Cone();

//---------------------------------------------------------------------------
//                          give additional data
//---------------------------------------------------------------------------

    /* Sets if the Cone prints verbose output.
     * The default value for the Cone is the global verbose.
     * returns the old value
     */
    bool setVerbose (bool v);

//---------------------------------------------------------------------------
//                           make computations
//---------------------------------------------------------------------------

    // return what was NOT computed
    // ConeProperties compute(ComputationMode mode = Mode::hilbertBasisSeries); //default: everything
    ConeProperties compute(ConeProperties ToCompute);
    // special case for up to 3 CPs
    ConeProperties compute(ConeProperty::Enum);
    ConeProperties compute(ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties compute(ConeProperty::Enum, ConeProperty::Enum, ConeProperty::Enum);

//---------------------------------------------------------------------------
//                         check what is computed
//---------------------------------------------------------------------------

    bool isComputed(ConeProperty::Enum prop) const;
    //returns true, when ALL properties in CheckComputed are computed
    bool isComputed(ConeProperties CheckComputed) const;

//---------------------------------------------------------------------------
//   get the results, these methods will start a computation if necessary
//   throws an NotComputableException if not succesful
//---------------------------------------------------------------------------

    // dimension and rank invariants
    size_t getEmbeddingDim() const { return dim; };   // is always known
    size_t getRank();                           // depends on ExtremeRays
    Number getIndex(); // depends on OriginalMonoidGenerators
    Number getInternalIndex(); // = getIndex()
    Number getUnitGroupIndex(); // ditto
    // only for inhomogeneous case:
    size_t getRecessionRank();
    long getAffineDim();
    size_t getModuleRank();

    const Matrix<Number>& getGeneratorsMatrix();
    const vector< vector<Number> >& getGenerators();
    size_t getNrGenerators();

    const Matrix<Number>& getExtremeRaysMatrix();
    const vector< vector<Number> >& getExtremeRays();
    size_t getNrExtremeRays();

    const Matrix<Number>& getVerticesOfPolyhedronMatrix();
    const vector< vector<Number> >& getVerticesOfPolyhedron();
    size_t getNrVerticesOfPolyhedron();

    const Matrix<Number>& getSupportHyperplanesMatrix();
    const vector< vector<Number> >& getSupportHyperplanes();
    size_t getNrSupportHyperplanes();
    
    const Matrix<Number>& getMaximalSubspaceMatrix();
    const vector< vector<Number> >& getMaximalSubspace();
    size_t getDimMaximalSubspace();

    // depends on the ConeProperty::s SupportHyperplanes and Sublattice
    map< InputType, vector< vector<Number> > > getConstraints();

    size_t getTriangulationSize();
    Number getTriangulationDetSum();

    // the actual grading is Grading/GradingDenom
    vector<Number> getGrading();
    Number getGradingDenom();

    vector<Number> getDehomogenization();
    
    bool inequalities_present;

    bool isPointed();
    bool isInhomogeneous();
    bool isDeg1ExtremeRays();

    const Matrix<Number>& getOriginalMonoidGeneratorsMatrix();
    const vector< vector<Number> >& getOriginalMonoidGenerators();
    size_t getNrOriginalMonoidGenerators();

    const Sublattice_Representation<Number>& getSublattice();
    // the following 2 methods give information about the last used triangulation
    // if no triangulation was computed so far they return false
    bool isTriangulationNested();
    bool isTriangulationPartial();
    const vector< pair<vector<key_t>, Number> >& getTriangulation();
    const vector< vector<bool> >& getOpenFacets();
    
    const renf_class* getRenf() const;

//---------------------------------------------------------------------------
//                          private part
//---------------------------------------------------------------------------

private:
    size_t dim;

    Sublattice_Representation<Number> BasisChange;  //always use compose_basis_change() !
    Sublattice_Representation<Number> BasisChangePointed; // to the pointed cone
    bool BC_set;
    bool verbose;
    ConeProperties is_Computed;
    // Matrix<Number> GeneratorsOfToricRing;
    Matrix<Number> OriginalMonoidGenerators;
    Matrix<Number> Generators;
    Matrix<Number> ExtremeRays;
    vector<bool> ExtremeRaysIndicator;
    Matrix<Number> VerticesOfPolyhedron;
    Matrix<Number> SupportHyperplanes;
    Matrix<Number> ExcludedFaces;
    Matrix<Number> PreComputedSupportHyperplanes;
    size_t TriangulationSize;
    Number TriangulationDetSum;
    bool triangulation_is_nested;
    bool triangulation_is_partial;
    vector< pair<vector<key_t>, Number> > Triangulation;
    vector<vector<bool> > OpenFacets;
    vector< pair<vector<key_t>, long> > InExData;
    // mpq_class multiplicity;
    vector<Number> WitnessNotIntegrallyClosed;
    Matrix<Number> HilbertBasis;
    Matrix<Number> BasisMaxSubspace;
    Matrix<Number> ModuleGeneratorsOverOriginalMonoid;
    Matrix<Number> Deg1Elements;

    vector<Number> Grading;
    vector<Number> Dehomogenization;
    Number GradingDenom;
    Number index;  // the internal index
    Number unit_group_index;

    bool pointed;
    bool inhomogeneous;

    int affine_dim; //dimension of polyhedron
    size_t recession_rank; // rank of recession monoid
    size_t module_rank; // for the inhomogeneous case
    Matrix<Number> ModuleGenerators;

    bool explicit_HilbertSeries;
    bool naked_dual;

    Matrix<Number> WeightsGrad;
    vector<bool> GradAbs;
    
    renf_class *Renf;

    bool no_lattice_restriction; // true if cine generators are known to be in the relevant lattice
    bool normalization; // true if input type normalization is used

    // if this is true we allow to change to a smaller integer type in the computation
    bool change_integer_type;
    

    void compose_basis_change(const Sublattice_Representation<Number>& SR); // composes SR

    // main input processing
    void process_multi_input(const map< InputType, vector< vector<Number> > >& multi_input_data);
    void prepare_input_lattice_ideal(map< InputType, vector< vector<Number> > >& multi_input_data);
    void prepare_input_constraints(const map< InputType, vector< vector<Number> > >& multi_input_data,
            Matrix<Number>& equations, Matrix<Number>& congruence, Matrix<Number>& Inequalities);
    void prepare_input_generators(map< InputType, vector< vector<Number> > >& multi_input_data,
                     Matrix<Number>& LatticeGenerators);
    void homogenize_input(map< InputType, vector< vector<Number> > >& multi_input_data);
    void check_precomputed_support_hyperplanes();

    
    void setWeights ();
    void setDehomogenization (const vector<Number>& lf);

    void checkDehomogenization();
    void check_vanishing_of_grading_and_dehom();
    void process_lattice_data(const Matrix<Number>& LatticeGenerators, Matrix<Number>& Congruences, Matrix<Number>& Equations);
    
    ConeProperties recursive_compute(ConeProperties ToCompute);

    Matrix<Number> prepare_input_type_2(const vector< vector<Number> >& Input);
    Matrix<Number> prepare_input_type_3(const vector< vector<Number> >& Input);
    void prepare_input_type_4(Matrix<Number>& Inequalities);

    /* only used by the constructors */
    void initialize();

    template<typename NumberFC>
    void compute_inner(ConeProperties& ToCompute);

    /* compute the generators using the support hyperplanes */
    void compute_generators();
    template<typename NumberFC>
    void compute_generators_inner();

    /* compute method for the dual_mode, used in compute(mode) */
    void compute_dual(ConeProperties& ToCompute);
    template<typename NumberFC>
    void compute_dual_inner(ConeProperties& ToCompute);
    
    void set_implicit_dual_mode(ConeProperties& ToCompute);

    /* extract the data from Full_Cone, this may remove data from Full_Cone!*/
    template<typename NumberFC>
    void extract_data(Full_Cone<NumberFC>& FC);
    template<typename NumberFC>
    void extract_supphyps(Full_Cone<NumberFC>& FC);
    
    void extract_supphyps(Full_Cone<Number>& FC);


    /* set OriginalMonoidGenerators */
    void set_original_monoid_generators(const Matrix<Number>&);

    /* set ExtremeRays, in inhomogeneous case also VerticesOfPolyhedron */
    void set_extreme_rays(const vector<bool>&);

    /* If the Hilbert basis and the original monoid generators are computed,
     * use them to check whether the original monoid is integrally closed. */
    void check_integrally_closed();
    void compute_unit_group_index();
    /* try to find a witness for not integrally closed in the Hilbert basis */
    void find_witness();

    Number compute_primary_multiplicity();
    template<typename NumberFC>
    Number compute_primary_multiplicity_inner();
    
    void compute_integer_hull();
    void complete_sublattice_comp(ConeProperties& ToCompute); // completes the sublattice computations
    void complete_HilbertSeries_comp(ConeProperties& ToCompute);
    
    void set_renf(renf_class *GivenRenf);

};

// helpers

template<typename Number>
vector<vector<Number> > find_input_matrix(const map< InputType, vector< vector<Number> > >& multi_input_data,
                               const InputType type);

template<typename Number>
void insert_zero_column(vector< vector<Number> >& mat, size_t col);

template<typename Number>
void insert_column(vector< vector<Number> >& mat, size_t col, Number entry);

}  //end namespace libQnormaliz

#endif /* CONE_H_ */
