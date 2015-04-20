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

#include "libnormaliz.h"
#include "cone_property.h"
#include "sublattice_representation.h"
#include "matrix.h"
#include "HilbertSeries.h"

namespace libnormaliz {
using std::vector;
using std::map;
using std::pair;

template<typename Integer> class Full_Cone;
//template<typename Integer> class Matrix;

// type for simplex, short in contrast to class Simplex
template<typename Integer> struct SHORTSIMPLEX {
    vector<key_t> key;                // full key of simplex
    Integer height;                   // height of last vertex over opposite facet
    Integer vol;                      // volume if computed, 0 else
};

template<typename Integer> struct STANLEYDATA {
    vector<key_t> key;
    Matrix<Integer> offsets;
};

template<typename Integer>
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
    Cone(InputType type, const vector< vector<Integer> >& input_data);

    Cone(InputType type1, const vector< vector<Integer> >& input_data1,
         InputType type2, const vector< vector<Integer> >& input_data2);

    Cone(InputType type1, const vector< vector<Integer> >& input_data1,
         InputType type2, const vector< vector<Integer> >& input_data2,
         InputType type3, const vector< vector<Integer> >& input_data3);

    /* give multiple input */
    Cone(const map< InputType , vector< vector<Integer> > >& multi_input_data);

    /* DEPRECATED ** use version with type as first argument! */
    Cone(const vector< vector<Integer> >& input_data,
         InputType type = Type::integral_closure);


//---------------------------------------------------------------------------
//                          give additional data
//---------------------------------------------------------------------------

    /* Sets the linear form which is used to grade. 
     * It has to be an N-grading, i.e. all generators must have a value >=1.
     * If it is not, a subclass of NormalizException will be thrown at the 
     * time of detection which can be in this method or later!
     * It will delete all data from the cone that depend on the grading!
     */
    void setGrading (const vector<Integer>& lf);
    void setDehomogenization (const vector<Integer>& lf);


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
//          get the results, these methods do not start a computation
//---------------------------------------------------------------------------

    // dimension and rank invariants
    size_t getEmbeddingDim() const { return dim; };   // is always known
    size_t getRank() const;                           // depends on ExtremeRays
    // only for inhomogeneous case:
    size_t getRecessionRank() const;
    long getAffineDim() const;
    size_t getModuleRank() const;

    const Matrix<Integer>& getGeneratorsMatrix() const;
    const vector< vector<Integer> >& getGenerators() const;
    size_t getNrGenerators() const;

    const Matrix<Integer>& getExtremeRaysMatrix() const;
    const vector< vector<Integer> >& getExtremeRays() const;
    size_t getNrExtremeRays() const;

    const Matrix<Integer>& getVerticesOfPolyhedronMatrix() const;
    const vector< vector<Integer> >& getVerticesOfPolyhedron() const;
    size_t getNrVerticesOfPolyhedron() const;

    // The following constraints depend all on ConeProperty::SupportHyperplanes
    const Matrix<Integer>& getSupportHyperplanesMatrix() const;
    const vector< vector<Integer> >& getSupportHyperplanes() const;
    size_t getNrSupportHyperplanes() const;
    const Matrix<Integer>& getEquationsMatrix() const;
    const vector< vector<Integer> >& getEquations() const;
    const Matrix<Integer>& getCongruencesMatrix() const;
    const vector< vector<Integer> >& getCongruences() const;
    map< InputType, vector< vector<Integer> > > getConstraints() const;

    const Matrix<Integer>& getExcludedFacesMatrix() const;
    const vector< vector<Integer> >& getExcludedFaces() const;

    size_t getTriangulationSize() const;
    Integer getTriangulationDetSum() const;

    const Matrix<Integer>& getHilbertBasisMatrix() const;
    const vector< vector<Integer> >& getHilbertBasis() const;
    size_t getNrHilbertBasis() const;

    const Matrix<Integer>& getModuleGeneratorsMatrix() const;
    const vector< vector<Integer> >& getModuleGenerators() const;
    size_t getNrModuleGenerators() const;

    const Matrix<Integer>& getDeg1ElementsMatrix() const;
    const vector< vector<Integer> >& getDeg1Elements() const;
    size_t getNrDeg1Elements() const;

    // the actual grading is Grading/GradingDenom
    vector<Integer> getGrading() const;
    Integer getGradingDenom() const;

    vector<Integer> getDehomogenization() const;

    mpq_class getMultiplicity() const;

    bool isPointed() const;
    bool isInhomogeneous() const;
    bool isDeg1ExtremeRays() const;
    bool isDeg1HilbertBasis() const;
    bool isIntegrallyClosed() const;
    bool isReesPrimary() const;
    Integer getReesPrimaryMultiplicity() const;
    Integer getShift() const;
    const Matrix<Integer>& getGeneratorsOfToricRingMatrix() const;
    const vector< vector<Integer> >& getGeneratorsOfToricRing() const;
    const Matrix<Integer>& getOriginalMonoidGeneratorsMatrix() const;
    const vector< vector<Integer> >& getOriginalMonoidGenerators() const;
    // the following methods return const refs to avoid copying of big data objects
    const Sublattice_Representation<Integer>& getBasisChange() const;
    const HilbertSeries& getHilbertSeries() const; //general purpose object
    const vector< pair<vector<key_t>, Integer> >& getTriangulation() const;
    const vector< pair<vector<key_t>, long> >& getInclusionExclusionData() const;
    const list< STANLEYDATA<Integer> >& getStanleyDec() const;

//---------------------------------------------------------------------------
//                          private part
//---------------------------------------------------------------------------

private:
    size_t dim;

    Sublattice_Representation<Integer> BasisChange;  //always use compose_basis_change() !
    bool BC_set;
    ConeProperties is_Computed;
    Matrix<Integer> GeneratorsOfToricRing;
    Matrix<Integer> OriginalMonoidGenerators;
    Matrix<Integer> Generators;
    Matrix<Integer> ExtremeRays;
    vector<bool> ExtremeRaysIndicator;
    Matrix<Integer> VerticesOfPolyhedron;
    Matrix<Integer> SupportHyperplanes;
    Matrix<Integer> ExcludedFaces;
    size_t TriangulationSize;
    Integer TriangulationDetSum;
    vector< pair<vector<key_t>, Integer> > Triangulation;
    vector< pair<vector<key_t>, long> > InExData;
    list< STANLEYDATA<Integer> > StanleyDec;
    mpq_class multiplicity;
    Matrix<Integer> HilbertBasis;
    Matrix<Integer> Deg1Elements;
    HilbertSeries HSeries;
    vector<Integer> Grading;
    vector<Integer> Dehomogenization;
    Integer GradingDenom;
    bool pointed;
    bool inhomogeneous;
    bool deg1_extreme_rays;
    bool deg1_hilbert_basis;
    bool integrally_closed;
    bool rees_primary;
    Integer ReesPrimaryMultiplicity;
    Integer shift; // needed in the inhomogeneous case to make degrees positive
    int affine_dim; //dimension of polyhedron
    size_t recession_rank; // rank of recession monoid
    size_t module_rank; // for the inhomogeneous case
    Matrix<Integer> ModuleGenerators;
    
    Matrix<Integer> WeightsGrad;
    vector<bool> GradAbs;

    void compose_basis_change(const Sublattice_Representation<Integer>& SR); // composes SR

    // main input processing
    void process_multi_input(const map< InputType, vector< vector<Integer> > >& multi_input_data);
    void prepare_input_lattice_ideal(const map< InputType, vector< vector<Integer> > >& multi_input_data);
    void prepare_input_constraints(const map< InputType, vector< vector<Integer> > >& multi_input_data);
    void prepare_input_generators(const map< InputType, vector< vector<Integer> > >& multi_input_data);
    void homogenize_input(map< InputType, vector< vector<Integer> > >& multi_input_data);
    void check_trunc_nonneg(const vector< vector<Integer> >& input_gens);

    /* Progress input for subtypes */
    void prepare_input_type_0(const vector< vector<Integer> >& Input);
    void prepare_input_type_1(const vector< vector<Integer> >& Input);
    void prepare_input_type_2(const vector< vector<Integer> >& Input);
    void prepare_input_type_3(const vector< vector<Integer> >& Input);

    void prepare_input_type_456(const Matrix<Integer>& Congruences, const Matrix<Integer>& Equations, Matrix<Integer>& Inequalities);
    void prepare_input_type_45(const Matrix<Integer>& Equations, Matrix<Integer>& Inequalities);

    /* only used by the constructors */
    void initialize();

    /* compute the generators using the support hyperplanes */
    void compute_generators();

    /* compute method for the dual_mode, used in compute(mode) */
    ConeProperties compute_dual(ConeProperties ToCompute);

    /* extract the data from Full_Cone, this may remove data from Full_Cone!*/
    void extract_data(Full_Cone<Integer>& FC);

    /* set OriginalMonoidGenerators and Generators */
    void set_original_monoid_generators(const Matrix<Integer>&);

    /* set ExtremeRays, in inhomogeneous case also VerticesOfPolyhedron */
    void set_extreme_rays(const vector<bool>&);

    void check_integrally_closed();

    /* set this object to the zero cone */
    void set_zero_cone();
    
    Integer compute_primary_multiplicity();
    
};

// helpers

template<typename Integer>
vector<vector<Integer> > find_input_matrix(const map< InputType, vector< vector<Integer> > >& multi_input_data, 
                               const InputType type);

template<typename Integer>
void insert_zero_column(vector< vector<Integer> >& mat, size_t col);

template<typename Integer>
void insert_column(vector< vector<Integer> >& mat, size_t col, Integer entry);

}  //end namespace libnormaliz

#endif /* CONE_H_ */
