/*
 * Normaliz
 * Copyright (C) 2007-2019  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#ifndef LIBNORMALIZ_AUTOMORPHISM_H
#define LIBNORMALIZ_AUTOMORPHISM_H

#include <set>
#include <boost/dynamic_bitset.hpp>
#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h" 
#include "libnormaliz/nmz_nauty.h"
// #include "libnormaliz/HilbertSeries.h"

namespace libnormaliz {
using namespace std;

string quality_to_string(AutomParam::Quality quality);

template<typename Integer> class Cone;
template<typename Integer> class Full_Cone;
template<typename Integer> class Isomorphism_Classes;

template<typename Integer>
class AutomorphismGroup {
    
    template<typename> friend class Cone;
    template<typename> friend class Full_Cone;
    template<typename> friend class Isomorphism_Classes;
    
    Matrix<Integer> GensRef, SpecialGensRef, LinFormsRef, SpecialLinFormsRef;
    // the data defining the cone. Usially Gens = extreme rays, LinForms = support hyperplanes
    // SpecialGens: vectors to be left fixed
    // SpecialLinforms: grading, dehomogenization and possibly others
    
    Matrix<Integer> GensComp, LinFormsComp; //for computation
    
    bool addedComputationGens, addedComputationLinForms;
    
    vector<vector<key_t> > GenPerms; 
    vector<vector<key_t> > LinFormPerms;
    
    vector<vector<key_t> > ExtRaysPerms; // used in Cone
    vector<vector<key_t> > VerticesPerms;// ditto
    vector<vector<key_t> > SuppHypsPerms; // ditto
    
    vector<vector<key_t> > GenOrbits;
    vector<vector<key_t> > LinFormOrbits;
    
    vector<vector<key_t> > ExtRaysOrbits; // used in Cone
    vector<vector<key_t> > VerticesOrbits;// ditto
    vector<vector<key_t> > SuppHypsOrbits; // ditto
    
    vector<key_t> CanLabellingGens;
    
    vector<Matrix<Integer> > LinMaps;
    
    mpz_class order;
    
    size_t nr_special_gens, nr_special_linforms;

    set<AutomParam::Goals> is_Computed;
    set<AutomParam::Quality> Qualities;
    AutomParam::Method method;
    
    bool make_linear_maps_primal(const Matrix<Integer>& GivenGens,const vector<vector<key_t> >& ComputedGenPerms);
    void gen_data_via_lin_maps();
    void linform_data_via_lin_maps();
    void linform_data_via_incidence();
    void reset();
    
    void set_basic_gens_and_lin_forms(const Matrix<Integer>& ExtRays,  const Matrix<Integer>& SpecialGens,
                                const Matrix<Integer>& SuppHyps, const Matrix<Integer>& SpecialLinForms);
    
    bool compute_inner(const AutomParam::Quality& desired_quality, const bool force_gens_x_linforms=false);
    bool compute_integral();
    bool compute_polytopal(const AutomParam::Quality& desired_quality);
    void dualize();
    void swap_data_from_dual(AutomorphismGroup<Integer> Dual);
    void swap_data_from(AutomorphismGroup<Integer> Copy);

    
public:
    
    BinaryMatrix CanType; // see nauty
    
    const Matrix<Integer>& getGens() const;
    const Matrix<Integer>& getLinForms() const;
    const Matrix<Integer>& getSpecialLinForms() const;
    
    mpz_class getOrder() const;
    const vector<vector<key_t> >& getExtremeRaysPerms() const;
    const vector<vector<key_t> >& getVerticesPerms() const;
    const vector<vector<key_t> >& getSupportHyperplanesPerms() const;
    const vector<vector<key_t> >& getExtremeRaysOrbits() const;
    const vector<vector<key_t> >& getVerticesOrbits() const;
    const vector<vector<key_t> >& getSupportHyperplanesOrbits() const;
    const vector<Matrix<Integer> >& getLinMaps() const;
    const vector<key_t>& getCanLabellingGens() const;
    
    set<AutomParam::Quality> getQualities() const;
    AutomParam::Method getMethod() const;
    bool Is_Computed(AutomParam::Goals goal) const;
    string getQualitiesString() const;
    
    list<vector<Integer> > orbit_primal(const vector<Integer>& v) const;
    void add_images_to_orbit(const vector<Integer>& v,set<vector<Integer> >& orbit) const;
    
    BinaryMatrix getCanType();
    
    // bool compute(const AutomParam::Quality& desired_quality, const set<AutomParam::Goals>& ToCompute); // not yet implemented
    
    bool compute(const AutomParam::Quality& desired_quality, const bool force_gens_x_linforms=false);

    AutomorphismGroup();
    
    AutomorphismGroup(const Matrix<Integer>& ExtRays, const Matrix<Integer>& SupHyps, 
                      const Matrix<Integer>& SpecialLinForms);
    
    AutomorphismGroup(const Matrix<Integer>& ExtRays, const Matrix<Integer>& SpecialGens, 
                                              const Matrix<Integer>& SuppHyps, const Matrix<Integer>& SpecialLinForms);
    
    void addComputationGens(const Matrix<Integer>& GivenGens);
    void addComputationLinForms(const Matrix<Integer>& GivenLinearForms);
    
}; // end class

template<typename Integer> class Isomorphism_Classes;

template<typename Integer>
class IsoType {
    
    template<typename> friend class Isomorphism_Classes;

    size_t rank;
    Matrix<Integer> ExtremeRays;
    size_t nrExtremeRays;
    // Matrix<Integer> SupportHyperplanes;
    size_t nrSupportHyperplanes;
    Matrix<Integer> HilbertBasis; // without extreme rays
    // size_t nrHilbertBasis; // with extreme rays, but not used
    vector<Integer> Grading;
    vector<Integer> Truncation;
    // HilbertSeries HilbertSer;
    mpq_class Multiplicity;
    bool needs_Hilbert_basis;

    // For the coordinate transformation to the canonical basis
    vector<key_t> CanLabellingGens;
    Matrix<Integer> CanTransform;
    Integer CanDenom;
    vector<key_t> CanBasisKey;
    
    BinaryMatrix CanType;
    IsoType(); // constructs a dummy object

public:
    
    bool isOfType(const Full_Cone<Integer>& C) const;
    // bool isOfType(Cone<Integer>& C) const;
    
    IsoType(const Full_Cone<Integer>& C,bool& success); // success indicates whether a class could be created
    // IsoType(Cone<Integer>& C, bool slim=true);
    
    // size_t getRank();
    // Matrix<Integer> getExtremeRays() const;
    // Matrix<Integer> getSupportHyperplanes() const;
    const Matrix<Integer>& getHilbertBasis() const;
    // vector<Integer> getGrading() const;
    // vector<Integer> getTruncation() const;
    // HilbertSeries getHilbertSeries() const;
    mpq_class getMultiplicity() const;
    const Matrix<Integer>& getCanTransform() const;
    Integer getCanDenom() const;
    BinaryMatrix getCanType();
};

template<typename Integer>
class Isomorphism_Classes {
    
    template<typename> friend class Cone;
    template<typename> friend class Full_Cone;
    
    list<IsoType<Integer> >  Classes;
    
public:
    
    Isomorphism_Classes();
    const IsoType<Integer>& find_type(Full_Cone<Integer>& C, bool& found) const;
    const IsoType<Integer>& find_type(Cone<Integer>& C, bool& found) const;
    void add_type(Full_Cone<Integer>& C, bool& success);
    void add_type(Cone<Integer>& C);
};

vector<vector<key_t> > convert_to_orbits(const vector<key_t>& raw_orbits);

vector<vector<key_t> > cycle_decomposition(vector<key_t> perm, bool with_fixed_points=false);

void pretty_print_cycle_dec(const vector<vector<key_t> >& dec, ostream& out);

vector<vector<key_t> > keys(const list<boost::dynamic_bitset<> >& Partition);

list<boost::dynamic_bitset<> > partition(size_t n, const vector<vector<key_t> >& Orbits);

list<boost::dynamic_bitset<> > join_partitions(const list<boost::dynamic_bitset<> >& P1,
                                               const list<boost::dynamic_bitset<> >& P2);

vector<vector<key_t> > orbits(const vector<vector<key_t> >& Perms, size_t N);

} // namespace

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
