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

#ifndef LIBNORMALIZ_DESCENT_H_
#define LIBNORMALIZ_DESCENT_H_

#include <vector>
#include <set>
#include <list>
#include <map>

#include <libnormaliz/general.h>
#include <libnormaliz/matrix.h>
#include <libnormaliz/sublattice_representation.h>
#include "libnormaliz/dynamic_bitset.h"

namespace libnormaliz {
using std::map;
using std::pair;
using std::vector;

template <typename Integer>
class DescentSystem;

template <typename Integer>
class DescentFace {
   public:
    bool dead; // to be skipped in next round.
    // size_t dim; // cone dimension of the face
    mpq_class coeff;
    // bool facets_computed;
    // bool multiplicity_computed;
    // bool simplicial;
    size_t tree_size;  // the number of paths in the tree from top to to this face
    // dynamic_bitset own_facets; // own_facets[i]==true <==> SuppHyps[i] contains this face
    
    dynamic_bitset FacetsOfFace; // an indicator picking for each facet F of *this a facet of the cone
                                 // cutting out F from *this <== a minimal subset of global supphyps
                                 // cutting out ther facets of this face
                                 
    vector<dynamic_bitset> FacetOrbits; // orbits of FacetsOfFace under automorphism group of face
                                        //as subsets of FacetsOfFace    
    DescentFace();
    // DescentFace(const size_t dim_given, const dynamic_bitset& facets_given);
    
    void compute(DescentSystem<Integer>& FF,
                                   size_t dim,
                                   const dynamic_bitset& facets_cutting_out,
                                   size_t mother_tree_size
                );
    
    void find_sublattice(Matrix<Integer>& Gens_this, Sublattice_Representation<Integer>& Sublatt_this, 
                         bool& sub_latt_computed, vector<key_t> mother_key, size_t dim, 
                         const Matrix<Integer>& FF_Gens);
    
    void find_facets(map<dynamic_bitset, dynamic_bitset>& FacetInds, map<dynamic_bitset, key_t>& CutOutBy,
                                       map<dynamic_bitset, vector<key_t> >& SimpKeys, map<dynamic_bitset, vector<bool> >& SimpInds,
                     
                                       const bool ind_better_than_keys,                                       
                                       const DescentSystem<Integer>& FF, const vector<key_t>& mother_key, 
                                       const dynamic_bitset& facets_cutting_mother_out, size_t dim);
    void find_optimal_vertex(key_t& m_ind,
                   const DescentSystem<Integer>& FF, const map<dynamic_bitset, dynamic_bitset>& FacetInds, const vector<key_t>& mother_key);
    
    void make_simplicial_facet(map<dynamic_bitset, vector<key_t> >& SimpKeys, map<dynamic_bitset, vector<bool> >& SimpInds,
                                                 map<dynamic_bitset, key_t>& CutOutBy, const bool ind_better_than_keys, 
                                                 const DescentSystem<Integer>& FF, const vector<key_t>& mother_key,
                                                 // map<dynamic_bitset, dynamic_bitset>& FacetInds,
                                                 const dynamic_bitset& facet_ind, vector<key_t> facet_key);
    
    void find_facets_from_FacetsOfFace(map<dynamic_bitset, dynamic_bitset>& FacetInds, map<dynamic_bitset, key_t>& CutOutBy,
                                       map<dynamic_bitset, vector<key_t> >& SimpKeys, map<dynamic_bitset, vector<bool> >& SimpInds,
                                       vector<key_t>& Orbit,
                                       
                                       const bool ind_better_than_keys,                                       
                                       const DescentSystem<Integer>& FF, const vector<key_t>& mother_key, 
                                       const dynamic_bitset& facets_cutting_mother_out, size_t dim);
};

template <typename Integer>
class DescentSystem {
   public:
    bool verbose;

    Matrix<Integer> Gens;
    Matrix<Integer> SuppHyps;
    vector<Integer> Grading;
    vector<Integer> GradGens;
    vector<mpz_class> GradGens_mpz;

    bool SimplePolytope;
    bool exploit_automorphisms;
    
    size_t dim;
    size_t nr_supphyps;
    size_t nr_gens;

    size_t descent_steps;
    size_t nr_simplicial;
    size_t tree_size;
    size_t system_size;

    vector<dynamic_bitset> SuppHypInd;

    map<dynamic_bitset, DescentFace<Integer> > OldFaces;
    map<dynamic_bitset, DescentFace<Integer> > NewFaces;
    map< IsoType<Integer>, DescentFace<Integer>* , IsoType_compare<Integer> > Isos; // associate faces to isomorphism classes

    vector<size_t> OldNrFacetsContainingGen;
    vector<size_t> NewNrFacetsContainingGen;

    mpq_class multiplicity;

    DescentSystem(const Matrix<Integer>& Gens, const Matrix<Integer>& SuppHyps, const vector<Integer>& Grading);
    DescentSystem();
    void compute();
    void collect_old_faces_in_iso_classes();
    bool set_verbose(bool onoff);
    void setExploitAutoms(bool on);
    mpq_class getMultiplicity();
};

}  // namespace libnormaliz

#endif /* DESCENT_H_ */
