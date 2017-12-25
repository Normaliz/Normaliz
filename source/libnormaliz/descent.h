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

#ifndef DESCENT_H_
#define DESCENT_H_

#include <vector>
#include <set>
#include <list>
#include <map>
#include <boost/dynamic_bitset.hpp>

namespace libnormaliz {
using std::vector;
using std::map;
using std::pair;

template<typename Integer> class DescentSystem;

template<typename Integer> class DescentFace {
    
public:
    
    size_t dim; // cone dimension of the face
    mpq_class multiplicity;
    // bool facets_computed;
    // bool multiplicity_computed;
    bool simplicial;
    size_t tree_size;
    boost::dynamic_bitset<> own_facets; // own_facets[i]==true <==> SuppHyps[i] contains this face
    
    
    libnormaliz::key_t selected_gen; // the generator of C selected for descent
    vector<boost::dynamic_bitset<> > opposite_facets; // facets opposite to the selected generator,
                                                       // identified by the SuppsHyps containing them    
    DescentFace();    
    DescentFace(const size_t dim_given, const boost::dynamic_bitset<>& facets_given);
    vector<boost::dynamic_bitset<> >& compute(DescentSystem<Integer>& FF); 
    // returns the list of facets for descent from *this (identified by mother)
    
    void compute_multiplicity(DescentSystem<Integer>& FF);    
};

template<typename Integer> class DescentSystem {
    
public:
    
    bool verbose;
    
    Matrix<Integer> Gens;
    Matrix<Integer> SuppHyps;
    vector<Integer> Grading;
    vector<Integer> GradGens;
    
    size_t dim;
    size_t nr_supphyps;
    size_t nr_gens;
    
    size_t descent_steps;
    size_t nr_simplicial;
    
    vector<boost::dynamic_bitset<> > SuppHypInd;
    
    vector<map<boost::dynamic_bitset<>, DescentFace<Integer> > > Faces;
    
    mpq_class multiplicity;
    
    DescentSystem(const Matrix<Integer>& Gens, const Matrix<Integer>& SuppHyps, const vector<Integer>& Grading);
    void build();
    void compute_multiplicities();
    bool set_verbose(bool onoff);
    mpq_class getMultiplicity();
    void compute();
};



} // end namespace

#endif /* DESCENT_H_ */