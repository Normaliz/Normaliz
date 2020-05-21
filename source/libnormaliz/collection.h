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

#ifndef LIBNORMALIZ_COLLECTION_H_
#define LIBNORMALIZ_COLLECTION_H_

#include <vector>
#include <map>
#include <set>
#include <string>
#include <utility>  // for pair

#include <libnormaliz/cone.h>

namespace libnormaliz {
    
template <typename Integer>
class ConeCollection;
    
template <typename Integer>
class MiniCone {
    
    friend class ConeCollection<Integer>;
    
    vector<key_t> GenKeys;
    bool is_simplex;
    bool dead;

    Matrix <Integer> Genereators;
    Matrix<Integer> SupportHyperplanes;
    Matrix<Integer> HilbertBasis;
    Integer multiplicity;
    
    ConeCollection<Integer>* Collection;
    
    // Cone<Integer> make_cone();
    list<MiniCone<Integer> > refine(const key_t key); 
    
    MiniCone<Integer>(const vector<key_t> GKeys, const Integer& mult, ConeCollection<Integer>& Coll);
        
};

template <typename Integer>
class ConeCollection {
    
    friend class MiniCone<Integer>;
    
    list<MiniCone<Integer> > Members;
    Matrix<Integer> Generators;
    set<vector<Integer> > AllRays;
    
    bool is_fan;
    bool is_triangulation;
    
    void refine(const key_t key);
    void add_minicone(const vector<key_t> GKeys, const Integer& multiplicity);

public:
    void insert_all_gens();
    void make_unimodular();    
    vector<pair<vector<key_t>, Integer> > getKeysAndMult() const;    
    ConeCollection(Cone<Integer> C, bool from_triangulation);
    
};

} // namespace

#endif /* LIBNORMALIZ_COLLECTION_H_ */
