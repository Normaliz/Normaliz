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

#include <stdlib.h>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

#include "libnormaliz/cone.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/map_operations.h"
#include "libnormaliz/convert.h"
#include "libnormaliz/my_omp.h"
#include "libnormaliz/collection.h"

namespace libnormaliz {
using namespace std;

template <typename Integer>
ConeCollection<Integer>::ConeCollection(const Cone<Integer> C, bool from_triangulation){
    
    assert(from_triangulation); // nfor the time being
    
    is_fan = true;
    is_triangulation = true;
    
    C.compute(ConeProperty::Generators, ConeProperty::Triangulation);
    
    C.BasisChangePointed.to_sublattice(Generators,C.Generators);
    
    for(auto& S: C.getTriangulation()){
        add_minicone(S.first, S.second);
    }    
}

template <typename Integer>
void ConeCollection<Integer>::add_minicone(const vector<key_t> GKeys, const Integer& multiplicity){
    
    MiniCone<Integer> MC(GKeys,multiplicity);
    MC.is_simplex = is_triangulation;
    Members.push_back(MC);
    for(auto& k:GKeys)
        AllRays.insert(GKeys[k]);
    
    
    
}

template <typename Integer>
MiniCone<Integer>::MiniCone(const vector<key_t> GKeys, const Integer& mult){
    GensKey = keys;
    multiplicity = mult;
    dead = false;    
}

template <typename Integer>
list<MiniCone<Integer> > MiniCone<Integer>::refine(const key_t key){
    
    if(SupportHyperplanes.nr_of_rows() == 0){
        Integer dummy;
        Collection->Generators(GensKey,SupportHyperplanes,dummy,false);
    }
    
    list<MiniCone<Integer> > Refinement;
    vector<key_t> opposite_facets;
    
   for(size_t i=0; i<SupportHyperplanes.nr_of_rows(); ++i){
       Integer test = v_scalar_product(Collection->Generators[key],SupportHyperplanes[i]);
        if(test < 0)
            return Refinement;
        if(test = 0)
            continue;
        opposite_facets.push_back(i);
   }
   
   if(opposite_facets.size() ==1) // not contained in this minicone or extreme ray of it
       return Refinement;
    
    for(size_t j=0; j<opposite_facets.size(); ++j){
        vector<key_t> NewGKey = GensKey;
        NewGKey[opposite_facets[j]] = key;
        sort(NewGKey.begin(),NewGKey.end());
        Integer new_mult = Collection->Generators.submatrix(NewGKey).vol();
        Refinement.push_back(MiniCone(NewGKey,new_mult));
    }
    
    dead = true; // will be replaced by refinement
    
    return Refinement;
}

template <typename Integer>
void ConeCollection<Integer>::refine(const key_t key){
    
    list<MiniCone<Integer> > NewMinis;
    for(auto& T: Members){
        T.refine(key).splice(NewMinis.begin(), NewMinis);        
    }
    AllRays.insert(Generators[key]);
    
    for(auto T = Members.bgin(); T!=Members.end();){
        if(T.dead)
            T = Members.erase(T);
        else
            ++T;        
    }
    Members.splice(Members.end(), NewMinis);
}

template <typename Integer>
void ConeCollection<Integer>::insert_all_gens(){
    
    for(size_t i = 0; i < Generators.size(); ++i)
        refine(i);    
}

template <typename Integer>
void ConeCollection<Integer>::make_unimodular(){
    
    while(true){        
        list<vector<Integer> > AllHilbs;
        
        for(auto& T: Members){
            if(T.multiplicity == 1)
                continue;
            Full_Cone<Integer> FC(Generators.submatrix(T.GenKeys));
            FC.do_hilbert_basis = true;
            FC.compute();
            AllHilbs.splice(AllHilbs.end(), FC.Hilbert_Basis);        
        }
        
        if(AllHilbs.empty())
            break;
        
        AllHilbs.sort();
        AllHilbs.unique();
        
        for(auto H = AllHilbs.begin(); H != AllHilbs.end();){
            if(AllRays.find(H) != AllRays.end())
                H = AllHilbs.erase(H);
            else
                ++H;           
        }
        
        if(AllHilbs.empty())
            break;
        
        key_t key = Generators.size();
        for(auto& H: AllHilbs){
            Generators.append(H);
            refine(key);            
        }    
    }
}

template <typename Integer>
vector<pair<vector<key_t>, Integer> >  ConeCollection<Integer>::getKeysAndMult() const{
    
    vector<pair<vector<key_t>, Integer> > KeysAndMult;
    for(auto& T: Members)
        KeysAndMult.push_back(make_pair(T.GensKey, T.multiplicity));
    return KeysAndMult;    
}


} // namespace
