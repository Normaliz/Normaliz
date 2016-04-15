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

//---------------------------------------------------------------------------

#include "libnormaliz/matrix.h"
#include "libnormaliz/nmz_nauty.h"

namespace libnormaliz {
using namespace std;

template<typename Integer>
vector<vector<key_t> > Automorphism_Group<Integer>::getGenPerms() const{
    return GenPerms;
}

template<typename Integer>
vector<vector<key_t> > Automorphism_Group<Integer>::getLinFormPerms() const{
    return LinFormPerms;
}

template<typename Integer>
vector<vector<key_t> > Automorphism_Group<Integer>::getGenOrbits() const{
    return GenOrbits;
}

template<typename Integer>
vector<vector<key_t> > Automorphism_Group<Integer>::getLinFormOrbits() const{
    return GenPerms;
}

template<typename Integer>
vector<Matrix<Integer> > Automorphism_Group<Integer>::getLinMaps() const{
    return GenPerms;
}

template<typename Integer>
bool Automorphism_Group<Integer>::isFromAmbeientSpace() const{
    return from_ambient_space;
}

template<typename Integer>
bool Automorphism_Group<Integer>::isLinMapsComputed() const{
    return LinMaps_computed;
}

template<typename Integer>
bool Automorphism_Group<Integer>::isGraded() const{
    return graded;
}

template<typename Integer>
bool Automorphism_Group<Integer>::isInhomogeneous() const{
    return inhomogeneous;
}

template<typename Integer>
void Automorphism_Group<Integer>::setFromAmbeientSpace(bool on_off){
    from_ambient_space=on_off;
}

template<typename Integer>
void Automorphism_Group<Integer>::setGraded(bool on_off){
    graded=on_off;
}

template<typename Integer>
void Automorphism_Group<Integer>::setInhomogeneous(bool on_off){
    inhomogeneous=on_off;
}

template<typename Integer>
void Automorphism_Group<Integer>::makeLinMaps(){
}

template<typename Integer>
Automorphism_Group<Integer>::Automorphism_Group(){
    LinMaps_computed=false;
    from_ambient_space=false;
    graded=false;
    inhomogeneous=false;
}

template<typename Integer>
void Automorphism_Group<Integer>::compute(const Matrix<Integer>& GivenGens,const Matrix<Integer>& GivenLinForms){

    Gens=GivenGens;
    LinForms=GivenLinForms;
    vector<vector<long> > result=compute_automs(Gens,LinForms,order);
    size_t nr_automs=result.size()/2-2;
    for(size_t i=0;i<nr_automs;++i){
        vector<key_t> dummy(result[0].size());
        for(size_t j=0;j<dummy.size();++j)
            dummy[j]=result[i][j];
        GenPerms.push_back(dummy);
        vector<key_t> dummy_too(result[nr_automs].size());
        for(size_t j=0;j<dummy_too.size();++j)
            dummy_too[j]=result[i][j];
        LinFormPerms.push_back(dummy_too);
    }
    GenOrbits=convert_to_orbits(result[result.size()-2]);
    // cout << GenOrbits;
    LinFormOrbits=convert_to_orbits(result[result.size()-1]);
    // cout << LinFormOrbits;
}

vector<vector<key_t> > convert_to_orbits(const vector<long>& raw_orbits){
    
    vector<key_t> key(raw_orbits.size());
    vector<vector<key_t> > orbits;
    for(key_t i=0;i<raw_orbits.size();++i){
        if(raw_orbits[i]==(long) i){
            orbits.push_back(vector<key_t>(1,i));
            key[i]=orbits.size()-1;
        }
        else{
            orbits[key[raw_orbits[i]]].push_back(i);
        }
    }
    return orbits;    
}
    
template<typename Integer>
vector<vector<long> > compute_automs(const Matrix<Integer>& Gens, const Matrix<Integer>& LinForms, mpz_class& group_order){

    vector<vector<long> > Automs=compute_automs_by_nauty(Gens.get_elements(), LinForms.get_elements(), group_order);
    return Automs;
}

} // namespace
