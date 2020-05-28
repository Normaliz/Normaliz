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
ConeCollection<Integer>::ConeCollection(){
    is_initialized = false;    
}

template <typename Integer>
void ConeCollection<Integer>::initialize_minicones(const vector<pair<vector<key_t>, Integer> >& Triangulation){
    
    is_fan = true;
    is_triangulation = true;
    
    vector<key_t> GKeys;    
    Members.resize(1);
     
    for(auto& S: Triangulation){
        add_minicone(0,0,S.first, S.second);
        for(auto& g:S.first){
            assert(g < Generators.nr_of_rows());
            AllRays.insert(Generators[g]);
        }
    }
    
    is_initialized = true;
}

template <typename Integer>
void ConeCollection<Integer>::set_up(const Matrix<Integer>& Gens, const vector<pair<vector<key_t>, Integer> >& Triangulation){
    
    Generators = Gens;
    initialize_minicones(Triangulation);
}

/*
template <typename Integer>
ConeCollection<Integer>::ConeCollection(Cone<Integer> C, bool from_triangulation){
    
    assert(from_triangulation); // nfor the time being
    
    is_fan = true;
    is_triangulation = true;
    
    C.compute(ConeProperty::Generators, ConeProperty::Triangulation);
    
    Generators = C.BasisChangePointed.to_sublattice(C.Generators);
    
    vector<key_t> GKeys;
    
    Members.resize(1);
     
    for(auto& S: C.getTriangulation()){
        add_minicone(0,0,S.first, S.second);
    }
    
    //print();    
    // cout << "Ende Konstruktor D0D0D --- " << endl;
    
    verbose = C.verbose;
    
    is_initialized = true;
}*/



template <typename Integer>
void ConeCollection<Integer>::add_extra_generators(const Matrix<Integer>& NewGens){
    
    assert(is_initialized);
    
    if(verbose)
        verboseOutput() << "Inserting " << NewGens.nr_of_rows() << " new generators" << endl;
    
    for(size_t i = 0; i< NewGens.nr_of_rows(); ++i){
        if(AllRays.find(NewGens[i]) == AllRays.end()){
            Generators.append(NewGens[i]);
            refine(Generators. nr_of_rows()-1);
        }
    }    
}

template <typename Integer>
void ConeCollection<Integer>::add_minicone(const int level, const key_t mother, const vector<key_t>& GKeys, const Integer& multiplicity){
    
    MiniCone<Integer> MC(GKeys,multiplicity,*this);
    MC.is_simplex = is_triangulation;
    MC.level = level;
    // cout << "level " << level << " " << Members.size() << endl;
    MC.my_place = Members[level].size();
    Members[level].push_back(MC);
    if(level >0)
        Members[level-1][mother].Daughters.push_back(MC.my_place);
    /* for(auto& k:GKeys){
        AllRays.insert(Generators[k]);
    } */
    // print();
    return;   
}

template <typename Integer>
MiniCone<Integer>::MiniCone(const vector<key_t> GKeys, const Integer& mult, ConeCollection<Integer>& Coll){
    GenKeys = GKeys;
    multiplicity = mult;
    Collection = &(Coll);
    // dead = false;
}

template <typename Integer>
void ConeCollection<Integer>::print() const{
    cout << "================= Number of levels " << Members.size() << endl;
    for(size_t k=0; k< Members.size();++k){
        cout << "Level " << k << " Size " << Members[k].size()<< endl;
        cout << "-------------" << endl;
        for(size_t i=0; i< Members[k].size(); ++i)
            Members[k][i].print();
    }
    cout << "=======================================" << endl;
    
}

template <typename Integer>
void MiniCone<Integer>::print() const{

 cout<< "***** Mini " << level << " " << my_place << " Gens " << GenKeys;
 cout << "mult " << multiplicity << " daughters " << Daughters;
 cout << "----------------------" << endl;    
}

template <typename Integer>
void MiniCone<Integer>::refine(const key_t key){
    
    // cout << "################################### refining minocone " << level << " " << my_place << endl;
    
    bool has_daughters = (Daughters.size() > 0);

    if(SupportHyperplanes.nr_of_rows() == 0){
        Integer dummy;
        // cout << "************ " << GenKeys.size() << " " << Collection->Generators.nr_of_columns() << endl;
        Collection->Generators.simplex_data(GenKeys,SupportHyperplanes,dummy,false);
    }

    // cout << "SuppHyps " << endl;
    // SupportHyperplanes.pretty_print(cout);
    // cout << "-----------" << endl;
    // cout << "key " << key << " VVV " << Collection->Generators[key];

    vector<key_t> opposite_facets;

    for(size_t i=0; i<SupportHyperplanes.nr_of_rows(); ++i){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        Integer test = v_scalar_product(Collection->Generators[key],SupportHyperplanes[i]);
        if(test < 0)
            return;
        if(test == 0)
            continue;
        opposite_facets.push_back(i);
    }
    
    if(opposite_facets.size() ==1) // not contained in this minicone or extreme ray of it
        return;

    if(has_daughters){
        for(auto& d:Daughters){
            // cout << "Calling " << level + 1 << " " << d << endl;
            Collection->Members[level+1][d].refine(key);
        }
        return;
    }

    // cout << "opposite facets " << opposite_facets;

    for(size_t j=0; j<opposite_facets.size(); ++j){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        vector<key_t> NewGKey = GenKeys;
        NewGKey[opposite_facets[j]]= key;
        sort(NewGKey.begin(),NewGKey.end());
        Integer new_mult = Collection->Generators.submatrix(NewGKey).vol();
        // cout << "Mother " << my_place << endl;
        Collection->add_minicone(level+1,my_place,NewGKey,new_mult);
    }

    // cout << "ref " << Refinement.size() <<endl;

    //dead = true; // will be replaced by refinement

    return;
}

template <typename Integer>
void ConeCollection<Integer>::refine(const key_t key){
    
    if(AllRays.find(Generators[key]) != AllRays.end())
        return;
    
    // cout << "+++++++++++++++++++++++++++++++++++++++ Refine with vector " << key << endl;
    
    if(! Members[Members.size()-1].empty()){
        Members.resize(Members.size()+1);
        if(verbose)
            verboseOutput() << "Adding new level to tree structure" << endl;
    }

    for(size_t i=0; i< Members[0].size(); ++i){
        // cout << "RRRRRR " << i << " KKKK " << key << endl;
        Members[0][i].refine(key);
    }
    
    AllRays.insert(Generators[key]);
    
    /* list<MiniCone<Integer> > NewMinis;
    for(auto& T: Members){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        NewMinis.splice(NewMinis.end(), T.refine(key));        
    }
    AllRays.insert(Generators[key]);
    
    for(auto T = Members.begin(); T!=Members.end();){
        if(T->dead)
            T = Members.erase(T);
        else
            ++T;        
    }
    Members.splice(Members.end(), NewMinis);*/
}

template <typename Integer>
void ConeCollection<Integer>::insert_all_gens(){
    
    if(verbose)
        verboseOutput() << "Inserting " << Generators.nr_of_rows() << " given generators" << endl;
    
    size_t nr_inserted = 0;
    
    for(size_t i = 0; i < Generators.nr_of_rows(); ++i){
        // cout << "iiii " << i << endl;
                    nr_inserted++;
            
        if(verbose && nr_inserted % 100000 == 0)
            verboseOutput() << nr_inserted << " vectors inserted" << endl;
        
        refine(i);
    }
}

template <typename Integer>
void ConeCollection<Integer>::make_unimodular(){
    
    int omp_start_level = omp_get_level();
    
    while(true){        
        list< pair<vector<Integer>, pair<key_t,key_t> > > AllHilbs;
        vector< list< pair<vector<Integer>, pair<key_t,key_t> > > > Hilbs_thread(omp_get_max_threads());

        if(verbose){
                verboseOutput() << "Computing Hilbert bases of simplicial cones" << endl;
        }
        
        size_t nr_hilb_comp = 0;
        
        for(key_t k=0; k< Members.size(); ++k){
            
            bool skip_remaining = false;
            std::exception_ptr tmp_exception;            

#pragma omp parallel
            {
            int tn;
            if (omp_get_level() == omp_start_level)
                tn = 0;
            else
                tn = omp_get_ancestor_thread_num(omp_start_level + 1);
            
#pragma omp for
            for(key_t i = 0; i< Members[k].size(); ++i){
                
                if(skip_remaining)
                    continue;
                
                try{
            
                INTERRUPT_COMPUTATION_BY_EXCEPTION
                
                // cout << "Keys " << T.GenKeys;
                // cout << "mult " << T.multiplicity << endl;
                if(Members[k][i].multiplicity == 1) // already unimodular
                    continue;
                if(!Members[k][i].Daughters.empty()) // already subdivided
                    continue;
                Full_Cone<Integer> FC(Generators.submatrix(Members[k][i].GenKeys));
                FC.do_Hilbert_basis = true;
                FC.compute();
                
#pragma omp atomic
                nr_hilb_comp++;
                
                if(verbose && nr_hilb_comp % 50000 == 0){
#pragma omp critical(VERBOSE)
                    verboseOutput() << nr_hilb_comp << " Hilbert bases, computed" << endl;                    
                }
 
                // remove extreme rays -- can perhaps be done more efiiciently
                for(auto H = FC.Hilbert_Basis.begin(); H != FC.Hilbert_Basis.end();){
                    if(AllRays.find(*H) != AllRays.end())
                        H = FC.Hilbert_Basis.erase(H);
                    else
                        ++H;
                }              
                
                for(auto H = FC.Hilbert_Basis.begin(); H != FC.Hilbert_Basis.end();++H){
                        Hilbs_thread[tn].push_back(make_pair(*H, make_pair(k,i)));
                }
                
                } catch (const std::exception&) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
#pragma omp flush(skip_remaining)
                }
            } // i
        }// parallel
        
        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);
        
        } // k
        
        for(int i=0; i< omp_get_max_threads(); ++i)
            AllHilbs.splice(AllHilbs.end(), Hilbs_thread[i]);
        
        // cout << "AllHilbs " << endl;
        // for(auto& H: AllHilbs)
        //    cout << H;
        
        if(AllHilbs.empty())
            break;
        
        AllHilbs.sort();
        
        if(verbose)
            verboseOutput() << "Inserting " << AllHilbs.size() << " Hilbert bais elements of  simplices" << endl; 
        
        size_t nr_inserted = 0;
        
        vector<Integer> last_inserted;
        key_t key = Generators.nr_of_rows(); // to make gcc happy
        for(auto& H: AllHilbs){
            
            INTERRUPT_COMPUTATION_BY_EXCEPTION
                            
            if(! Members[Members.size()-1].empty()){
                Members.resize(Members.size()+1);
                if(verbose)
                    verboseOutput() << "Adding new level to tree structure" << endl;
            }
            if(H.first != last_inserted){
                last_inserted = H.first;
                key = Generators.nr_of_rows();
                Generators.append(H.first);
            }
            Members[H.second.first][H.second.second].refine(key);
            nr_inserted++;
            
            if(verbose && nr_inserted % 100000 == 0)
                verboseOutput() << nr_inserted << " vectors inserted" << endl;
        }
        
        for(auto& H: AllHilbs){
            AllRays.insert(H.first);
        }
        //cout << "Ende Durchgang " << endl;
        // for(auto& M:Members)
        // cout << "MMMM " << M.multiplicity << M.GenKeys;
    }
    // cout << "Ende unimod ------------------------------- " << endl;
    // for(auto& M:Members)
    //    cout << "MMMM " << M.multiplicity << M.GenKeys;
}
template <typename Integer>
void ConeCollection<Integer>::flatten(){
    
    // print();
    size_t tree_depth = 0;
    for(size_t k = 0; k< Members.size(); ++k){
        if(Members[k].size() >0)
            tree_depth++;
        for(key_t i = 0; i< Members[k].size(); ++i){
            // cout << "Out " << k << " " << i << " " << Members[k][i].Daughters.size() << endl;
            if(Members[k][i].Daughters.size() == 0)
                KeysAndMult.push_back(make_pair(Members[k][i].GenKeys, Members[k][i].multiplicity));
        }
    }
    if(verbose)
        verboseOutput() << "Tree depth " << tree_depth << ", Number of subcones " << KeysAndMult.size() 
                    << ", Number of generetors " << Generators.nr_of_rows() << endl;    
}

template <typename Integer>
const vector<pair<vector<key_t>, Integer> >&  ConeCollection<Integer>::getKeysAndMult() const{

    return KeysAndMult;    
}

template <typename Integer>
const Matrix<Integer>&  ConeCollection<Integer>::getGenerators() const{
    
    /*Matrix<Integer> Copy = Generators;
    Copy.remove_duplicate_and_zero_rows();
    cout << "Gen " << Generators.nr_of_rows() << " Copy " << Copy.nr_of_rows() << endl;*/
    return Generators;
}
    


} // namespace
