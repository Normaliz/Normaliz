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

#include "libnormaliz/descent.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/my_omp.h"
#include "libnormaliz/sublattice_representation.h"

namespace libnormaliz {

template<typename Integer>
DescentFace<Integer>::DescentFace(){

    simplicial=false;
    coeff=0;
    tree_size=0;
}

template<typename Integer>
DescentFace<Integer>::DescentFace( const size_t dim_given, const boost::dynamic_bitset<>& facets_given){
    
    dim=dim_given;
    simplicial=false;
    own_facets=facets_given;
    tree_size=0;
    coeff=0;
}

template<typename Integer>
DescentSystem<Integer>::DescentSystem(const Matrix<Integer>& Gens_given, const Matrix<Integer>& SuppHyps_given, const vector<Integer>& Grading_given){

    descent_steps=0;
    nr_simplicial=0;
    
    Gens=Gens_given;
    SuppHyps=SuppHyps_given;
    Grading=Grading_given;
    
    nr_gens=Gens.nr_of_rows();
    nr_supphyps=SuppHyps.nr_of_rows();    
    dim=Gens.nr_of_columns();
    
    GradGens.resize(nr_gens);
    for(size_t i=0;i<nr_gens;++i)
        GradGens[i]=v_scalar_product(Grading,Gens[i]);
    
    multiplicity=0;

    SuppHypInd.resize(nr_supphyps);
    for(size_t i=0;i<nr_supphyps;++i){
        SuppHypInd[i].resize(nr_gens);
        for(size_t j=0;j<nr_gens;++j)
            if(v_scalar_product(SuppHyps[i],Gens[j])==0)
                SuppHypInd[i][j]=true;        
    }
}

template<typename Integer>
void  DescentFace<Integer>::compute(DescentSystem<Integer>& FF){
    
    size_t nr_supphyps=FF.nr_supphyps;
    size_t nr_gens=FF.nr_gens;
    
    size_t d=dim;    
    
    boost::dynamic_bitset<> GensInd(nr_gens);
    GensInd.set();    
    for(size_t i=0;i<nr_supphyps;++i){ // find Gens in this
        if(own_facets[i]==true) 
            GensInd= GensInd & FF.SuppHypInd[i];            
    }
    
    vector<libnormaliz::key_t> mother_key; // contains indices of Gens of *this
    for(size_t i=0;i<nr_gens;++i)
        if(GensInd[i])
            mother_key.push_back(i);
        
    Matrix<Integer> Gens_this=FF.Gens.submatrix(mother_key);
    Sublattice_Representation<Integer> Sublatt_this(Gens_this,true); // must take the saturation
        
    if(mother_key.size()==dim){ // *this is simplicial{
        simplicial=true;
        Matrix<Integer> Embedded_Gens=Sublatt_this.to_sublattice(Gens_this);
        Integer det=Embedded_Gens.vol();
        mpz_class mpz_det=convertTo<mpz_class>(det);
        mpq_class multiplicity=mpz_det;
        for(size_t i=0;i<Gens_this.nr_of_rows();++i)
            multiplicity/=convertTo<mpz_class>(FF.GradGens[mother_key[i]]);
        // cout << "M " << multiplicity << " C " << coeff << endl;
        #pragma omp critical(ADD_MULT)
        FF.multiplicity+=multiplicity*coeff;
        // tree_size=1;
        FF.nr_simplicial++;
        return;

    }   
    
    // Now we find the facets of *this.

    boost::dynamic_bitset<> facet_ind(nr_gens); // lists Gens
    map<boost::dynamic_bitset<>, boost::dynamic_bitset<> > FacetInds;
    map<boost::dynamic_bitset<>, key_t > CutOutBy; // the facet citting it out
    // entry is (facet_ind,indicator(SuppHyps))

    for(size_t i=0;i<nr_supphyps;++i){
        if(own_facets[i]==true) // contains *this
            continue;
        
        // we can identify rhge facet(*this) uniquely only via the Gens in it
        vector<libnormaliz::key_t> facet_key;
        for(size_t k=0;k<mother_key.size();++k){
            if(FF.SuppHypInd[i][mother_key[k]]==true)
                facet_key.push_back(mother_key[k]);            
        }
        if(facet_key.size() < d-1) // can't be a facet(*this)
            continue;
        
        // now we emake facet_ind
        facet_ind.reset();
        for(size_t i=0;i<facet_key.size();++i)
            facet_ind[facet_key[i]]=true;
        
        // next we check whether we have the intersection already
        if(FacetInds.find(facet_ind)!=FacetInds.end()){ // already found, we need it only once
            FacetInds[facet_ind][i]=true;// but we must add SuppHyps[i] to the facets(C) containing the current facet(*this)
            continue;            
        }

        // Must check the dimension
        if(FF.Gens.rank_submatrix(facet_key)<d-1) // dimension is too small
            continue; 

        // now we have a new facet
        FacetInds[facet_ind]=own_facets;
        FacetInds[facet_ind][i]=true;  //plus the facet cutting out facet_ind
        CutOutBy[facet_ind]=i;
    }
    
    // At this point we know the facets of *this.
    // The map FacetKeys assigns the set of containing SuppHyps to the facet_ind(Gens).
    // The set of containing SuppHyps is a unique signature as well.
    
    // Now we want to find the generator with the lrast number opf opposite facets(*this)    
    vector<size_t> count_in_facets(nr_gens);
    for(size_t i=0;i<mother_key.size();++i){
        size_t k=mother_key[i];
        for(auto F=FacetInds.begin();F!=FacetInds.end();++F)
            if((F->first)[k]==true)
                count_in_facets[k]++;        
    }
        
    size_t m=count_in_facets[0]; // we must have at least one facet (actually 3, since dim 2 is simplicial)
    libnormaliz::key_t m_ind=0;
        
    for(size_t i=1;i<count_in_facets.size();++i)
        if(count_in_facets[i]>m){
            m=count_in_facets[i];
            m_ind=i;
        }
        
    selected_gen=m_ind; // this is the selected generator (minimal number of opposite facets)
    vector<Integer> embedded_selected_gen=Sublatt_this.to_sublattice(FF.Gens[m_ind]);
    
    // now we must find the facets opposite to thge selected generator
    
    auto G=FacetInds.begin();    
    for(;G!=FacetInds.end();++G){
       if((G->first)[m_ind]==false){ // is opposite
            opposite_facets.push_back(G->second);
            vector<Integer> embedded_supphyp=Sublatt_this.to_sublattice_dual(FF.SuppHyps[CutOutBy[G->first]]);
            Integer ht=v_scalar_product(embedded_selected_gen,embedded_supphyp);
            heights.push_back(ht);
       }       
    }
    // cout << "H " << heights;
}

template<typename Integer>
void DescentSystem<Integer>::compute(){
    
    const size_t ReportBound=400;
    
    boost::dynamic_bitset<> empty(nr_supphyps);
    DescentFace<Integer> top(dim,empty);
    OldFaces[empty]=top;
    OldFaces[empty].coeff=1;
    long d=(long) dim;
    
    while(!OldFaces.empty()){
        
        if(verbose)
            verboseOutput() << "Descent from dim " << d << ", size " << OldFaces.size() << endl;
        
        auto F=OldFaces.begin();
        
        size_t nr_F=OldFaces.size();
        size_t kkpos=0;
        bool skip_remaining=false;
        
        const long VERBOSE_STEPS = 50;
        long step_x_size = nr_F-VERBOSE_STEPS;
        size_t total=nr_F;
        
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif
        #pragma omp parallel for firstprivate(kkpos,F) schedule(dynamic)
        for(size_t kk=0;kk< nr_F;++kk){
            
            if(skip_remaining)
                continue;
            
            if(verbose && nr_F>=ReportBound){
                #pragma omp critical(VERBOSE)
                while ((long)(kk*VERBOSE_STEPS) >= step_x_size) {
                    step_x_size += total;
                    verboseOutput() << "." <<flush;
                }
            }
            
#ifndef NCATCH
        try {
#endif
 
            INTERRUPT_COMPUTATION_BY_EXCEPTION
            
            for(;kk > kkpos; kkpos++, F++) ;
            for(;kk < kkpos; kkpos--, F--) ;
            
            F->second.compute(*this);
            if(F->second.simplicial)
                continue;
            
            auto G=(F->second).opposite_facets.begin();
            mpz_class deg_mpz=convertTo<mpz_class>(GradGens[(F->second).selected_gen]);
            mpq_class divided_coeff=(F->second).coeff/deg_mpz;
            #pragma omp critical(INSERT)
            {
            size_t j=0;
            for(;G!=(F->second).opposite_facets.end();++G){
                if(NewFaces.find(*G)==NewFaces.end())
                    NewFaces[*G]=DescentFace<Integer>(d-1,*G);
               NewFaces[*G].coeff+=divided_coeff*convertTo<mpz_class>((F->second).heights[j]);
               ++j;
            }
            }
            
#ifndef NCATCH
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
                #pragma omp flush(skip_remaining)
            }
#endif
        } // parallel for kk
        
        if(verbose && nr_F>=ReportBound)
            verboseOutput() << endl;
        
#ifndef NCATCH
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
        
        OldFaces=NewFaces;
        NewFaces.clear();
        d--;
        
    }    // while        
        
}
    


template<typename Integer>
bool DescentSystem<Integer>::set_verbose(bool onoff){
    bool old_verbose=verbose;
    verbose=onoff;
    return old_verbose;
}

template<typename Integer>
mpq_class DescentSystem<Integer>::getMultiplicity(){
    cout << "mult " << multiplicity << endl;
    return multiplicity;
}


template class DescentFace<long>;
template class DescentFace<long long>;
template class DescentFace<mpz_class>;

template class DescentSystem<long>;
template class DescentSystem<long long>;
template class DescentSystem<mpz_class>;

} // namespace

