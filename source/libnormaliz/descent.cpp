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

    // facets_computed=false;
    // multiplicity_computed=false;
    simplicial=false;
    tree_size=0;
}

template<typename Integer>
DescentFace<Integer>::DescentFace( const size_t dim_given, const boost::dynamic_bitset<>& facets_given){
    
    dim=dim_given;
    // facets_computed=false;
    // multiplicity_computed=false;
    simplicial=false;
    own_facets=facets_given;
    tree_size=0;
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
vector<boost::dynamic_bitset<> >&  DescentFace<Integer>::compute(DescentSystem<Integer>& FF){
    
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
        
    if(mother_key.size()==dim){ // *this is simplicial{
        simplicial=true;
        return opposite_facets;
    }   
    
    // Now we find the facets of *this.

    boost::dynamic_bitset<> facet_ind(nr_gens); // lists Gens
    map<boost::dynamic_bitset<>, boost::dynamic_bitset<> > FacetInds; 
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
    
    // now we must find the facets opposite to thge selected generator
    
    auto G=FacetInds.begin();    
    for(;G!=FacetInds.end();++G){
       if((G->first)[m_ind]==false)
            opposite_facets.push_back(G->second);   
    }
    
    return opposite_facets;
}

template<typename Integer>
void DescentSystem<Integer>::build(){
    
    const size_t ReportBound=400;
    
    Faces.resize(dim+1);
    boost::dynamic_bitset<> empty(nr_supphyps);
    DescentFace<Integer> top(dim,empty);
    Faces[dim][empty]=(top);
    for(long d=(long) dim;d>=1;--d){
        if(verbose)
            verboseOutput() << "Start descent from dim " << d << ", size " << Faces[d].size() << endl;
        
        auto F=Faces[d].begin();
        
        size_t nr_F=Faces[d].size();
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
            
            vector<boost::dynamic_bitset<> >& facets_of_F=F->second.compute(*this);
            auto G=facets_of_F.begin();
            #pragma omp critical(INSERT)
            {
            for(;G!=facets_of_F.end();++G)
                if(Faces[d-1].find(*G)==Faces[d-1].end())
                    Faces[d-1][*G]=DescentFace<Integer>(d-1,*G);
            }
            
#ifndef NCATCH
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
                #pragma omp flush(skip_remaining)
            }
#endif
        }
        
        if(verbose && nr_F>=ReportBound)
            verboseOutput() << endl;
        
#ifndef NCATCH
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
        
    }    
}

template<typename Integer>
void DescentFace<Integer>::compute_multiplicity(DescentSystem<Integer>& FF){
    
    size_t nr_supphyps=FF.nr_supphyps;
    size_t nr_gens=FF.nr_gens;
    
    // reconstruct the generators of *this   
    boost::dynamic_bitset<> gens_in_G(nr_gens);
    gens_in_G.set();
    for(size_t i=0;i<nr_supphyps;++i)
        if(own_facets[i]==true)
            gens_in_G &= FF.SuppHypInd[i];
        
    vector<libnormaliz::key_t> GensKey;
    for(size_t i=0;i<nr_gens;++i)
        if(gens_in_G[i]==true)
            GensKey.push_back(i);
        
    Matrix<Integer> Gens_G=FF.Gens.submatrix(GensKey);
    Sublattice_Representation<Integer> Sublatt_G(Gens_G,true); // must take the saturation
    
    if(simplicial){
        Matrix<Integer> Embedded_Gens=Sublatt_G.to_sublattice(Gens_G);
        Integer det=Embedded_Gens.vol();
        mpz_class mpz_det=convertTo<mpz_class>(det);
        multiplicity=mpz_det;
        for(size_t i=0;i<Gens_G.nr_of_rows();++i)
            multiplicity/=convertTo<mpz_class>(FF.GradGens[GensKey[i]]);
        tree_size=1;
        FF.nr_simplicial++;
        return;
    }
    
    
    vector<Integer> embedded_selected_gen=Sublatt_G.to_sublattice(FF.Gens[selected_gen]);
    mpz_class deg_mpz=convertTo<mpz_class>(FF.GradGens[selected_gen]);   

    multiplicity=0;
    
    // now we go over the facets opposite to the selected generator
    
    for(size_t i=0;i<opposite_facets.size();++i){
        // find SuppHyp defining this opposite facet
        size_t j;
        for(j=0;j<nr_supphyps;++j)
            if(own_facets[j]==false && opposite_facets[i][j]==true)
                break; // it is SuppHyps[j]

        vector<Integer> embedded_supphyp=Sublatt_G.to_sublattice_dual(FF.SuppHyps[j]);
        Integer ht=v_scalar_product(embedded_selected_gen,embedded_supphyp);

        mpq_class mult_facet=FF.Faces[dim-1][opposite_facets[i]].multiplicity;
        
        multiplicity+=mult_facet*convertTo<mpz_class>(ht)/deg_mpz;
        tree_size+=FF.Faces[dim-1][opposite_facets[i]].tree_size;
    }
    
    FF.descent_steps+=opposite_facets.size();
    
    if(dim==FF.dim){
        FF.multiplicity=multiplicity;
        if(verbose){
            verboseOutput() << "Mult " << multiplicity << endl;
            verboseOutput() << "Mult (float) " << mpq_to_nmz_float(multiplicity) << endl;
            verboseOutput() << "Full tree size " << tree_size << endl;
            verboseOutput() << "Number of descent steps " << FF.descent_steps << endl;
            verboseOutput() << "Number of simplicial Faces " << FF.nr_simplicial << endl;
            size_t total=0;
            for(size_t d=0;d<=FF.dim;++d)
                total+=FF.Faces[d].size();
            verboseOutput() << "Total number of faces " << total << endl;
        }
    }
}

template<typename Integer>
void DescentSystem<Integer>::compute_multiplicities(){
    
    const size_t ReportBound=1000;

    for(size_t d=0;d<=dim;++d){
        if(verbose)
            verboseOutput() << "compute multiplicity in dim " << d << ", size " << Faces[d].size() << endl;
        
        auto F=Faces[d].begin();
        size_t nr_F=Faces[d].size();
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
            
            if (skip_remaining) continue;
            
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
            
            F->second.compute_multiplicity(*this);
            
#ifndef NCATCH
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
                #pragma omp flush(skip_remaining)
            }
#endif
            
        }
        
        if(verbose && nr_F>=ReportBound)
            verboseOutput() << endl;
        
#ifndef NCATCH
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
    }    
}

template<typename Integer>
bool DescentSystem<Integer>::set_verbose(bool onoff){
    bool old_verbose=verbose;
    verbose=onoff;
    return old_verbose;
}

template<typename Integer>
mpq_class DescentSystem<Integer>::getMultiplicity(){
    return multiplicity;
}

template<typename Integer>
void DescentSystem<Integer>::compute(){
    build();
    compute_multiplicities();
}

template class DescentFace<long>;
template class DescentFace<long long>;
template class DescentFace<mpz_class>;

template class DescentSystem<long>;
template class DescentSystem<long long>;
template class DescentSystem<mpz_class>;

} // namespace

