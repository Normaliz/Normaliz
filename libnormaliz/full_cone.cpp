/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 */

//---------------------------------------------------------------------------

#include <stdlib.h>
#include <set>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <time.h>
#include <deque>

#include "full_cone.h"
#include "vector_operations.h"
#include "lineare_transformation.h"
#include "list_operations.h"
#include "my_omp.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//private
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::add_hyperplane(const size_t& new_generator, const FACETDATA & positive,const FACETDATA & negative){
    size_t k;
    
    // NEW: indgen is the index of the generator being inserted    
    
    // vector<Integer> hyperplane(hyp_size,0); // initialized with 0
    
    FACETDATA NewFacet; NewFacet.Hyp.resize(dim); NewFacet.GenInHyp.resize(nr_gen);
    
    Integer used_for_tests;
    if (test_arithmetic_overflow) {  // does arithmetic tests
        for (k = 0; k <dim; k++) {
            NewFacet.Hyp[k]=positive.ValNewGen*negative.Hyp[k]-negative.ValNewGen*positive.Hyp[k];
            used_for_tests =(positive.ValNewGen%overflow_test_modulus)*(negative.Hyp[k]%overflow_test_modulus)-(negative.ValNewGen%overflow_test_modulus)*(positive.Hyp[k]%overflow_test_modulus);
            if (((NewFacet.Hyp[k]-used_for_tests) % overflow_test_modulus)!=0) {
                errorOutput()<<"Arithmetic failure in Full_Cone::add_hyperplane. Possible arithmetic overflow.\n";
                throw ArithmeticException();
            }
        }
    }
    else  {                      // no arithmetic tests
        for (k = 0; k <dim; k++) {
            NewFacet.Hyp[k]=positive.ValNewGen*negative.Hyp[k]-negative.ValNewGen*positive.Hyp[k];
        }
    }
    v_make_prime(NewFacet.Hyp);
    NewFacet.ValNewGen=0; // not really needed, only for completeness
    
    NewFacet.GenInHyp=positive.GenInHyp & negative.GenInHyp; // new hyperplane contains old gen iff both pos and neg do
    NewFacet.GenInHyp.set(new_generator);  // new hyperplane contains new generator
    
    if(parallel_in_pyramid){
        #pragma omp critical(HYPERPLANE)
        Facets.push_back(NewFacet);
    }
    else
       Facets.push_back(NewFacet);
}


//---------------------------------------------------------------------------


template<typename Integer>
void Full_Cone<Integer>::find_new_facets(const size_t& new_generator){

    //to see if possible to replace the function .end with constant iterator since push-back is performed.

    // NEW: new_generator is the index of the generator being inserted

    size_t i,k,nr_zero_i;
    size_t subfacet_dim=dim-2; // NEW dimension of subfacet
    size_t facet_dim=dim-1; // NEW dimension of facet
    
    const bool tv_verbose =  false; // true && !is_pyramid; //verbose && Support_Hyperplanes.size()>10000; //verbose in this method call
    
        
    // preparing the computations
    deque <FACETDATA*> Pos_Simp,Pos_Non_Simp;
    deque <FACETDATA*> Neg_Simp,Neg_Non_Simp;
    deque <FACETDATA*> Neutral_Simp, Neutral_Non_Simp;
    
    boost::dynamic_bitset<> Zero_Positive(nr_gen),Zero_Negative(nr_gen);

    bool simplex;
    
    if (tv_verbose) verboseOutput()<<"transform_values: create SZ,Z,PZ,P,NS,N"<<endl<<flush;
    size_t ipos=0;
    
    typename list<FACETDATA>::iterator ii = Facets.begin();
    
    size_t listsize=Facets.size();

    for (size_t kk=0; kk<listsize; ++kk) {
        for(;kk > ipos; ++ipos, ++ii) ;
        for(;kk < ipos; --ipos, --ii) ;
        
        simplex=true;
        nr_zero_i=0;
        for(size_t j=0;j<nr_gen;j++){
            if(ii->GenInHyp.test(j))
                nr_zero_i++;
            if(nr_zero_i>facet_dim){
                simplex=false;
                break;
            }
        }
        
        if(ii->ValNewGen>0)
            Zero_Positive|=ii->GenInHyp;
        else if(ii->ValNewGen<0)
            Zero_Negative|=ii->GenInHyp;       
            
        if (ii->ValNewGen==0) {
            ii->GenInHyp.set(new_generator);  // Must be set explicitly !!
            if (simplex) {
                Neutral_Simp.push_back(&(*ii));
            }   else {
                Neutral_Non_Simp.push_back(&(*ii));
            }
        }
        else if (ii->ValNewGen>0) {
            if (simplex) {
                Pos_Simp.push_back(&(*ii));
            } else {
                Pos_Non_Simp.push_back(&(*ii));
            }
        } 
        else if (ii->ValNewGen<0) {
            if (simplex) {
                Neg_Simp.push_back(&(*ii));
            } else {
                Neg_Non_Simp.push_back(&(*ii));
            }
        }
    }
    
    boost::dynamic_bitset<> Zero_PN(nr_gen);
    Zero_PN=Zero_Positive & Zero_Negative;
    
    size_t nr_PosSimp  = Pos_Simp.size();
    size_t nr_PosNonSimp = Pos_Non_Simp.size();
    size_t nr_NegSimp  = Neg_Simp.size();
    size_t nr_NegNonSimp = Neg_Non_Simp.size();
    size_t nr_NeuSimp  = Neutral_Simp.size();
    size_t nr_NeuNonSimp = Neutral_Non_Simp.size();
    
    bool ranktest=((nr_PosNonSimp+nr_NegNonSimp+nr_NeuNonSimp>dim*dim*dim/6));  

    if (tv_verbose) verboseOutput()<<"PS "<<nr_PosSimp<<" P "<<nr_PosNonSimp<<" NS "<<nr_NegSimp<<" N "<<nr_NegNonSimp<<" ZS "<<nr_NeuSimp<<" Z "<<nr_NeuNonSimp<<endl<<flush;

    if (tv_verbose) verboseOutput()<<"transform_values: create lst of pairs <subfacet,facet> of NS"<<endl<<flush;
    
    vector< list<pair < boost::dynamic_bitset<>, int> > > Neg_Subfacet_Multi(omp_get_max_threads()) ;

    boost::dynamic_bitset<> zero_i(nr_gen);
    boost::dynamic_bitset<> subfacet(nr_gen);

    #pragma omp parallel for firstprivate(zero_i,subfacet) private(k,nr_zero_i) schedule(dynamic)
    for (i=0; i<nr_NegSimp;i++){
        zero_i=Zero_PN & Neg_Simp[i]->GenInHyp;
        
        nr_zero_i=0;
        for(size_t j=0;j<nr_gen;j++){
            if(zero_i.test(j))
                nr_zero_i++;
            if(nr_zero_i>subfacet_dim){
                break;
            }
        }

        if(nr_zero_i==subfacet_dim) // NEW This case treated separately
            Neg_Subfacet_Multi[omp_get_thread_num()].push_back(pair <boost::dynamic_bitset<>, int> (zero_i,i));
            
        else{       
            for (k =0; k<nr_gen; k++) {  
                if(zero_i.test(k)) {              
                    subfacet=zero_i;
                    subfacet.reset(k);  // remove k-th element from facet to obtain subfacet
                    Neg_Subfacet_Multi[omp_get_thread_num()].push_back(pair <boost::dynamic_bitset<>, int> (subfacet,i));
                }
            }
        }
    }
    
    list < pair < boost::dynamic_bitset<>, int> > Neg_Subfacet_Multi_United;
    for(int i=0;i<omp_get_max_threads();++i)
        Neg_Subfacet_Multi_United.splice(Neg_Subfacet_Multi_United.begin(),Neg_Subfacet_Multi[i]);       
    Neg_Subfacet_Multi_United.sort();


    if (tv_verbose) verboseOutput()<<"transform_values: discard double subfacets of NS, list size "<< Neg_Subfacet_Multi_United.size() <<endl<<flush;

    list< pair < boost::dynamic_bitset<>, int > >::iterator jj;
    list< pair < boost::dynamic_bitset<>, int > >::iterator del;
    jj =Neg_Subfacet_Multi_United.begin();                               // remove negative subfecets shared
    while (jj!= Neg_Subfacet_Multi_United.end()) {                       // by two neg simpl facets
        del=jj++;
        if (jj!=Neg_Subfacet_Multi_United.end() && (*jj).first==(*del).first) {   //delete since is the intersection of two negative simplicies
            Neg_Subfacet_Multi_United.erase(del);
            del=jj++;
            Neg_Subfacet_Multi_United.erase(del);
        }
    }

    size_t nr_NegSubfMult = Neg_Subfacet_Multi_United.size();
    if (tv_verbose) verboseOutput()<<"transform_values: after discarding double subfacets of NS list size "<<nr_NegSubfMult<<endl<<flush;
    
    map < boost::dynamic_bitset<>, int > Neg_Subfacet;
    size_t nr_NegSubf=0;
    
    #pragma omp parallel private(jj)
    {
    size_t i,j,k,nr_zero_i;
    boost::dynamic_bitset<> subfacet(dim-2);
    jj = Neg_Subfacet_Multi_United.begin();
    size_t jjpos=0;

    map < boost::dynamic_bitset<>, int > ::iterator last_inserted=Neg_Subfacet.begin(); // used to speedup insertion into the new map
    bool found;
    #pragma omp for schedule(dynamic)
    for (size_t j=0; j<nr_NegSubfMult; ++j) {             // remove negative subfacets shared
        for(;j > jjpos; ++jjpos, ++jj) ;                // by non-simpl neg or neutral facets 
        for(;j < jjpos; --jjpos, --jj) ;

        subfacet=(*jj).first;
        found=false; 
        for (i = 0; i <nr_NeuSimp; i++) {
            found=subfacet.is_subset_of(Neutral_Simp[i]->GenInHyp);
            if(found)
                break;
        }
        if (!found) {
            for (i = 0; i <nr_NeuNonSimp; i++) {
                found=subfacet.is_subset_of(Neutral_Non_Simp[i]->GenInHyp);
                if(found)
                    break;                    
                }
            if(!found) {
                for (i = 0; i <nr_NegNonSimp; i++) {
                    found=subfacet.is_subset_of(Neg_Non_Simp[i]->GenInHyp);
                    if(found)
                        break; 
                }
            }
        }
        if (!found) {
            if(parallel_in_pyramid){
                #pragma omp critical(NEGATIVE_SUBFACET)
                {last_inserted=Neg_Subfacet.insert(last_inserted,*jj);}
            } else {
                last_inserted=Neg_Subfacet.insert(last_inserted,*jj);
            }
        }
    }
    
    #pragma omp single
    {nr_NegSubf=Neg_Subfacet.size();}
    
    #pragma omp single
    if (tv_verbose) {
        verboseOutput()<<"transform_values: reduced map size "<< nr_NegSubf <<endl<<flush;
    } 
    #pragma omp single nowait
    {Neg_Subfacet_Multi_United.clear();}

    #pragma omp single
    if (tv_verbose) {
        verboseOutput()<<"transform_values: PS vs NS"<<endl<<flush;
    }
    
    boost::dynamic_bitset<> zero_i(nr_gen);
    map <boost::dynamic_bitset<>, int> ::iterator jj_map;

    #pragma omp for schedule(dynamic) nowait           // Now matching positive and negative (sub)facets
    for (i =0; i<nr_PosSimp; i++){ //Positive Simp vs.Negative Simp
        zero_i=Pos_Simp[i]->GenInHyp & Zero_PN;
        nr_zero_i=0;
        for(size_t m=0;m<nr_gen;m++){
            if(zero_i.test(m))
                nr_zero_i++;
            if(nr_zero_i>subfacet_dim){
                break;
            }
        }
        
        if (nr_zero_i==subfacet_dim) {                 // NEW slight change in logic. Positive simpl facet shared at most
            jj_map=Neg_Subfacet.find(zero_i);           // one subfacet with negative simpl facet
            if (jj_map!=Neg_Subfacet.end()) {
                add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Simp[(*jj_map).second]);
                (*jj_map).second = -1;  // block subfacet in further searches
            }
        }
        if (nr_zero_i==facet_dim){    // now there could be more such subfacets. We make all and search them.      
            for (k =0; k<nr_gen; k++) {  // BOOST ROUTINE
                if(zero_i.test(k)) { 
                    subfacet=zero_i;
                    subfacet.reset(k);  // remove k-th element from facet to obtain subfacet
                    jj_map=Neg_Subfacet.find(subfacet);
                    if (jj_map!=Neg_Subfacet.end()) {
                        add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Simp[(*jj_map).second]);
                        (*jj_map).second = -1;
                    }
                }
            }
        }
    }

    #pragma omp single
    if (tv_verbose) {
        verboseOutput()<<"transform_values: NS vs P"<<endl << flush;
    }


    jj_map = Neg_Subfacet.begin();
    jjpos=0;
    #pragma omp for schedule(dynamic) nowait
    for (size_t j=0; j<nr_NegSubf; ++j) {
        for( ; j > jjpos; ++jjpos, ++jj_map) ;
        for( ; j < jjpos; --jjpos, --jj_map) ;

        if ( (*jj_map).second != -1 ) {  // skip used subfacets
            for (i = 0; i <nr_PosNonSimp; i++) {
                if(jj_map->first.is_subset_of(Pos_Non_Simp[i]->GenInHyp)){
                    add_hyperplane(new_generator,*Pos_Non_Simp[i],*Neg_Simp[(*jj_map).second]);
                    break;
                }
            }
        }
    }
    
    #pragma omp single nowait
    if (tv_verbose) {
        verboseOutput()<<"transform_values: PS vs N"<<endl << flush;
    }

    vector<key_t> key(nr_gen);
    size_t nr_missing;
    bool common_subfacet;
    #pragma omp for schedule(dynamic) nowait
    for (size_t i =0; i<nr_PosSimp; i++){ //Positive Simp vs.Negative Non Simp
        zero_i=Zero_PN & Pos_Simp[i]->GenInHyp;
        nr_zero_i=0;
        for(j=0;j<nr_gen && nr_zero_i<=facet_dim;j++)
            if(zero_i.test(j)){
                key[nr_zero_i]=j;
                nr_zero_i++;
            }        
        if(nr_zero_i>=subfacet_dim) {
            for (j=0; j<nr_NegNonSimp; j++){ // search negative facet with common subfacet
                nr_missing=0; 
                common_subfacet=true;               
                for(k=0;k<nr_zero_i;k++) {
                    if(!Neg_Non_Simp[j]->GenInHyp.test(key[k])) {
                        nr_missing++;
                        if(nr_missing==2 || nr_zero_i==subfacet_dim) {
                            common_subfacet=false;
                            break;
                        }
                    }
                 }
                    
                 if(common_subfacet){                 
                    add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Non_Simp[j]);
                    if(nr_zero_i==subfacet_dim) // only one subfacet can lie in negative hyperplane
                        break;
                 }
            }           
        }            
    } // PS vs N



    #pragma omp single nowait
    if (tv_verbose) {
        verboseOutput()<<"transform_values: P vs N"<<endl<<flush;
    }
    
    list<FACETDATA*> AllNonSimpHyp;
    typename list<FACETDATA*>::iterator a;
    if(!ranktest){
        for(i=0;i<nr_PosNonSimp;++i)
            AllNonSimpHyp.push_back(&(*Pos_Non_Simp[i]));
        for(i=0;i<nr_NegNonSimp;++i)
            AllNonSimpHyp.push_back(&(*Neg_Non_Simp[i]));
        for(i=0;i<nr_NeuNonSimp;++i)
            AllNonSimpHyp.push_back(&(*Neutral_Non_Simp[i]));
    }  
    // cout << "ranktest " << ranktest << endl;   
    
    
    bool exactly_two;
    FACETDATA *hp_i, *hp_j, *hp_t; // pointers to current hyperplanes
    
    size_t missing_bound, nr_common_zero;
    boost::dynamic_bitset<> common_zero(nr_gen);
    vector<key_t> common_key(nr_gen);
    
    #pragma omp for schedule(dynamic) nowait
    for (size_t i =0; i<nr_PosNonSimp; i++){ //Positive Non Simp vs.Negative Non Simp

        hp_i=Pos_Non_Simp[i];
        zero_i=Zero_PN & hp_i->GenInHyp;
        nr_zero_i=0;
        for(j=0;j<nr_gen;j++)
            if(zero_i.test(j)){
                key[nr_zero_i]=j;
                nr_zero_i++;
            }
            
        if (nr_zero_i>=subfacet_dim) {
        
            missing_bound=nr_zero_i-subfacet_dim; // at most this number of generators can be missing
                                                  // to have a chance for common subfacet
                for (j=0; j<nr_NegNonSimp; j++){
                hp_j=Neg_Non_Simp[j];
                
                nr_missing=0; 
                nr_common_zero=0;
                common_subfacet=true;               
                for(k=0;k<nr_zero_i;k++) {
                    if(!hp_j->GenInHyp.test(key[k])) {
                        nr_missing++;
                        if(nr_missing>missing_bound || nr_zero_i==subfacet_dim) {
                            common_subfacet=false;
                            break;
                        }
                    }
                    else {
                        common_key[nr_common_zero]=key[k];
                        nr_common_zero++;
                    }
                 }
                 

                if(common_subfacet){//intersection of *i and *j may be a subfacet
                    common_zero=zero_i & hp_j->GenInHyp;
                    exactly_two=true;

                    if (ranktest) {
                        Matrix<Integer> Test(nr_common_zero,dim);
                        for (k = 0; k < nr_common_zero; k++)
                            Test.write(k,Generators.read(common_key[k])); 

                        if (Test.rank_destructive()<subfacet_dim) {
                            exactly_two=false;
                        }
                    } // ranktest
                    else{                 // now the comparison test
                        for (a=AllNonSimpHyp.begin();a!=AllNonSimpHyp.end();++a){
                            hp_t=*a;
                            if ((hp_t!=hp_i) && (hp_t!=hp_j) && common_zero.is_subset_of(hp_t->GenInHyp)) {                                
                                exactly_two=false;
                                AllNonSimpHyp.splice(AllNonSimpHyp.begin(),AllNonSimpHyp,a);
                                break;
                            }
                        }                       
                    } // else
                    if (exactly_two) {  //intersection of i and j is a subfacet
                        add_hyperplane(new_generator,*hp_i,*hp_j);
                    }
                }
            }
        }
    }
    } //END parallel

    //removing the negative hyperplanes
    // now done in build_cone

    if (tv_verbose) verboseOutput()<<"transform_values: done"<<endl;
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::extend_triangulation(const size_t& new_generator){

    size_t listsize = Facets.size();
    vector<typename list<FACETDATA>::iterator> visible;
    visible.reserve(listsize);
    typename list<FACETDATA>::iterator i = Facets.begin();

    // #pragma omp critical(VERBOSE)
    // verboseOutput() << "L " << pyr_level << " H " << listsize << " T " << TriangulationSize << endl << flush;
    
    for (; i!=Facets.end(); ++i) 
        if (i->ValNewGen < 0) // visible facet
            visible.push_back(i);

    listsize = visible.size();
    // cout << "Pyr Level " << pyr_level << " Visible " << listsize <<  " Triang " << TriangulationSize << endl;


    typename list<SHORTSIMPLEX>::iterator oldTriBack = --Triangulation.end();
    #pragma omp parallel private(i)
    {
    size_t k,l;
    bool one_not_in_i, not_in_facet;
    size_t not_in_i=0;
    size_t facet_dim=dim-1;
    size_t nr_in_i=0;
    list<SHORTSIMPLEX> Triangulation_kk;
    
    typename list<SHORTSIMPLEX>::iterator j;
    
    vector<key_t> key(dim);
    
    #pragma omp for schedule(dynamic)
    for (size_t kk=0; kk<listsize; ++kk) {
        i=visible[kk];
        
        nr_in_i=0;
        for(size_t m=0;m<nr_gen;m++){
            if(i->GenInHyp.test(m))
                nr_in_i++;
            if(nr_in_i>facet_dim){
                break;
            }
        }

        if (nr_in_i==facet_dim){  // simplicial
            l=0;
            for (k = 0; k <nr_gen; k++) {
                if (i->GenInHyp[k]==1) {
                    key[l]=k;
                    l++;
                }
            }
            key[dim-1]=new_generator;
 
            if(parallel_in_pyramid) {
                #pragma omp critical(TRIANG) // critical only on top level
                store_key(key,-i->ValNewGen,Triangulation);
            } else {
                store_key(key,-i->ValNewGen,Triangulation);
            }
            continue;
        }
        
        size_t irrelevant_vertices=0;
        for(size_t vertex=0;vertex<VertInTri.size();++vertex){
        
            if(i->GenInHyp[VertInTri[vertex]]==0) // lead vertex not in hyperplane
                continue;
                
            if(irrelevant_vertices<dim-2){
                ++irrelevant_vertices;
                continue;
            }       
        
            j=TriSectionFirst[vertex];
            bool done=false;
            for(;!done;j++)
            {
              done=(j==TriSectionLast[vertex]);
              key=j->key;
              one_not_in_i=false;  // true indicates that one gen of simplex is not in hyperplane
              not_in_facet=false;  // true indicates that a second gen of simplex is not in hyperplane
              for(k=0;k<dim;k++){            
                 if ( !i->GenInHyp.test(key[k])) {
                     if(one_not_in_i){
                         not_in_facet=true;
                         break;
                     }
                     one_not_in_i=true;
                     not_in_i=k;
                  }
              }
              
              if(not_in_facet) // simplex does not share facet with hyperplane
                 continue;
              
              key[not_in_i]=new_generator;              
              store_key(key,-i->ValNewGen,Triangulation_kk);
                       
            } // j
            
        } // for vertex

        if(true) { // parallel_in_pyramid
            #pragma omp critical(TRIANG)
                Triangulation.splice(Triangulation.end(),Triangulation_kk);
            }
        else 
            Triangulation.splice(Triangulation.end(),Triangulation_kk);
    } // for kk

    } // parallel

    VertInTri.push_back(new_generator);
    TriSectionFirst.push_back(++oldTriBack);
    TriSectionLast.push_back(--Triangulation.end());    
    
    // cout << " Aus extend Tri " << omp_get_level() <<  endl;
} 

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::store_key(const vector<key_t>& key, const Integer& height,
                                list<SHORTSIMPLEX>& Triangulation){

    SHORTSIMPLEX newsimplex;
    newsimplex.key=key;
    newsimplex.height=height;
    
    #pragma omp atomic
    TriangulationSize++;
            
    
    if(keep_triangulation){
        Triangulation.push_back(newsimplex);
        return;  
    }
    
    bool Simpl_available=true;
    int tn;
    if(omp_get_level()==0)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(1);
    typename list<SHORTSIMPLEX>::iterator F;

    if(Top_Cone->FS[tn].empty()){
        #pragma omp critical(FREESIMPL)
        {
        if(Top_Cone->FreeSimpl.empty())
            Simpl_available=false;
        else{
            F=Top_Cone->FreeSimpl.begin();  // take 100 simplices from FreeSimpl
            size_t q; for(q=0;q<1000;++q){            // or what you can get
                if(F==Top_Cone->FreeSimpl.end())
                    break;
                ++F;
            }
        
            if(q<1000)
                Top_Cone->FS[tn].splice(Top_Cone->FS[tn].begin(),
                    Top_Cone->FreeSimpl);
            else
                Top_Cone->FS[tn].splice(Top_Cone->FS[tn].begin(),
                              Top_Cone->FreeSimpl,Top_Cone->FreeSimpl.begin(),++F);
        } // else
        } // critical
    } // if empty
          

    if(Simpl_available){
        Triangulation.splice(Triangulation.end(),Top_Cone->FS[tn],
                        Top_Cone->FS[tn].begin());
        Triangulation.back()=newsimplex;
    }
    else
        Triangulation.push_back(newsimplex);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::process_pyramids(const size_t new_generator,const bool recursive){
          
    // cout << "In PP" << endl;
    
    size_t start_level=omp_get_level(); // allows us to check that we are on level 0
                                        // outside the loop and can therefore call evaluation

    vector<key_t> Pyramid_key;
    Pyramid_key.reserve(nr_gen);
    boost::dynamic_bitset<> in_Pyramid(nr_gen); 
     

    typename list< FACETDATA >::iterator l;
    size_t i,lpos, listsize=Facets.size();
    //use deque here to have independent entries
    deque<bool> done(listsize,false);
    
    size_t nr_done=0;
    size_t store_level;
    if(recursion_allowed)
        store_level=0;
    else
        store_level=pyr_level+1;
        
    bool skip_remaining_tri,skip_remaining_pyr;

    do{    

    lpos=0;
    skip_remaining_tri=skip_remaining_pyr=false;
    l=Facets.begin();

    // if the next loop is de-parallelized, then
    // you MUST change "false" to "true" in process_pyramid, Pyramid.parallel_in_pyramid=false;
    // (or commnent it out)
    
    // #pragma omp parallel for private(i) firstprivate(lpos,l,Pyramid_key,in_Pyramid) schedule(dynamic) 
    for (size_t k=0; k<listsize; k++) {
    
        if(skip_remaining_tri || skip_remaining_pyr )
            continue;
            
        for(;k > lpos; lpos++, l++) ;
        for(;k < lpos; lpos--, l--) ;
                
        if(done[lpos])
            continue;
            
        done[lpos]=true;
        
        #pragma omp atomic 
        nr_done++;
            
        if(l->ValNewGen>=0 ||(!recursive && do_partial_triangulation && l->ValNewGen>=-1)){
            continue;
        }
            
        boost::dynamic_bitset<> in_Pyramid(nr_gen,false);           
        Pyramid_key.push_back(new_generator);
        in_Pyramid.set(new_generator);
        for(i=0;i<nr_gen;i++){
            if(in_triang[i] && v_scalar_product(l->Hyp,Generators.read(i))==0){ // from the second extension on the incidence data
                Pyramid_key.push_back(i);                                      // in Facets is no longer up-to-date
                in_Pyramid.set(i);
            }
        }
        
        
        process_pyramid(Pyramid_key, in_Pyramid, new_generator, recursive);
        Pyramid_key.clear();
        
        if(check_evaluation_buffer_size() && start_level==0  && nr_done < listsize){  // we interrupt parallel execution if it is really parallel
                                                           //  to keep the triangulation buffer under control
            skip_remaining_tri=true;
        }
        
        if(Top_Cone->nrPyramids[store_level] > 200000 && start_level==0  && nr_done < listsize){  // we interrupt parallel execution if it is really parallel
                                                           //  to keep the triangulation buffer under control
            skip_remaining_pyr=true;                      // CHOOSE SAME VALUE IN evaluate_stored_pyramids
        }
            
    } // end parallel for k
    
    if(skip_remaining_tri)
    {
        // verboseOutput() << nr_done << " of " << listsize << " pyramids done. " << omp_get_level() << endl;
        Top_Cone->evaluate_triangulation();
    }
    if(skip_remaining_pyr){
            if (verbose)
                verboseOutput() << ">>>>>>>>>>  descending to level " << store_level << endl;
                Top_Cone->evaluate_stored_pyramids(store_level);
    }    
    
    } while(skip_remaining_tri || skip_remaining_pyr);
    
    //verboseOutput() << nr_done << " of " << listsize << " pyramids done." << endl;
        
    if(check_evaluation_buffer()){
          Top_Cone->evaluate_triangulation();
    } 
    
    // cout << "Aus PP " << omp_get_level() << endl;

}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::process_pyramid(const vector<key_t> Pyramid_key, const boost::dynamic_bitset<> in_Pyramid, 
                          const size_t new_generator,const bool recursive){
    
    #pragma omp atomic
    Top_Cone->totalNrPyr++;
    
    list<vector<Integer> > NewFacets;
    
    if(Pyramid_key.size()==dim){  // simplicial pyramid done here
        #pragma omp atomic
        Top_Cone->nrSimplicialPyr++;
        Simplex<Integer> S(Pyramid_key, Generators);
        if(recursive){ // the facets may be facets of the mother cone and if recursive==true must be given back
            Matrix<Integer> H=S.read_support_hyperplanes();
            for (size_t i=0; i<dim;i++)
                NewFacets.push_back(H.read(i));
        }    
        if(do_triangulation || (do_partial_triangulation && S.read_volume()>1)){
            if(parallel_in_pyramid) {
                #pragma omp critical(TRIANG) // critical only on top level
                store_key(Pyramid_key,S.read_volume(),Triangulation);  
            } else {
                store_key(Pyramid_key,S.read_volume(),Triangulation);
            }  
        }
    }
    else {  // non-simplicial
        if(recursive){
            Full_Cone<Integer> Pyramid(*this,Pyramid_key);
            Pyramid.pyr_level=0; // value is in principle irrelevant
            Pyramid.do_all_hyperplanes=true;

/*          AT PRESENT WE DO NOT LIMIT THE RECURSION FOR PARTIAL TRIANGULATION HERE
A REASONABLE SOLUTION WOULD REQUIRE BOOKKEEPING OF THE NATURAL LEVEL OF THE PYRAMIDS
SEE ALSO evaluate_stored_pyramids
            
            // the next line has the effect that we fully triangulate the pyramid in order to avoid
            // recursion down to simplices in case of partial triangulation
            Pyramid.do_triangulation= (do_partial_triangulation && pyr_level >=1) || do_triangulation;
            
            // cout << "In py " << omp_get_level() << " " << omp_get_thread_num() << " tn " << Pyramid.thread_num << endl;
            
            if(Pyramid.do_triangulation)
                Pyramid.do_partial_triangulation=false;
*/

            Pyramid.build_cone(); // build and evaluate pyramid
            // Pyramid.parallel_in_pyramid=false;
            NewFacets.splice(NewFacets.begin(),Pyramid.Support_Hyperplanes);
       }
       else{ // if recursive==false we only store the pyramid
           vector<key_t> key_wrt_top(Pyramid_key.size());
           for(size_t i=0;i<Pyramid_key.size();i++)
                key_wrt_top[i]=Top_Key[Pyramid_key[i]];
           #pragma omp critical(STOREPYRAMIDS)
           {
           if(recursion_allowed){    
                Top_Cone->Pyramids[0].push_back(key_wrt_top); // if we come from top cone
                Top_Cone->nrPyramids[0]++;// Pyramids go level 0
           }
           else{                                                         
                Top_Cone->Pyramids[pyr_level+1].push_back(key_wrt_top);
                Top_Cone->nrPyramids[pyr_level+1]++;           
           }
           }
       }
    }   

    // now we give support hyperplanes back to the mother cone if only if
    // they are not done on the "mother" level if and only if recursive==true
    

    if(recursive){
        size_t i;                        
        typename list<vector<Integer> >::iterator pyr_hyp = NewFacets.begin();
        bool new_global_hyp;
        FACETDATA NewFacet;
        Integer test;   
        for(;pyr_hyp!= NewFacets.end();pyr_hyp++){
            if(v_scalar_product(Generators.read(new_generator),*pyr_hyp)>0)
                continue;
            new_global_hyp=true;
            for(i=0;i<nr_gen;i++){
                if(in_Pyramid.test(i) || !in_triang[i])
                    continue;
                test=v_scalar_product(Generators.read(i),*pyr_hyp);
                if(test<=0){
                    new_global_hyp=false;
                    break;
                }
            
            }
            if(new_global_hyp){
                NewFacet.Hyp=*pyr_hyp;                
                if(parallel_in_pyramid){
                    #pragma omp critical(HYPERPLANE) 
                    Facets.push_back(NewFacet);
                } else {
                    Facets.push_back(NewFacet);
                }
            }
        }
    }
    NewFacets.clear();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_and_evaluate_start_simplex(){

    size_t i,j;
    Integer factor;

    
    Simplex<Integer> S = find_start_simplex();
    vector<key_t> key=S.read_key();   // generators indexed from 0
        
    for (i = 0; i < dim; i++) {
        in_triang[key[i]]=true;
        if (ht1_triangulation && isComputed(ConeProperty::LinearForm))
            ht1_triangulation = (gen_degrees[i] == 1);
    }
       
    Matrix<Integer> H=S.read_support_hyperplanes();
    for (i = 0; i <dim; i++) {
        FACETDATA NewFacet; NewFacet.Hyp.resize(dim); NewFacet.GenInHyp.resize(nr_gen);
        NewFacet.Hyp=H.read(i);
        for(j=0;j < dim;j++)
            if(j!=i)
                NewFacet.GenInHyp.set(key[j]);
        NewFacet.ValNewGen=-1;         // must be taken negative since opposite facet
        Facets.push_back(NewFacet);    // was visible before adding this vertex
    }
    
    if(!is_pyramid){
        //define Order_Vector, decides which facets of the simplices are excluded
        Order_Vector = vector<Integer>(dim,0);
        Matrix<Integer> G=S.read_generators();
        //srand(12345);
        for(i=0;i<dim;i++){
            factor=(unsigned long)(2*(rand()%(2*dim))+3);
            for(j=0;j<dim;j++)
                Order_Vector[j]+=factor*G.read(i,j);        
        }
    }

    //the volume is an upper bound for the height
    if(do_triangulation || (do_partial_triangulation && S.read_volume()>1))
    {
        store_key(key,S.read_volume(),Triangulation);  // height understood positive
    }
    
    if(do_triangulation){ // we must prepare the sections of the triangulation
        for(i=0;i<dim;i++)
        {
            VertInTri.push_back(key[i]);
            TriSectionFirst.push_back(Triangulation.begin());
            TriSectionLast.push_back(Triangulation.begin());
        }
    }
    
}



//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_stored_pyramids(const size_t level){

    assert(omp_get_level()==0);

    if(Pyramids[level].empty())
        return;
    Pyramids.resize(level+2); // provide space for a new generation
    nrPyramids.resize(level+2);
    nrPyramids[level+1]=0;

    size_t nr_done=0;
    size_t nr_pyramids=nrPyramids[level];
    vector<short> Done(nr_pyramids,0);
    if (verbose) {
        verboseOutput() << "************************************************" << endl;
        verboseOutput() << "Evaluating " << nr_pyramids << " pyramids on level " << level << endl;
        verboseOutput() << "************************************************" << endl;
    }
    typename list<vector<key_t> >::iterator p;
    size_t ppos;
    bool skip_remaining_tri,skip_remaining_pyr;

    do
    {

       p = Pyramids[level].begin();
       ppos=0;
       skip_remaining_tri=false;
       skip_remaining_pyr=false;
    
       #pragma omp parallel for firstprivate(p,ppos) schedule(dynamic) 
       for(size_t i=0; i<nr_pyramids; i++){
       
           if(skip_remaining_tri || skip_remaining_pyr)
                continue;
                
           for(; i > ppos; ++ppos, ++p) ;
           for(; i < ppos; --ppos, --p) ;
           
           if(Done[i])
               continue;
           Done[i]=1;
           
           #pragma omp atomic
           nr_done++;
           
           Full_Cone<Integer> Pyramid(*this,*p);
           Pyramid.recursion_allowed=false; // ABSOLUTELY NECESSARY HERE
           Pyramid.pyr_level=level;
           Pyramid.do_all_hyperplanes=false;
           Pyramid.parallel_in_pyramid=false;
           if(level>=2 && do_partial_triangulation){ // limits the descent of do_partial_triangulation
               Pyramid.do_triangulation=true;
               Pyramid.do_partial_triangulation=false;
           }
           Pyramid.build_cone();
           if(check_evaluation_buffer_size() && nr_done < nr_pyramids)  // we interrupt parallel execution if it is really parallel
                skip_remaining_tri=true;                         //  to keep the triangulation buffer under control
                
            if(nrPyramids[level+1]>200000 && nr_done < nr_pyramids) // CHOOSE SAME VALUE IN process_pytamids
                 skip_remaining_pyr=true;
       }
       
        if(skip_remaining_tri){
            if (verbose)
                verboseOutput() << nr_done << " of " << nr_pyramids << 
                    " pyramids done on level " << level << ", ";
            Top_Cone->evaluate_triangulation();
        }

        if(skip_remaining_pyr){
            if (verbose)
                verboseOutput() << ">>>>>>>>>>>>>>> " << nr_done << " of " << nr_pyramids << 
                            " pyramids done, descending to level " << level+1 << endl;
            evaluate_stored_pyramids(level+1);
       }
    
     }while(skip_remaining_tri || skip_remaining_pyr);
     
     if (verbose)
         verboseOutput() << nr_done << " of " << nr_pyramids << " pyramids on level "<< level << " done!"<<endl;
     if(check_evaluation_buffer())
     {
        Top_Cone->evaluate_triangulation();
     }
     
     Pyramids[level].clear();
     nrPyramids[level]=0;
     evaluate_stored_pyramids(level+1);  
}

//---------------------------------------------------------------------------

/* builds the cone successively by inserting generators, computes all essential data
except global reduction */
template<typename Integer>
void Full_Cone<Integer>::build_cone() {
    if(dim>0){            //correction needed to include the 0 cone;
    if (verbose && !is_pyramid) {
        verboseOutput()<<"\n************************************************************\n";
        verboseOutput()<<"starting primal algorithm ";
        if (do_partial_triangulation) verboseOutput()<<"with partial triangulation ";
        if (do_triangulation) {
            verboseOutput()<<"with full triangulation ";
            if (!keep_triangulation) verboseOutput()<<"and direct evaluation ";
        }
        if (!do_triangulation && !do_partial_triangulation) verboseOutput()<<"(only support hyperplanes) ";
        verboseOutput()<<"..."<<endl;
    }
    size_t i;
    
    bool supphyp_recursion=false;
    bool tri_recursion=false;

    // DECIDE WHETHER TO USE RECURSION
    long long RecBoundSuppHyp = dim*dim*dim;
    RecBoundSuppHyp *= RecBoundSuppHyp * 10; //dim^6 * 10
    int bound_div = nr_gen-dim+1;
    if(bound_div > 3* (int) dim) bound_div = 3*dim;
    RecBoundSuppHyp /= bound_div;

    size_t RecBoundTriang = 1000000;  // 1 Mio      5000000; // 5Mio
    size_t EvalBoundTriang = 2500000; // 2.5 Mio --- should coincide with the value in check_evaluation_buffer_size()


    // if(!is_pyramid) cout << "RecBoundSuppHyp = "<<RecBoundSuppHyp<<endl;

    find_and_evaluate_start_simplex();
    
    Integer scalar_product;
    bool new_generator;
    size_t last_to_be_inserted; // good to know in case of do_all_hyperplanes==false
    last_to_be_inserted=nr_gen-1; 
    for(int j=nr_gen-1;j>=0;--j){
        if(isComputed(ConeProperty::ExtremeRays)){
            if(!in_triang[j] && Extreme_Rays[j]){
                last_to_be_inserted=j;
                break;
            }
        }
        else
            if(!in_triang[j]){
                last_to_be_inserted=j;
                break;
            }
    }

    for (i=0;i<nr_gen;++i) {
    
        if(in_triang[i] || (isComputed(ConeProperty::ExtremeRays) && !Extreme_Rays[i]))
            continue;

        new_generator=false;

        typename list< FACETDATA >::iterator l=Facets.begin();

        long long nr_pos=0; long long nr_neg=0;
        vector<Integer> L;
        size_t old_nr_supp_hyps=Facets.size();                
        
        size_t lpos=0;
        #pragma omp parallel for private(L,scalar_product) firstprivate(lpos,l) reduction(+: nr_pos, nr_neg) schedule(dynamic)
        for (size_t k=0; k<old_nr_supp_hyps; k++) {
            for(;k > lpos; lpos++, l++) ;
            for(;k < lpos; lpos--, l--) ;

            L=Generators.read(i);
            scalar_product=v_scalar_product(L,(*l).Hyp);            
            // l->ValPrevGen=l->ValNewGen;  // last new generator is now previous generator
            l->ValNewGen=scalar_product;
            if (scalar_product<0) {
                new_generator=true;
                nr_neg++;
            }
            if (scalar_product>0) {
                nr_pos++;
            }
        }  //end parallel for
        
        if(!new_generator)
            continue;

        // the i-th generator is used in the triangulation
        in_triang[i]=true;
        if (ht1_triangulation && isComputed(ConeProperty::LinearForm))
            ht1_triangulation = (gen_degrees[i] == 1);
        
            
        // Magic Bounds to decide whether to use pyramids
        if( supphyp_recursion || (recursion_allowed && nr_neg*nr_pos>RecBoundSuppHyp)){  // go to pyramids because of supphyps
             if(check_evaluation_buffer()){
                // cout << "Evaluation Build Mitte" << endl;
                    Top_Cone->evaluate_triangulation();
            }
            // cout << "In SuppHyp Rec" << endl;   
            supphyp_recursion=true;
            process_pyramids(i,true); //recursive
            // cout << "Zurück " << endl; 
        }
        else{
            if( tri_recursion || (do_triangulation // 
                         && (nr_neg*TriangulationSize > RecBoundTriang 
                                || 3*omp_get_max_threads()*TriangulationSize>EvalBoundTriang ))){ // go to pyramids because of triangulation
                if(check_evaluation_buffer()){
                    Top_Cone->evaluate_triangulation();
                }
                // if(!is_pyramid)
                // cout << nr_neg << " " <<TriangulationSize << " " << omp_get_max_threads() << " " << nr_neg*TriangulationSize << " " <<  3*omp_get_max_threads()*TriangulationSize << endl;
                tri_recursion=true;
                process_pyramids(i,false); //non-recursive
            }
            else{  // no pyramids necesary
                if(do_partial_triangulation)
                    process_pyramids(i,false); // non-recursive
                if(do_triangulation)
                    extend_triangulation(i);
            }

            if(do_all_hyperplanes || i!=last_to_be_inserted)
                find_new_facets(i);
        }
        
        // removing the negative hyperplanes if necessary
        if(do_all_hyperplanes || i!=last_to_be_inserted){
            l=Facets.begin();
            for (size_t j=0; j<old_nr_supp_hyps;j++){
                if (l->ValNewGen<0) 
                    l=Facets.erase(l);
            else 
                l++;
            }
        }

        if (verbose && !is_pyramid) {
            verboseOutput() << "generator="<< i+1 <<" and "<<Facets.size()<<" hyperplanes... ";
            if(keep_triangulation)
                verboseOutput() << TriangulationSize << " simplices ";
            verboseOutput()<<endl;
        }
    }

    if(do_all_hyperplanes){
        typename list<FACETDATA>::const_iterator IHV=Facets.begin();
        for(;IHV!=Facets.end();IHV++){
            Support_Hyperplanes.push_back(IHV->Hyp);
        }
    }
    

    } // end if (dim>0)
    
    Facets.clear();
    is_Computed.set(ConeProperty::SupportHyperplanes);
    
    if(!is_pyramid)
    {
        evaluate_stored_pyramids(0);                    
    }

    transfer_triangulation_to_top(); // transfer remaining simplices to top
    if(check_evaluation_buffer()){
        // cout << "Evaluating in build_cone at end, pyr level " << pyr_level << endl;
        // cout << "Evaluation Build Ende " << is_pyramid << endl;
        Top_Cone->evaluate_triangulation();
    }
    
    if(!is_pyramid && !keep_triangulation) // force evaluation of remaining simplices
        // cout << " Top cone evaluation in build_cone" << endl;
        evaluate_triangulation();          // on top level

    if(!is_pyramid && keep_triangulation)  // in this case triangulation now complete
        is_Computed.set(ConeProperty::Triangulation);  // and stored 
               
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::extreme_rays_and_ht1_check() {
    check_pointed();
    if(!pointed) return;
    compute_extreme_rays();
    ht1_check();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::set_degrees() {
    if(gen_degrees.size()==0 && isComputed(ConeProperty::LinearForm)) // now we set the degrees
    {
        gen_degrees.resize(nr_gen);
        vector<Integer> gen_degrees_Integer=Generators.MxV(Linear_Form);
        for (size_t i=0; i<nr_gen; i++) {
            if (gen_degrees_Integer[i] < 1) {
                errorOutput() << "Grading gives non-positive value " << gen_degrees_Integer[i] << " for generator " << i+1 << "." << endl;
                throw BadInputException();
            }
            gen_degrees[i] = explicit_cast_to_long(gen_degrees_Integer[i]);
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::sort_gens_by_degree() {
    if(gen_degrees.size()==0 || ht1_extreme_rays)
        return;
    
    list<vector<Integer> > genList;
    vector<Integer> v(dim+3);
    vector<Integer> w(dim);
    size_t i,j;
    
    for(i=0;i<nr_gen;i++){
        v[0]=gen_degrees[i];
        v[1]=i;                // keep the input order as far as possible
        w=Generators[i];
        for(j=0;j<dim;j++)
            v[j+2]=w[j];
        v[dim+2]=0;
        if(Extreme_Rays[i]) // after sorting we must recover the extreme rays
            v[dim+2]=1;
        genList.push_back(v);
    }
    genList.sort();
    
    i=0;
    typename list<vector<Integer> >::iterator g=genList.begin();
    for(;g!=genList.end();++g){
        v=*g;
        gen_degrees[i]=explicit_cast_to_long<Integer>(v[0]);
        Extreme_Rays[i]=false;
        if(v[dim+2]>0)
            Extreme_Rays[i]=true;
        for(j=0;j<dim;j++)
            w[j]=v[j+2];
        Generators[i]=w;
        i++;
    }
    
    cout << "Degrees after sort" << endl;
    for(size_t k=0;k<nr_gen;k++)
        cout << gen_degrees[k] << " " ;
    cout << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_support_hyperplanes(){
    if(isComputed(ConeProperty::SupportHyperplanes))
        return;
    bool save_tri      = do_triangulation;
    bool save_part_tri = do_partial_triangulation;
    bool save_HB       = do_Hilbert_basis;
    bool save_ht1_el   = do_ht1_elements;
    bool save_h_vector = do_h_vector;
    do_triangulation         = false;
    do_partial_triangulation = false;
    do_Hilbert_basis         = false; 
    do_ht1_elements          = false; 
    do_h_vector              = false;
    build_cone();
    do_triangulation         = save_tri;
    do_partial_triangulation = save_part_tri;
    do_Hilbert_basis         = save_HB;
    do_ht1_elements          = save_ht1_el;
    do_h_vector              = save_h_vector;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::check_evaluation_buffer(){

    return(omp_get_level()==0 && check_evaluation_buffer_size());
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::check_evaluation_buffer_size(){

    size_t EvalBoundTriang = 2500000; // 2.5 Mio --- should coincide with the value in build_cone()

    return(!Top_Cone->keep_triangulation && 
               Top_Cone->TriangulationSize > EvalBoundTriang);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::transfer_triangulation_to_top(){  // NEW EVA

    size_t i;

    // cout << "Pyr level " << pyr_level << endl;
    
    if(!is_pyramid) {  // we are in top cone
        if(check_evaluation_buffer()){
            evaluate_triangulation();
        }
        return;      // no transfer necessary
    }

    // now we are in a pyramid

    // cout << "In pyramid " << endl;
  
    typename list<SHORTSIMPLEX>::iterator pyr_simp=Triangulation.begin();
    for(;pyr_simp!=Triangulation.end();pyr_simp++)
        for(i=0;i<dim;i++)
            pyr_simp->key[i]=Top_Key[pyr_simp->key[i]];

    // cout << "Keys transferred " << endl;
    #pragma omp critical(TRIANG)
    {
        Top_Cone->Triangulation.splice(Top_Cone->Triangulation.end(),Triangulation);
        Top_Cone->TriangulationSize+=TriangulationSize;
    }
    TriangulationSize  =  0;

    // cout << "Done." << endl;
  
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_triangulation(){

    assert(omp_get_level()==0);
   
    if(TriangulationSize>0)
    {
    const long VERBOSE_STEPS = 50;
    long step_x_size = TriangulationSize-VERBOSE_STEPS;
    if (verbose) {
        verboseOutput() << "evaluating "<<TriangulationSize<<" simplices" <<endl;
        /* verboseOutput() << "---------+---------+---------+---------+---------+"
                        << " (one | per 2%)" << endl;*/
    }
    
    totalNrSimplices+=TriangulationSize;

    #pragma omp parallel 
    {
        typename list<SHORTSIMPLEX>::iterator s = Triangulation.begin();
        size_t spos=0;
        SimplexEvaluator<Integer> simp(*this);
        #pragma omp for schedule(dynamic) 
        for(size_t i=0; i<TriangulationSize; i++){
            for(; i > spos; ++spos, ++s) ;
            for(; i < spos; --spos, --s) ;

            s->height = simp.evaluate(s->key,s->height);
            if(keep_triangulation)
                sort(s->key.begin(),s->key.end());
            if (verbose) {
                #pragma omp critical(VERBOSE)
                while ((long)(i*VERBOSE_STEPS) >= step_x_size) {
                    step_x_size += TriangulationSize;
                    verboseOutput() << "|" <<flush;
                }
            }
        }
        // get accumulated data from the SimplexEvaluator
        #pragma omp critical(MULTIPLICITY)
        multiplicity += simp.getMultiplicitySum(); 
        if (do_h_vector) {
            #pragma omp critical(HSERIES)
            Hilbert_Series += simp.getHilbertSeriesSum();
        }
    }

    if (do_partial_triangulation) {
        if (verbose) {
            cout<<endl<<"ht1: "<<Ht1_Elements.size();
            cout<<"   cand: "<<Candidates.size()<<flush;
        }
        Ht1_Elements.sort();
        Ht1_Elements.unique();
        Candidates.sort();
        Candidates.unique();
    }


    if (verbose)
    {
        verboseOutput()  << endl  << totalNrSimplices << " simplices";
        if(do_Hilbert_basis)
            verboseOutput() << ", " << Candidates.size() << " HB candidates";
        if(do_ht1_elements)
            verboseOutput() << ", " << Ht1_Elements.size()<< " ht1 vectors";
        verboseOutput() << " accumulated." << endl;
    }
    
    if(!keep_triangulation){
        // Triangulation.clear();
        #pragma omp critical(FREESIMPL)
        FreeSimpl.splice(FreeSimpl.begin(),Triangulation);       
        TriangulationSize=0;
    }
    
} // TriangulationSize

}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm(){

    // if keep_triangulation==false we must first find a grading if it is needed
    if (!keep_triangulation && !isComputed(ConeProperty::LinearForm)
        && (do_triangulation || do_ht1_elements || do_h_vector)) {
        ht1_check();
        if(!ht1_generated && !isComputed(ConeProperty::ExtremeRays)) {
            if (verbose) {
                verboseOutput() << "Cannot find grading s.t. all generators have degree 1! Computing Extreme rays first:" << endl;
            }
            compute_support_hyperplanes();
            extreme_rays_and_ht1_check();
            if(!pointed) return;

            Support_Hyperplanes.clear();  // will be computed again by build_cone
            is_Computed.reset(ConeProperty::SupportHyperplanes);
            for(size_t i=0;i<nr_gen;i++)
                in_triang[i]=false;
        }
    }
    set_degrees();
    sort_gens_by_degree();

    /***** Main Work is done in build_cone() *****/
    build_cone();  // evaluates if keep_triangulation==false
    /***** Main Work is done in build_cone() *****/
    
    cout << "TotPyr "<< totalNrPyr << endl;
    cout << "Uni "<< Unimod << " Ht1NonUni " << Ht1NonUni << " NonDecided " << NonDecided << " TotNonDec " << NonDecidedHyp<< endl;

    extreme_rays_and_ht1_check();
    if(!pointed) return;

    if (keep_triangulation) {
        evaluate_triangulation();
    }
    FreeSimpl.clear();
    
    if (do_triangulation)
        is_Computed.set(ConeProperty::TriangulationSize,true);
    if (isComputed(ConeProperty::LinearForm) && do_triangulation)
        is_Computed.set(ConeProperty::Multiplicity,true);
        
    if (do_Hilbert_basis) {
        global_reduction();
        is_Computed.set(ConeProperty::HilbertBasis,true);
        check_integrally_closed();
    }
    
    if (isComputed(ConeProperty::LinearForm) && do_Hilbert_basis) {
        select_ht1_elements();
        check_ht1_hilbert_basis();
    }
    if (do_ht1_elements) {
        for(size_t i=0;i<nr_gen;i++)
            if(v_scalar_product(Linear_Form,Generators.read(i))==1)
                Ht1_Elements.push_front(Generators.read(i));
        Ht1_Elements.sort();
        Ht1_Elements.unique();
        is_Computed.set(ConeProperty::Ht1Elements,true);
    }
    if (do_h_vector) {
        Hilbert_Series.simplify();
        is_Computed.set(ConeProperty::HilbertSeries);
    }

}

   
//---------------------------------------------------------------------------
// Normaliz modes (public)
//---------------------------------------------------------------------------

// pure dualization
template<typename Integer>
void Full_Cone<Integer>::dualize_cone() {  
    compute_support_hyperplanes();
    reset_tasks();
}

// -s
template<typename Integer>
void Full_Cone<Integer>::support_hyperplanes() {
    // recursion_allowed=true; 
    compute_support_hyperplanes();
    extreme_rays_and_ht1_check();
    reset_tasks();
}

// -v
template<typename Integer>
void Full_Cone<Integer>::support_hyperplanes_triangulation() {
    do_triangulation=true;
    keep_triangulation=true;
    primal_algorithm();
    reset_tasks();
}


// -V
template<typename Integer>
void Full_Cone<Integer>::support_hyperplanes_triangulation_pyramid() {   
    do_triangulation=true; 
    primal_algorithm();
    reset_tasks();
}

//-n with -a
template<typename Integer>
void Full_Cone<Integer>::triangulation_hilbert_basis() {
    do_Hilbert_basis=true;
    do_triangulation=true;
    keep_triangulation=true;
    primal_algorithm();
    reset_tasks();
}

//-n without -a
template<typename Integer>
void Full_Cone<Integer>::multiplicity_hilbert_basis() {
    do_Hilbert_basis=true;
    do_triangulation=true;
    // keep_triangulation=true;
    primal_algorithm();
    reset_tasks();
}

// -N
template<typename Integer>
void Full_Cone<Integer>::hilbert_basis() {
    do_Hilbert_basis=true;
    do_partial_triangulation=true;
    primal_algorithm();
    reset_tasks();
}

// -h
template<typename Integer>
void Full_Cone<Integer>::hilbert_basis_polynomial() {
    do_Hilbert_basis=true;
    do_h_vector=true;
    do_triangulation=true;
    keep_triangulation=true;
    primal_algorithm();
    reset_tasks();    
}

// -H
template<typename Integer>
void Full_Cone<Integer>::hilbert_basis_polynomial_pyramid() {
    do_Hilbert_basis=true;
    do_h_vector=true;
    do_triangulation=true;
    primal_algorithm();
    reset_tasks();    
}

// -p
template<typename Integer>
void Full_Cone<Integer>::hilbert_polynomial() {
    do_ht1_elements=true;
    do_h_vector=true;
    do_triangulation=true;
    keep_triangulation=true;
    primal_algorithm();
    reset_tasks();
}

// -P
template<typename Integer>
void Full_Cone<Integer>::hilbert_polynomial_pyramid() {
    do_ht1_elements=true;
    do_h_vector=true;
    do_triangulation=true;
    primal_algorithm();
    reset_tasks();
}

// -1
template<typename Integer>
void Full_Cone<Integer>::ht1_elements() {
    do_ht1_elements=true;
    do_partial_triangulation=true;
    primal_algorithm();
    reset_tasks();
}

template<typename Integer>
void Full_Cone<Integer>::dual_mode() {
    Support_Hyperplanes.sort();
    Support_Hyperplanes.unique();
    Support_Hyperplanes.remove(vector<Integer>(dim,0));

    if(dim>0) {            //correction needed to include the 0 cone;
        ht1_check();
        if (isComputed(ConeProperty::LinearForm)) {
            if (verbose) { 
                verboseOutput() << "Find height 1 elements" << endl;
            }
            select_ht1_elements();
        }
    } else {
        ht1_extreme_rays = ht1_generated = true;
        Linear_Form=vector<Integer>(dim);
        is_Computed.set(ConeProperty::IsHt1ExtremeRays);
        is_Computed.set(ConeProperty::IsHt1Generated);
        is_Computed.set(ConeProperty::LinearForm);
    }
    if (isComputed(ConeProperty::LinearForm)) check_ht1_hilbert_basis();
    check_integrally_closed();
}

//---------------------------------------------------------------------------
// Checks and auxiliary algorithms
//---------------------------------------------------------------------------


template<typename Integer>
Simplex<Integer> Full_Cone<Integer>::find_start_simplex() const {
    if (isComputed(ConeProperty::ExtremeRays)) {
        vector<key_t> marked_extreme_rays(0);
        for (size_t i=0; i<nr_gen; i++) {
            if (Extreme_Rays[i])
                marked_extreme_rays.push_back(i);
        }
        vector<key_t> key_extreme = Generators.submatrix(Extreme_Rays).max_rank_submatrix_lex(dim);
        assert(key_extreme.size() == dim);
        vector<key_t> key(dim);
        for (key_t i=0; i<dim; i++) {
            key[i] = marked_extreme_rays[key_extreme[i]];
        }
        return Simplex<Integer>(key, Generators);
    } 
    else {
    // assert(Generators.rank()>=dim); 
        return Simplex<Integer>(Generators);
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::select_matrix_from_list(const list<vector<Integer> >& S,
                                   vector<size_t>& selection){

    sort(selection.begin(),selection.end());
    assert(selection.back()<S.size());
    size_t i=0,j=0;
    size_t k=selection.size();
    Matrix<Integer> M(selection.size(),S.front().size());
    typename list<vector<Integer> >::const_iterator ll=S.begin();
    for(;ll!=S.end()&&i<k;++ll){
        if(j==selection[i]){
            M[i]=*ll;
            i++;
        }
        j++;
    }
    return M;
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_extreme_rays(){

    if (isComputed(ConeProperty::ExtremeRays))
        return;
    assert(isComputed(ConeProperty::SupportHyperplanes));

    if(dim*Support_Hyperplanes.size() < nr_gen)
         compute_extreme_rays_rank();
    else
         compute_extreme_rays_compare();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_extreme_rays_rank(){

    
    size_t i,j;
    typename list<vector<Integer> >::iterator s;
    vector<size_t> gen_in_hyperplanes;
    gen_in_hyperplanes.reserve(Support_Hyperplanes.size());
    Matrix<Integer> M;
    
    for(i=0;i<nr_gen;++i){
        Extreme_Rays[i]=false;
        j=0;
        gen_in_hyperplanes.clear();
        for(s=Support_Hyperplanes.begin();s!=Support_Hyperplanes.end();++s){
            if(v_scalar_product(Generators[i],*s)==0)
                gen_in_hyperplanes.push_back(j);
            j++;
        }
        if(gen_in_hyperplanes.size()< dim-1)
            continue;
        M=select_matrix_from_list(Support_Hyperplanes,gen_in_hyperplanes);
        if(M.rank_destructive()>=dim-1)
            Extreme_Rays[i]=true;   
    }

    is_Computed.set(ConeProperty::ExtremeRays);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_extreme_rays_compare(){

    size_t i,j,k,l,t;
    // Matrix<Integer> SH=getSupportHyperplanes().transpose();
    // Matrix<Integer> Val=Generators.multiplication(SH);
    size_t nc=Support_Hyperplanes.size();
    
    vector<vector<bool> > Val(nr_gen);
    for (i=0;i<nr_gen;++i)
       Val[i].resize(nc);
        
    // Attention: in this routine Val[i][j]==0, i.e. false, indicates that
    // the i-th generator is contained in the j-th support hyperplane
    
    vector<key_t> Zero(nc);
    vector<key_t> nr_zeroes(nr_gen);
    typename list<vector<Integer> >::iterator s;

    for (i = 0; i <nr_gen; i++) {
        k=0;
        Extreme_Rays[i]=true;
        s=Support_Hyperplanes.begin();
        for (j = 0; j <nc; ++j,++s) {
            if (v_scalar_product(Generators[i],*s)==0) {
                k++;
                Val[i][j]=false;                
            }
            else
                Val[i][j]=true;  
        }
        nr_zeroes[i]=k;
        if (k<dim-1||k==nc)  // not contained in enough facets or in all (0 as generator)
            Extreme_Rays[i]=false;
    }

    for (i = 0; i <nr_gen; i++) {
        if(!Extreme_Rays[i])  // already known to be non-extreme
            continue;

        k=0;
        for (j = 0; j <nc; j++) {
            if (Val[i][j]==false) {
                Zero[k]=j;
                k++;
            }
        }

        for (j = 0; j <nr_gen; j++) {
            if (i!=j && Extreme_Rays[j]                // not compare with itself or a known nonextreme ray
                     && nr_zeroes[i]<nr_zeroes[j]) {   // or something whose zeroes cannot be a superset
                l=0;
                for (t = 0; t < nr_zeroes[i]; t++) {
                    if (Val[j][Zero[t]]==false)
                        l++;
                    if (l>=nr_zeroes[i]) {
                        Extreme_Rays[i]=false;
                        break;
                    }
                }
            }
        }
    }

    is_Computed.set(ConeProperty::ExtremeRays);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::select_ht1_elements() {

    typename list<vector<Integer> >::iterator h = Hilbert_Basis.begin();
    for(;h!=Hilbert_Basis.end();h++)
        if(v_scalar_product(Linear_Form,*h)==1)
            Ht1_Elements.push_back(*h);
    is_Computed.set(ConeProperty::Ht1Elements,true);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_pointed() {
    assert(isComputed(ConeProperty::SupportHyperplanes));
    if (isComputed(ConeProperty::IsPointed))
        return;
    Matrix<Integer> SH = getSupportHyperplanes();
    pointed = (SH.rank_destructive() == dim);
    is_Computed.set(ConeProperty::IsPointed);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::ht1_check() {
    if (!isComputed(ConeProperty::LinearForm)          // we still need it and
     && !isComputed(ConeProperty::IsHt1ExtremeRays)) { // we have not tried it
        if (isComputed(ConeProperty::ExtremeRays)) {
            Matrix<Integer> Extreme=Generators.submatrix(Extreme_Rays);
            Linear_Form = Extreme.homogeneous(ht1_extreme_rays);
            is_Computed.set(ConeProperty::IsHt1ExtremeRays);
            if (ht1_extreme_rays) {
                is_Computed.set(ConeProperty::LinearForm);
            }
        } else // extreme rays not known
        if (!isComputed(ConeProperty::IsHt1Generated)) {
            Linear_Form = Generators.homogeneous(ht1_generated);
            is_Computed.set(ConeProperty::IsHt1Generated);
            if (ht1_generated) {
                ht1_extreme_rays=true;
                is_Computed.set(ConeProperty::IsHt1ExtremeRays);
                is_Computed.set(ConeProperty::LinearForm);
            }
        }
    }

    //now we hopefully have a grading

    if (!isComputed(ConeProperty::LinearForm)) {
        if (isComputed(ConeProperty::ExtremeRays)) {
            // there is no hope to find a grading later
            ht1_generated = false;
            ht1_extreme_rays = false;
            is_Computed.set(ConeProperty::IsHt1ExtremeRays);
            is_Computed.set(ConeProperty::IsHt1Generated);
            if (do_ht1_elements || do_h_vector) {
                errorOutput() << "No grading specified and cannot find one. "
                              << "Disabling some computations!" << endl;
                do_ht1_elements = false;
                do_h_vector = false;
            }
        }
        return; // we are done
    }
    
    set_degrees();
        
    if (!isComputed(ConeProperty::IsHt1Generated)) {
        ht1_generated = true;
        for (size_t i = 0; i < nr_gen; i++) {
            if (gen_degrees[i] != 1) {
                ht1_generated = false;
                break;
            }
        }
        is_Computed.set(ConeProperty::IsHt1Generated);
        if (ht1_generated) {
            ht1_extreme_rays = true;
            is_Computed.set(ConeProperty::IsHt1ExtremeRays);
        }
    }
    if (!isComputed(ConeProperty::IsHt1ExtremeRays)
      && isComputed(ConeProperty::ExtremeRays)) {
        ht1_extreme_rays = true;
        for (size_t i = 0; i < nr_gen; i++) {
            if (Extreme_Rays[i] && gen_degrees[i] != 1) {
                ht1_extreme_rays = false;
                break;
            }
        }
        is_Computed.set(ConeProperty::IsHt1ExtremeRays);
    }
}

template<typename Integer>
void Full_Cone<Integer>::check_ht1_hilbert_basis() {
    if (isComputed(ConeProperty::IsHt1HilbertBasis))
        return;

    if ( !isComputed(ConeProperty::LinearForm) || !isComputed(ConeProperty::HilbertBasis)) {
        errorOutput() << "WARNING: unsatisfied preconditions in check_ht1_hilbert_basis()!" <<endl;
        return;
    }
    
    if (isComputed(ConeProperty::Ht1Elements)) {
        ht1_hilbert_basis = (Ht1_Elements.size() == Hilbert_Basis.size());
    } else {
        ht1_hilbert_basis = true;
        typename list< vector<Integer> >::iterator h;
        for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
            if (v_scalar_product((*h),Linear_Form)!=1) {
                ht1_hilbert_basis = false;
                break;
            }
        }
    }
    is_Computed.set(ConeProperty::IsHt1HilbertBasis);
}

template<typename Integer>
void Full_Cone<Integer>::check_integrally_closed() {
    if (isComputed(ConeProperty::IsIntegrallyClosed))
        return;

    if ( !isComputed(ConeProperty::HilbertBasis)) {
        errorOutput() << "WARNING: unsatisfied preconditions in check_integrally_closed()!" <<endl;
        return;
    }
    integrally_closed = false;
    if (Hilbert_Basis.size() <= nr_gen) {
        integrally_closed = true;
        typename list< vector<Integer> >::iterator h;
        for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
            integrally_closed = false;
            for (size_t i=0; i< nr_gen; i++) {
                if ((*h) == Generators.read(i)) {
                    integrally_closed = true;
                    break;
                }
            }
            if (!integrally_closed) {
                break;
            }
        }
    }
    is_Computed.set(ConeProperty::IsIntegrallyClosed);
}

//---------------------------------------------------------------------------
// Global reduction
//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::is_reducible(list< vector<Integer>* >& Irred, const vector< Integer >& new_element){
    size_t i;
    size_t s=Support_Hyperplanes.size();
    // new_element can be longer than dim (it has one extra entry for the norm)
    // the scalar product function just takes the first dim entries
    vector <Integer> scalar_product=l_multiplication(Support_Hyperplanes,new_element);
    typename list< vector<Integer>* >::iterator j;
    vector<Integer> *reducer;
    for (j =Irred.begin(); j != Irred.end(); j++) {
        reducer=(*j);
        for (i = 0; i < s; i++) {
            if ((*reducer)[i]>scalar_product[i]){
                break;
            }
        }
        if (i==s) {
            //found a "reducer" and move it to the front
            Irred.push_front(*j);
            Irred.erase(j);
            return true;
        }
    }
    return false;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::global_reduction() {
    Integer norm;
    
    list <vector<Integer> > HB;
    typename list <vector<Integer> >::iterator c;
    
    for (size_t i = 0; i <nr_gen; i++) {
        Candidates.push_front(Generators.read(i));
    }
/*    if(verbose) verboseOutput()<<"sorting the candidates... "<<flush;
    Candidates.sort();
    if(verbose) verboseOutput()<<"make them unique... "<<flush;
    Candidates.unique();
    if(verbose) verboseOutput()<<"done."<<endl;
*/  // Duplicates are avoided or removed earlier
    if (nr_gen == dim) { // cone is simplicial, therefore no global reduction is necessary
        if (verbose) {
            verboseOutput()<<"Cone is simplicial, no global reduction necessary."<<endl;
        }
        Hilbert_Basis.splice(Hilbert_Basis.end(), Candidates);
        return;
    }
    

    vector<Integer> degree_function=compute_degree_function();

    c = Candidates.begin();
    size_t cpos = 0;
    size_t csize=Candidates.size();
    
    if(verbose) {
        verboseOutput()<<"computing the degrees of the candidates... "<<flush;
    }
    //go over candidates: do single scalar product and save it at the end of the candidate
    //for (c = Candidates.begin(); c != Candidates.end(); c++) 
    vector<Integer> scalar_product;
    for (size_t j=0; j<csize; ++j) {
        for(;j > cpos; ++cpos, ++c) ;
        for(;j < cpos; --cpos, --c) ;

        norm=v_scalar_product(degree_function,(*c));
        c->push_back(norm);
    }
    if(verbose) {
        verboseOutput()<<"sorting the list... "<<endl;
    }
    Candidates.sort(compare_last<Integer>);
    if (verbose) {
        verboseOutput()<< csize <<" candidate vectors sorted."<<endl;
    }
    
    // do global reduction
    list< vector<Integer> > HBtmp;
    Integer norm_crit;
    while ( !Candidates.empty() ) {
        //use norm criterion to find irreducible elements
        c=Candidates.begin();
        norm_crit=(*c)[dim]*2;  //candidates with smaller norm are irreducible
        if ( Candidates.back()[dim] < norm_crit) { //all candidates are irreducible
            if (verbose) {
                verboseOutput()<<Hilbert_Basis.size()+Candidates.size();
                verboseOutput()<<" Hilbert Basis elements of degree <= "<<norm_crit-1<<"; done"<<endl;
            }
            for (; c!=Candidates.end(); ++c) {
                c->pop_back();
            }
            Hilbert_Basis.splice(Hilbert_Basis.end(), Candidates);
            break;
        }
        while ( (*c)[dim] < norm_crit ) { //can't go over the end because of the previous if
            // remove norm
            c->pop_back();
            // push the scalar products to the reducer list
            HBtmp.push_back(l_multiplication(Support_Hyperplanes, *c));
            // and the candidate itself to the Hilbert basis
            Hilbert_Basis.splice(Hilbert_Basis.end(), Candidates, c++);
        }
        csize = Candidates.size();
        if (verbose) {
            verboseOutput()<<Hilbert_Basis.size()<< " Hilbert Basis elements of degree <= "<<norm_crit-1<<"; "<<csize<<" candidates left"<<endl;
        }

        // reduce candidates against HBtmp
        // fill pointer list
        list < vector <Integer>* >  HBpointers;  // used to put "reducer" to the front
        c = HBtmp.begin();
        while (c != HBtmp.end()) {
            HBpointers.push_back(&(*(c++)));
        }

        long VERBOSE_STEPS = 50;      //print | for 2%
        if (verbose && csize>50000) { //print | for 1000 candidates
            VERBOSE_STEPS=csize/1000;
        }
        long step_x_size = csize-VERBOSE_STEPS;
        long counter = 0;
        long steps_done = 0;
        if (verbose) {
            verboseOutput() << "---------+---------+---------+---------+---------+";
            if (VERBOSE_STEPS == 50) {
                verboseOutput() << " (one | per 2%)" << endl;
            } else { 
                verboseOutput() << " (one | per 1000 candidates)" << endl;
            }
        }


        #pragma omp parallel private(c,cpos) firstprivate(HBpointers)
        {
        
        c=Candidates.begin();
        cpos=0;
        #pragma omp for schedule(dynamic)
        for (size_t k=0; k<csize; ++k) {
            for(;k > cpos; ++cpos, ++c) ;
            for(;k < cpos; --cpos, --c) ;
            
            if ( is_reducible(HBpointers, *c) ) {
                (*c)[dim]=-1; //mark as reducible
            }

            if (verbose) {
                #pragma omp critical(VERBOSE)
                {
                counter++;

                while (counter*VERBOSE_STEPS >= step_x_size) {
                    steps_done++;
                    step_x_size += csize;
                    verboseOutput() << "|" <<flush;
                    if(VERBOSE_STEPS > 50 && steps_done%50 == 0) {
                        verboseOutput() << "  " << (steps_done) << "000" << endl;
                    }
                }
                } //end critical(VERBOSE)
            }
        } //end for
        } //end parallel
        if (verbose) verboseOutput() << endl;

        // delete reducible candidates
        c = Candidates.begin();
        while (c != Candidates.end()) {
            if ((*c)[dim]==-1) {
                c = Candidates.erase(c);
            } else {
                ++c;
            }
        }
        HBtmp.clear();
    }

    if (verbose) {
        verboseOutput()<<Hilbert_Basis.size()<< " Hilbert Basis elements"<<endl;
    }
}


//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Full_Cone<Integer>::compute_degree_function() const {
    if(verbose) {
        verboseOutput()<<"computing degree function... "<<flush;
    }
    size_t i;  
    vector<Integer> degree_function(dim,0);
    if(isComputed(ConeProperty::LinearForm)){ //use Linear_From if we have one
        for (i=0; i<dim; i++) {
            degree_function[i] = Linear_Form[i];
        }
        if(verbose) {
            verboseOutput()<<"using given or homogenous linear form."<<endl;
        }
    } else { // add hyperplanes to get a degree function
        typename list< vector<Integer> >::const_iterator h;
        for (h=Support_Hyperplanes.begin(); h!=Support_Hyperplanes.end(); ++h) {
            for (i=0; i<dim; i++) {
                degree_function[i]+=(*h)[i];
            }
        } 
        v_make_prime(degree_function); //TODO maybe not needed
        if(verbose) {
            verboseOutput()<<"done."<<endl;
        }
    }
    return degree_function;
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Full_Cone<Integer>::primary_multiplicity() const{
    size_t i,j,k;
    Integer primary_multiplicity=0;
    vector <key_t> key,new_key(dim-1);
    Matrix<Integer> Projection(nr_gen,dim-1);
    for (i = 0; i < nr_gen; i++) {
        for (j = 0; j < dim-1; j++) {
            Projection.write(i,j,Generators.read(i,j));
        }
    }
    typename list< vector<Integer> >::const_iterator h;
    typename list<SHORTSIMPLEX>::const_iterator t;
    for (h =Support_Hyperplanes.begin(); h != Support_Hyperplanes.end(); ++h){
        if ((*h)[dim-1]!=0) {
            for (t =Triangulation.begin(); t!=Triangulation.end(); ++t){
                key=t->key;
                for (i = 0; i <dim; i++) {
                    k=0;
                    for (j = 0; j < dim; j++) {
                        if (j!=i && Generators.read(key[j],dim-1)==1) {
                            if (v_scalar_product(Generators.read(key[j]),(*h))==0) {
                                k++;
                            }
                        }
                        if (k==dim-1) {
                            for (j = 0; j <i; j++) {
                                new_key[j]=key[j];
                            }
                            for (j = i; j <dim-1; j++) {
                                new_key[j]=key[j+1];
                            }
                            // add the volume of the projected simplex
                            primary_multiplicity +=
                              Projection.submatrix(new_key).vol_destructive();
                        }
                    }
                }
            }
        }
    }
    return primary_multiplicity;
}
//---------------------------------------------------------------------------
// Constructors
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::reset_tasks(){
    do_triangulation = false;
    do_partial_triangulation = false;
    do_Hilbert_basis = false;
    do_ht1_elements = false;
    keep_triangulation = false;
    do_h_vector=false;
    is_pyramid = false;
    
     nrSimplicialPyr=0;
     totalNrPyr=0;
}

//---------------------------------------------------------------------------

template<typename Integer>
Full_Cone<Integer>::Full_Cone(Matrix<Integer> M){ // constructor of the top cone
    dim=M.nr_of_columns();
    if (dim!=M.rank()) {
        error_msg("error: Matrix with rank = number of columns needed in the constructor of the object Full_Cone<Integer>.\nProbable reason: Cone not full dimensional (<=> dual cone not pointed)!");
        throw NormalizException();
    }
    Generators = M;
    nr_gen=Generators.nr_of_rows();
    if (nr_gen != static_cast<size_t>(static_cast<key_t>(nr_gen))) {
        error_msg("To many generators to fit in range of key_t!");
        throw NormalizException();
    }
    //make the generators coprime and remove 0 rows
    vector<Integer> gcds = Generators.make_prime();
    vector<key_t> key=v_non_zero_pos(gcds);
    if (key.size() < nr_gen) {
        Generators=Generators.submatrix(key);
        nr_gen=Generators.nr_of_rows();
    }
    multiplicity = 0;
    is_Computed =  bitset<ConeProperty::EnumSize>();  //initialized to false
    is_Computed.set(ConeProperty::Generators);
    pointed = false;
    ht1_extreme_rays = false;
    ht1_generated = false;
    ht1_hilbert_basis = false;
    integrally_closed = false;
    
    reset_tasks();
    
    Extreme_Rays = vector<bool>(nr_gen,false);
    in_triang = vector<bool> (nr_gen,false);
    ht1_triangulation = true;
    if(dim==0){            //correction needed to include the 0 cone;
        multiplicity = 1;
        Hilbert_Series.add(vector<num_t>(1,1),vector<denom_t>());
        is_Computed.set(ConeProperty::HilbertSeries);
        is_Computed.set(ConeProperty::Triangulation);
    }
    pyr_level=-1;
    Top_Cone=this;
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=i;
    totalNrSimplices=0;
    TriangulationSize=0;
    
    FS.resize(omp_get_max_threads());
    
    Pyramids.resize(1);  // prepare storage for pyramids
    nrPyramids.resize(1);
    nrPyramids[0]=0;
    recursion_allowed=true;
    
    do_all_hyperplanes=true;
    parallel_in_pyramid=true;
}

//---------------------------------------------------------------------------

template<typename Integer>
Full_Cone<Integer>::Full_Cone(const Cone_Dual_Mode<Integer> &C) {

    dim = C.dim;
    Generators = C.get_generators();
    nr_gen = Generators.nr_of_rows();

    multiplicity = 0;
    is_Computed =  bitset<ConeProperty::EnumSize>();  //initialized to false
    is_Computed.set(ConeProperty::Generators);
    pointed = true;
    is_Computed.set(ConeProperty::IsPointed);
    ht1_extreme_rays = false;
    ht1_generated = false;
    ht1_triangulation = false;
    ht1_hilbert_basis = false;
    integrally_closed = false;
    
    reset_tasks();
    
    Extreme_Rays = vector<bool>(nr_gen,true); //all generators are extreme rays
    is_Computed.set(ConeProperty::ExtremeRays);
    Matrix<Integer> SH = C.SupportHyperplanes;
    for (size_t i=0; i < SH.nr_of_rows(); i++) {
        Support_Hyperplanes.push_back(SH.read(i));
    }
    is_Computed.set(ConeProperty::SupportHyperplanes);
    in_triang = vector<bool>(nr_gen,false);
    Hilbert_Basis = C.Hilbert_Basis;
    is_Computed.set(ConeProperty::HilbertBasis);
    if(dim==0){            //correction needed to include the 0 cone;
        multiplicity = 1;
        Hilbert_Series.add(vector<num_t>(1,1),vector<denom_t>());
        is_Computed.set(ConeProperty::HilbertSeries);
    }
    pyr_level=-1;
    Top_Cone=this;
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=i;
    totalNrSimplices=0;
    TriangulationSize=0;
    
    do_all_hyperplanes=true;
}
//---------------------------------------------------------------------------

/* constructor for pyramids */
template<typename Integer>
Full_Cone<Integer>::Full_Cone(Full_Cone<Integer>& C, const vector<key_t>& Key) {

    Generators = C.Generators.submatrix(Key);
    dim = Generators.nr_of_columns();
    nr_gen = Generators.nr_of_rows();
    
    Top_Cone=C.Top_Cone; // relate to top cone
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=C.Top_Key[Key[i]];
  
    multiplicity = 0;
    
    Extreme_Rays = vector<bool>(nr_gen,false);
    is_Computed.set(ConeProperty::ExtremeRays, C.isComputed(ConeProperty::ExtremeRays));
    if(isComputed(ConeProperty::ExtremeRays))
        for(size_t i=0;i<nr_gen;i++)
            Extreme_Rays[i]=C.Extreme_Rays[Key[i]];
    in_triang = vector<bool> (nr_gen,false);
    ht1_triangulation = true;
    
    Linear_Form=C.Linear_Form;
    is_Computed.set(ConeProperty::LinearForm, C.isComputed(ConeProperty::LinearForm));
    Order_Vector=C.Order_Vector;
    
    do_triangulation=C.do_triangulation;
    do_partial_triangulation=C.do_partial_triangulation;
    do_ht1_elements=C.do_ht1_elements;
    do_h_vector=C.do_h_vector;
    do_Hilbert_basis=C.do_Hilbert_basis;
    keep_triangulation=C.keep_triangulation;
    is_pyramid=true;
    
    // pyr_level set by the calling routine
    
    totalNrSimplices=0;
    if(C.gen_degrees.size()>0){ // now we copy the degrees
    	gen_degrees.resize(nr_gen);
        for (size_t i=0; i<nr_gen; i++) {
            gen_degrees[i] = C.gen_degrees[Key[i]];
        }
    }
    TriangulationSize=0;
    
    recursion_allowed=C.recursion_allowed;
    // do_all_hyperplanes to be set by calling routine
    parallel_in_pyramid=true; // this is the safe option
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::isComputed(ConeProperty::Enum prop) const{
    return is_Computed.test(prop);
}

//---------------------------------------------------------------------------
// Data access
//---------------------------------------------------------------------------

template<typename Integer>
size_t Full_Cone<Integer>::getDimension()const{
    return dim;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Full_Cone<Integer>::getNrGenerators()const{
    return nr_gen;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::isPointed()const{
    return pointed;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::isHt1ExtremeRays() const{
    return ht1_extreme_rays;
}

template<typename Integer>
bool Full_Cone<Integer>::isHt1HilbertBasis() const{
    return ht1_hilbert_basis;
}

template<typename Integer>
bool Full_Cone<Integer>::isIntegrallyClosed() const{
    return integrally_closed;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Full_Cone<Integer>::getLinearForm() const{
    return Linear_Form;
}

//---------------------------------------------------------------------------

template<typename Integer>
mpq_class Full_Cone<Integer>::getMultiplicity()const{
    return multiplicity;
}

//---------------------------------------------------------------------------

template<typename Integer>
const Matrix<Integer>& Full_Cone<Integer>::getGenerators()const{
    return Generators;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<bool> Full_Cone<Integer>::getExtremeRays()const{
    return Extreme_Rays;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getSupportHyperplanes()const{
    size_t s= Support_Hyperplanes.size();
    Matrix<Integer> M(s,dim);
    size_t i=0;
    typename list< vector<Integer> >::const_iterator l;
    for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++) {
        M.write(i,(*l));
        i++;
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::getTriangulation(list< vector<key_t> >& Triang, list<Integer>& TriangVol) const {
    Triang.clear();
    TriangVol.clear();
    vector<key_t> key(dim);
    typename list<SHORTSIMPLEX>::const_iterator l;
    for (l =Triangulation.begin(); l != Triangulation.end(); l++) {
        key=l->key;
        Triang.push_back(key);
        TriangVol.push_back(l->height);
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getHilbertBasis()const{
    size_t s= Hilbert_Basis.size();
    Matrix<Integer> M(s,dim);
    size_t i=0;
    typename list< vector<Integer> >::const_iterator l;
    for (l =Hilbert_Basis.begin(); l != Hilbert_Basis.end(); l++) {
        M.write(i,(*l));
        i++;
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getHt1Elements()const{
    size_t s= Ht1_Elements.size();
    Matrix<Integer> M(s,dim);
    size_t i=0;
    typename list< vector<Integer> >::const_iterator l;
    for (l =Ht1_Elements.begin(); l != Ht1_Elements.end(); l++) {
        M.write(i,(*l));
        i++;
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::error_msg(string s) const{
    errorOutput() <<"\nFull Cone "<< s<<"\n";
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::print()const{
    verboseOutput()<<"\ndim="<<dim<<".\n";
    verboseOutput()<<"\nnr_gen="<<nr_gen<<".\n";
    verboseOutput()<<"\nhyp_size="<<hyp_size<<".\n";
    verboseOutput()<<"\nHomogeneous is "<<ht1_generated<<".\n";
    verboseOutput()<<"\nLinear_Form is:\n";
    v_read(Linear_Form);
    verboseOutput()<<"\nMultiplicity is "<<multiplicity<<".\n";
    verboseOutput()<<"\nGenerators are:\n";
    Generators.read();
    verboseOutput()<<"\nExtreme_rays are:\n";
    v_read(Extreme_Rays);
    verboseOutput()<<"\nSupport Hyperplanes are:\n";
    l_read(Support_Hyperplanes);
    verboseOutput()<<"\nTriangulation is:\n";
    l_read(Triangulation);
    verboseOutput()<<"\nHilbert basis is:\n";
    l_read(Hilbert_Basis);
    verboseOutput()<<"\nHt1 elements are:\n";
    l_read(Ht1_Elements);
    verboseOutput()<<"\nHilbert Series  is:\n";
    verboseOutput()<<Hilbert_Series;
}

} //end namespace


