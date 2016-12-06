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

#include <stdlib.h>
#include <set>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <time.h>
#include <deque>

#include "libQnormaliz/Qfull_cone.h"
#include "libQnormaliz/Qcone_helper.h"
#include "libQnormaliz/Qvector_operations.h"
#include "libQnormaliz/Qlist_operations.h"
#include "libQnormaliz/Qmap_operations.h"
#include "libQnormaliz/Qmy_omp.h"
#include "libQnormaliz/Qinteger.h"
// #include "libQnormaliz/Qsublattice_representation.h"
// #include "libQnormaliz/Qoffload_handler.h"

//---------------------------------------------------------------------------

const size_t RecBoundTriang=1000000;   //  if number(supphyps)*size(triang) > RecBoundTriang
                                       // we pass to (non-recirsive) pyramids

const size_t EvalBoundTriang=2500000; // if more than EvalBoundTriang simplices have been stored
                               // evaluation is started (whenever possible)

const size_t EvalBoundPyr=200000;   // the same for stored pyramids of level > 0

const size_t EvalBoundLevel0Pyr=200000; // 1000000;   // the same for stored level 0 pyramids

// const size_t EvalBoundRecPyr=200000;   // the same for stored RECURSIVE pyramids

// const size_t IntermedRedBoundHB=2000000;  // bound for number of HB elements before 
                                              // intermediate reduction is called
                                              
const int largePyramidFactor=20;  // pyramid is large if largePyramidFactor*Comparisons[Pyramid_key.size()-dim] > old_nr_supp_hyps

const int SuppHypRecursionFactor=100; // pyramids for supphyps formed if Pos*Neg > this factor*dim^4

const size_t RAM_Size=1000000000; // we assume that there is at least 1 GB of RAM

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//private
//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::check_simpliciality_hyperplane(const FACETDATA& hyp) const{
    size_t nr_gen_in_hyp=0;
    for(size_t i=0; i<nr_gen;++i)
        if(in_triang[i]&& hyp.GenInHyp.test(i))
            nr_gen_in_hyp++;
    if((hyp.simplicial &&  nr_gen_in_hyp!=dim-2) || (!hyp.simplicial &&  nr_gen_in_hyp==dim-2)){
        // NOTE: in_triang set at END of main loop in build_cone
        cout << "Simplicial " << hyp.simplicial << " dim " << dim << " gen_in_hyp " << nr_gen_in_hyp << endl;
        assert(false);
    }
}

template<typename Number>
void Full_Cone<Number>::set_simplicial(FACETDATA& hyp){
    size_t nr_gen_in_hyp=0;
    for(size_t i=0; i<nr_gen;++i)
        if(in_triang[i]&& hyp.GenInHyp.test(i))
            nr_gen_in_hyp++;
    hyp.simplicial=(nr_gen_in_hyp==dim-2);
}

template<typename Number>
void Full_Cone<Number>::number_hyperplane(FACETDATA& hyp, const size_t born_at, const size_t mother){
// add identifying number, the birth day and the number of mother 

    hyp.Mother=mother;
    hyp.BornAt=born_at;
    if(!multithreaded_pyramid){
        hyp.Ident=HypCounter[0];
        HypCounter[0]++;
        return;
    }
    
    int tn;
    if(omp_get_level()==0)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(1);
    hyp.Ident=HypCounter[tn];
    HypCounter[tn]+=omp_get_max_threads();
    
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::is_hyperplane_included(FACETDATA& hyp) {
    if (!is_pyramid) { // in the topcone we always have ov_sp > 0
        return true;
    }
    //check if it would be an excluded hyperplane
    Number ov_sp = v_scalar_product(hyp.Hyp,Order_Vector);
    if (ov_sp > 0) {
        return true;
    } else if (ov_sp == 0) {
        for (size_t i=0; i<dim; i++) {
            if (hyp.Hyp[i]>0) {
                return true;
            } else if (hyp.Hyp[i]<0) {
                return false;
            }
        }
    }
    return false;
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::add_hyperplane(const size_t& new_generator, const FACETDATA & positive,const FACETDATA & negative,
                            list<FACETDATA>& NewHyps, bool known_to_be_simplicial){
// adds a new hyperplane found in find_new_facets to this cone (restricted to generators processed)

    size_t k;
    
    FACETDATA NewFacet; NewFacet.Hyp.resize(dim); NewFacet.GenInHyp.resize(nr_gen);        
    
    for (k = 0; k <dim; k++) {
        NewFacet.Hyp[k]=positive.ValNewGen*negative.Hyp[k]-negative.ValNewGen*positive.Hyp[k];
        if(!check_range(NewFacet.Hyp[k]))
            break;    
    }
    
    /* cout << "==========================================" << endl;
    cout << NewFacet.Hyp;
    cout << "==========================================" << endl; */
    
    v_simplify(NewFacet.Hyp);
    
    NewFacet.ValNewGen=0;    
    NewFacet.GenInHyp=positive.GenInHyp & negative.GenInHyp; // new hyperplane contains old gen iff both pos and neg do
    if(known_to_be_simplicial){
        NewFacet.simplicial=true;
        check_simpliciality_hyperplane(NewFacet);
    }
    else
        set_simplicial(NewFacet);
    NewFacet.GenInHyp.set(new_generator);  // new hyperplane contains new generator
    number_hyperplane(NewFacet,nrGensInCone,positive.Ident);
    
    NewHyps.push_back(NewFacet);
}


//---------------------------------------------------------------------------


template<typename Number>
void Full_Cone<Number>::find_new_facets(const size_t& new_generator){
// our Fourier-Motzkin implementation
// the special treatment of simplicial facets was inserted because of line shellings.
// At present these are not computed.

    //to see if possible to replace the function .end with constant iterator since push-back is performed.

    // for dimension 0 and 1 F-M is never necessary and can lead to problems
    // when using dim-2
    if (dim <= 1)
        return;

    // NEW: new_generator is the index of the generator being inserted

    size_t i,k,nr_zero_i;
    size_t subfacet_dim=dim-2; // NEW dimension of subfacet
    size_t facet_dim=dim-1; // NEW dimension of facet
    
    const bool tv_verbose = false; //verbose && !is_pyramid; // && Support_Hyperplanes.nr_of_rows()>10000; //verbose in this method call
    
        
    // preparing the computations, the various types of facets are sorted into the deques
    deque <FACETDATA*> Pos_Simp,Pos_Non_Simp;
    deque <FACETDATA*> Neg_Simp,Neg_Non_Simp;
    deque <FACETDATA*> Neutral_Simp, Neutral_Non_Simp;
    
    boost::dynamic_bitset<> Zero_Positive(nr_gen),Zero_Negative(nr_gen); // here we collect the vertices that lie in a
                                        // postive resp. negative hyperplane

    bool simplex;
    
    if (tv_verbose) verboseOutput()<<"transform_values:"<<flush;
    
    typename list<FACETDATA>::iterator ii = Facets.begin();
    
    for (; ii != Facets.end(); ++ii) {
        // simplex=true;
        // nr_zero_i=0;
        simplex=ii->simplicial; // at present simplicial, will become nonsimplicial if neutral
        /* for (size_t j=0; j<nr_gen; j++){
            if (ii->GenInHyp.test(j)) {
                if (++nr_zero_i > facet_dim) {
                    simplex=false;
                    break;
                }
            }
        }*/
        
        if (ii->ValNewGen==0) {
            ii->GenInHyp.set(new_generator);  // Must be set explicitly !!
            ii->simplicial=false;  // simpliciality definitly gone with the new generator
            if (simplex) {
                Neutral_Simp.push_back(&(*ii)); // simplicial without the new generator
            }   else {
                Neutral_Non_Simp.push_back(&(*ii)); // nonsim¸plicial already without the new generator
            }
        }
        else if (ii->ValNewGen>0) {
            Zero_Positive |= ii->GenInHyp;
            if (simplex) {
                Pos_Simp.push_back(&(*ii));
            } else {
                Pos_Non_Simp.push_back(&(*ii));
            }
        } 
        else if (ii->ValNewGen<0) {
            Zero_Negative |= ii->GenInHyp;
            if (simplex) {
                Neg_Simp.push_back(&(*ii));
            } else {
                Neg_Non_Simp.push_back(&(*ii));
            }
        }
    }
    
    // TO DO: Negativliste mit Zero_Positive verfeinern, also die aussondern, die nicht genug positive Erz enthalten
    // Eventuell sogar Rang-Test einbauen.
    // Letzteres k√∂nnte man auch bei den positiven machen, bevor sie verarbeitet werden
    
    boost::dynamic_bitset<> Zero_PN(nr_gen);
    Zero_PN = Zero_Positive & Zero_Negative;
    
    size_t nr_PosSimp  = Pos_Simp.size();
    size_t nr_PosNonSimp = Pos_Non_Simp.size();
    size_t nr_NegSimp  = Neg_Simp.size();
    size_t nr_NegNonSimp = Neg_Non_Simp.size();
    size_t nr_NeuSimp  = Neutral_Simp.size();
    size_t nr_NeuNonSimp = Neutral_Non_Simp.size();
    
    if (tv_verbose) verboseOutput()<<" PS "<<nr_PosSimp<<", P "<<nr_PosNonSimp<<", NS "<<nr_NegSimp<<", N "<<nr_NegNonSimp<<", ZS "<<nr_NeuSimp<<", Z "<<nr_NeuNonSimp<<endl;

    if (tv_verbose) verboseOutput()<<"transform_values: subfacet of NS: "<<flush;
    
    vector< list<pair < boost::dynamic_bitset<>, int> > > Neg_Subfacet_Multi(omp_get_max_threads()) ;

    boost::dynamic_bitset<> zero_i, subfacet;

    // This parallel region cannot throw a NormalizException
    #pragma omp parallel for private(zero_i,subfacet,k,nr_zero_i)
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
            
        if(nr_zero_i==facet_dim){
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


    if (tv_verbose) verboseOutput()<<Neg_Subfacet_Multi_United.size() << ", " << flush;

    list< pair < boost::dynamic_bitset<>, int > >::iterator jj;
    list< pair < boost::dynamic_bitset<>, int > >::iterator del;
    jj =Neg_Subfacet_Multi_United.begin();           // remove negative subfacets shared
    while (jj!= Neg_Subfacet_Multi_United.end()) {   // by two neg simpl facets
        del=jj++;
        if (jj!=Neg_Subfacet_Multi_United.end() && (*jj).first==(*del).first) {   //delete since is the intersection of two negative simplicies
            Neg_Subfacet_Multi_United.erase(del);
            del=jj++;
            Neg_Subfacet_Multi_United.erase(del);
        }
    }

    size_t nr_NegSubfMult = Neg_Subfacet_Multi_United.size();
    if (tv_verbose) verboseOutput() << nr_NegSubfMult << ", " << flush;
    
    vector<list<FACETDATA> > NewHypsSimp(nr_PosSimp);
    vector<list<FACETDATA> > NewHypsNonSimp(nr_PosNonSimp);

    map < boost::dynamic_bitset<>, int > Neg_Subfacet;
    size_t nr_NegSubf=0;
    
    // size_t NrMatches=0, NrCSF=0, NrRank=0, NrComp=0, NrNewF=0;
    
    /* deque<bool> Indi(nr_NegNonSimp);
    for(size_t j=0;j<nr_NegNonSimp;++j)
        Indi[j]=false; */
        
    if(multithreaded_pyramid){
        #pragma omp atomic
        nrTotalComparisons+=nr_NegNonSimp*nr_PosNonSimp;
    }
    else{
        nrTotalComparisons+=nr_NegNonSimp*nr_PosNonSimp; 
    } 

    
//=====================================================================
// parallel from here

    bool skip_remaining = false;
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif

    #pragma omp parallel private(jj)
    {
    size_t i,j,k,nr_zero_i;
    boost::dynamic_bitset<> subfacet(dim-2);
    jj = Neg_Subfacet_Multi_United.begin();
    size_t jjpos=0;
    int tn = omp_get_ancestor_thread_num(1);

    bool found;
    // This for region cannot throw a NormalizException
    #pragma omp for schedule(dynamic)
    for (size_t j=0; j<nr_NegSubfMult; ++j) {  // remove negative subfacets shared
        for(;j > jjpos; ++jjpos, ++jj) ;       // by non-simpl neg or neutral facets 
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
        if (found) {
            jj->second=-1;
        }
    }
    
    #pragma omp single
    { //remove elements that where found in the previous loop
    jj = Neg_Subfacet_Multi_United.begin();
    map < boost::dynamic_bitset<>, int > ::iterator last_inserted=Neg_Subfacet.begin(); // used to speedup insertion into the new map
    for (; jj!= Neg_Subfacet_Multi_United.end(); ++jj) {
        if ((*jj).second != -1) {
            last_inserted = Neg_Subfacet.insert(last_inserted,*jj);
        }
    }
    nr_NegSubf=Neg_Subfacet.size();
    }
    
    #pragma omp single nowait
    {Neg_Subfacet_Multi_United.clear();}

    
    boost::dynamic_bitset<> zero_i(nr_gen);
    map <boost::dynamic_bitset<>, int> ::iterator jj_map;

    
    #pragma omp single nowait
    if (tv_verbose) {
        verboseOutput() << "PS vs NS and PS vs N , " << flush;
    }

    vector<key_t> key(nr_gen);
    size_t nr_missing;
    bool common_subfacet;
    // we cannot use nowait here because of the way we handle exceptions in this loop
    #pragma omp for schedule(dynamic) //nowait
    for (size_t i =0; i<nr_PosSimp; i++){

        if (skip_remaining) continue;
#ifndef NCATCH
        try {
#endif
        zero_i=Zero_PN & Pos_Simp[i]->GenInHyp;
        nr_zero_i=0;
        for(j=0;j<nr_gen && nr_zero_i<=facet_dim;j++)
            if(zero_i.test(j)){
                key[nr_zero_i]=j;
                nr_zero_i++;
            } 
            
        if(nr_zero_i<subfacet_dim)
            continue;
            
        // first PS vs NS
        
        if (nr_zero_i==subfacet_dim) {                 // NEW slight change in logic. Positive simpl facet shared at most
            jj_map=Neg_Subfacet.find(zero_i);           // one subfacet with negative simpl facet
            if (jj_map!=Neg_Subfacet.end()) {
                add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Simp[(*jj_map).second],NewHypsSimp[i],true);
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
                        add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Simp[(*jj_map).second],NewHypsSimp[i],true);
                        (*jj_map).second = -1;
                        // Indi[j]=true;
                    }
                }
            }
        }

        // now PS vs N

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
               add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Non_Simp[j],NewHypsSimp[i],true);
               if(nr_zero_i==subfacet_dim) // only one subfacet can lie in negative hyperplane
                   break;
            }
       }
#ifndef NCATCH
       } catch(const std::exception& ) {
           tmp_exception = std::current_exception();
           skip_remaining = true;
           #pragma omp flush(skip_remaining)
       }
#endif

    } // PS vs NS and PS vs N

    if (!skip_remaining) {
    #pragma omp single nowait
    if (tv_verbose) {
        verboseOutput() << "P vs NS and P vs N" << endl;
    }

    list<FACETDATA*> AllNonSimpHyp;
    typename list<FACETDATA*>::iterator a;

    for(i=0;i<nr_PosNonSimp;++i)
        AllNonSimpHyp.push_back(&(*Pos_Non_Simp[i]));
    for(i=0;i<nr_NegNonSimp;++i)
        AllNonSimpHyp.push_back(&(*Neg_Non_Simp[i]));
    for(i=0;i<nr_NeuNonSimp;++i)
        AllNonSimpHyp.push_back(&(*Neutral_Non_Simp[i])); 
    size_t nr_NonSimp = nr_PosNonSimp+nr_NegNonSimp+nr_NeuNonSimp;
   
    bool ranktest;
    FACETDATA *hp_i, *hp_j, *hp_t; // pointers to current hyperplanes
    
    size_t missing_bound, nr_common_zero;
    boost::dynamic_bitset<> common_zero(nr_gen);
    vector<key_t> common_key;
    common_key.reserve(nr_gen);
    vector<int> key_start(nrGensInCone);
    
    #pragma omp for schedule(dynamic) // nowait
    for (size_t i =0; i<nr_PosNonSimp; i++){ //Positive Non Simp vs.Negative Simp and Non Simp

        if (skip_remaining) continue;

#ifndef NCATCH
        try {
#endif
        jj_map = Neg_Subfacet.begin();       // First the Simp
        for (j=0; j<nr_NegSubf; ++j,++jj_map) {
            if ( (*jj_map).second != -1 ) {  // skip used subfacets
                if(jj_map->first.is_subset_of(Pos_Non_Simp[i]->GenInHyp)){
                    add_hyperplane(new_generator,*Pos_Non_Simp[i],*Neg_Simp[(*jj_map).second],NewHypsNonSimp[i],true);
                    (*jj_map).second = -1; // has now been used
                }
            }
        }
        
        // Now the NonSimp

        hp_i=Pos_Non_Simp[i];
        zero_i=Zero_PN & hp_i->GenInHyp; // these are the potential vertices in an intersection
        nr_zero_i=0;
        int last_existing=-1;
        for(size_t jj=0;jj<nrGensInCone;jj++) // we make a "key" of the potential vertices in the intersection
        {
            j=GensInCone[jj];
            if(zero_i.test(j)){
                key[nr_zero_i]=j;
                for(size_t kk= last_existing+1;kk<=jj;kk++)  // used in the extension test
                    key_start[kk]=nr_zero_i;                 // to find out from which generator on both have existed
                nr_zero_i++;
                last_existing= jj;
            }
        }
        if(last_existing< (int)nrGensInCone-1)
            for(size_t kk=last_existing+1;kk<nrGensInCone;kk++)
                key_start[kk]=nr_zero_i;
                
        if (nr_zero_i<subfacet_dim) 
            continue;
        
        // now nr_zero_i is the number of vertices in hp_i that have a chance to lie in a negative facet
        // and key contains the indices
        
       missing_bound=nr_zero_i-subfacet_dim; // at most this number of generators can be missing
                                             // to have a chance for common subfacet                                            
       
       for (j=0; j<nr_NegNonSimp; j++){
    
        
           hp_j=Neg_Non_Simp[j];
           
           if(hp_i->Ident==hp_j->Mother || hp_j->Ident==hp_i->Mother){   // mother and daughter coming together
               add_hyperplane(new_generator,*hp_i,*hp_j,NewHypsNonSimp[i],false);  // their intersection is a subfacet
               continue;                                                           // simplicial set in add_hyperplane
           } 
           
           
           bool extension_test=hp_i->BornAt==hp_j->BornAt || (hp_i->BornAt<hp_j->BornAt && hp_j->Mother!=0)
                                                          || (hp_j->BornAt<hp_i->BornAt && hp_i->Mother!=0);
                                                          
           // extension_test=false;
                                                          
           size_t both_existing_from=key_start[max(hp_i->BornAt,hp_j->BornAt)];
                      
           nr_missing=0; 
           nr_common_zero=0;
           common_key.clear();
           size_t second_loop_bound=nr_zero_i;
           common_subfacet=true;
           
           // We use the following criterion:
           // if the two facets are not mother and daughter (taken care of already), then
           // they cannot have intersected in a subfacet at the time when the second was born.
           // In other words: they can only intersect in a subfacet now, if at least one common vertex
           // has been added after the birth of the younger one.
           // this is indicated by "extended".
           
           if(extension_test){
               bool extended=false;
               second_loop_bound=both_existing_from;  // fisrt we find the common vertices inserted from the step
                                                      // where both facets existed the first time
               for(k=both_existing_from;k<nr_zero_i;k++){
                   if(!hp_j->GenInHyp.test(key[k])) {
                       nr_missing++;
                       if(nr_missing>missing_bound) {
                           common_subfacet=false;
                           break;
                       }
                   }
                   else {
                       extended=true;  // in this case they have a common vertex added after their common existence
                       common_key.push_back(key[k]);
                       nr_common_zero++;
                   }
               }

               if(!extended || !common_subfacet) // 
                   continue;
           }
                    
           
           for(k=0;k<second_loop_bound;k++) {  // now the remaining 
               if(!hp_j->GenInHyp.test(key[k])) {
                   nr_missing++;
                   if(nr_missing>missing_bound) {
                       common_subfacet=false;
                       break;
                   }
               }
               else {
                   common_key.push_back(key[k]);
                   nr_common_zero++;
               }
            }
            
           if(!common_subfacet)
                continue;
           /* #pragma omp atomic
           NrCSF++;*/
           
           if(using_GMP<Number>())           
                ranktest = (nr_NonSimp > 10*dim*dim*nr_common_zero/3); // in this case the rank computation takes longer
           else
               ranktest = (nr_NonSimp > dim*dim*nr_common_zero/3);

           if(ranktest) {
               
            // cout << "Rangtest" << endl;
           
           /* #pragma omp atomic
            NrRank++; */
            
               Matrix<Number>& Test = Top_Cone->RankTest[tn];
               if (Test.rank_submatrix(Generators,common_key)<subfacet_dim) {
                   common_subfacet=false;
               }
           } // ranktest
           else{                 // now the comparison test
           
           /* #pragma omp atomic
            NrComp++; */
            
               common_zero = zero_i & hp_j->GenInHyp;
               for (a=AllNonSimpHyp.begin();a!=AllNonSimpHyp.end();++a){
                   hp_t=*a;
                   if ((hp_t!=hp_i) && (hp_t!=hp_j) && common_zero.is_subset_of(hp_t->GenInHyp)) {                                
                       common_subfacet=false;
                       AllNonSimpHyp.splice(AllNonSimpHyp.begin(),AllNonSimpHyp,a); // for the "darwinistic" mewthod
                       break;
                   }
               }                       
           } // else
           if (common_subfacet) {  //intersection of i and j is a subfacet
               add_hyperplane(new_generator,*hp_i,*hp_j,NewHypsNonSimp[i],false); //simplicial set in add_hyperplane
               /* #pragma omp atomic
                NrNewF++; */
                // Indi[j]=true;
           }
        }
#ifndef NCATCH
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
        }
#endif
    } // end for
    } // end !skip_remaining
    } //END parallel
    
#ifndef NCATCH
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
//=====================================================================
// parallel until here


    /* if(!is_pyramid)
      cout << "Matches " << NrMatches << " pot. common subf " << NrCSF << " rank test " <<  NrRank << " comp test "
        << NrComp << " neww hyps " << NrNewF << endl; */


    for(i=0;i<nr_PosSimp;i++)
        Facets.splice(Facets.end(),NewHypsSimp[i]);

    for(i=0;i<nr_PosNonSimp;i++)
        Facets.splice(Facets.end(),NewHypsNonSimp[i]);

    //removing the negative hyperplanes
    // now done in build_cone

    if (tv_verbose) verboseOutput()<<"transform_values: done"<<endl;
    
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::extend_triangulation(const size_t& new_generator){
// extends the triangulation of this cone by including new_generator
// simplicial facets save us from searching the "brother" in the existing triangulation
// to which the new simplex gets attached

    size_t listsize =old_nr_supp_hyps; // Facets.size();
    vector<typename list<FACETDATA>::iterator> visible;
    visible.reserve(listsize);
    typename list<FACETDATA>::iterator i = Facets.begin();

    listsize=0;
    for (; i!=Facets.end(); ++i) 
        if (i->ValNewGen < 0){ // visible facet
            visible.push_back(i);
            listsize++;
        }

#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif

    typename list< SHORTSIMPLEX<Number> >::iterator oldTriBack = --TriangulationBuffer.end();
    #pragma omp parallel private(i)
    {
    size_t k,l;
    bool one_not_in_i, not_in_facet;
    size_t not_in_i=0;
    // size_t facet_dim=dim-1;
    // size_t nr_in_i=0;

    list< SHORTSIMPLEX<Number> > Triangulation_kk;
    typename list< SHORTSIMPLEX<Number> >::iterator j;
    
    vector<key_t> key(dim);
    
    // if we only want a partial triangulation but came here because of a deep level
    // mark if this part of the triangulation has not to be evaluated
    bool skip_eval = false;

    #pragma omp for schedule(dynamic)
    for (size_t kk=0; kk<listsize; ++kk) {

#ifndef NCATCH
    try {
#endif
        i=visible[kk];
        
        /* nr_in_i=0;
        for(size_t m=0;m<nr_gen;m++){
            if(i->GenInHyp.test(m))
                nr_in_i++;
            if(nr_in_i>facet_dim){
                break;
            }
        }*/
        
        skip_eval = Top_Cone->do_partial_triangulation && i->ValNewGen == -1
                    && is_hyperplane_included(*i);

        if (i->simplicial){  // simplicial
            l=0;
            for (k = 0; k <nr_gen; k++) {
                if (i->GenInHyp[k]==1) {
                    key[l]=k;
                    l++;
                }
            }
            key[dim-1]=new_generator;
 
           if (skip_eval)
                store_key(key,0,0,Triangulation_kk);
            else
                store_key(key,-i->ValNewGen,0,Triangulation_kk);
            continue;
        } // end simplicial
        
        size_t irrelevant_vertices=0;
        for(size_t vertex=0;vertex<nrGensInCone;++vertex){
        
            if(i->GenInHyp[GensInCone[vertex]]==0) // lead vertex not in hyperplane
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
              if (skip_eval)
                  store_key(key,0,j->vol,Triangulation_kk);
              else
                  store_key(key,-i->ValNewGen,j->vol,Triangulation_kk);
                       
            } // j
            
        } // for vertex

#ifndef NCATCH
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
        }
#endif

    } // omp for kk

    if (multithreaded_pyramid) {
        #pragma omp critical(TRIANG)
        TriangulationBuffer.splice(TriangulationBuffer.end(),Triangulation_kk);
    } else
        TriangulationBuffer.splice(TriangulationBuffer.end(),Triangulation_kk);

    } // parallel

#ifndef NCATCH
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif

    // GensInCone.push_back(new_generator); // now in extend_cone
    TriSectionFirst.push_back(++oldTriBack);
    TriSectionLast.push_back(--TriangulationBuffer.end());
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::store_key(const vector<key_t>& key, const Number& height,
            const Number& mother_vol, list< SHORTSIMPLEX<Number> >& Triangulation){
// stores a simplex given by key and height in Triangulation
// mother_vol is the volume of the simplex to which the new one is attached

    SHORTSIMPLEX<Number> newsimplex;
    newsimplex.key=key;
    newsimplex.height=height;
    newsimplex.vol=0;
    
    if(multithreaded_pyramid){
        #pragma omp atomic
        TriangulationBufferSize++;
    }
    else {
        TriangulationBufferSize++;
    }
    int tn;
    if(omp_get_level()==0)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(1);
    
    if (height == 0) Top_Cone->triangulation_is_partial = true;
    
    if (keep_triangulation){
        Triangulation.push_back(newsimplex);
        return;  
    }
    
    bool Simpl_available=true;

    typename list< SHORTSIMPLEX<Number> >::iterator F;

    if(Top_Cone->FS[tn].empty()){
        if (Top_Cone->FreeSimpl.empty()) {
            Simpl_available=false;
        } else {
            #pragma omp critical(FREESIMPL)
            {
            if (Top_Cone->FreeSimpl.empty()) {
                Simpl_available=false;
            } else {
                // take 1000 simplices from FreeSimpl or what you can get
                F = Top_Cone->FreeSimpl.begin();
                size_t q;
                for (q = 0; q < 1000; ++q, ++F) {
                    if (F == Top_Cone->FreeSimpl.end())
                        break;
                }

                if(q<1000)
                    Top_Cone->FS[tn].splice(Top_Cone->FS[tn].begin(),
                        Top_Cone->FreeSimpl);
                else
                    Top_Cone->FS[tn].splice(Top_Cone->FS[tn].begin(),
                                  Top_Cone->FreeSimpl,Top_Cone->FreeSimpl.begin(),F);
            } // if empty global (critical)
            } // critical
        } // if empty global
    } // if empty thread

    if (Simpl_available) {
        Triangulation.splice(Triangulation.end(),Top_Cone->FS[tn],
                        Top_Cone->FS[tn].begin());
        Triangulation.back() = newsimplex;
    } else {
        Triangulation.push_back(newsimplex);
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::process_pyramids(const size_t new_generator,const bool recursive){

    /*

    We distinguish two types of pyramids:

    (i) recursive pyramids that give their support hyperplanes back to the mother.
    (ii) independent pyramids that are not linked to the mother.

    The parameter "recursive" indicates whether the pyramids that will be created
    in process_pyramid(s) are of type (i) or (ii).

    Every pyramid can create subpyramids of both types (not the case in version 2.8 - 2.10).

    Whether "this" is of type (i) or (ii) is indicated by do_all_hyperplanes.

    The creation of (sub)pyramids of type (i) can be blocked by setting recursion_allowed=false.
    (Not done in this version.)

    is_pyramid==false for the top_cone and ==true else.

    multithreaded_pyramid indicates whether parallelization takes place within the
    computation of a pyramid or whether it is computed in a single thread (defined in build_cone).

    Recursie pyramids are processed immediately after creation (as in 2.8). However, there
    are two exceptions:

    (a) In order to avoid very long waiting times for the computation of the "large" ones,
    these are treated as follows: the support hyperplanes of "this" coming from their bases
    (as negative hyperplanes of "this") are computed by matching them with the
    positive hyperplanes of "this". This Fourier-Motzkin step is much more
    efficient if a pyramid is large. For triangulation a large recursive
    pyramid is then stored as a pyramid of type (ii).

    (b) If "this" is processed in a parallelized loop calling process_pyramids, then
    the loop in process_pyramids cannot be interrupted for the evaluation of simplices. As a
    consequence an extremely long lst of simplices could arise if many small subpyramids of "this"
    are created in process_pyramids. In order to prevent this dangeous effect, small recursive
    subpyramids are stored for later triangulation if the simplex buffer has reached its
    size bound.

    Pyramids of type (ii) are stpred in Pyramids. The store_level of the created pyramids is 0
    for all pyramids created (possibly recursively) from the top cone. Pyramids created
    in evaluate_stored_pyramids get the store level for their subpyramids in that routine and
    transfer it to their recursive daughters. (correction March 4, 2015).

    Note: the top cone has pyr_level=-1. The pyr_level has no algorithmic relevance
    at present, but it shows the depth of the pyramid recursion at which the pyramid has been
    created.

    */


    size_t start_level=omp_get_level(); // allows us to check that we are on level 0
                                        // outside the loop and can therefore call evaluation
                                        // in order to empty the buffers
    vector<key_t> Pyramid_key;
    Pyramid_key.reserve(nr_gen);
    bool skip_triang; // make hyperplanes but skip triangulation (recursive pyramids only)

    deque<bool> done(old_nr_supp_hyps,false);
    bool skip_remaining;
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif
    typename list< FACETDATA >::iterator hyp;
    size_t nr_done=0;

    do{  // repeats processing until all hyperplanes have been processed

    hyp=Facets.begin();
    size_t hyppos=0;
    skip_remaining = false;
    
    const long VERBOSE_STEPS = 50;
    long step_x_size = old_nr_supp_hyps-VERBOSE_STEPS;
    const size_t RepBound=10000;

    #pragma omp parallel for private(skip_triang) firstprivate(hyppos,hyp,Pyramid_key) schedule(dynamic) reduction(+: nr_done)
    for (size_t kk=0; kk<old_nr_supp_hyps; ++kk) {

        if (skip_remaining) continue;
        
        if(verbose && old_nr_supp_hyps>=RepBound){
            #pragma omp critical(VERBOSE)
            while ((long)(kk*VERBOSE_STEPS) >= step_x_size) {
                step_x_size += old_nr_supp_hyps;
                verboseOutput() << "." <<flush;
            }
        }
        
#ifndef NCATCH
        try {
#endif
            for(;kk > hyppos; hyppos++, hyp++) ;
            for(;kk < hyppos; hyppos--, hyp--) ;

            if(done[hyppos])
                continue;

            done[hyppos]=true;

            nr_done++;

            if (hyp->ValNewGen == 0){                   // MUST BE SET HERE
                hyp->GenInHyp.set(new_generator);
                if(recursive) hyp->simplicial=false;                  // in the recursive case
            }

            if (hyp->ValNewGen >= 0) // facet not visible
                continue;

            skip_triang = false;
            if (Top_Cone->do_partial_triangulation && hyp->ValNewGen>=-1) { //ht1 criterion
                skip_triang = is_hyperplane_included(*hyp);
                if (skip_triang) {
                    Top_Cone->triangulation_is_partial = true;
                    if (!recursive) {
                        continue;
                    }
                }
            }

            Pyramid_key.clear(); // make data of new pyramid
            Pyramid_key.push_back(new_generator);
            for(size_t i=0;i<nr_gen;i++){
                if(in_triang[i] && hyp->GenInHyp.test(i)) {
                    Pyramid_key.push_back(i);
                }
            }

            // now we can store the new pyramid at the right place (or finish the simplicial ones)
            if (recursive && skip_triang) { // mark as "do not triangulate"
                process_pyramid(Pyramid_key, new_generator,store_level,0, recursive,hyp,start_level);
            } else { //default
                process_pyramid(Pyramid_key, new_generator,store_level,-hyp->ValNewGen, recursive,hyp,start_level);
            }
            // interrupt parallel execution if it is really parallel
            // to keep the triangulationand pyramid buffers under control
            if (start_level==0) {
                if (check_evaluation_buffer_size() || Top_Cone->check_pyr_buffer(store_level)) {
                    skip_remaining = true;
                }
            }

#ifndef NCATCH
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
        }
#endif
    } // end parallel loop over hyperplanes

#ifndef NCATCH
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
    
    if (start_level==0 && check_evaluation_buffer_size()) {
        Top_Cone->evaluate_triangulation();
    }
    
    if (start_level==0 && Top_Cone->check_pyr_buffer(store_level)) {
        Top_Cone->evaluate_stored_pyramids(store_level);
    }
    
    if (verbose && old_nr_supp_hyps>=RepBound)
        verboseOutput() << endl;

    } while (nr_done < old_nr_supp_hyps);
    
    
    evaluate_large_rec_pyramids(new_generator);

}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::process_pyramid(const vector<key_t>& Pyramid_key,
                          const size_t new_generator,const size_t store_level, Number height, const bool recursive,
                          typename list< FACETDATA >::iterator hyp, size_t start_level){
// processes simplicial pyramids directly, stores other pyramids into their depots

    #pragma omp atomic
    Top_Cone->totalNrPyr++;

    if(Pyramid_key.size()==dim){  // simplicial pyramid completely done here
        #pragma omp atomic        // only for saving memory
        Top_Cone->nrSimplicialPyr++;
        if(recursive){ // the facets may be facets of the mother cone and if recursive==true must be given back
            Matrix<Number> H(dim,dim);
            Number dummy_vol;
            Generators.simplex_data(Pyramid_key,H, dummy_vol,false);
            list<FACETDATA> NewFacets;
            FACETDATA NewFacet;
            NewFacet.GenInHyp.resize(nr_gen);
            for (size_t i=0; i<dim;i++) {
                NewFacet.Hyp = H[i];
                NewFacet.GenInHyp.set();
                NewFacet.GenInHyp.reset(i);
                NewFacet.simplicial=true;
                NewFacets.push_back(NewFacet);
            }
            select_supphyps_from(NewFacets,new_generator,Pyramid_key); // takes itself care of multithreaded_pyramid
        }
        if (height != 0 && (do_triangulation || do_partial_triangulation)) {
            if(multithreaded_pyramid) {
#ifndef NCATCH
                std::exception_ptr tmp_exception;
#endif
                #pragma omp critical(TRIANG)
                {
#ifndef NCATCH
                try{
#endif
                    store_key(Pyramid_key,height,0,TriangulationBuffer);
                    nrTotalComparisons+=dim*dim/2;
#ifndef NCATCH
                } catch(const std::exception& ) {
                    tmp_exception = std::current_exception();
                }
#endif
                } // end critical
#ifndef NCATCH
                if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
            } else {
                store_key(Pyramid_key,height,0,TriangulationBuffer);
                nrTotalComparisons+=dim*dim/2;
            }
        }
    }
    else {  // non-simplicial
    
        bool large=(largePyramidFactor*Comparisons[Pyramid_key.size()-dim] > old_nr_supp_hyps); // Pyramid_key.size()>largePyramidFactor*dim;
        
        if (!recursive || (large && (do_triangulation || do_partial_triangulation) && height!=0) ) {  // must also store for triangulation if recursive and large
            vector<key_t> key_wrt_top(Pyramid_key.size());
            for(size_t i=0;i<Pyramid_key.size();i++)
                key_wrt_top[i]=Top_Key[Pyramid_key[i]];
            #pragma omp critical(STOREPYRAMIDS)
            {
            //      cout << "store_level " << store_level << " large " << large << " pyr level " << pyr_level << endl;
            Top_Cone->Pyramids[store_level].push_back(key_wrt_top);
            Top_Cone->nrPyramids[store_level]++;
            } // critical
            if(!recursive)    // in this case we need only store for future triangulation, and that has been done
                return;
        }
        // now we are in the recursive case and must compute support hyperplanes of the subpyramid
        if(large){  // large recursive pyramid
            if(multithreaded_pyramid){
                #pragma omp critical(LARGERECPYRS)
                LargeRecPyrs.push_back(*hyp);  // LargeRecPyrs are kept and evaluated locally
            }
            else
                LargeRecPyrs.push_back(*hyp);
            return; // done with the large recusive pyramids
        }

        // only recursive small ones left

        Full_Cone<Number> Pyramid(*this,Pyramid_key);
        Pyramid.Mother = this;
        Pyramid.Mother_Key = Pyramid_key;    // need these data to give back supphyps
        Pyramid.apex=new_generator;
        if (height == 0) { //indicates "do not triangulate"
            Pyramid.do_triangulation = false;
            Pyramid.do_partial_triangulation = false;
            // Pyramid.do_Hilbert_basis = false;
            // Pyramid.do_deg1_elements=false;
        }

        bool store_for_triangulation=(store_level!=0) //loop in process_pyramids cannot be interrupted
            && (Pyramid.do_triangulation || Pyramid.do_partial_triangulation) // we must (partially) triangulate
            && (start_level!=0 && Top_Cone->TriangulationBufferSize > 2*EvalBoundTriang); // evaluation buffer already full  // EvalBoundTriang

        if (store_for_triangulation) {
            vector<key_t> key_wrt_top(Pyramid_key.size());
            for(size_t i=0;i<Pyramid_key.size();i++)
                key_wrt_top[i]=Top_Key[Pyramid_key[i]];
            #pragma omp critical(STOREPYRAMIDS)
            {
            Top_Cone->Pyramids[store_level].push_back(key_wrt_top);
            Top_Cone->nrPyramids[store_level]++;
            } // critical
            // Now we must suppress immediate triangulation
            Pyramid.do_triangulation = false;
            Pyramid.do_partial_triangulation = false;
            // Pyramid.do_Hilbert_basis = false;
            // Pyramid.do_deg1_elements=false;
        }

        Pyramid.build_cone();

        if(multithreaded_pyramid){
            #pragma omp atomic
            nrTotalComparisons+=Pyramid.nrTotalComparisons;
        } else
            nrTotalComparisons+=Pyramid.nrTotalComparisons;
    }  // else non-simplicial
}


//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::find_and_evaluate_start_simplex(){

    size_t i,j;
    Number factor;

    
    /* Simplex<Number> S = find_start_simplex();
    vector<key_t> key=S.read_key();   // generators indexed from 0 */
    
    vector<key_t> key=find_start_simplex();
    assert(key.size()==dim); // safety heck
    if(verbose){
        verboseOutput() << "Start simplex ";
        for(size_t i=0;i<key.size();++i)
            verboseOutput() <<  key[i]+1 << " ";
        verboseOutput() << endl;
    }
    Matrix<Number> H(dim,dim);
    Number vol;
    Generators.simplex_data(key,H,vol,do_partial_triangulation || do_triangulation);
    
    // H.pretty_print(cout);
    
        
    for (i = 0; i < dim; i++) {
        in_triang[key[i]]=true;
        GensInCone.push_back(key[i]);
        if (deg1_triangulation && isComputed(ConeProperty::Grading))
            deg1_triangulation = (gen_degrees[key[i]] == 1);
    }
    
    nrGensInCone=dim;
    
    nrTotalComparisons=dim*dim/2;
    Comparisons.push_back(nrTotalComparisons);
       
    for (i = 0; i <dim; i++) {
        FACETDATA NewFacet; NewFacet.GenInHyp.resize(nr_gen);
        NewFacet.Hyp=H[i];
        NewFacet.simplicial=true; // indeed, the start simplex is simplicial
        for(j=0;j < dim;j++)
            if(j!=i)
                NewFacet.GenInHyp.set(key[j]);
        NewFacet.ValNewGen=-1;         // must be taken negative since opposite facet
        number_hyperplane(NewFacet,0,0); // created with gen 0
        Facets.push_back(NewFacet);    // was visible before adding this vertex
    }
    
    if(!is_pyramid){
        //define Order_Vector, decides which facets of the simplices are excluded
        Order_Vector = vector<Number>(dim,0);
        // Matrix<Number> G=S.read_generators();
        for(i=0;i<dim;i++){
            factor=(unsigned long) (1+i%10);  // (2*(rand()%(2*dim))+3);
            for(j=0;j<dim;j++)
                Order_Vector[j]+=factor*Generators[key[i]][j];        
        }
    }

    //the volume is an upper bound for the height
    if(do_triangulation || (do_partial_triangulation && vol>1))
    {
        store_key(key,vol,1,TriangulationBuffer);
        if(do_only_multiplicity) {
            #pragma omp atomic
            TotDet++;
        }
    } else if (do_partial_triangulation) {
        triangulation_is_partial = true;
    }
    
    if(do_triangulation){ // we must prepare the sections of the triangulation
        for(i=0;i<dim;i++)
        {
            // GensInCone.push_back(key[i]); // now done in first loop since always needed
            TriSectionFirst.push_back(TriangulationBuffer.begin());
            TriSectionLast.push_back(TriangulationBuffer.begin());
        }
    }
    
}


//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::select_supphyps_from(const list<FACETDATA>& NewFacets, 
                    const size_t new_generator, const vector<key_t>& Pyramid_key){
// the mother cone (=this) selects supphyps from the list NewFacets supplied by the daughter
// the daughter provides the necessary information via the parameters

    size_t i;
    boost::dynamic_bitset<> in_Pyr(nr_gen);
    for (i=0; i<Pyramid_key.size(); i++) {
        in_Pyr.set(Pyramid_key[i]);
    }
    // the new generator is always the first in the pyramid
    assert(Pyramid_key[0] == new_generator);


    typename list<FACETDATA>::const_iterator pyr_hyp = NewFacets.begin();
    bool new_global_hyp;
    FACETDATA NewFacet;
    NewFacet.GenInHyp.resize(nr_gen);
    Number test;
    for (; pyr_hyp!=NewFacets.end(); ++pyr_hyp) {
        if(!pyr_hyp->GenInHyp.test(0)) // new gen not in hyp
            continue;
        new_global_hyp=true;
        for (i=0; i<nr_gen; ++i){
            if(in_Pyr.test(i) || !in_triang[i])
                continue;
            test=v_scalar_product(Generators[i],pyr_hyp->Hyp);
            if(test<=0){
                new_global_hyp=false;
                break;
            }

        }
        if(new_global_hyp){
            NewFacet.Hyp=pyr_hyp->Hyp;
            NewFacet.GenInHyp.reset();
            // size_t gens_in_facet=0;
            for (i=0; i<Pyramid_key.size(); ++i) {
                if (pyr_hyp->GenInHyp.test(i) && in_triang[Pyramid_key[i]]) {
                    NewFacet.GenInHyp.set(Pyramid_key[i]);
                    // gens_in_facet++;
                }
            }
            /* for (i=0; i<nr_gen; ++i) {
                if (NewFacet.GenInHyp.test(i) && in_triang[i]) {
                    gens_in_facet++;
                }
            }*/
            // gens_in_facet++; // Note: new generator not yet in in_triang
            NewFacet.GenInHyp.set(new_generator);
            NewFacet.simplicial=pyr_hyp->simplicial; // (gens_in_facet==dim-1); 
            check_simpliciality_hyperplane(NewFacet);
            number_hyperplane(NewFacet,nrGensInCone,0); //mother unknown
            if(multithreaded_pyramid){
                #pragma omp critical(GIVEBACKHYPS) 
                Facets.push_back(NewFacet);
            } else {
                Facets.push_back(NewFacet);
            }
        }
    }
}

//---------------------------------------------------------------------------
template<typename Number>
void Full_Cone<Number>::match_neg_hyp_with_pos_hyps(const FACETDATA& hyp, size_t new_generator,list<FACETDATA*>& PosHyps, boost::dynamic_bitset<>& Zero_P){

    size_t missing_bound, nr_common_zero;
    boost::dynamic_bitset<> common_zero(nr_gen);
    vector<key_t> common_key;
    common_key.reserve(nr_gen);
    vector<key_t> key(nr_gen);
    bool common_subfacet;
    list<FACETDATA> NewHyp;
    size_t subfacet_dim=dim-2;
    size_t nr_missing;
    typename list<FACETDATA*>::iterator a;
    list<FACETDATA> NewHyps;
    Matrix<Number> Test(0,dim);
    
    boost::dynamic_bitset<> zero_hyp=hyp.GenInHyp & Zero_P;  // we intersect with the set of gens in positive hyps
    
    size_t nr_zero_hyp=0;
    vector<int> key_start(nrGensInCone);
    size_t j;
    int last_existing=-1;
    for(size_t jj=0;jj<nrGensInCone;jj++)
    {
        j=GensInCone[jj];
        if(zero_hyp.test(j)){
            key[nr_zero_hyp]=j;
            for(size_t kk= last_existing+1;kk<=jj;kk++)
                key_start[kk]=nr_zero_hyp;
            nr_zero_hyp++;
            last_existing= jj;
        }
    }
    if(last_existing< (int)nrGensInCone-1)
        for(size_t kk=last_existing+1;kk<nrGensInCone;kk++)
            key_start[kk]=nr_zero_hyp;
            
    if (nr_zero_hyp<dim-2) 
        return;
    
    int tn = omp_get_ancestor_thread_num(1);
    missing_bound=nr_zero_hyp-subfacet_dim; // at most this number of generators can be missing
                                          // to have a chance for common subfacet
                                          
    typename list< FACETDATA*>::iterator hp_j_iterator=PosHyps.begin();
    
    FACETDATA* hp_j;

    for (;hp_j_iterator!=PosHyps.end();++hp_j_iterator){ //match hyp with the given Pos
        hp_j=*hp_j_iterator;


       if(hyp.Ident==hp_j->Mother || hp_j->Ident==hyp.Mother){   // mother and daughter coming together
                                            // their intersection is a subfacet
            add_hyperplane(new_generator,*hp_j,hyp,NewHyps,false);    // simplicial set in add_hyperplane
            continue;           
       }
       
       
       bool extension_test=hyp.BornAt==hp_j->BornAt || (hyp.BornAt<hp_j->BornAt && hp_j->Mother!=0)
                                                      || (hp_j->BornAt<hyp.BornAt && hyp.Mother!=0);
                                                      
       size_t both_existing_from=key_start[max(hyp.BornAt,hp_j->BornAt)];
                  
       nr_missing=0; 
       nr_common_zero=0;
       common_key.clear();
       size_t second_loop_bound=nr_zero_hyp;
       common_subfacet=true;  
       
       if(extension_test){
           bool extended=false;
           second_loop_bound=both_existing_from;
           for(size_t k=both_existing_from;k<nr_zero_hyp;k++){
               if(!hp_j->GenInHyp.test(key[k])) {
                   nr_missing++;
                   if(nr_missing>missing_bound) {
                       common_subfacet=false;
                       break;
                   }
               }
               else {
                   extended=true;
                   common_key.push_back(key[k]);
                   nr_common_zero++;
               }
           }

           if(!extended || !common_subfacet) // 
               continue;
       }
                
       for(size_t k=0;k<second_loop_bound;k++) {
           if(!hp_j->GenInHyp.test(key[k])) {
               nr_missing++;
               if(nr_missing>missing_bound) {
                   common_subfacet=false;
                   break;
               }
           }
           else {
               common_key.push_back(key[k]);
               nr_common_zero++;
           }
        }
        
       if(!common_subfacet)
            continue;
       
       assert(nr_common_zero >=subfacet_dim);
            
        // only rank test since we have many supphyps anyway
        if (!hp_j->simplicial){
            Matrix<Number>& Test = Top_Cone->RankTest[tn];
            if(Test.rank_submatrix(Generators,common_key)<subfacet_dim)
                common_subfacet=false;     // don't make a hyperplane
        }
        
        if(common_subfacet)
            add_hyperplane(new_generator,*hp_j,hyp,NewHyps,false);  // simplicial set in add_hyperplane
    } // for           

    if(multithreaded_pyramid)
        #pragma omp critical(GIVEBACKHYPS)
        Facets.splice(Facets.end(),NewHyps);
    else
        Facets.splice(Facets.end(),NewHyps);

}

//---------------------------------------------------------------------------
template<typename Number>
void Full_Cone<Number>::collect_pos_supphyps(list<FACETDATA*>& PosHyps, boost::dynamic_bitset<>& Zero_P, size_t& nr_pos){
           
    // positive facets are collected in a list
    
    typename list<FACETDATA>::iterator ii = Facets.begin();
    nr_pos=0;
    
    for (size_t ij=0; ij< old_nr_supp_hyps; ++ij, ++ii)
        if (ii->ValNewGen>0) {
            Zero_P |= ii->GenInHyp;
            PosHyps.push_back(&(*ii));
            nr_pos++;
        }
}

//---------------------------------------------------------------------------
template<typename Number>
void Full_Cone<Number>::evaluate_large_rec_pyramids(size_t new_generator){
    
    size_t nrLargeRecPyrs=LargeRecPyrs.size();
    if(nrLargeRecPyrs==0)
        return;
        
    if(verbose)
        verboseOutput() << "large pyramids " << nrLargeRecPyrs << endl;
    
    list<FACETDATA*> PosHyps;
    boost::dynamic_bitset<> Zero_P(nr_gen);
    size_t nr_pos;
    collect_pos_supphyps(PosHyps,Zero_P,nr_pos);
    
    nrTotalComparisons+=nr_pos*nrLargeRecPyrs;
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif
    
    const long VERBOSE_STEPS = 50;
    long step_x_size = nrLargeRecPyrs-VERBOSE_STEPS;
    const size_t RepBound=100;
    
    #pragma omp parallel
    {
    size_t ppos=0;
    typename list<FACETDATA>::iterator p=LargeRecPyrs.begin(); 
    
    #pragma omp for schedule(dynamic) 
    for(size_t i=0; i<nrLargeRecPyrs; i++){
        for(; i > ppos; ++ppos, ++p) ;
        for(; i < ppos; --ppos, --p) {};

        if(verbose && nrLargeRecPyrs>=RepBound){
            #pragma omp critical(VERBOSE)
            while ((long)(i*VERBOSE_STEPS) >= step_x_size) {
                step_x_size += nrLargeRecPyrs;
                verboseOutput() << "." <<flush;
            }
        }
        
#ifndef NCATCH
        try {
#endif
            match_neg_hyp_with_pos_hyps(*p,new_generator,PosHyps,Zero_P);
#ifndef NCATCH
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
        }
#endif
    }
    } // parallel
#ifndef NCATCH
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
    
    if(verbose && nrLargeRecPyrs>=RepBound)
        verboseOutput() << endl;

    LargeRecPyrs.clear();
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::check_pyr_buffer(const size_t level){
    if(level==0)
        return(nrPyramids[0] > EvalBoundLevel0Pyr);
    else
        return(nrPyramids[level] > EvalBoundPyr);
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::evaluate_stored_pyramids(const size_t level){
// evaluates the stored non-recursive pyramids

    assert(!omp_in_parallel());

    if(Pyramids[level].empty())
        return;
    if (Pyramids.size() < level+2) {
        Pyramids.resize(level+2);      // provide space for a new generation
        nrPyramids.resize(level+2, 0);
    }

    size_t eval_down_to = 0;

#ifdef NMZ_MIC_OFFLOAD
#ifndef __MIC__
    // only on host and if offload is available
    if (level == 0 && nrPyramids[0] > EvalBoundLevel0Pyr) {
        eval_down_to = EvalBoundLevel0Pyr;
    }
#endif
#endif

    vector<char> Done(nrPyramids[level],0);
    if (verbose) {
        verboseOutput() << "**************************************************" << endl;

        for (size_t l=0; l<=level; ++l) {
            if (nrPyramids[l]>0) {
                verboseOutput() << "level " << l << " pyramids remaining: "
                                << nrPyramids[l] << endl;
            }
        }
        verboseOutput() << "**************************************************" << endl;
    }
    typename list<vector<key_t> >::iterator p;
    size_t ppos;
    bool skip_remaining;
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif

    while (nrPyramids[level] > eval_down_to) {

       p = Pyramids[level].begin();
       ppos=0;
       skip_remaining = false;
    
       #pragma omp parallel for firstprivate(p,ppos) schedule(dynamic) 
       for(size_t i=0; i<nrPyramids[level]; i++){
           if (skip_remaining)
                continue;
           for(; i > ppos; ++ppos, ++p) ;
           for(; i < ppos; --ppos, --p) ;
           
           if(Done[i])
               continue;
           Done[i]=1;

#ifndef NCATCH
           try {
#endif
               Full_Cone<Number> Pyramid(*this,*p);
               // Pyramid.recursion_allowed=false;
               Pyramid.do_all_hyperplanes=false;
               if (level>=2 && do_partial_triangulation){ // limits the descent of do_partial_triangulation
                   Pyramid.do_triangulation=true;
                   Pyramid.do_partial_triangulation=false;
               }
               Pyramid.store_level=level+1;
               Pyramid.build_cone();
               if (check_evaluation_buffer_size() || Top_Cone->check_pyr_buffer(level+1)) {
                   // interrupt parallel execution to keep the buffer under control
                   skip_remaining = true;
               }
#ifndef NCATCH
           } catch(const std::exception& ) {
               tmp_exception = std::current_exception();
               skip_remaining = true;
               #pragma omp flush(skip_remaining)
           }
#endif
        } //end parallel for
#ifndef NCATCH
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif

        // remove done pyramids
        p = Pyramids[level].begin();
        for(size_t i=0; p != Pyramids[level].end(); i++){
            if (Done[i]) {
                p=Pyramids[level].erase(p);
                nrPyramids[level]--;
                Done[i]=0;
            } else {
                ++p;
            }
        }

        if (check_evaluation_buffer_size()) {
            if (verbose)
                verboseOutput() << nrPyramids[level] <<
                    " pyramids remaining on level " << level << ", ";
            Top_Cone->evaluate_triangulation();
        }

        if (Top_Cone->check_pyr_buffer(level+1)) {
            evaluate_stored_pyramids(level+1);
        }
    
    } //end while (nrPyramids[level] > 0)
     
    if (verbose) {
        verboseOutput() << "**************************************************" << endl;
        verboseOutput() << "all pyramids on level "<< level << " done!"<<endl;
        if (nrPyramids[level+1] == 0) {
            for (size_t l=0; l<=level; ++l) {
                if (nrPyramids[l]>0) {
                    verboseOutput() << "level " << l << " pyramids remaining: "
                                    << nrPyramids[l] << endl;
                }
            }
            verboseOutput() << "**************************************************" << endl;
        }
    }
    if(check_evaluation_buffer())
    {
        Top_Cone->evaluate_triangulation();
    }
     
    evaluate_stored_pyramids(level+1);
}
    


//---------------------------------------------------------------------------

/* builds the cone successively by inserting generators */
template<typename Number>
void Full_Cone<Number>::build_cone() {
    // if(dim>0){            //correction needed to include the 0 cone;
    
    // cout << "Pyr " << pyr_level << endl;

    long long RecBoundSuppHyp = dim*dim;
    RecBoundSuppHyp *= RecBoundSuppHyp*SuppHypRecursionFactor; //dim^4 * 3000
    
    tri_recursion=false; 
    
    multithreaded_pyramid=(omp_get_level()==0);
    
    if(!use_existing_facets){
        if(multithreaded_pyramid){
            HypCounter.resize(omp_get_max_threads());
            for(size_t i=0;i<HypCounter.size();++i)
                HypCounter[i]=i+1;
        } else{
            HypCounter.resize(1);
            HypCounter[0]=1;    
        }
        
        find_and_evaluate_start_simplex();
    }
    
    size_t last_to_be_inserted; // good to know in case of do_all_hyperplanes==false
    last_to_be_inserted=nr_gen-1;  // because we don't need to compute support hyperplanes in this case 
    for(int j=nr_gen-1;j>=0;--j){
        if(!in_triang[j]){
            last_to_be_inserted=j;
            break;
        }
    } // last_to_be_inserted now determined
    
    bool is_new_generator;
    typename list< FACETDATA >::iterator l;


    for (size_t i=start_from;i<nr_gen;++i) { 
    
        start_from=i;
    
        if (in_triang[i])
            continue;
            
        if(do_triangulation && TriangulationBufferSize > 2*RecBoundTriang) // emermergency brake
            tri_recursion=true;               // to switch off production of simplices in favor
                                              // of non-recursive pyramids
        Number scalar_product;                                              
        is_new_generator=false;
        l=Facets.begin();
        old_nr_supp_hyps=Facets.size(); // Facets will be xtended in the loop 

        long long nr_pos=0, nr_neg=0;
        long long nr_neg_simp=0, nr_pos_simp=0;
        vector<Number> L;           
#ifndef NCATCH
        std::exception_ptr tmp_exception;
#endif
        
        size_t lpos=0;
        #pragma omp parallel for private(L,scalar_product) firstprivate(lpos,l) reduction(+: nr_pos, nr_neg)
        for (size_t k=0; k<old_nr_supp_hyps; k++) {
#ifndef NCATCH
            try {
#endif
                for(;k > lpos; lpos++, l++) ;
                for(;k < lpos; lpos--, l--) ;

                L=Generators[i];
                scalar_product=v_scalar_product(L,(*l).Hyp);
                l->ValNewGen=scalar_product;
                if (scalar_product<0) {
                    is_new_generator=true;
                    nr_neg++;
                    if(l->simplicial)
                        #pragma omp atomic
                        nr_neg_simp++;
                }
                if (scalar_product>0) {
                    nr_pos++;
                    if(l->simplicial)
                        #pragma omp atomic
                        nr_pos_simp++;
                }
#ifndef NCATCH
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
            }
#endif
        }  //end parallel for
#ifndef NCATCH
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif

        if(!is_new_generator)
            continue;

        // the i-th generator is used in the triangulation
        // in_triang[i]=true; // now at end of loop
        if (deg1_triangulation && isComputed(ConeProperty::Grading))
            deg1_triangulation = (gen_degrees[i] == 1);
        
        /* if(!is_pyramid && verbose ) 
            verboseOutput() << "Neg " << nr_neg << " Pos " << nr_pos << " NegSimp " <<nr_neg_simp << " PosSimp " <<nr_pos_simp << endl;*/
        // First we test whether to go to recursive pyramids because of too many supphyps
        if (recursion_allowed && nr_neg*nr_pos-(nr_neg_simp*nr_pos_simp) > RecBoundSuppHyp) {  // use pyramids because of supphyps
            /* if(!is_pyramid && verbose )
                verboseOutput() << "Building pyramids" << endl; */
            if (do_triangulation)
                tri_recursion = true; // We can not go back to classical triangulation
            if(check_evaluation_buffer()){
                Top_Cone->evaluate_triangulation();
            }

            process_pyramids(i,true); //recursive
            lastGen=i;
            nextGen=i+1; 
        }
        else{ // now we check whether to go to pyramids because of the size of triangulation
              // once we have done so, we must stay with it
            if( tri_recursion || (do_triangulation 
                && (nr_neg*TriangulationBufferSize > RecBoundTriang
                    || 3*omp_get_max_threads()*TriangulationBufferSize>EvalBoundTriang ))){ // go to pyramids because of triangulation
                if(check_evaluation_buffer()){
                    Top_Cone->evaluate_triangulation();
                }
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
                if (l->ValNewGen<0) {
                    l=Facets.erase(l);
                }
                else
                    ++l;
            }
        }
        
        GensInCone.push_back(i);
        nrGensInCone++;
        
        Comparisons.push_back(nrTotalComparisons);
        
        if(verbose) {
            verboseOutput() << "gen="<< i+1 <<", ";
            if (do_all_hyperplanes || i!=last_to_be_inserted) {
                verboseOutput() << Facets.size()<<" hyp";
            } else {
                verboseOutput() << Support_Hyperplanes.nr_of_rows()<<" hyp";
            }
            if(nrPyramids[0]>0)
                verboseOutput() << ", " << nrPyramids[0] << " pyr"; 
            if(do_triangulation||do_partial_triangulation)
                verboseOutput() << ", " << TriangulationBufferSize << " simpl";
            verboseOutput()<< endl;
        }
        
        in_triang[i]=true;
        
    }  // loop over i
    
    start_from=nr_gen;
    
    if (is_pyramid && do_all_hyperplanes)  // must give supphyps back to mother
        Mother->select_supphyps_from(Facets, apex, Mother_Key);
    
    // transfer Facets --> SupportHyperplanes
    if (do_all_hyperplanes) {
        nrSupport_Hyperplanes = Facets.size();
        Support_Hyperplanes = Matrix<Number>(nrSupport_Hyperplanes,0);
        typename list<FACETDATA>::iterator IHV=Facets.begin();
        for (size_t i=0; i<nrSupport_Hyperplanes; ++i, ++IHV) {
            swap(Support_Hyperplanes[i],IHV->Hyp);
        }
        is_Computed.set(ConeProperty::SupportHyperplanes);
    } 
    Support_Hyperplanes.set_nr_of_columns(dim);
   
    
    if(do_extreme_rays && do_all_hyperplanes)
        compute_extreme_rays(true);
    
    transfer_triangulation_to_top(); // transfer remaining simplices to top
    if(check_evaluation_buffer()){
        Top_Cone->evaluate_triangulation();
    }  

    // } // end if (dim>0)
    
    Facets.clear(); 

}

template<typename Number>
void Full_Cone<Number>::start_message() {
    
       if (verbose) {
        verboseOutput()<<"************************************************************"<<endl;
        verboseOutput()<<"starting primal algorithm ";
        if (do_partial_triangulation) verboseOutput()<<"with partial triangulation ";
        if (do_triangulation) {
            verboseOutput()<<"with full triangulation ";
        }
        if (!do_triangulation && !do_partial_triangulation) verboseOutput()<<"(only support hyperplanes) ";
        verboseOutput()<<"..."<<endl;
    }   
}

template<typename Number>
void Full_Cone<Number>::end_message() {
    
       if (verbose) {
        verboseOutput() << "------------------------------------------------------------"<<endl;
    }   
}


//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::build_top_cone() {
    
    if(dim==0)
        return;
 
    // if( ( !do_bottom_dec || deg1_generated || dim==1 || (!do_triangulation && !do_partial_triangulation))) {        
        build_cone();
    // }

    evaluate_stored_pyramids(0);  // force evaluation of remaining pyramids
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::check_evaluation_buffer(){

    return(omp_get_level()==0 && check_evaluation_buffer_size());
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::check_evaluation_buffer_size(){

    return(!Top_Cone->keep_triangulation && 
               Top_Cone->TriangulationBufferSize > EvalBoundTriang);
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::transfer_triangulation_to_top(){  // NEW EVA

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
    int tn = 0;
    if (omp_in_parallel())
        tn = omp_get_ancestor_thread_num(1);
  
    typename list< SHORTSIMPLEX<Number> >::iterator pyr_simp=TriangulationBuffer.begin();
    while (pyr_simp!=TriangulationBuffer.end()) {
        if (pyr_simp->height == 0) { // it was marked to be skipped
            Top_Cone->FS[tn].splice(Top_Cone->FS[tn].end(), TriangulationBuffer, pyr_simp++);
            --TriangulationBufferSize;
        } else {
            for (i=0; i<dim; i++)  // adjust key to topcone generators
                pyr_simp->key[i]=Top_Key[pyr_simp->key[i]];
            ++pyr_simp;
        }
    }

    // cout << "Keys transferred " << endl;
    #pragma omp critical(TRIANG)
    {
        Top_Cone->TriangulationBuffer.splice(Top_Cone->TriangulationBuffer.end(),TriangulationBuffer);
        Top_Cone->TriangulationBufferSize += TriangulationBufferSize;
    }
    TriangulationBufferSize = 0;
  
}

//---------------------------------------------------------------------------
template<typename Number>
void Full_Cone<Number>::get_supphyps_from_copy(bool from_scratch){

    if(isComputed(ConeProperty::SupportHyperplanes)) // we have them already
        return;
    
    Full_Cone copy((*this).Generators);
    copy.verbose=verbose;
    if(!from_scratch){
        copy.start_from=start_from;
        copy.use_existing_facets=true;
        copy.keep_order=true;
        copy.HypCounter=HypCounter;
        copy.Extreme_Rays_Ind=Extreme_Rays_Ind;
        copy.in_triang=in_triang;
        copy.old_nr_supp_hyps=old_nr_supp_hyps;
        if(isComputed(ConeProperty::ExtremeRays))
            copy.is_Computed.set(ConeProperty::ExtremeRays);
        copy.GensInCone=GensInCone;
        copy.nrGensInCone=nrGensInCone;
        copy.Comparisons=Comparisons;
        if(!Comparisons.empty())
            copy.nrTotalComparisons=Comparisons[Comparisons.size()-1];
        
        typename list< FACETDATA >::const_iterator l=Facets.begin();
        
        for(size_t i=0;i<old_nr_supp_hyps;++i){
            copy.Facets.push_back(*l);
            ++l;
        }
    }
    
    copy.dualize_cone();
    
    std::swap(Support_Hyperplanes,copy.Support_Hyperplanes);
    nrSupport_Hyperplanes = copy.nrSupport_Hyperplanes;
    is_Computed.set(ConeProperty::SupportHyperplanes);
    do_all_hyperplanes = false;
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::evaluate_triangulation(){

    assert(omp_get_level()==0);
    
    if (TriangulationBufferSize == 0)
        return;
    
    if (keep_triangulation) {
        Triangulation.splice(Triangulation.end(),TriangulationBuffer);
    } else {
        // #pragma omp critical(FREESIMPL)
        FreeSimpl.splice(FreeSimpl.begin(),TriangulationBuffer);
    }
    TriangulationBufferSize=0;

}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::primal_algorithm(){

    primal_algorithm_initialize();

    /***** Main Work is done in build_top_cone() *****/
    build_top_cone();  // evaluates if keep_triangulation==false
    /***** Main Work is done in build_top_cone() *****/

    check_pointed();
    if(!pointed){
        throw NonpointedException();
    }

    primal_algorithm_finalize();
    primal_algorithm_set_computed();
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::primal_algorithm_initialize() {

}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::primal_algorithm_finalize() {

    if (keep_triangulation) {
        is_Computed.set(ConeProperty::Triangulation);
    }
    if (do_cone_dec) {
        is_Computed.set(ConeProperty::ConeDecomposition);
    }

    evaluate_triangulation();
    FreeSimpl.clear();
    
    if(verbose) {
        verboseOutput() << "Total number of pyramids = "<< totalNrPyr << ", among them simplicial " << nrSimplicialPyr << endl;
        // cout << "Uni "<< Unimod << " Ht1NonUni " << Ht1NonUni << " NonDecided " << NonDecided << " TotNonDec " << NonDecidedHyp<< endl;
    }

}


//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::primal_algorithm_set_computed() {

    extreme_rays_and_deg1_check();
    if(!pointed){
        throw NonpointedException();
    }

    if (do_triangulation || do_partial_triangulation) {
        is_Computed.set(ConeProperty::TriangulationSize,true);
    }    
}

   
//---------------------------------------------------------------------------
// Normaliz modes (public)
//---------------------------------------------------------------------------

// check the do_* bools, they must be set in advance
// this method (de)activate them according to dependencies between them
template<typename Number>
void Full_Cone<Number>::do_vars_check(bool with_default) {

    do_extreme_rays=true; // we always want to do this if compute() is called

    /* if (do_default_mode && with_default) {
        do_Hilbert_basis = true;
        do_h_vector = true;
        if(!inhomogeneous)
            do_class_group=true;
    }
    */
    
    if (do_integrally_closed) {
        if (do_Hilbert_basis) {
            do_integrally_closed = false; // don't interrupt the computation
        } else {
            do_Hilbert_basis = true;
        }
    }

    // activate implications
    if (do_module_gens_intcl) do_Hilbert_basis= true;
    if (do_module_gens_intcl) use_bottom_points= false;
    //if (do_hsop)            do_Hilbert_basis = true;
    if (do_Stanley_dec)     keep_triangulation = true;
    if (do_cone_dec)        keep_triangulation = true;
    if (keep_triangulation) do_determinants = true;
    if (do_multiplicity)    do_determinants = true;
    if ((do_multiplicity || do_h_vector) && inhomogeneous)    do_module_rank = true;
    if (do_determinants)    do_triangulation = true;
    if (do_h_vector && (with_default || explicit_h_vector))        do_triangulation = true;
    if (do_deg1_elements)   do_partial_triangulation = true;
    if (do_Hilbert_basis)   do_partial_triangulation = true;
    // activate 
    do_only_multiplicity = do_determinants;
    stop_after_cone_dec = true;
    if(do_cone_dec)          do_only_multiplicity=false;
        
    if (do_Stanley_dec || do_h_vector || do_deg1_elements 
                     || do_Hilbert_basis) {
        do_only_multiplicity = false;
        stop_after_cone_dec = false;
        do_evaluation = true;
    }
    if (do_determinants)    do_evaluation = true;

    // deactivate
    if (do_triangulation)   do_partial_triangulation = false;
    if (do_Hilbert_basis)   do_deg1_elements = false; // they will be extracted later
}


// general purpose compute method
// do_* bools must be set in advance, this method does sanity checks for it
// if no bool is set it does support hyperplanes and extreme rays
template<typename Number>
void Full_Cone<Number>::compute() {
    
    if(dim==0){
        set_zero_cone();
        return;
    }
    

    do_vars_check(false);
    explicit_full_triang=do_triangulation; // to distinguish it from do_triangulation via default mode
    if(do_default_mode)
        do_vars_check(true);

    start_message();
    
    if(Support_Hyperplanes.nr_of_rows()==0 && !do_Hilbert_basis && !do_h_vector && !do_multiplicity && !do_deg1_elements
        && !do_Stanley_dec && !do_triangulation && !do_determinants)
        assert(Generators.max_rank_submatrix_lex().size() == dim);

    minimize_support_hyperplanes(); // if they are given
    if (inhomogeneous)
        set_levels();

    if ((!do_triangulation && !do_partial_triangulation)
            || (Grading.size()>0 && !isComputed(ConeProperty::Grading))){
            // in the second case there are only two possibilities:
            // either nonpointed or bad grading
        do_triangulation=false;
        do_partial_triangulation=false;
        support_hyperplanes();
    }
    else{
        if(isComputed(ConeProperty::IsPointed) && !pointed){
            end_message();
            return;
        }

        sort_gens_by_degree(true);

        primal_algorithm();        
    }  
    
    end_message();
}




// -s
template<typename Number>
void Full_Cone<Number>::support_hyperplanes() { 
    if(!isComputed(ConeProperty::SupportHyperplanes)){
        sort_gens_by_degree(false); // we do not want to triangulate here
        build_top_cone();           
    }
    extreme_rays_and_deg1_check();
    if(inhomogeneous){
        find_level0_dim();
    }
}

//---------------------------------------------------------------------------
// Checks and auxiliary algorithms
//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::extreme_rays_and_deg1_check() {
    check_pointed();
    if(!pointed){
        throw NonpointedException();
    }
    //cout << "Generators" << endl;
    //Generators.pretty_print(cout);
    //cout << "SupportHyperplanes" << endl;
    //Support_Hyperplanes.pretty_print(cout);
    compute_extreme_rays();
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::find_level0_dim(){

    if(!isComputed(ConeProperty::Generators)){
        throw FatalException("Missing Generators.");
    }
    
    Matrix<Number> Help(nr_gen,dim);
    for(size_t i=0; i<nr_gen;++i)
        if(gen_levels[i]==0)
            Help[i]=Generators[i];
        
    ProjToLevel0Quot=Help.kernel();
    
    level0_dim=dim-ProjToLevel0Quot.nr_of_rows();
    is_Computed.set(ConeProperty::RecessionRank);
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::set_levels() {
    if(inhomogeneous && Truncation.size()!=dim){
        throw FatalException("Truncation not defined in inhomogeneous case.");
    }    
    
    // cout <<"trunc " << Truncation;

    if(gen_levels.size()!=nr_gen) // now we compute the levels
    {
        gen_levels.resize(nr_gen);
        vector<Number> gen_levels_Number=Generators.MxV(Truncation);
        for (size_t i=0; i<nr_gen; i++) {
            if (gen_levels_Number[i] < 0) {
                throw FatalException("Truncation gives non-positive value "
                        + toString(gen_levels_Number[i]) + " for generator "
                        + toString(i+1) + ".");
            }
            convert(gen_levels[i], gen_levels_Number[i]);
            // cout << "Gen " << Generators[i];
            // cout << "level " << gen_levels[i] << endl << "----------------------" << endl;
        }
    }
    
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::sort_gens_by_degree(bool triangulate) {
    // if(deg1_extreme_rays)  // gen_degrees.size()==0 || 
    // return;
    
    if(keep_order)
        return;
    
    Matrix<Number> Weights(0,dim);
    vector<bool> absolute;
    if(triangulate){
            Weights.append(vector<Number>(dim,1));
            absolute.push_back(true);
    }
    
    vector<key_t> perm=Generators.perm_by_weights(Weights,absolute);
    Generators.order_rows_by_perm(perm);
    order_by_perm(Extreme_Rays_Ind,perm);
    if(inhomogeneous && gen_levels.size()==nr_gen)
        order_by_perm(gen_levels,perm);
    compose_perm_gens(perm);
    
    if (verbose) {
        if(triangulate){
            if(isComputed(ConeProperty::Grading)){
                verboseOutput() <<"Generators sorted by degree and lexicographically" << endl;
                verboseOutput() << "Generators per degree:" << endl;
                verboseOutput() << count_in_map<long,long>(gen_degrees);
            }
            else
                verboseOutput() << "Generators sorted by 1-norm and lexicographically" << endl;
        }
        else{
            verboseOutput() << "Generators sorted lexicographically" << endl;
        }
    }
    keep_order=true;
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::compose_perm_gens(const vector<key_t>& perm) {
    order_by_perm(PermGens,perm);
}

//---------------------------------------------------------------------------

// an alternative to compute() for the basic tasks that need no triangulation
template<typename Number>
void Full_Cone<Number>::dualize_cone(bool print_message){
    
    if(dim==0){
        set_zero_cone();
        return;
    }

    // DO NOT CALL do_vars_check!!

    bool save_tri      = do_triangulation;
    bool save_part_tri = do_partial_triangulation;
    do_triangulation         = false;
    do_partial_triangulation = false;
    
    if(print_message) start_message();
    
    sort_gens_by_degree(false);
    
    if(!isComputed(ConeProperty::SupportHyperplanes))
        build_top_cone();
    
    if(do_pointed)
        check_pointed();

    if(do_extreme_rays) // in case we have known the support hyperplanes
        compute_extreme_rays();

    do_triangulation         = save_tri;
    do_partial_triangulation = save_part_tri;
    if(print_message) end_message();
}

//---------------------------------------------------------------------------

template<typename Number>
vector<key_t> Full_Cone<Number>::find_start_simplex() const {
        return Generators.max_rank_submatrix_lex();
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Full_Cone<Number>::select_matrix_from_list(const list<vector<Number> >& S,
                                   vector<size_t>& selection){

    sort(selection.begin(),selection.end());
    assert(selection.back()<S.size());
    size_t i=0,j=0;
    size_t k=selection.size();
    Matrix<Number> M(selection.size(),S.front().size());
    typename list<vector<Number> >::const_iterator ll=S.begin();
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

template<typename Number>

void Full_Cone<Number>::minimize_support_hyperplanes(){
    if(Support_Hyperplanes.nr_of_rows() == 0)
        return;
    if(isComputed(ConeProperty::SupportHyperplanes)){
        nrSupport_Hyperplanes=Support_Hyperplanes.nr_of_rows();
        return;
    }
    if (verbose) {
        verboseOutput() << "Minimize the given set of support hyperplanes by "
                        << "computing the extreme rays of the dual cone" << endl;
    }
    Full_Cone<Number> Dual(Support_Hyperplanes);
    Dual.verbose=verbose;
    Dual.Support_Hyperplanes = Generators;
    Dual.is_Computed.set(ConeProperty::SupportHyperplanes);
    Dual.compute_extreme_rays();
    Support_Hyperplanes = Dual.Generators.submatrix(Dual.Extreme_Rays_Ind); //only essential hyperplanes
    is_Computed.set(ConeProperty::SupportHyperplanes);
    nrSupport_Hyperplanes=Support_Hyperplanes.nr_of_rows();
    do_all_hyperplanes=false;
}
    

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::compute_extreme_rays(bool use_facets){

    if (isComputed(ConeProperty::ExtremeRays))
        return;
    // when we do approximation, we do not have the correct hyperplanes
    // and cannot compute the extreme rays
    if (is_approximation)
        return;
    assert(isComputed(ConeProperty::SupportHyperplanes));
    
    check_pointed();
    if(!pointed){
        throw NonpointedException();
    }

    if(dim*Support_Hyperplanes.nr_of_rows() < nr_gen) {
         compute_extreme_rays_rank(use_facets);
    } else {
         compute_extreme_rays_compare(use_facets);
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::compute_extreme_rays_rank(bool use_facets){

    if (verbose) verboseOutput() << "Select extreme rays via rank ... " << flush;

    size_t i;
    vector<key_t> gen_in_hyperplanes;
    gen_in_hyperplanes.reserve(Support_Hyperplanes.nr_of_rows());
    Matrix<Number> M(Support_Hyperplanes.nr_of_rows(),dim);

    deque<bool> Ext(nr_gen,false);
    #pragma omp parallel for firstprivate(gen_in_hyperplanes,M)
    for(i=0;i<nr_gen;++i){
//        if (isComputed(ConeProperty::Triangulation) && !in_triang[i])
//            continue;
        gen_in_hyperplanes.clear();
        if(use_facets){
            typename list<FACETDATA>::const_iterator IHV=Facets.begin();            
            for (size_t j=0; j<Support_Hyperplanes.nr_of_rows(); ++j, ++IHV){
                if(IHV->GenInHyp.test(i))
                    gen_in_hyperplanes.push_back(j);
            }            
        }
        else{
            for (size_t j=0; j<Support_Hyperplanes.nr_of_rows(); ++j){
            if(v_scalar_product(Generators[i],Support_Hyperplanes[j])==0)
                gen_in_hyperplanes.push_back(j);
            }
        }
        if (gen_in_hyperplanes.size() < dim-1)
            continue;
        if (M.rank_submatrix(Support_Hyperplanes,gen_in_hyperplanes) >= dim-1)
            Ext[i]=true;   
    }
    for(i=0; i<nr_gen;++i)
        Extreme_Rays_Ind[i]=Ext[i];

    is_Computed.set(ConeProperty::ExtremeRays);
    if (verbose) verboseOutput() << "done." << endl;
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::compute_extreme_rays_compare(bool use_facets){

    if (verbose) verboseOutput() << "Select extreme rays via comparison ... " << flush;

    size_t i,j,k;
    // Matrix<Number> SH=getSupportHyperplanes().transpose();
    // Matrix<Number> Val=Generators.multiplication(SH);
    size_t nc=Support_Hyperplanes.nr_of_rows();
    
    vector<vector<bool> > Val(nr_gen);
    for (i=0;i<nr_gen;++i)
       Val[i].resize(nc);
        
    // In this routine Val[i][j]==1, i.e. true, indicates that
    // the i-th generator is contained in the j-th support hyperplane
    
    vector<key_t> Zero(nc);
    vector<key_t> nr_ones(nr_gen);

    for (i = 0; i <nr_gen; i++) {
        k=0;
        Extreme_Rays_Ind[i]=true;
        if(use_facets){
            typename list<FACETDATA>::const_iterator IHV=Facets.begin();            
            for (j=0; j<Support_Hyperplanes.nr_of_rows(); ++j, ++IHV){
                if(IHV->GenInHyp.test(i)){
                    k++;
                    Val[i][j]=true;                
                }
                else
                Val[i][j]=false;  
            }          
        }
        else{
            for (j = 0; j <nc; ++j) {
                if (v_scalar_product(Generators[i],Support_Hyperplanes[j])==0) {
                    k++;
                    Val[i][j]=true;                
                }
                else
                    Val[i][j]=false;  
            }
        }
        nr_ones[i]=k;
        if (k<dim-1||k==nc)  // not contained in enough facets or in all (0 as generator)
            Extreme_Rays_Ind[i]=false;
    }
    
    maximal_subsets(Val,Extreme_Rays_Ind);    

    is_Computed.set(ConeProperty::ExtremeRays);
    if (verbose) verboseOutput() << "done." << endl;
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::contains(const vector<Number>& v) {
    for (size_t i=0; i<Support_Hyperplanes.nr_of_rows(); ++i)
        if (v_scalar_product(Support_Hyperplanes[i],v) < 0)
            return false;
    return true;
}
//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::contains(const Full_Cone& C) {
    for(size_t i=0;i<C.nr_gen;++i)
        if(!contains(C.Generators[i])){
            cerr << "Missing generator " << C.Generators[i] << endl;
            return(false);
    }
    return(true);
}
//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::check_pointed() {
    if (isComputed(ConeProperty::IsPointed))
        return;
    assert(isComputed(ConeProperty::SupportHyperplanes));
    if (verbose) verboseOutput() << "Checking pointedness ... " << flush;

    pointed = (Support_Hyperplanes.max_rank_submatrix_lex().size() == dim);
    is_Computed.set(ConeProperty::IsPointed);
    if (verbose) verboseOutput() << "done." << endl;
}


//---------------------------------------------------------------------------
template<typename Number>
void Full_Cone<Number>::disable_grading_dep_comp() {

    if (do_multiplicity || do_deg1_elements || do_h_vector) {
        if (do_default_mode) {
            // if (verbose)
            //    verboseOutput() << "No grading specified and cannot find one. "
            //                    << "Disabling some computations!" << endl;
            do_deg1_elements = false;
            do_h_vector = false;
            if(!explicit_full_triang){
                do_triangulation=false;
                do_partial_triangulation=true;
            }
        } else {
            throw BadInputException("No grading specified and cannot find one. Cannot compute some requested properties!");
        }
    }
}

//---------------------------------------------------------------------------

/* computes a degree function, s.t. every generator has value >0 */
template<typename Number>
vector<Number> Full_Cone<Number>::compute_degree_function() const {
    size_t i;  
    vector<Number> degree_function(dim,0);
    // add hyperplanes to get a degree function
        if(verbose) {
            verboseOutput()<<"computing degree function... "<<flush;
        }
        size_t h;
        for (h=0; h<Support_Hyperplanes.nr_of_rows(); ++h) {
            for (i=0; i<dim; i++) {
                degree_function[i] += Support_Hyperplanes.get_elem(h,i);
            }
        }
        v_simplify(degree_function);
        if(verbose) {
            verboseOutput()<<"done."<<endl;
        }
    return degree_function;
}

//---------------------------------------------------------------------------
// Constructors
//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::reset_tasks(){
    do_triangulation = false;
    do_partial_triangulation = false;
    do_determinants = false;
    do_multiplicity=false;
    do_integrally_closed = false;
    do_Hilbert_basis = false;
    do_deg1_elements = false;
    keep_triangulation = false;
    do_Stanley_dec=false;
    do_h_vector=false;
    do_hsop = false;
    do_excluded_faces=false;
    do_approximation=false;
    do_default_mode=false;
    do_class_group = false;
    do_module_gens_intcl = false;
    do_module_rank = false;
    do_cone_dec=false;
    stop_after_cone_dec=false;
    
    do_extreme_rays=false;
    do_pointed=false;
    
    do_evaluation = false;
    do_only_multiplicity=false;

    use_bottom_points = true;

    nrSimplicialPyr=0;
    totalNrPyr=0;
    is_pyramid = false;
    triangulation_is_nested = false;
    triangulation_is_partial = false;
}


//---------------------------------------------------------------------------

template<typename Number>
Full_Cone<Number>::Full_Cone(const Matrix<Number>& M, bool do_make_prime){ // constructor of the top cone
// do_make_prime left for syntax reasons, irrelevantant

    dim=M.nr_of_columns();
    if(dim>0)
        Generators=M;
    // M.pretty_print(cout);
    // assert(M.row_echelon()== dim); rank check now done later 
    
    /*index=1;                      // not used at present
    for(size_t i=0;i<dim;++i)
        index*=M[i][i];
    index=Iabs(index); */

    //make the generators coprime, remove 0 rows and duplicates
    has_generator_with_common_divisor = false; // irrelevant

    Generators.simplify_rows();

    Generators.remove_duplicate_and_zero_rows();
    nr_gen = Generators.nr_of_rows();

    if (nr_gen != static_cast<size_t>(static_cast<key_t>(nr_gen))) {
        throw FatalException("Too many generators to fit in range of key_t!");
    }
    
    multiplicity = 0;
    is_Computed = bitset<ConeProperty::EnumSize>();  //initialized to false
    is_Computed.set(ConeProperty::Generators);
    pointed = false;
    is_simplicial = nr_gen == dim;
    deg1_extreme_rays = false;
    deg1_generated = false;
    deg1_generated_computed = false;
    deg1_hilbert_basis = false;
    
    reset_tasks();
    
    Extreme_Rays_Ind = vector<bool>(nr_gen,false);
    in_triang = vector<bool> (nr_gen,false);
    deg1_triangulation = true;
    if(dim==0){            //correction needed to include the 0 cone;
        is_Computed.set(ConeProperty::Triangulation);
    }
    pyr_level=-1;
    Top_Cone=this;
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=i;
    totalNrSimplices=0;
    TriangulationBufferSize=0;

    detSum = 0;
    shift = 0;
    
    FS.resize(omp_get_max_threads());
    
    Pyramids.resize(20);  // prepare storage for pyramids
    nrPyramids.resize(20,0);
      
    recursion_allowed=true;
    
    do_all_hyperplanes=true;
    // multithreaded_pyramid=true; now in build_cone where it is defined dynamically

    
    nextGen=0;
    store_level=0;
    
    Comparisons.reserve(nr_gen);
    nrTotalComparisons=0;

    inhomogeneous=false;
    
    level0_dim=dim; // must always be defined
    
    use_existing_facets=false;
    start_from=0;
    old_nr_supp_hyps=0;
    
    verbose=false;
    
    RankTest = vector< Matrix<Number> >(omp_get_max_threads(), Matrix<Number>(0,dim));
    
    do_bottom_dec=false;
    suppress_bottom_dec=false;
    keep_order=false;

    approx_level = 1;
    is_approximation=false;
    
    PermGens.resize(nr_gen);
    for(size_t i=0;i<nr_gen;++i)
        PermGens[i]=i;
}

//---------------------------------------------------------------------------

/* constructor for pyramids */
template<typename Number>
Full_Cone<Number>::Full_Cone(Full_Cone<Number>& C, const vector<key_t>& Key) {

    Generators = C.Generators.submatrix(Key);
    dim = Generators.nr_of_columns();
    nr_gen = Generators.nr_of_rows();
    has_generator_with_common_divisor = C.has_generator_with_common_divisor;
    is_simplicial = nr_gen == dim;
    
    Top_Cone=C.Top_Cone; // relate to top cone
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=C.Top_Key[Key[i]];
  
    multiplicity = 0;
    
    Extreme_Rays_Ind = vector<bool>(nr_gen,false);
    is_Computed.set(ConeProperty::ExtremeRays, C.isComputed(ConeProperty::ExtremeRays));
    if(isComputed(ConeProperty::ExtremeRays))
        for(size_t i=0;i<nr_gen;i++)
            Extreme_Rays_Ind[i]=C.Extreme_Rays_Ind[Key[i]];
    in_triang = vector<bool> (nr_gen,false);
    deg1_triangulation = true;

    // not used in a pyramid, but set precaution
    deg1_extreme_rays = false;
    deg1_generated = false;
    deg1_generated_computed = false;
    deg1_hilbert_basis = false;
    
    Grading=C.Grading;
    is_Computed.set(ConeProperty::Grading, C.isComputed(ConeProperty::Grading));
    Order_Vector=C.Order_Vector;

    do_extreme_rays=false;
    do_triangulation=C.do_triangulation;
    do_partial_triangulation=C.do_partial_triangulation;
    do_determinants=C.do_determinants;
    do_multiplicity=C.do_multiplicity;
    do_deg1_elements=C.do_deg1_elements;
    do_h_vector=C.do_h_vector;
    do_Hilbert_basis=C.do_Hilbert_basis;
    keep_triangulation=C.keep_triangulation;
    do_only_multiplicity=C.do_only_multiplicity;
    do_evaluation=C.do_evaluation;
    do_Stanley_dec=C.do_Stanley_dec;
    inhomogeneous=C.inhomogeneous;   // at present not used in proper pyramids
    is_pyramid=true;
    
    pyr_level=C.pyr_level+1;

    totalNrSimplices=0;
    detSum = 0;
    shift = C.shift;
    if(C.gen_degrees.size()>0){ // now we copy the degrees
        gen_degrees.resize(nr_gen);
        for (size_t i=0; i<nr_gen; i++) {
            gen_degrees[i] = C.gen_degrees[Key[i]];
        }
    }
    if(C.gen_levels.size()>0){ // now we copy the levels
        gen_levels.resize(nr_gen);
        for (size_t i=0; i<nr_gen; i++) {
            gen_levels[i] = C.gen_levels[Key[i]];
        }
    }
    TriangulationBufferSize=0;

    
    recursion_allowed=C.recursion_allowed; // must be reset if necessary 
    do_all_hyperplanes=true; //  must be reset for non-recursive pyramids
    // multithreaded_pyramid=false; // SEE ABOVE
    
    nextGen=0;
    store_level = C.store_level;
    
    Comparisons.reserve(nr_gen);
    nrTotalComparisons=0;
    
    level0_dim = C.level0_dim; // must always be defined
    
    use_existing_facets=false;
    start_from=0;
    old_nr_supp_hyps=0;
    verbose=false;
    
    approx_level = C.approx_level;
    is_approximation = C.is_approximation;
    
    do_bottom_dec=false;
    suppress_bottom_dec=false;
    keep_order=true;
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::set_zero_cone() {
    
    assert(dim==0);
    
    if (verbose) {
        verboseOutput() << "Zero cone detected!" << endl;
    }
    
    // The basis change already is transforming to zero.
    is_Computed.set(ConeProperty::Sublattice);
    is_Computed.set(ConeProperty::Generators);
    is_Computed.set(ConeProperty::ExtremeRays);
    Support_Hyperplanes=Matrix<Number> (0);
    is_Computed.set(ConeProperty::SupportHyperplanes);    
    totalNrSimplices = 0;
    is_Computed.set(ConeProperty::TriangulationSize);    
    detSum = 0;
    is_Computed.set(ConeProperty::Triangulation);
    
    pointed = true;
    is_Computed.set(ConeProperty::IsPointed);
    
    deg1_extreme_rays = true;
    is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
    
    if (inhomogeneous) {  // empty set of solutions
        is_Computed.set(ConeProperty::VerticesOfPolyhedron);        
        module_rank = 0;
        is_Computed.set(ConeProperty::ModuleRank);
        is_Computed.set(ConeProperty::ModuleGenerators);             
        level0_dim=0;
        is_Computed.set(ConeProperty::RecessionRank);
    }
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::isComputed(ConeProperty::Enum prop) const{
    return is_Computed.test(prop);
}

//---------------------------------------------------------------------------
// Data access
//---------------------------------------------------------------------------

template<typename Number>
size_t Full_Cone<Number>::getDimension()const{
    return dim;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Full_Cone<Number>::getNrGenerators()const{
    return nr_gen;
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::isPointed()const{
    return pointed;
}

//---------------------------------------------------------------------------

template<typename Number>
bool Full_Cone<Number>::isDeg1ExtremeRays() const{
    return deg1_extreme_rays;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Full_Cone<Number>::getModuleRank()const{
    return module_rank;
}


//---------------------------------------------------------------------------

template<typename Number>
const Matrix<Number>& Full_Cone<Number>::getGenerators()const{
    return Generators;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<bool> Full_Cone<Number>::getExtremeRays()const{
    return Extreme_Rays_Ind;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Full_Cone<Number>::getSupportHyperplanes()const{
    return Support_Hyperplanes;
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::error_msg(string s) const{
    errorOutput() <<"\nFull Cone "<< s<<"\n";
}

//---------------------------------------------------------------------------

template<typename Number>
void Full_Cone<Number>::print()const{
    verboseOutput()<<"\ndim="<<dim<<".\n";
    verboseOutput()<<"\nnr_gen="<<nr_gen<<".\n";
    // verboseOutput()<<"\nhyp_size="<<hyp_size<<".\n";
    verboseOutput()<<"\nGrading is:\n";
    verboseOutput()<< Grading;
    verboseOutput()<<"\nMultiplicity is "<<multiplicity<<".\n";
    verboseOutput()<<"\nGenerators are:\n";
    Generators.pretty_print(verboseOutput());
    verboseOutput()<<"\nExtreme_rays are:\n";
    verboseOutput()<< Extreme_Rays_Ind;
    verboseOutput()<<"\nSupport Hyperplanes are:\n";
    Support_Hyperplanes.pretty_print(verboseOutput());
    verboseOutput()<<"\nHilbert basis is:\n";
    verboseOutput()<< Hilbert_Basis;
    verboseOutput()<<"\nDeg1 elements are:\n";
    verboseOutput()<< Deg1_Elements;
}

} //end namespace
