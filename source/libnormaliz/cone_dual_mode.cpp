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
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <string>
#include <algorithm>

#include "cone_dual_mode.h"
#include "vector_operations.h"
#include "lineare_transformation.h"
#include "list_operations.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//private
//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::splice_them_sort(CandidateList< Integer>& Total, vector<CandidateList< Integer> >& Parts){

    CandidateList<Integer> New;
    for(int i=0;i<omp_get_max_threads();i++)
        Total.Candidates.splice(New.Candidates.end(),Parts[i].Candidates);
    New.sort_by_val();
    Total.merge_by_val(New);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::set_generation_0(CandidateList< Integer>& L){

    typename list<Candidate<Integer> >::iterator c;
    for(c=L.Candidates.begin();c!=L.Candidates.end();++c)
        c->generation=0;
}


//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::select_HB(CandidateList<Integer>& Cand, size_t guaranteed_HB_deg, 
                    CandidateList<Integer>& Irred, bool all_irreducible){

    if(all_irreducible){
        Irred.merge_by_val(Cand);
        return;
    }

    typename list<Candidate<Integer> >::iterator h;
    for(h=Cand.Candidates.begin(); h!=Cand.Candidates.end();){
        if(h->old_tot_deg<=guaranteed_HB_deg){
            Irred.Candidates.splice(Irred.Candidates.end(),Cand.Candidates,h++);
        }
        else{
            ++h;
        }
    }
    Irred.auto_reduce_sorted();  // necessary since the guaranteed HB degree only determines 
                          // in which degrees we can already decide whether an element belongs to the HB            
}


//---------------------------------------------------------------------------


//public
//---------------------------------------------------------------------------

template<typename Integer>
Cone_Dual_Mode<Integer>::Cone_Dual_Mode(const Matrix<Integer>& M){
    dim=M.nr_of_columns();
    if (dim!=M.rank()) {
        errorOutput()<<"Cone_Dual_Mode error: constraints do not define pointed cone!"<<endl;
        // M.pretty_print(errorOutput());
        throw BadInputException();
    }
    SupportHyperplanes = M;
    // support hyperplanes are already coprime (except for truncation/grading)
    // so just remove 0 rows
    SupportHyperplanes.remove_zero_rows();
    nr_sh = SupportHyperplanes.nr_of_rows();
    // hyp_size = dim + nr_sh;
    first_pointed = true;
    Intermediate_HB.dual=true;

    if (nr_sh != static_cast<size_t>(static_cast<key_t>(nr_sh))) {
        errorOutput()<<"Too many support hyperplanes to fit in range of key_t!"<<endl;
        throw FatalException();
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Cone_Dual_Mode<Integer>::get_support_hyperplanes() const {
    return SupportHyperplanes;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Cone_Dual_Mode<Integer>::get_generators()const{
    return Generators;
}

template<typename Integer>
vector<bool> Cone_Dual_Mode<Integer>::get_extreme_rays() const{
    return ExtremeRays;

}


size_t counter=0,counter1=0;

//---------------------------------------------------------------------------

// In the inhomogeneous case or when only degree 1 elements are to be found,
// we truncate the Hilbert basis at level 1. The level is the ordinaryl degree in
// for degree 1 elements and the degree of the homogenizing variable 
// in the inhomogeneous case.
//
// As soon as there are no positive or neutral (with respect to the current hyperplane)
// elements in the current Hilbert basis and truncate==true, new elements can only 
// be produced as sums of positive irreds of level 1 and negative irreds of level 0.
// In particular no new negative elements can be produced, and the only type of
// reduction on the positive side is the elimination of duplicates.
//
// If there are no elements on level 0 at all, then new elements cannot be produced anymore,
// and the production of new elements can be skipped.

template<typename Integer>
void Cone_Dual_Mode<Integer>::cut_with_halfspace_hilbert_basis(const size_t& hyp_counter, 
         const bool lifting, vector<Integer>& old_lin_subspace_half, bool pointed){
    if (verbose==true) {
        verboseOutput()<<"=======================================" << endl;
        verboseOutput()<<"cut with halfspace "<<hyp_counter+1 <<" ..."<<endl;
    }
    
    truncate=inhomogeneous || do_only_Deg1_Elements;
    
    size_t i;
    int sign;

    CandidateList<Integer> Positive_Irred(true),Negative_Irred(true),Neutral_Irred(true);
    Integer orientation, scalar_product,diff,factor;
    vector <Integer> hyperplane=SupportHyperplanes[hyp_counter]; // the current hyperplane dividing the old cone
    typename list<Candidate<Integer> >::iterator h;

    if (lifting==true) {
        orientation=v_scalar_product<Integer>(hyperplane,old_lin_subspace_half);
        if(orientation<0){
            orientation=-orientation;
            v_scalar_multiplication<Integer>(old_lin_subspace_half,-1); //transforming into the generator of the positive half of the old max lin subsapce
        }
        // from now on orientation > 0 
        
        for (h = Intermediate_HB.Candidates.begin(); h != Intermediate_HB.Candidates.end(); ++h) { //reduction  modulo  the generators of the two halves of the old max lin subspace
            scalar_product=v_scalar_product(hyperplane,h->cand); //  allows us to declare "old" HB candiadtes as irreducible
            sign=1;                                                               
            if (scalar_product<0) {
                scalar_product=-scalar_product;
                sign=-1;
            }
            factor=scalar_product/orientation;  // we reduce all elements by the generator of the halfspace
            for (i = 0; i < dim; i++) {
                h->cand[i]=h->cand[i]-sign*factor*old_lin_subspace_half[i];
            }
        }
        
        //adding the generators of the halves of the old max lin subspaces to the the "positive" and the "negative" generators
        // ABSOLUTELY NECESSARY since we need a monoid system of generators of the full "old" cone

        Candidate<Integer> halfspace_gen_as_cand(old_lin_subspace_half,nr_sh);
        halfspace_gen_as_cand.generation=0;
        halfspace_gen_as_cand.mother=0;
        halfspace_gen_as_cand.old_tot_deg=0;
        (halfspace_gen_as_cand.values)[hyp_counter]=orientation; // value under the new linear form
        halfspace_gen_as_cand.sort_deg=explicit_cast_to_long(orientation);
        assert(orientation!=0);
        Positive_Irred.push_back(halfspace_gen_as_cand);
        v_scalar_multiplication<Integer>(halfspace_gen_as_cand.cand,-1);    
        Negative_Irred.push_back(halfspace_gen_as_cand);
        
    } //end lifting
    
    long gen0_mindeg;  // minimal degree of a generator
    if(lifting)
        gen0_mindeg=0;  // sort_deg has already been set > 0 for half_space_gen
    else
        gen0_mindeg=Intermediate_HB.Candidates.begin()->sort_deg;
    typename list<Candidate<Integer> >::const_iterator hh;
    for(hh=Intermediate_HB.Candidates.begin();hh!=Intermediate_HB.Candidates.end();++hh)
        if(hh->sort_deg < gen0_mindeg)
            gen0_mindeg=hh->sort_deg;
        
    bool gen1_pos=false, gen1_neg=false;    
    bool no_pos_in_level0=pointed;
    bool all_positice_level=pointed;
    for (h = Intermediate_HB.Candidates.begin(); h != Intermediate_HB.Candidates.end(); ++h) { //dividing into negative and positive
        Integer new_val=v_scalar_product<Integer>(hyperplane,h->cand);
        long new_val_long=explicit_cast_to_long(new_val);
        h->generation=1;
        h->reducible=false;
        h->mother=0;
        h->old_tot_deg=h->sort_deg;
        if (new_val>0) {
            gen1_pos=true;
            h->values[hyp_counter]=new_val;
            h->sort_deg+=new_val_long;
            Positive_Irred.Candidates.push_back(*h); // could be spliced
            if(h->values[0]==0){
                no_pos_in_level0=false;
                all_positice_level=false;
            }
        }
        if (new_val<0) {
            gen1_neg=true;
            h->values[hyp_counter]=-new_val;
            h->sort_deg+=-new_val_long;
            Negative_Irred.Candidates.push_back(*h);
            if(h->values[0]==0){
                all_positice_level=false;
            }
        }
        if (new_val==0) {
            Neutral_Irred.Candidates.push_back(*h);
            if(h->values[0]==0){
                no_pos_in_level0=false;
                all_positice_level=false;
            }
        }       
    }
   

    if((truncate && (no_pos_in_level0 && !all_positice_level))){
        if(verbose){
            verboseOutput() << "Eliminating negative generators of level > 0" << endl;
        }
        for (h = Negative_Irred.Candidates.begin(); h != Negative_Irred.Candidates.end();){
            if(h->values[0]>0)
                h=Negative_Irred.Candidates.erase(h);
            else
                ++h;
        }
    }
    


    #pragma omp parallel
    {

        #pragma omp single nowait
        {
        check_range(Negative_Irred);
        Negative_Irred.sort_by_val();
        Negative_Irred.last_hyp=hyp_counter;
        }

        #pragma omp single nowait
        {
        check_range(Positive_Irred);
        Positive_Irred.sort_by_val();
        Positive_Irred.last_hyp=hyp_counter;
        }

        #pragma omp single nowait
        Neutral_Irred.sort_by_val();
        Neutral_Irred.last_hyp=hyp_counter;
    }
    
    CandidateList<Integer> New_Positive_Irred(true),New_Negative_Irred(true),New_Neutral_Irred(true);
    New_Negative_Irred.last_hyp=hyp_counter;
    New_Positive_Irred.last_hyp=hyp_counter;
    New_Neutral_Irred.last_hyp=hyp_counter;
    
    CandidateList<Integer> Positive_Depot(true),Negative_Depot(true),Neutral_Depot(true);
    
    vector<CandidateList<Integer> > New_Positive_thread(omp_get_max_threads()),
                      New_Negative_thread(omp_get_max_threads()),
                      New_Neutral_thread(omp_get_max_threads());
    for(long i=0;i<omp_get_max_threads();++i){
        New_Positive_thread[i].dual=true;
        New_Negative_thread[i].dual=true;   
        New_Neutral_thread[i].dual=true;
    }

    typename list<Candidate<Integer> >::const_iterator n,p;
    typename list<Candidate<Integer> >::iterator c;
    
    bool not_done;
    if(lifting)
        not_done=gen1_pos || gen1_neg;
    else
        not_done=gen1_pos && gen1_neg;
    
    while(not_done && !(truncate && all_positice_level)) {

        //generating new elements

        size_t psize=Positive_Irred.size();

        if (verbose) {
            size_t nsize=Negative_Irred.size();
            size_t zsize=Neutral_Irred.size();
            // if (psize*nsize>1000000)
                verboseOutput()<<"Positive: "<<psize<<"  Negative: "<<nsize<<"  Neutral: "<<zsize<<endl;
        }
        
        #pragma omp parallel private(p,n,diff)
        {
        Candidate<Integer> new_candidate(dim,nr_sh);
        new_candidate.generation=1;  // the new generation
        size_t ppos=0;
        p = Positive_Irred.Candidates.begin();
        #pragma omp for schedule(dynamic)
        for(i = 0; i<psize; ++i){
            for(;i > ppos; ++ppos, ++p) ;
            for(;i < ppos; --ppos, --p) ;

            for (n = Negative_Irred.Candidates.begin(); n != Negative_Irred.Candidates.end(); ++n){

                if(truncate && p->values[0]+n->values[0] >=2) // in the inhomogeneous case we truncate at level 1
                    continue;
                
                Integer neg_val=n->values[hyp_counter];
                Integer pos_val=p->values[hyp_counter];

                if (p->generation==0 && n->generation==0)
                    continue; // two "old" candidates have been paired already
                
                if ( (p->mother!=0 && p->mother<=neg_val)|| (n->mother!=0 && n->mother<=pos_val) ){  
                    #pragma omp atomic     // sum would be irreducible by mother + the vector on the opposite side
                    counter1++;
                    continue;
                }
                
                #pragma omp atomic
                counter++;
                
                new_candidate.old_tot_deg=p->old_tot_deg+n->old_tot_deg;
                diff=pos_val-neg_val;
                v_add_result(new_candidate.values,p->values,n->values);   // new_candidate=v_add

                if (diff>0) {
                    new_candidate.values[hyp_counter]=diff;
                    new_candidate.sort_deg=p->sort_deg+n->sort_deg-2*explicit_cast_to_long(neg_val);
                    if(!(truncate && no_pos_in_level0) && (Positive_Irred.is_reducible_last_hyp(new_candidate) ||
                                Neutral_Irred.is_reducible(new_candidate)))
                        continue;
                    v_add_result(new_candidate.cand,p->cand,n->cand);
                    new_candidate.mother=pos_val;                    
                    New_Positive_thread[omp_get_thread_num()].push_back(new_candidate);
                }
                if (diff<0) {
                    if(truncate && no_pos_in_level0) // don't need new negative elements anymore
                        continue;
                    new_candidate.values[hyp_counter]=-diff;
                    new_candidate.sort_deg=p->sort_deg+n->sort_deg-2*explicit_cast_to_long(pos_val);
                    if(Negative_Irred.is_reducible_last_hyp(new_candidate)) {
                        continue;
                    }
                    if(Neutral_Irred.is_reducible(new_candidate)) {
                        continue;
                    }
                    v_add_result(new_candidate.cand,p->cand,n->cand);
                    new_candidate.mother=neg_val;;
                    New_Negative_thread[omp_get_thread_num()].push_back(new_candidate);
                }
                if (diff==0) {
                    new_candidate.values[hyp_counter]=0;
                    new_candidate.sort_deg=p->sort_deg+n->sort_deg-2*explicit_cast_to_long(pos_val);  //pos_val==neg_val
                    if(!(truncate && no_pos_in_level0) && Neutral_Irred.is_reducible(new_candidate)) {
                        continue;
                    }
                    v_add_result(new_candidate.cand,p->cand,n->cand);
                    new_candidate.mother=0; // irrelevant
                    New_Neutral_thread[omp_get_thread_num()].push_back(new_candidate);
                }
            }
        } //end generation of new elements

        } //END PARALLEL


        splice_them_sort(Neutral_Depot,New_Neutral_thread); // sort by sort_deg and values   
        
        splice_them_sort(Positive_Depot,New_Positive_thread);

        splice_them_sort(Negative_Depot,New_Negative_thread);
        
        if(Positive_Depot.empty() && Negative_Depot.empty())
            not_done=false;
            
        // Attention: the element with smallest old_tot_gen need not be the first in the list which is ordered by sort_deg
        size_t gen1_mindeg=0;  // minimal old_tot_deg of a new element used for generation
        bool first=true;
        for(p = Positive_Depot.Candidates.begin();p!=Positive_Depot.Candidates.end();++p){
            if(first){
                first=false;
                gen1_mindeg=p->old_tot_deg;
            }
            if(p->old_tot_deg<gen1_mindeg)
                gen1_mindeg=p->old_tot_deg;
        }
        
        for(p = Negative_Depot.Candidates.begin();p!=Negative_Depot.Candidates.end();++p){
            if(first){
                first=false;
                gen1_mindeg=p->old_tot_deg;
            }
            if(p->old_tot_deg<gen1_mindeg)
                gen1_mindeg=p->old_tot_deg;

        }
        
        size_t min_deg_new=gen0_mindeg+gen1_mindeg;
        if(not_done)
            assert(min_deg_new>0);

        size_t all_known_deg=min_deg_new-1;
        size_t guaranteed_HB_deg=2*all_known_deg+1;  // the degree up to which we can decide whether an element belongs to the HB
        
        if(not_done){
            select_HB(Neutral_Depot,guaranteed_HB_deg,New_Neutral_Irred,truncate && no_pos_in_level0);         
        }
        else{
            Neutral_Depot.auto_reduce();                    // in this case new elements will not be produced anymore
            Neutral_Irred.merge_by_val(Neutral_Depot);      // and there is nothing to do for positive or negative elements
            Neutral_Irred.unique_vectors();                 // but the remaining neutral elements must be auto-reduced.             
       }
       
        // There is a subtle point in the 3 merge operations that follow.
        // In order to have the generations correct, it is important that an
        // element of the Irred lits precede equivalent elements of the New_Irred.
        // This is guaranteed in C++ STL list merge (see Josutis, 2nd ed., p.423)

       if (!New_Neutral_Irred.empty()) {
            if(!(truncate && no_pos_in_level0)){
                Positive_Depot.reduce_by(New_Neutral_Irred);                
                Neutral_Depot.reduce_by(New_Neutral_Irred);
            }
            Negative_Depot.reduce_by(New_Neutral_Irred);
            Neutral_Irred.merge_by_val(New_Neutral_Irred); 
            Neutral_Irred.unique_vectors();
        }                                                 
       
        select_HB(Positive_Depot,guaranteed_HB_deg,New_Positive_Irred,truncate && no_pos_in_level0);

        select_HB(Negative_Depot,guaranteed_HB_deg,New_Negative_Irred,truncate && no_pos_in_level0);
                                                                 
        set_generation_0(Positive_Irred);
        if (!New_Positive_Irred.empty()) {
            if(!(truncate && no_pos_in_level0))
                Positive_Depot.reduce_by(New_Positive_Irred);
            check_range(New_Positive_Irred);  // check for danger of overflow
            Positive_Irred.merge_by_val(New_Positive_Irred);
            Positive_Irred.unique_vectors();
        }

        set_generation_0(Negative_Irred);
        if (!New_Negative_Irred.empty()) {
            Negative_Depot.reduce_by(New_Negative_Irred);
            check_range(New_Negative_Irred);
            Negative_Irred.merge_by_val(New_Negative_Irred);
            Negative_Irred.unique_vectors();
        }

    }  // while(not_done)

    Intermediate_HB.clear();
    Intermediate_HB.Candidates.splice(Intermediate_HB.Candidates.begin(),Positive_Irred.Candidates);
    Intermediate_HB.Candidates.splice(Intermediate_HB.Candidates.end(),Neutral_Irred.Candidates);    
    Intermediate_HB.sort_by_val();       

    if (verbose) {
        verboseOutput()<<"Hilbert basis size="<<Intermediate_HB.size()<<endl;
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Cone_Dual_Mode<Integer>::cut_with_halfspace(const size_t& hyp_counter, const Matrix<Integer>& Basis_Max_Subspace){
    size_t i,j,rank_subspace=Basis_Max_Subspace.nr_of_rows();

    vector <Integer> scalar_product,hyperplane=SupportHyperplanes[hyp_counter],old_lin_subspace_half;
    bool lifting=false;
    Matrix<Integer> New_Basis_Max_Subspace=Basis_Max_Subspace;
    if (rank_subspace!=0) {
        scalar_product=Basis_Max_Subspace.MxV(hyperplane);
        for (i = 0; i <rank_subspace; i++)
            if (scalar_product[i]!=0)
                break;
        if (i!=rank_subspace) {    // the new hyperplane does not contain the intersection of the previous hyperplanes
                                   // so we must intersect the new hyperplane and Max_Subspace
            lifting=true;
            //computing new maximal subspace
            Matrix<Integer> M(1,rank_subspace); // this is the restriction of the new linear form to Max_Subspace
            M[0]=scalar_product;
 
            Lineare_Transformation<Integer> LT=Transformation(M);
            Matrix<Integer> Lifted_Basis_Factor_Space_over_Ker_and_Ker=LT.get_right();
            // the coordinate transfprmation yields a splitting of Max_Subspace into the direct sum of the kernel
            // of the linear form (columns 1^,..) and a complementary 1-dimensional space (column 0)
            // First we dualize from columns to rows:
            Lifted_Basis_Factor_Space_over_Ker_and_Ker=Lifted_Basis_Factor_Space_over_Ker_and_Ker.transpose();

            // Now we must embed the subspaces of Max_Subspace into the full ambient space
            // First the new maximal subspace
            Matrix<Integer>  Ker(rank_subspace-1,rank_subspace);
            for (j = 0; j < rank_subspace-1; j++) {
                Ker[j]= Lifted_Basis_Factor_Space_over_Ker_and_Ker[j+1];
            }
            New_Basis_Max_Subspace=Ker.multiplication(Basis_Max_Subspace);
                        // and then the complementary 1-dim space
            // old_lin_subspace_half refers to the fact that the complementary space is subdivided into
            // two halfspaces generated by old_lin_subspace_half and -old_lin_subspace_half (taken care of in cut_with_halfspace_hilbert_basis
            old_lin_subspace_half=Basis_Max_Subspace.VxM(Lifted_Basis_Factor_Space_over_Ker_and_Ker[0]);
        }
    }
    bool pointed=(Basis_Max_Subspace.nr_of_rows()==0);

    cut_with_halfspace_hilbert_basis(hyp_counter, lifting,old_lin_subspace_half,pointed);

    return New_Basis_Max_Subspace;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::hilbert_basis_dual(){
    if(dim>0){            //correction needed to include the 0 cone;
        if (verbose==true) {
            verboseOutput()<<"\n************************************************************\n";
            verboseOutput()<<"computing Hilbert basis ..."<<endl;
        }
        
        if(Generators.nr_of_rows()!=ExtremeRays.size()){
            errorOutput() << "Mismatch of extreme rays and generators in cone dual mode. THIS SHOULD NOT HAPPEN." << endl;
            throw FatalException(); 
        }
        
        size_t hyp_counter;      // current hyperplane
        Matrix<Integer> Basis_Max_Subspace(dim);      //identity matrix
        for (hyp_counter = 0; hyp_counter < nr_sh; hyp_counter++) {
            Basis_Max_Subspace=cut_with_halfspace(hyp_counter,Basis_Max_Subspace);
        }
        
        if(ExtremeRays.size()==0){  // no precomputed generators
            extreme_rays_rank();
            relevant_support_hyperplanes();
            ExtremeRayList.clear();
            
        }
        else{  // must produce the relevant support hyperplanes from the generators
               // since the Hilbert basis may have been truncated
            vector<Integer> test(SupportHyperplanes.nr_of_rows());
            vector<key_t> key;
            vector <key_t> relevant_sh;
            size_t realdim=Generators.rank();
            for(key_t h=0;h<SupportHyperplanes.nr_of_rows();++h){
                key.clear();
                vector<Integer> test=Generators.MxV(SupportHyperplanes[h]);
                for(key_t i=0;i<test.size();++i)
                    if(test[i]==0)
                        key.push_back(i);
                if (key.size() >= realdim-1 && Generators.submatrix(key).rank() >= realdim-1)
                    relevant_sh.push_back(h);
            }    
            SupportHyperplanes = SupportHyperplanes.submatrix(relevant_sh);
        }
            
        if(verbose)
            verboseOutput() << "matches = " << counter << endl << "avoided = " << counter1 << endl;
            
        Intermediate_HB.extract(Hilbert_Basis);
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::extreme_rays_rank(){
    if (verbose) {
        verboseOutput() << "Find extreme rays" << endl;
    }
    
    typename list < Candidate <Integer> >::iterator c;
    vector <key_t> zero_list;
    size_t i,k;
    for (c=Intermediate_HB.Candidates.begin(); c!=Intermediate_HB.Candidates.end(); ++c){
        zero_list.clear();
        for (i = 0; i < nr_sh; i++) {
            if(c->values[i]==0) {
                zero_list.push_back(i);
            }
        }
        k=zero_list.size();
        if (k>=dim-1) {

            Matrix<Integer> Test=SupportHyperplanes.submatrix(zero_list);
            if (Test.rank_destructive()>=dim-1) {
                ExtremeRayList.push_back(&(*c));
            }
        }
    }
    size_t s = ExtremeRayList.size();
    Generators = Matrix<Integer>(s,dim);
   
    typename  list< Candidate<Integer>* >::const_iterator l;
    for (i=0, l=ExtremeRayList.begin(); l != ExtremeRayList.end(); ++l, ++i) {
        Generators[i]= (*l)->cand;
    }
    ExtremeRays=vector<bool>(s,true);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::relevant_support_hyperplanes(){
    if (verbose) {
        verboseOutput() << "Find relevant support hyperplanes" << endl;
    }
    list <key_t> zero_list;
    typename list<Candidate<Integer>* >::iterator gen_it;
    vector <key_t> relevant_sh;
    relevant_sh.reserve(nr_sh);
    size_t i,k;
    
    size_t realdim = Generators.rank();

    for (i = 0; i < nr_sh; ++i) {
        Matrix<Integer> Test(0,dim);
        k = 0;
        for (gen_it = ExtremeRayList.begin(); gen_it != ExtremeRayList.end(); ++gen_it) {
            if ((*gen_it)->values[i]==0) {
                Test.append((*gen_it)->cand);
                k++;
            }
        }
        if (k >= realdim-1 && Test.rank_destructive()>=realdim-1) {
            relevant_sh.push_back(i);
        }
    }
    SupportHyperplanes = SupportHyperplanes.submatrix(relevant_sh);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::to_sublattice(const Sublattice_Representation<Integer>& SR) {
    assert(SR.get_dim() == dim);

    dim = SR.get_rank();
    // hyp_size = dim+nr_sh;
    SupportHyperplanes = SR.to_sublattice_dual(SupportHyperplanes);
    typename list<vector<Integer> >::iterator it;
    vector<Integer> tmp;
    
    Generators = SR.to_sublattice(Generators);

    for (it = Hilbert_Basis.begin(); it != Hilbert_Basis.end(); ) {
        tmp = SR.to_sublattice(*it);
        it = Hilbert_Basis.erase(it);
        Hilbert_Basis.insert(it,tmp);
    }
}

} //end namespace libnormaliz
