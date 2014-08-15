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
void Cone_Dual_Mode<Integer>::splice_them(CandidateList< Integer>& Total, vector<CandidateList< Integer> >& Parts){

    for(int i=0;i<omp_get_max_threads();i++)
        Total.Candidates.splice(Total.Candidates.end(),Parts[i].Candidates);
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

// PARTLY OBSOLETE because if new data structure

// The vectors in GeneratorList are composed as follows:
// [0]: current total degree (sum of linear forms applied to vector so far)
// [1.. nr_sh]: values under the linear forms
// [nr_sh+1]: generation (in the sense of creation): we match vectors of which at least one has gen >=1
//            and the sum gets gen = 2.
// [nr_sh+2]: value of "mother" under the current linear form (given to newly produced vector)
// last dim coordinates: original coordinates of vector
//
// [nr_sh+1] is used as the degree function as soon as a pointed cone is to be split
// into two halfs. [0] is the total degree in each of the two halves, and not additive 
// in the whole cone.//
//
// In the inhomogeneous case or when only degree 1 elements are to be found,
// we truncate the Hilbert basis at level 1. The level is the ordinaryl degree in
// for degree 1 elements and the degree of the homogenizing variable 
// in the inhomogeneous case.
//
// As soon as there are no positive or neutral (with respect to the current hyperplane)
// elements in the cirrent Hilbert basis and truncate==true, new elements can only 
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
        verboseOutput()<<"cut with halfspace "<<hyp_counter<<" ..."<<endl;
    }
    truncate=inhomogeneous || do_only_Deg1_Elements;
    
    size_t i;
    int sign;
    bool not_done;
    CandidateList<Integer> Positive_Irred(true),Negative_Irred(true),Neutral_Irred(true);
    Integer orientation, scalar_product,diff,factor;
    vector <Integer> hyperplane=SupportHyperplanes[hyp_counter-1]; // the current hyperplane dividing the old cone
    typename list<Candidate<Integer> >::iterator h;

    if (lifting==true) {
        orientation=v_scalar_product<Integer>(hyperplane,old_lin_subspace_half);
        if(orientation<0){
            orientation=-orientation;
            v_scalar_multiplication<Integer>(old_lin_subspace_half,-1); //transforming into the generator of the positive half of the old max lin subsapce
        }
        // from now on orientation > 0 (actually =1)
        
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
        (halfspace_gen_as_cand.cand)[hyp_counter-1]=orientation; // value under the new linear form
        halfspace_gen_as_cand.sort_deg=explicit_cast_to_long(orientation);
        assert(orientation!=0);
        Positive_Irred.Candidates.push_back(halfspace_gen_as_cand);
        v_scalar_multiplication<Integer>(halfspace_gen_as_cand.cand,-1);    
        Negative_Irred.Candidates.push_back(halfspace_gen_as_cand);
        
    } //end lifting
    
    bool no_pos_in_level0=pointed;
    bool all_positice_level=pointed;
    for (h = Intermediate_HB.Candidates.begin(); h != Intermediate_HB.Candidates.end(); ++h) { //dividing into negative and positive
        Integer new_val=v_scalar_product<Integer>(hyperplane,h->cand);
        long new_val_long=explicit_cast_to_long(new_val);
        h->generation=1;
        h->reducible=false;
        h->mother=0;
        h->values.resize(hyp_counter);
        if (new_val>0) {
            h->values[hyp_counter-1]=new_val;
            h->sort_deg+=new_val_long;
            Positive_Irred.Candidates.push_back(*h); // could be spliced
            if(h->values[0]==0){
                no_pos_in_level0=false;
                all_positice_level=false;
            }
        }
        if (new_val<0) {
            h->values[hyp_counter-1]=-new_val;
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
    
    
    // cout << "no_pos_in_level0 = " << no_pos_in_level0 <<endl;
    // cout << "all_positive_level = " << all_positice_level <<endl;

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
        Negative_Irred.sort_it();
        }

        #pragma omp single nowait
        {
        check_range(Positive_Irred);
        Positive_Irred.sort_it();
        }

        #pragma omp single nowait
        Neutral_Irred.sort_it();
    }
    

    //long int counter=0;
    CandidateList<Integer> New_Positive(true),New_Negative(true),New_Neutral(true);
    vector<CandidateList<Integer> > New_Positive_thread(omp_get_max_threads()),
                      New_Negative_thread(omp_get_max_threads()),
                      New_Neutral_thread(omp_get_max_threads());
    for(size_t i=0;i<omp_get_max_threads();++i){
        New_Positive_thread[i].dual=true;
        New_Negative_thread[i].dual=true;   
        New_Neutral_thread[i].dual=true;
    }

    typename list<Candidate<Integer> >::const_iterator n,p;
    typename list<Candidate<Integer> >::iterator c;
    not_done=true;
    while(not_done && !(truncate && all_positice_level)) {
        not_done=false;
        New_Positive.clear();
        New_Negative.clear();
        New_Neutral.clear();

        //generating new elements

        size_t psize=Positive_Irred.size();

        if (verbose) {
            size_t nsize=Negative_Irred.size();
            size_t zsize=Neutral_Irred.size();
            if (psize*nsize>1000000)
                verboseOutput()<<"Positive: "<<psize<<"  Negative: "<<nsize<<"  Neutral: "<<zsize<<endl;
        }
        
        #pragma omp parallel private(p,n,diff)
        {
        Candidate<Integer> new_candidate(nr_sh);
        size_t ppos=0;
        p = Positive_Irred.Candidates.begin();
        #pragma omp for schedule(dynamic)
        for(i = 0; i<psize; ++i){
            for(;i > ppos; ++ppos, ++p) ;
            for(;i < ppos; --ppos, --p) ;

            for (n = Negative_Irred.Candidates.begin(); n != Negative_Irred.Candidates.end(); ++n){

                // cout << "In Schleife " << endl;


                if(truncate && p->values[0]+n->values[0] >=2) // in the inhomogeneous case we truncate at level 1
                    continue;
                    
                assert(p->generation<=1);
                assert(n->generation<=1);
                
                Integer neg_val=n->values[hyp_counter-1];
                Integer pos_val=p->values[hyp_counter-1];

                if (p->generation==0 && n->generation==0)
                    continue; // two "old" candidates have been paired already
                
                if ( (p->mother!=0 && p->mother<=neg_val)|| (n->mother!=0 && n->mother<=pos_val) ){  
                    #pragma omp atomic     // sum would be irreducible by mother + the vector on the opposite side
                    counter1++;
                    continue;
                }
                
                #pragma omp atomic
                counter++;
                diff=pos_val-neg_val;
                v_add_result(new_candidate.values,p->values,n->values);   // new_candidate=v_add

                if (diff>0) {
                    new_candidate.values[hyp_counter-1]=diff;
                    new_candidate.sort_deg-=2*explicit_cast_to_long(neg_val);
                    if(!(truncate && no_pos_in_level0) && (Positive_Irred.is_reducible(new_candidate) ||
                                Neutral_Irred.is_reducible(new_candidate)))
                        continue;
                    v_add_result(new_candidate.cand,p->cand,n->cand);
                    new_candidate.generation=2;  // the new generation
                    new_candidate.mother=pos_val;                    
                    New_Positive_thread[omp_get_thread_num()].push_back(new_candidate);
                }
                if (diff<0) {
                    if(truncate && no_pos_in_level0) // don't need new negative elements anymore
                        continue;
                    new_candidate.values[hyp_counter-1]=-diff;
                    new_candidate.sort_deg-=2*explicit_cast_to_long(pos_val);
                    if(Negative_Irred.is_reducible(new_candidate)) {
                        continue;
                    }
                    if(Neutral_Irred.is_reducible(new_candidate)) {
                        continue;
                    }
                    v_add_result(new_candidate.cand,p->cand,n->cand);
                    new_candidate.generation=2;
                    new_candidate.mother=neg_val;;
                    New_Negative_thread[omp_get_thread_num()].push_back(new_candidate);
                }
                if (diff==0) {
                    new_candidate.values[hyp_counter-1]=0;
                    new_candidate.sort_deg-=2*explicit_cast_to_long(pos_val);  //pos_val==neg_val
                    if(!(truncate && no_pos_in_level0) && Neutral_Irred.is_reducible(new_candidate)) {
                        continue;
                    }
                    v_add_result(new_candidate.cand,p->cand,n->cand);
                    new_candidate.mother=0;
                    new_candidate.generation=0;
                    New_Neutral_thread[omp_get_thread_num()].push_back(new_candidate);
                }
            }
        } //end generation of new elements

        } //END PARALLEL


        splice_them(New_Neutral,New_Neutral_thread);
        New_Neutral.unique_auto_reduce(truncate && no_pos_in_level0);
        
        splice_them(New_Positive,New_Positive_thread);
        New_Positive.unique_auto_reduce(truncate && no_pos_in_level0);

        splice_them(New_Negative,New_Negative_thread);
        New_Negative.unique_auto_reduce(truncate && no_pos_in_level0);

        if (!New_Neutral.empty()) {

            if(!(truncate && no_pos_in_level0)){
                New_Positive.reduce_by(New_Neutral);
                Positive_Irred.reduce_by(New_Neutral);
                
                Neutral_Irred.reduce_by(New_Neutral);
            }

            New_Negative.reduce_by(New_Neutral);
            Negative_Irred.reduce_by(New_Neutral);

            Neutral_Irred.merge(New_Neutral);
        }

        if (!New_Positive.empty()) {
            not_done=true;
            if(!(truncate && no_pos_in_level0))
                Positive_Irred.reduce_by(New_Positive);
            check_range(New_Positive);  // check for danger of overflow
            Positive_Irred.merge(New_Positive);
        }
        #pragma omp single nowait
        if (!New_Negative.empty()) {
            not_done=true;
            Negative_Irred.reduce_by(New_Negative);
            check_range(New_Negative);
            Negative_Irred.merge(New_Negative);
        }


        // adjust generation

        #pragma omp parallel
        {
        #pragma omp single nowait
        for (c = Positive_Irred.Candidates.begin(); c != Positive_Irred.Candidates.end(); ++c){
            if(c->generation>0) {
                c->generation--;
            }
        }
        #pragma omp single nowait
        for (typename list < Candidate <Integer> >::iterator c2 = Negative_Irred.Candidates.begin(); c2 != Negative_Irred.Candidates.end(); ++c2){
            if(c->generation>0) {
                c->generation--;
            }
        }        
        } // END PARALLEL

//      verboseOutput()<<not_done;
    }


    Intermediate_HB.clear();
    Intermediate_HB.Candidates.splice(Intermediate_HB.Candidates.begin(),Positive_Irred.Candidates);
    Intermediate_HB.Candidates.splice(Intermediate_HB.Candidates.end(),Neutral_Irred.Candidates);

    //still possible to have double elements in the Hilbert basis, coming from different generation    
    Intermediate_HB.unique_vectors();

    // (Matrix<Integer>(Intermediate_HB)).print(cout);

    if (verbose) {
        verboseOutput()<<"Hilbert basis size="<<Intermediate_HB.size()<<endl;
    }


}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Cone_Dual_Mode<Integer>::cut_with_halfspace(const size_t& hyp_counter, const Matrix<Integer>& Basis_Max_Subspace){
    size_t i,j,rank_subspace=Basis_Max_Subspace.nr_of_rows();
    // cout << "Dim Unterraum vorher" << Basis_Max_Subspace.nr_of_rows() << endl;
    vector <Integer> scalar_product,hyperplane=SupportHyperplanes[hyp_counter-1],old_lin_subspace_half;
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
void Cone_Dual_Mode<Integer>::extreme_rays_rank(){
    if (verbose) {
        verboseOutput() << "Find extreme rays" << endl;
    }
    
    typename list < vector <Integer> >::iterator c;
    list <key_t> zero_list;
    size_t i,j,k;
    for (c=Hilbert_Basis.begin(); c!=Hilbert_Basis.end(); ++c){
        zero_list.clear();
        for (i = 0; i < nr_sh; i++) {
            if ((*c)[i+1]==0) {
                zero_list.push_back(i);
            }
        }
        k=zero_list.size();
        if (k>=dim-1) {
            vector <key_t> zero_vector(k);
            for (j = 0; j < k; j++) {
                zero_vector[j]=zero_list.front();
                zero_list.pop_front();
            }
            Matrix<Integer> Test=SupportHyperplanes.submatrix(zero_vector);
            if (Test.rank()>=dim-1) {
                GeneratorList.push_back((*c));
            }
        }
    }
    size_t s = GeneratorList.size();
    Generators = Matrix<Integer>(s,dim);
   
    typename  list< vector<Integer> >::const_iterator l;
    for (i=0, l=GeneratorList.begin(); l != GeneratorList.end(); ++l, ++i) {
        Generators[i]= v_cut_front(*l, dim);
    ExtremeRays=vector<bool>(s,true);
    }
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
        for (hyp_counter = 1; hyp_counter <= nr_sh; hyp_counter++) {
            Basis_Max_Subspace=cut_with_halfspace(hyp_counter,Basis_Max_Subspace);
        }
        if(ExtremeRays.size()==0){  // no precomputed generators
            extreme_rays_rank();
            relevant_support_hyperplanes();
            GeneratorList.clear();
            
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
            
        l_cut_front(Hilbert_Basis,dim);
        if(verbose)
            verboseOutput() << "matches = " << counter << endl << "avoided = " << counter1 << endl;
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Cone_Dual_Mode<Integer>::relevant_support_hyperplanes(){
    if (verbose) {
        verboseOutput() << "Find relevant support hyperplanes" << endl;
    }
    list <key_t> zero_list;
    typename list<vector<Integer> >::iterator gen_it;
    vector <key_t> relevant_sh;
    relevant_sh.reserve(nr_sh);
    size_t i,k;
    
    size_t realdim = Generators.rank();

    for (i = 0; i < nr_sh; ++i) {
        Matrix<Integer> Test(0,dim);
        k = 0;
        for (gen_it = GeneratorList.begin(); gen_it != GeneratorList.end(); ++gen_it) {
            if ((*gen_it)[i+1]==0) {
                Test.append( v_cut_front(*gen_it,dim) );
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
