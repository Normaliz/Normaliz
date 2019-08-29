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

//---------------------------------------------------------------------------

#include <stdlib.h>
#include <set>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <time.h>
#include <deque>
#include <math.h>

#include "libnormaliz/cone.h"
#include "libnormaliz/full_cone.h"
#include "libnormaliz/project_and_lift.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/list_operations.h"
#include "libnormaliz/map_operations.h"
#include "libnormaliz/integer.h"
#include "libnormaliz/sublattice_representation.h"
#include "libnormaliz/offload_handler.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

clock_t pyrtime;

const size_t EvalBoundTriang=2500000; // if more than EvalBoundTriang simplices have been stored
                               // evaluation is started (whenever possible)

const size_t EvalBoundPyr=200000;   // the same for stored pyramids of level > 0

const size_t EvalBoundLevel0Pyr=200000; // 1000000;   // the same for stored level 0 pyramids
                                              
const int largePyramidFactor=20;  // pyramid is large if largePyramidFactor*Comparisons[Pyramid_key.size()-dim] > old_nr_supp_hyps

const int SuppHypRecursionFactor=320000; // pyramids for supphyps formed if Pos*Neg > ...

const size_t RAM_Size=1000000000; // we assume that there is at least 1 GB of RAM

const long GMP_time_factor=10; // factor by which GMP arithmetic differs from long long
const long renf_time_factor=20; //N the same for renf
const long renf_time_factor_pyr=5; //used for control of pyramid building without triangulation

const long ticks_norm_quot=155; // approximately the quotient of the ticks row/cont in A553 with GMP

template<typename Integer>
void Full_Cone<Integer>::compute_automorphisms( size_t nr_special_gens){
    
     if(!do_automorphisms || isComputed(ConeProperty::Automorphisms)){
        return;
    }
    
    bool only_from_god_father=false;
    if(do_integrally_closed && descent_level>0) // we can only work with automprphisms induced by God_Father
        only_from_god_father=true;

    if(quality_of_automorphisms!=AutomParam::ambient)
        get_supphyps_from_copy(true); // of course only if they haven't been computed
        extreme_rays_and_deg1_check(); // ditto

        if(!isComputed(ConeProperty::SupportHyperplanes) || !isComputed(ConeProperty::ExtremeRays)){
            throw FatalException("Trying to compute austomorphism group without sufficient data! THIS SHOULD NOT HAPPEN!");
    }
    
    if(!inhomogeneous && quality_of_automorphisms==AutomParam::rational && !isComputed(ConeProperty::Grading))
        throw BadInputException("Rational austomorphism group only computable for polytopes");
    
    if(verbose)
        verboseOutput() << "Computing automorphism group" << endl;
    
    Matrix<Integer> SpecialLinForms(0,dim);
    if(inhomogeneous){
        SpecialLinForms.append(Truncation);
    }
    if(isComputed(ConeProperty::Grading) && Grading.size()>0){
        SpecialLinForms.append(Grading);
    }
    
    if(quality_of_automorphisms!=AutomParam::ambient)    
        Automs=AutomorphismGroup<Integer>(Generators.submatrix(Extreme_Rays_Ind),
                   Support_Hyperplanes,SpecialLinForms);
    else
        Automs=AutomorphismGroup<Integer>(Generators,Support_Hyperplanes,SpecialLinForms);
        
    
    bool success=Automs.compute(quality_of_automorphisms);

    if(!success){
        if(only_from_god_father){
            if(verbose)
                verboseOutput() << "Coputation of automorphism group from extreme rays failed" << endl;
            return;
        }
        if(verbose)
            verboseOutput() << "Coputation of integral automorphism group from extreme rays failed, using Hilbert basis" << endl;
        if(!isComputed(ConeProperty::HilbertBasis)){
            if(verbose)
                verboseOutput() << "Must compute Hilbert basis first, making copy" << endl;
            Full_Cone<Integer> Copy(Generators);
            Copy.do_Hilbert_basis=true;
            Copy.keep_order=true;
            Copy.verbose=verbose;
            Copy.Support_Hyperplanes=Support_Hyperplanes;
            Copy.nrSupport_Hyperplanes=nrSupport_Hyperplanes;
            Copy.is_Computed.set(ConeProperty::SupportHyperplanes);
            Copy.Extreme_Rays_Ind=Extreme_Rays_Ind;
            Copy.is_Computed.set(ConeProperty::ExtremeRays);
            Copy.compute();
            if(Copy.isComputed(ConeProperty::HilbertBasis)){
                Hilbert_Basis.clear();
                Hilbert_Basis.splice(Hilbert_Basis.begin(),Copy.Hilbert_Basis);
                is_Computed.set(ConeProperty::HilbertBasis);
                do_Hilbert_basis=false;
            }
            // do_Hilbert_basis=true; <-- makes no sense            
        }
        
            Automs=AutomorphismGroup<Integer>(Generators.submatrix(Extreme_Rays_Ind),
                   Support_Hyperplanes,SpecialLinForms);
            
            Automs.addComputationGens(Matrix<Integer>(Hilbert_Basis));        
            success=Automs.compute(AutomParam::integral); 
    }
    assert(success==true);
    if(only_from_god_father){
        if(!check_extension_to_god_father())
            return;        
    }
    is_Computed.set(ConeProperty::Automorphisms);
    if(verbose)
        verboseOutput() << Automs.getQualitiesString() << "automorphism group of order " << Automs.getOrder() << "  done" << endl;

}

template<>
void Full_Cone<renf_elem_class>::compute_automorphisms( size_t nr_special_gens){
    
     if(!do_automorphisms || isComputed(ConeProperty::Automorphisms)){
        return;
    }
    
    
    get_supphyps_from_copy(true); // of course only if they haven't been computed
    extreme_rays_and_deg1_check(); // ditto

    if(!isComputed(ConeProperty::SupportHyperplanes) || !isComputed(ConeProperty::ExtremeRays)){
        throw FatalException("Trying to compute austomorphism group without sufficient data! THIS SHOULD NOT HAPPEN!");
            return;
    }

    if(verbose)
        verboseOutput() << "Computing automorphism group" << endl;
    
    Matrix<renf_elem_class> HelpGen=Generators.submatrix(Extreme_Rays_Ind);
    vector<renf_elem_class> HelpGrading;
    if(!inhomogeneous){
        if(!isComputed(ConeProperty::Grading))
            throw NotComputableException("For automorphisms of algebraic polyhedra input must define a polytope");
        HelpGrading=Grading;
    }
    else{
        HelpGrading=Truncation;
    }
    
    /*for(size_t i=0;i<HelpGen.nr_of_rows();++i){ // norm the extreme rays to vertices of polytope
        renf_elem_class test=v_scalar_product(HelpGen[i],HelpGrading);
        if(test==0)
            throw NotComputableException("For automorphisms of algebraic polyhedra input must defime a polytope!");
        v_scalar_division(HelpGen[i],test);       
    }*/

    Matrix<renf_elem_class> SpecialLinForms(0,dim);
    if(HelpGrading.size()>0)
        SpecialLinForms.append(HelpGrading);
    
    Automs=AutomorphismGroup<renf_elem_class>(HelpGen,Support_Hyperplanes,SpecialLinForms);  
    Automs.compute(AutomParam::algebraic);

    is_Computed.set(ConeProperty::Automorphisms);
    if(verbose)
        verboseOutput() << Automs.getQualitiesString() << "automorphism group of order " << Automs.getOrder() << "  done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_facet(const FACETDATA<Integer>& Fac, const size_t& new_generator) const {

    for(size_t jj=0;jj<nr_gen;++jj)
        if(in_triang[jj] && v_scalar_product(Fac.Hyp,Generators[jj])<0){
            cerr << "Hyp negative on generator " << jj << endl;
            assert(false);
        }
        
    vector<key_t> FacetKey;
    for(size_t jj=0;jj<nr_gen;++jj){
        if(in_triang[jj] || jj==new_generator){
            if(Fac.GenInHyp[jj])
                FacetKey.push_back(jj);
        }
        else{
            if(Fac.GenInHyp[jj]){
                cerr << "in_triang error generator " << jj << endl;
                assert(false);
            }
        }
    }
    
    if(Generators.rank_submatrix(FacetKey)<dim-1){
        cerr << "Redundant hyperplane" << endl;
        assert(false);
    }
    
    bool correct=true;
    for(size_t jj=0;jj<nr_gen;++jj){
        /*if(in_triang[jj])
            cerr << jj << endl;*/
        if(in_triang[jj] && Fac.GenInHyp[jj] && v_scalar_product(Fac.Hyp,Generators[jj])!=0){
            cerr << "Damned " << " Index " << jj << endl;
            correct=false;
        }
        if(in_triang[jj] && !Fac.GenInHyp[jj] && v_scalar_product(Fac.Hyp,Generators[jj])==0){
            cerr << "Damned 2" << " Index " << jj << endl;
            correct=false;
        }
    }
    if(!correct){
        cerr << "--------------- ";
        if(is_pyramid) 
            cerr << "pyr";
        cerr << endl;
        assert(false);
    }
    
}
//---------------------------------------------------------------------------


template<typename Integer>
double Full_Cone<Integer>::rank_time() {
    
    size_t nr_tests=50;
    /*
    if(using_GMP<Integer>())
        nr_tests/=GMP_time_factor;
    if(using_renf<Integer>())
        nr_tests/=renf_time_factor;*/
    size_t nr_selected=min(3*dim, nr_gen);
    
    clock_t cl;
    cl= clock();

    #pragma omp parallel
    {
    Matrix<Integer> Test(0,dim);
    #pragma omp for
    for(size_t kk=0;kk<omp_get_max_threads();++kk){
        for(size_t i=0;i<nr_tests;++i){
            vector<key_t> test_key;

            for(size_t j=0;j<nr_selected;++j)
                test_key.push_back(rand() % nr_gen);
            
            Test.rank_submatrix(Generators, test_key);
        }
    }
    }// parallel
    
    cl=clock()-cl;
    
    ticks_rank_per_row=cl;
    ticks_rank_per_row/=(nr_tests*nr_selected*omp_get_max_threads());
    ticks_rank_per_row/= 0.033*(double) omp_get_max_threads()+0.96; // correction taking into account the effort of parallelization

    if(verbose)
        verboseOutput() << "Per row " << ticks_rank_per_row << " ticks " <<endl;
    
    return ticks_rank_per_row;
}

template<typename Integer>
double Full_Cone<Integer>::cmp_time() {
    
    vector<list<boost::dynamic_bitset<> > > Facets_0_1(omp_get_max_threads());

    auto Fac=Facets.begin();
    for(size_t i=0;i<old_nr_supp_hyps;++i,++Fac){
        if(Fac->simplicial)
            continue;            
        Facets_0_1[0].push_back(Fac->GenInHyp);      
    }
    for(size_t i=1;i<omp_get_max_threads();++i)
        Facets_0_1[i]=Facets_0_1[0];
    
    clock_t cl;
    cl= clock();

    #pragma omp parallel
    {
    #pragma omp for
    for(size_t i=0;i<omp_get_max_threads();++i){  
            for(auto p=Facets_0_1[i].begin();p!=Facets_0_1[i].end();++p){
                /*bool contained=*/Facets.begin()->GenInHyp.is_subset_of(*p)
                    && (*p)!=(*Facets_0_1[i].begin())
                    && (*p)!=(*Facets_0_1[i].end()); 
            }
    }
    }
    
    cl=clock()-cl;
    
    ticks_comp_per_supphyp=cl;
    ticks_comp_per_supphyp/=(old_nr_supp_hyps*omp_get_max_threads());
    
    if(verbose)
        verboseOutput() << "Per comparison " << ticks_comp_per_supphyp << " ticks " << endl;
    
    return ticks_comp_per_supphyp;
}
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::set_zero_cone() {
    
    assert(dim==0);
    
    if (verbose) {
        verboseOutput() << "Zero cone detected!" << endl;
    }
    
    // The basis change already is transforming to zero.
    is_Computed.set(ConeProperty::Sublattice);
    is_Computed.set(ConeProperty::Generators);
    is_Computed.set(ConeProperty::ExtremeRays);
    Support_Hyperplanes=Matrix<Integer> (0);
    is_Computed.set(ConeProperty::SupportHyperplanes);    
    totalNrSimplices = 0;
    is_Computed.set(ConeProperty::TriangulationSize);    
    detSum = 0;
    is_Computed.set(ConeProperty::TriangulationDetSum);
    is_Computed.set(ConeProperty::Triangulation);
    is_Computed.set(ConeProperty::StanleyDec);
    multiplicity = 1;
    is_Computed.set(ConeProperty::Multiplicity);
    is_Computed.set(ConeProperty::HilbertBasis);
    if(!inhomogeneous)
        is_Computed.set(ConeProperty::Deg1Elements);
    
    Hilbert_Series = HilbertSeries(vector<num_t>(1,1),vector<denom_t>()); // 1/1
    is_Computed.set(ConeProperty::HilbertSeries);
    
    if (!is_Computed.test(ConeProperty::Grading)) {
        Grading = vector<Integer>(dim);
        // GradingDenom = 1;
        is_Computed.set(ConeProperty::Grading);
    }
    
    pointed = true;
    is_Computed.set(ConeProperty::IsPointed);
    
    deg1_extreme_rays = true;
    is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
    
    deg1_hilbert_basis = true;
    is_Computed.set(ConeProperty::IsDeg1HilbertBasis);
    
    if (inhomogeneous) {  // empty set of solutions
        is_Computed.set(ConeProperty::VerticesOfPolyhedron);        
        module_rank = 0;
        is_Computed.set(ConeProperty::ModuleRank);
        is_Computed.set(ConeProperty::ModuleGenerators);             
        level0_dim=0;
        is_Computed.set(ConeProperty::RecessionRank);
    }
    
    if (!inhomogeneous) {
        ClassGroup.resize(1,0);
        is_Computed.set(ConeProperty::ClassGroup);
    }
    
    if (inhomogeneous || ExcludedFaces.nr_of_rows() != 0) {
        multiplicity = 0;
        is_Computed.set(ConeProperty::Multiplicity);        
        Hilbert_Series.reset(); // 0/1
        is_Computed.set(ConeProperty::HilbertSeries);        
    }
}

template<>
void Full_Cone<renf_elem_class>::set_zero_cone() {
    
    assert(dim==0);
    
    if (verbose) {
        verboseOutput() << "Zero cone detected!" << endl;
    }
    
    // The basis change already is transforming to zero.
    is_Computed.set(ConeProperty::Sublattice);
    is_Computed.set(ConeProperty::Generators);
    is_Computed.set(ConeProperty::ExtremeRays);
    Support_Hyperplanes=Matrix<renf_elem_class> (0);
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

//===========================================================

template<typename Integer>
void Full_Cone<Integer>::check_simpliciality_hyperplane(const FACETDATA<Integer>& hyp) const{
    size_t nr_gen_in_hyp=0;
    for(size_t i=0; i<nr_gen;++i)
        if(in_triang[i]&& hyp.GenInHyp.test(i))
            nr_gen_in_hyp++;
    if((hyp.simplicial &&  nr_gen_in_hyp!=dim-2) || (!hyp.simplicial &&  nr_gen_in_hyp==dim-2)){
        // NOTE: in_triang set at END of main loop in build_cone
        errorOutput() << "Simplicial " << hyp.simplicial << " dim " << dim << " gen_in_hyp " << nr_gen_in_hyp << endl;
        assert(false);
    }
}

template<typename Integer>
void Full_Cone<Integer>::set_simplicial(FACETDATA<Integer>& hyp){
    size_t nr_gen_in_hyp=0;
    for(size_t i=0; i<nr_gen;++i)
        if(in_triang[i]&& hyp.GenInHyp.test(i))
            nr_gen_in_hyp++;
    hyp.simplicial=(nr_gen_in_hyp==dim-2);
}

template<typename Integer>
void Full_Cone<Integer>::number_hyperplane(FACETDATA<Integer>& hyp, const size_t born_at, const size_t mother){
// add identifying number, the birth day and the number of mother
    
    if(don_t_add_hyperplanes)
        return;

    hyp.Mother=mother;
    hyp.BornAt=born_at;
    if(!multithreaded_pyramid){
        hyp.Ident=HypCounter[0];
        HypCounter[0]++;
        return;
    }
    
    int tn;
    if(omp_get_level()==omp_start_level)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(omp_start_level+1);
    hyp.Ident=HypCounter[tn];
    HypCounter[tn]+=omp_get_max_threads();
    assert(HypCounter[tn]%omp_get_max_threads() == (tn+1)%omp_get_max_threads());
    
}

//---------------------------------------------------------------------------

// used to decide if a hyperplane has the order vector on the positive side
// plus lex criterion
template<typename Integer>
bool Full_Cone<Integer>::is_hyperplane_included(FACETDATA<Integer>& hyp) {
    if (!is_pyramid) { // in the topcone we always have ov_sp > 0
        return true;
    }
    //check if it would be an excluded hyperplane
    Integer ov_sp = v_scalar_product(hyp.Hyp,Order_Vector);
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
// produces the linear combination needed for a Fourier-Motzkin step
template<typename Integer>
vector<Integer> Full_Cone<Integer>::FM_comb(const vector<Integer>& Pos, const Integer& PosVal, 
                                            const vector<Integer>& Neg, const Integer& NegVal, bool extract_gcd){
  size_t k;
  vector<Integer> NewFacet(dim);
  for (k = 0; k <dim; k++) {
        NewFacet[k]=PosVal*Neg[k]-NegVal*Pos[k];
        if(!check_range(NewFacet[k]))
            break;    
    }
    
    if(k==dim){
        if(extract_gcd)
            v_make_prime(NewFacet);
    }
    else{
        #pragma omp atomic
        GMP_hyp++;
        vector<mpz_class> mpz_neg(dim), mpz_pos(dim), mpz_sum(dim);
        convert(mpz_neg, Neg);
        convert(mpz_pos, Pos);
        mpz_class mpz_NV, mpz_PV;
        mpz_NV=convertTo<mpz_class>(NegVal);
        mpz_PV=convertTo<mpz_class>(PosVal);
        for (k = 0; k <dim; k++)
            mpz_sum[k]=mpz_PV*mpz_neg[k]-mpz_NV*mpz_pos[k];
        if(extract_gcd)
            v_make_prime(NewFacet);
        v_make_prime(mpz_sum);
        convert(NewFacet,mpz_sum);
    }
    
    return NewFacet;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::add_hyperplane(const size_t& new_generator, const FACETDATA<Integer> & positive,const FACETDATA<Integer> & negative,
                            list<FACETDATA<Integer>>& NewHyps, bool known_to_be_simplicial){
// adds a new hyperplane found in find_new_facets to this cone (restricted to generators processed)
    
    if(don_t_add_hyperplanes)
        return;

    size_t k;
    
    FACETDATA<Integer> NewFacet; NewFacet.Hyp.resize(dim); NewFacet.GenInHyp.resize(nr_gen); 
    NewFacet.is_positive_on_all_original_gens=false;
    NewFacet.is_negative_on_some_original_gen=false;
    
    for (k = 0; k <dim; k++) {
        NewFacet.Hyp[k]=positive.ValNewGen*negative.Hyp[k]-negative.ValNewGen*positive.Hyp[k];
        if(!check_range(NewFacet.Hyp[k]))
            break;    
    }
    
    if(k==dim)
        v_make_prime(NewFacet.Hyp);
    else{
        #pragma omp atomic
        GMP_hyp++;
        vector<mpz_class> mpz_neg(dim), mpz_pos(dim), mpz_sum(dim);
        convert(mpz_neg, negative.Hyp);
        convert(mpz_pos, positive.Hyp);
        for (k = 0; k <dim; k++)
            mpz_sum[k]=convertTo<mpz_class>(positive.ValNewGen)*mpz_neg[k]-convertTo<mpz_class>(negative.ValNewGen)*mpz_pos[k];
        v_make_prime(mpz_sum);
        convert(NewFacet.Hyp, mpz_sum);
    }
    
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
    
    // check_facet(NewFacet, new_generator);
    
    NewHyps.push_back(NewFacet);
}


//---------------------------------------------------------------------------


template<typename Integer>
void Full_Cone<Integer>::find_new_facets(const size_t& new_generator){
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
    deque <FACETDATA<Integer>*> Pos_Simp,Pos_Non_Simp;
    deque <FACETDATA<Integer>*> Neg_Simp,Neg_Non_Simp;
    deque <FACETDATA<Integer>*> Neutral_Simp, Neutral_Non_Simp;
    
    boost::dynamic_bitset<> Zero_Positive(nr_gen),Zero_Negative(nr_gen); // here we collect the vertices that lie in a
                                        // postive resp. negative hyperplane

    bool simplex;
    
    if (tv_verbose) verboseOutput()<<"transform_values:"<<flush;
    
    for (auto &facet : Facets) {
        // simplex=true;
        // nr_zero_i=0;
        simplex=facet.simplicial; // at present simplicial, will become nonsimplicial if neutral
        /* for (size_t j=0; j<nr_gen; j++){
            if (facet.GenInHyp.test(j)) {
                if (++nr_zero_i > facet_dim) {
                    simplex=false;
                    break;
                }
            }
        }*/
        
        if (facet.ValNewGen==0) {
            facet.GenInHyp.set(new_generator);  // Must be set explicitly !!
            facet.simplicial=false;  // simpliciality definitly gone with the new generator
            if (simplex) {
                Neutral_Simp.push_back(&facet); // simplicial without the new generator
            }   else {
                Neutral_Non_Simp.push_back(&facet); // nonsimplicial already without the new generator
            }
        }
        else if (facet.ValNewGen>0) {
            Zero_Positive |= facet.GenInHyp;
            if (simplex) {
                Pos_Simp.push_back(&facet);
            } else {
                Pos_Non_Simp.push_back(&facet);
            }
        } 
        else if (facet.ValNewGen<0) {
            Zero_Negative |= facet.GenInHyp;
            if (simplex) {
                Neg_Simp.push_back(&facet);
            } else {
                Neg_Non_Simp.push_back(&facet);
            }
        }
    }
    
    // TO DO: Negativliste mit Zero_Positive verfeinern, also die aussondern, die nicht genug positive Erz enthalten
    // Eventuell sogar Rang-Test einbauen.
    // Letzteres könnte man auch bei den positiven machen, bevor sie verarbeitet werden
    
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

    
    // remove negative subfacets shared by two neg simpl facets
    for (auto jj =Neg_Subfacet_Multi_United.begin(); jj!= Neg_Subfacet_Multi_United.end();) {
        auto del=jj++;
        if (jj!=Neg_Subfacet_Multi_United.end() && (*jj).first==(*del).first) {   //delete since is the intersection of two negative simplicies
            Neg_Subfacet_Multi_United.erase(del);
            del=jj++;
            Neg_Subfacet_Multi_United.erase(del);
        }
    }

    size_t nr_NegSubfMult = Neg_Subfacet_Multi_United.size();
    if (tv_verbose) verboseOutput() << nr_NegSubfMult << ", " << flush;
    
    vector<list<FACETDATA<Integer>> > NewHypsSimp(nr_PosSimp);
    vector<list<FACETDATA<Integer>> > NewHypsNonSimp(nr_PosNonSimp);

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
    std::exception_ptr tmp_exception;
    
    if((using_GMP<Integer>() || using_renf<Integer>()) && Generators_float.nr_of_rows()==0 ){
        bool potential_ranktest=false;
        if(using_GMP<Integer>())           
            potential_ranktest = (old_nr_supp_hyps > GMP_time_factor*dim*dim*dim/3); // in this case the rank computation takes longer
        if(using_renf<Integer>())           
            potential_ranktest = (old_nr_supp_hyps > renf_time_factor*dim*dim*dim/3);     
        if(potential_ranktest)
            convert(Generators_float,Generators);
    }

    #pragma omp parallel
    {
    size_t i,j,k,nr_zero_i;
    boost::dynamic_bitset<> subfacet(dim-2);
    auto jj = Neg_Subfacet_Multi_United.begin();
    size_t jjpos=0;
    int tn = omp_get_ancestor_thread_num(omp_start_level+1);


    if(nr_PosNonSimp >= nr_NegNonSimp/10){ // to prevent a desaster if there are very few positive facets,
        bool found;                      // but many negative ones and pyramid decomposition is not applied
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
    }

    #pragma omp single
    { //remove elements that where found in the previous loop
    auto last_inserted=Neg_Subfacet.begin(); // used to speedup insertion into the new map
    for (auto jj = Neg_Subfacet_Multi_United.begin(); jj!= Neg_Subfacet_Multi_United.end(); ++jj) {
        if ((*jj).second != -1) {
            last_inserted = Neg_Subfacet.insert(last_inserted,*jj);
        }
    }
    nr_NegSubf=Neg_Subfacet.size();
    }
    
    #pragma omp single nowait
    {Neg_Subfacet_Multi_United.clear();}

    
    boost::dynamic_bitset<> zero_i(nr_gen);
    
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
        try {

        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
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
            auto jj_map=Neg_Subfacet.find(zero_i);           // one subfacet with negative simpl facet
            if (jj_map!=Neg_Subfacet.end()) {
                add_hyperplane(new_generator,*Pos_Simp[i],*Neg_Simp[(*jj_map).second],NewHypsSimp[i],true);
                (*jj_map).second = -1;  // block subfacet in further searches
            }
        }
        if (nr_zero_i==facet_dim){    // now there could be more such subfacets. We make all and search them.      
            for (k =0; k<nr_gen; k++) {  // BOOST ROUTINE
                
                INTERRUPT_COMPUTATION_BY_EXCEPTION
                
                if(zero_i.test(k)) { 
                    subfacet=zero_i;
                    subfacet.reset(k);  // remove k-th element from facet to obtain subfacet
                    auto jj_map=Neg_Subfacet.find(subfacet);
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
           
           INTERRUPT_COMPUTATION_BY_EXCEPTION
            
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
       } catch(const std::exception& ) {
           tmp_exception = std::current_exception();
           skip_remaining = true;
           #pragma omp flush(skip_remaining)
       }

    } // PS vs NS and PS vs N

    if (!skip_remaining) {
    #pragma omp single nowait
    if (tv_verbose) {
        verboseOutput() << "P vs NS and P vs N" << endl;
    }

    /* list<FACETDATA<Integer>*> AllNonSimpHyp;
    typename list<FACETDATA<Integer>*>::iterator a;*/

    list<boost::dynamic_bitset<> > Facets_0_1_thread;
    for(i=0;i<nr_PosNonSimp;++i)
        Facets_0_1_thread.push_back(Pos_Non_Simp[i]->GenInHyp);
    for(i=0;i<nr_NegNonSimp;++i)
        Facets_0_1_thread.push_back(Neg_Non_Simp[i]->GenInHyp);
    for(i=0;i<nr_NeuNonSimp;++i)
        Facets_0_1_thread.push_back(Neutral_Non_Simp[i]->GenInHyp); 
    size_t nr_NonSimp = nr_PosNonSimp+nr_NegNonSimp+nr_NeuNonSimp;
    
    /*list<boost::dynamic_bitset<> > Facets_0_1_thread;
    auto Fac=Facets.begin();
    for(size_t i=0;i<old_nr_supp_hyps;++i){
        if(!Fac->simplicial)
            Facets_0_1_thread.push_back(Fac->GenInHyp);
        ++Fac;        
    }
    assert(Facets_0_1_thread.size()==nr_NonSimp);*/
    
    /* list<boost::dynamic_bitset<>* > AllNonSimpHyp;
    for(auto F01=Facets_0_1.begin(); F01!=Facets_0_1.end();++F01)
        AllNonSimpHyp.push_back(&(*F01)) ;   */     
   
    bool ranktest;
    FACETDATA<Integer> *hp_i, *hp_j; // pointers to current hyperplanes
    
    size_t missing_bound, nr_common_zero;
    boost::dynamic_bitset<> common_zero(nr_gen);
    vector<key_t> common_key;
    common_key.reserve(nr_gen);
    vector<int> key_start(nrGensInCone);
    
    #pragma omp for schedule(dynamic) // nowait
    for (size_t i =0; i<nr_PosNonSimp; i++){ //Positive Non Simp vs.Negative Simp and Non Simp

        if (skip_remaining) continue;

        try {
        INTERRUPT_COMPUTATION_BY_EXCEPTION

        auto jj_map = Neg_Subfacet.begin();       // First the Simp
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
           
            if (subfacet_dim <=2) {  //intersection of i and j is a subfacet
               add_hyperplane(new_generator,*hp_i,*hp_j,NewHypsNonSimp[i],false); //simplicial set in add_hyperplane
               /* #pragma omp atomic
                NrNewF++; */
                // Indi[j]=true;
                // cout << "Subfacet" << endl;
                continue;
           }
           
           /* #pragma omp atomic
           NrCSF++;*/
           
           // clock_t cl=clock();
           
            // a priori values
            if(using_GMP<Integer>())           
                    ranktest = (nr_NonSimp > GMP_time_factor*dim*dim*nr_common_zero/3); // in this case the rank computation takes longer
            else{
                    if(using_renf<Integer>())           
                        ranktest = (nr_NonSimp > renf_time_factor*dim*dim*nr_common_zero/3);
                    else            
                        ranktest = (nr_NonSimp > dim*dim*nr_common_zero/3);               
            }
            
            if(Generators_float.nr_of_rows()>0){
                Matrix<nmz_float>& Test_float = Top_Cone->RankTest_float[tn];
                if(Test_float.rank_submatrix(Generators_float,common_key)<subfacet_dim){
                    ranktest=false;
                }
            }

           if(ranktest) {  //                
               // cout  << "Rang" << endl;            
               Matrix<Integer>& Test = Top_Cone->RankTest[tn];
               if (Test.rank_submatrix(Generators,common_key)<subfacet_dim) {
                   common_subfacet=false;
               }
           } // ranktest
           else{                 // now the comparison test
           
           /* #pragma omp atomic
            NrComp++; */
               auto a=Facets_0_1_thread.begin();
               common_zero = zero_i & hp_j->GenInHyp;
               for (;a!=Facets_0_1_thread.end();++a){
                   if ( common_zero.is_subset_of(*a) && (*a!=hp_i->GenInHyp) && (*a!=hp_j->GenInHyp)) {                                
                       common_subfacet=false;
                       Facets_0_1_thread.splice(Facets_0_1_thread.begin(),Facets_0_1_thread,a); // for the "darwinistic" mewthod
                       break;
                   }
               }                       
           } // else
           
           // cout << nr_common_zero*ticks_rank_per_row << " " << nr_NonSimp*ticks_comp_per_supphyp << " " << clock()-cl << endl;
           
           if (common_subfacet) {  //intersection of i and j is a subfacet
               add_hyperplane(new_generator,*hp_i,*hp_j,NewHypsNonSimp[i],false); //simplicial set in add_hyperplane
               /* #pragma omp atomic
                NrNewF++; */
                // Indi[j]=true;
                // cout << "Subfacet" << endl;
           }
        }
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
        }
    } // end for
    } // end !skip_remaining
    } //END parallel
    
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
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

template<typename Integer>
void Full_Cone<Integer>::extend_triangulation(const size_t& new_generator){
// extends the triangulation of this cone by including new_generator
// simplicial facets save us from searching the "brother" in the existing triangulation
// to which the new simplex gets attached

    size_t listsize =old_nr_supp_hyps; // Facets.size();
    vector<typename list<FACETDATA<Integer>>::iterator> visible;
    visible.reserve(listsize);

    listsize=0;
    for (auto i = Facets.begin(); i!=Facets.end(); ++i) {
        if (i->ValNewGen < 0){ // visible facet
            visible.push_back(i);
            listsize++;
        }
    }

    std::exception_ptr tmp_exception;

    auto oldTriBack = --TriangulationBuffer.end();
    #pragma omp parallel
    {
    size_t k,l;
    bool one_not_in_i, not_in_facet;
    size_t not_in_i=0;
    // size_t facet_dim=dim-1;
    // size_t nr_in_i=0;

    list< SHORTSIMPLEX<Integer> > Triangulation_kk;
    
    vector<key_t> key(dim);
    
    // if we only want a partial triangulation but came here because of a deep level
    // mark if this part of the triangulation has not to be evaluated
    bool skip_eval = false;

    #pragma omp for schedule(dynamic)
    for (size_t kk=0; kk<listsize; ++kk) {

    try {
        INTERRUPT_COMPUTATION_BY_EXCEPTION

        auto i=visible[kk];
        
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
        
            auto j=TriSectionFirst[vertex];
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

        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
        }

    } // omp for kk

    if (multithreaded_pyramid) {
        #pragma omp critical(TRIANG)
        TriangulationBuffer.splice(TriangulationBuffer.end(),Triangulation_kk);
    } else
        TriangulationBuffer.splice(TriangulationBuffer.end(),Triangulation_kk);

    } // parallel

    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);

    // GensInCone.push_back(new_generator); // now in extend_cone
    TriSectionFirst.push_back(++oldTriBack);
    TriSectionLast.push_back(--TriangulationBuffer.end());
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::store_key(const vector<key_t>& key, const Integer& height,
            const Integer& mother_vol, list< SHORTSIMPLEX<Integer> >& Triangulation){
// stores a simplex given by key and height in Triangulation
// mother_vol is the volume of the simplex to which the new one is attached

    SHORTSIMPLEX<Integer> newsimplex;
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
    if(omp_get_level()==omp_start_level)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(omp_start_level+1);
    
    if (do_only_multiplicity) {
        // directly compute the volume
        if (mother_vol==1)
            newsimplex.vol = height;
        // the multiplicity is computed in SimplexEvaluator
        for(size_t i=0; i<dim; i++) // and needs the key in TopCone numbers
            newsimplex.key[i]=Top_Key[newsimplex.key[i]];

        if (keep_triangulation)
            sort(newsimplex.key.begin(),newsimplex.key.end());
        Top_Cone->SimplexEval[tn].evaluate(newsimplex);
        // restore the local generator numbering, needed in extend_triangulation
        newsimplex.key=key;
    }
    if (height == 0) Top_Cone->triangulation_is_partial = true;
    
    if (keep_triangulation){
        Triangulation.push_back(newsimplex);
        return;  
    }
    
    bool Simpl_available=true;

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
                auto F = Top_Cone->FreeSimpl.begin();
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

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::store_key(const vector<key_t>& key, const renf_elem_class& height,
            const renf_elem_class& mother_vol, list< SHORTSIMPLEX<renf_elem_class> >& Triangulation){
// stores a simplex given by key and height in Triangulation
// mother_vol is the volume of the simplex to which the new one is attached

    SHORTSIMPLEX<renf_elem_class> newsimplex;
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
                auto F = Top_Cone->FreeSimpl.begin();
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
#endif

//---------------------------------------------------------------------------

// We measure times for small and large pyramids in order to better control
// which pyramids should be declared large.
//
// THIS IS A CRITICAL SUBROUTINE because of evaluate_large_rec_pyramids.
// One must take care that it does not change the state of *this.
// don_t_add_hyperplanes should keep it under control.
// One must not only avoid adding hyperplanes, but also changing the 
// numbering scheme for hyperplanes.
//
template<typename Integer>
void Full_Cone<Integer>::small_vs_large(const size_t new_generator){
    
    IsLarge=vector<bool> (nr_gen,false);
    
    don_t_add_hyperplanes=true; // during time measurement the addition of hyperplanes is blocked
    
    int save_nr_threads=omp_get_max_threads(); // must block parallelization
    omp_set_num_threads(1);                    // in small pyramids
    
    nr_pyrs_timed=vector<size_t>(nr_gen);
    time_of_large_pyr=vector<clock_t> (nr_gen);
    time_of_small_pyr=vector<clock_t>(nr_gen);
    
    auto hyp=Facets.begin();
    vector<key_t> Pyramid_key;
    size_t start_level=omp_get_level();
    
    for (size_t kk=0; kk<old_nr_supp_hyps; ++kk, ++hyp) {
        
        if(kk%50 != 0) 
            continue;
        
        if (hyp->ValNewGen >= 0) // facet not visible
                continue;
        
        Pyramid_key.clear(); // make data of new pyramid
        Pyramid_key.push_back(new_generator);
        for(size_t i=0;i<nr_gen;i++){
            if(in_triang[i] && hyp->GenInHyp.test(i)) {
                Pyramid_key.push_back(i);
            }
        }
        
        bool large=(largePyramidFactor*Comparisons[Pyramid_key.size()-dim] > old_nr_supp_hyps); // a priori decision
        if(large)
            continue;
        
        if(nr_pyrs_timed[Pyramid_key.size()]>=5)
            continue;
        
        // we first treat it as small pyramid 
        clock_t cl=clock();
        process_pyramid(Pyramid_key, new_generator,store_level,0, true,hyp,start_level); // is recursive, 0 blocks triangulation
        time_of_small_pyr[Pyramid_key.size()]+=clock()-cl;
        nr_pyrs_timed[Pyramid_key.size()]++;
        //now as large pyramid
        LargeRecPyrs.push_back(*hyp);
    }
    
    take_time_of_large_pyr=true;
    bool save_verbose=verbose;
    verbose=false;
    evaluate_large_rec_pyramids(new_generator);
    verbose=save_verbose;
    take_time_of_large_pyr=false;
    

    
    /* for(size_t i=dim-1;i<nr_gen;++i)
        cout << i << " " << nr_pyrs_timed[i] << " " << time_of_small_pyr[i] << " " << time_of_large_pyr[i] << endl;
 
    for(size_t i=dim-1;i<nr_gen;++i){
        if(nr_pyrs_timed[i]!=0)
            if(time_of_small_pyr[i] > time_of_large_pyr[i])
                cout << i << " " << time_of_small_pyr[i] << " " << time_of_large_pyr[i] << endl;
        
    }*/
    
    int kk;
    for(kk=nr_gen-1;kk>=dim;--kk){
        if(time_of_small_pyr[kk]==0)
            continue;            
        if(time_of_small_pyr[kk]>time_of_large_pyr[kk])
            IsLarge[kk]=true;
        else
            break;
    }
    
    // cout << "First large " << kk+1 << endl;
    
    don_t_add_hyperplanes=false;
    omp_set_num_threads(save_nr_threads);
    
    assert(Facets.size()==old_nr_supp_hyps);
        
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::process_pyramids(const size_t new_generator,const bool recursive){

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
    
    if((using_renf<Integer>() || using_GMP<Integer>()) && Generators_float.nr_of_rows()==0){
        // Generators.pretty_print(cout);
        convert(Generators_float,Generators);
    }
    
    if( !is_pyramid && !time_measured && recursive){ // (using_GMP<Integer>() || using_renf<Integer>()) &&
        rank_time();
        cmp_time();
        /* ticks_quot=(ticks_rank_per_row/ticks_comp_per_supphyp)/ticks_norm_quot;
        if(verbose)
            verboseOutput() << "Normed quotient " << ticks_quot << endl;*/
        time_measured=true;
    }
    
    IsLarge.clear();
    
    if(using_renf<Integer>()  && recursive && !is_pyramid && (!do_partial_triangulation || do_triangulation)){
        if(ticks_rank_per_row>20.0)
            small_vs_large(new_generator);
    }

    size_t start_level=omp_get_level(); // allows us to check that we are on level 0
                                        // outside the loop and can therefore call evaluation
                                        // in order to empty the buffers
    vector<key_t> Pyramid_key;
    Pyramid_key.reserve(nr_gen);
    bool skip_triang; // make hyperplanes but skip triangulation (recursive pyramids only)

    deque<bool> done(old_nr_supp_hyps,false);
    bool skip_remaining;
    std::exception_ptr tmp_exception;
    size_t nr_done=0;

    do{  // repeats processing until all hyperplanes have been processed

    auto hyp=Facets.begin();
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
        
        try {
            for(;kk > hyppos; hyppos++, hyp++) ;
            for(;kk < hyppos; hyppos--, hyp--) ;
            
            INTERRUPT_COMPUTATION_BY_EXCEPTION

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

        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
        }
    } // end parallel loop over hyperplanes

    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);

    if (!omp_in_parallel())
        try_offload(0);
    
    if (start_level==0 && check_evaluation_buffer_size()) {
        Top_Cone->evaluate_triangulation();
    }
    
    if (start_level==0 && Top_Cone->check_pyr_buffer(store_level)) {
        Top_Cone->evaluate_stored_pyramids(store_level);
    }
    
    if (verbose && old_nr_supp_hyps>=RepBound)
        verboseOutput() << endl;

    } while (nr_done < old_nr_supp_hyps);
    
    // cout << float_comp << " " << wrong << " " << wrong_short << endl;
    
    // wrong_positive=0;
    // wrong_negative=0;
    // total_comp_large_pyr=0;
    
    evaluate_large_rec_pyramids(new_generator);
    
    // cout << total_comp_large_pyr << " " << wrong_positive << " " << wrong_negative << endl;

}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::process_pyramid(const vector<key_t>& Pyramid_key,
                          const size_t new_generator,const size_t store_level, Integer height, const bool recursive,
                          typename list< FACETDATA<Integer> >::iterator hyp, size_t start_level){
// processes simplicial pyramids directly, stores other pyramids into their depots

    #pragma omp atomic
    Top_Cone->totalNrPyr++;

    if(Pyramid_key.size()==dim){  // simplicial pyramid completely done here
        #pragma omp atomic        // only for saving memory
        Top_Cone->nrSimplicialPyr++;
        if(recursive){ // the facets may be facets of the mother cone and if recursive==true must be given back
            Matrix<Integer> H(dim,dim);
            Integer dummy_vol;
            Generators.simplex_data(Pyramid_key,H, dummy_vol,false);
            list<FACETDATA<Integer>> NewFacets;
            FACETDATA<Integer> NewFacet;
            NewFacet.GenInHyp.resize(nr_gen);
            for (size_t i=0; i<dim;i++) {
                NewFacet.Hyp = H[i];
                NewFacet.GenInHyp.set();
                NewFacet.GenInHyp.reset(i);
                NewFacet.simplicial=true;
                NewFacet.is_positive_on_all_original_gens=false;
                NewFacet.is_negative_on_some_original_gen=false;
                NewFacets.push_back(NewFacet);
            }
            vector<bool> Pyr_in_triang(dim,true);
            select_supphyps_from(NewFacets,new_generator,Pyramid_key,Pyr_in_triang); // takes itself care of multithreaded_pyramid
        }
        if (height != 0 && (do_triangulation || do_partial_triangulation)) {
            if(multithreaded_pyramid) {
                std::exception_ptr tmp_exception;
                #pragma omp critical(TRIANG)
                {
                try{
                    store_key(Pyramid_key,height,0,TriangulationBuffer);
                    nrTotalComparisons+=dim*dim/2;
                } catch(const std::exception& ) {
                    tmp_exception = std::current_exception();
                }
                } // end critical
                if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
            } else {
                store_key(Pyramid_key,height,0,TriangulationBuffer);
                nrTotalComparisons+=dim*dim/2;
            }
        }
    }
    else {  // non-simplicial

        bool large;
        
        if(IsLarge.size()==0){ // no measurement
            long large_factor=largePyramidFactor;
            if(time_measured){
                mpq_class large_factor_mpq(ticks_rank_per_row);
                mpz_class add=round(large_factor_mpq);
                large_factor+=convertTo<long>(add);
            }
            large=(large_factor*Comparisons[Pyramid_key.size()-dim] > old_nr_supp_hyps);
        }
        else{ // with measurement
            large=(largePyramidFactor*Comparisons[Pyramid_key.size()-dim] > old_nr_supp_hyps);
            large =large || IsLarge[Pyramid_key.size()];
        }
        
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
        
        // pyrtime=clock();

        Full_Cone<Integer> Pyramid(*this,Pyramid_key);
        Pyramid.Mother = this;
        Pyramid.Mother_Key = Pyramid_key;    // need these data to give back supphyps
        Pyramid.apex=new_generator;
        if (height == 0) { //indicates "do not triangulate"
            Pyramid.do_triangulation = false;
            Pyramid.do_partial_triangulation = false;
            Pyramid.do_Hilbert_basis = false;
            Pyramid.do_deg1_elements=false;
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
            Pyramid.do_Hilbert_basis = false;
            Pyramid.do_deg1_elements=false;
        }

        Pyramid.build_cone();
        
        // cout << "Pyramid ticks " << clock() - pyrtime << endl;

        if(multithreaded_pyramid){
            #pragma omp atomic
            nrTotalComparisons+=Pyramid.nrTotalComparisons;
        } else
            nrTotalComparisons+=Pyramid.nrTotalComparisons;
    }  // else non-simplicial
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_and_evaluate_start_simplex(){
    
    // it is absolutely necessary to chooses the start simplex as the lex smallest+
    // in order to be consistent with pyramid decomposition

    size_t i,j; 
    
    vector<key_t> key=find_start_simplex();
    assert(key.size()==dim); // safety heck
    if(verbose){
        verboseOutput() << "Start simplex ";
        for(unsigned int i : key)
            verboseOutput() <<  i+1 << " ";
        verboseOutput() << endl;
    }
    Matrix<Integer> H(dim,dim);
    Integer vol;
    Generators.simplex_data(key,H,vol,do_partial_triangulation || do_triangulation);
    
    assert(key.size()==dim); // safety heck
    
    // cout << "Nach First " << clock()-pyrtime << endl;
    
    // cout << "Nach LinAl " << clock()-pyrtime << endl;
    
    // H.pretty_print(cout);
    
        
    for (i = 0; i < dim; i++) {
        in_triang[key[i]]=true;
        GensInCone.push_back(key[i]);
        if (deg1_triangulation && isComputed(ConeProperty::Grading))
            deg1_triangulation = (gen_degrees[key[i]] == 1);
    }
    
    nrGensInCone=dim;
    
    nrTotalComparisons=dim*dim/2;
    //if(!time_measured){
        if(using_GMP<Integer>())
            nrTotalComparisons*=(GMP_time_factor/4); // because of the linear algebra involved in this routine
        if(using_renf<Integer>())           
            nrTotalComparisons*=(renf_time_factor/4);
    /*}
    else{
        if(using_GMP<Integer>())
            nrTotalComparisons*=(GMP_time_factor/4)*ticks_quot; 
        if(using_renf<Integer>())
            nrTotalComparisons*=(renf_time_factor/4)*ticks_quot; 
    }*/
    
    Comparisons.push_back(nrTotalComparisons);
       
    for (i = 0; i <dim; i++) {
        FACETDATA<Integer> NewFacet; NewFacet.GenInHyp.resize(nr_gen);
        NewFacet.is_positive_on_all_original_gens=false;
        NewFacet.is_negative_on_some_original_gen=false;
        NewFacet.Hyp=H[i];
        if(using_renf<Integer>())
            v_standardize(NewFacet.Hyp,Norm);
        NewFacet.simplicial=true; // indeed, the start simplex is simplicial
        for(j=0;j < dim;j++)
            if(j!=i)
                NewFacet.GenInHyp.set(key[j]);
        NewFacet.ValNewGen=-1;         // must be taken negative since opposite facet
        number_hyperplane(NewFacet,0,0); // created with gen 0
        Facets.push_back(NewFacet);    // was visible before adding this vertex
    }
    
    Integer factor;
    if(!is_pyramid){
        //define Order_Vector, decides which facets of the simplices are excluded
        Order_Vector = vector<Integer>(dim,0);
        // Matrix<Integer> G=S.read_generators();
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
    
    // cout << "Nach Start " << clock()-pyrtime << endl;
    
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::select_supphyps_from(const list<FACETDATA<Integer>>& NewFacets, 
                    const size_t new_generator, const vector<key_t>& Pyramid_key, const vector<bool>& Pyr_in_triang){
// the mother cone (=this) selects supphyps from the list NewFacets supplied by the daughter
// the daughter provides the necessary information via the parameters

    size_t i;
    boost::dynamic_bitset<> in_Pyr(nr_gen);
    for (i=0; i<Pyramid_key.size(); i++) {
        in_Pyr.set(Pyramid_key[i]);
    }
    // the new generator is always the first in the pyramid
    assert(Pyramid_key[0] == new_generator);

    bool new_global_hyp;
    FACETDATA<Integer> NewFacet;
    NewFacet.is_positive_on_all_original_gens=false;
    NewFacet.is_negative_on_some_original_gen=false;
    NewFacet.GenInHyp.resize(nr_gen);
    Integer test;
    for (const auto& pyr_hyp : NewFacets) {
        if(!pyr_hyp.GenInHyp.test(0)) // new gen not in hyp
            continue;
        new_global_hyp=true;
        for (i=0; i<nr_gen; ++i){
            if(in_Pyr.test(i) || !in_triang[i])
                continue;
            test=v_scalar_product(Generators[i],pyr_hyp.Hyp);
            if(test<=0){
                new_global_hyp=false;
                break;
            }

        }
        if(new_global_hyp){
            NewFacet.Hyp=pyr_hyp.Hyp;
            NewFacet.GenInHyp.reset();
            // size_t gens_in_facet=0;
            for (i=0; i<Pyramid_key.size(); ++i) {
                if(in_triang[Pyramid_key[i]]) // this satisfied in the standard setting where the pyramid key is strictly ascending after
                    assert(Pyr_in_triang[i]); // the first entry and the start simplex of the pyramid is lex smallest
                if (pyr_hyp.GenInHyp.test(i) && in_triang[Pyramid_key[i]]) {
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
            NewFacet.simplicial=pyr_hyp.simplicial; // (gens_in_facet==dim-1); 
            check_simpliciality_hyperplane(NewFacet);
            number_hyperplane(NewFacet,nrGensInCone,0); //mother unknown
            
            if(don_t_add_hyperplanes)
                continue;
            
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
template<typename Integer>
void Full_Cone<Integer>::match_neg_hyp_with_pos_hyps(const FACETDATA<Integer>& hyp, size_t new_generator,list<FACETDATA<Integer>*>& PosHyps, boost::dynamic_bitset<>& Zero_P, vector<list<boost::dynamic_bitset<> > >& Facets_0_1){

    size_t missing_bound, nr_common_zero;
    boost::dynamic_bitset<> common_zero(nr_gen);
    vector<key_t> common_key;
    common_key.reserve(nr_gen);
    vector<key_t> key(nr_gen);
    bool common_subfacet;
    list<FACETDATA<Integer>> NewHyp;
    size_t subfacet_dim=dim-2;
    size_t nr_missing;
    list<FACETDATA<Integer>> NewHyps;
    Matrix<Integer> Test(0,dim);
    
    int tn;
    if(omp_get_level()==omp_start_level)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(omp_start_level+1);
    
    boost::dynamic_bitset<> zero_hyp=hyp.GenInHyp & Zero_P;  // we intersect with the set of gens in positive hyps
    
    vector<int> key_start(nrGensInCone);
    size_t nr_zero_hyp=0;
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
    
    missing_bound=nr_zero_hyp-subfacet_dim; // at most this number of generators can be missing
                                          // to have a chance for common subfacet

    for (const auto& hp_j : PosHyps){ //match hyp with the given Pos
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
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
       boost::dynamic_bitset<> common_zero(nr_gen);
       
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
                   common_zero.set(key[k]);
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
               common_zero.set(key[k]);
               nr_common_zero++;
           }
        }
        
       if(!common_subfacet)
            continue;
       
       assert(nr_common_zero >=subfacet_dim);
       
        //bool negative_predicted=false, positive_predicted=false;
        
        if (!hp_j->simplicial){
            
            // #pragma omp atomic
            // total_comp_large_pyr++;
            
           bool ranktest;
           
           if (time_measured){
              ranktest=(ticks_rank_per_row*nr_common_zero < ticks_per_cand);
           }
           else{ // a priori values
                if(using_GMP<Integer>())           
                        ranktest = (old_nr_supp_hyps > GMP_time_factor*dim*dim*nr_common_zero/3); // in this case the rank computation takes longer
                else{
                        if(using_renf<Integer>())           
                            ranktest = (old_nr_supp_hyps > renf_time_factor*dim*dim*nr_common_zero/3);
                        else            
                            ranktest = (old_nr_supp_hyps > dim*dim*nr_common_zero/3);               
                } 
           }
           
           //if(!ranktest)
           //    cout << " No Rank " << endl;
           
            // ranktest=true;
            
            /* if(nr_common_zero==subfacet_dim)
                ranktest=false;*/
            
            /* clock_t cl;
            
           if(verbose && time_measured){
                    verboseOutput() << ticks_rank_per_row*nr_common_zero << " " << ticks_per_cand << endl;
                    cl=clock();                
            }*/
            
            // ranktest=true;
            
            if(ranktest && Generators_float.nr_of_rows()>0){
                Matrix<nmz_float>& Test_float = Top_Cone->RankTest_float[tn];
                if(Test_float.rank_submatrix(Generators_float,common_key)<subfacet_dim){
                    ranktest=false;
                    //negative_predicted=true;
                }
                else{
                    //positive_predicted=true;
                }
            }
            
            if(ranktest){
                // cout << "Rank" << endl;
                Matrix<Integer>& Test = Top_Cone->RankTest[tn];
                if(Test.rank_submatrix(Generators,common_key)<subfacet_dim)
                    common_subfacet=false;     // don't make a hyperplane
            }
            else{                 // now the comparison test
                // cout << "Compare" << endl;
                for (auto hp_t=Facets_0_1[tn].begin();hp_t!=Facets_0_1[tn].end();++hp_t){
                    if (common_zero.is_subset_of(*hp_t) && (*hp_t!=hyp.GenInHyp) && (*hp_t!=hp_j->GenInHyp) ) {
                        Facets_0_1[tn].splice(Facets_0_1[tn].begin(),Facets_0_1[tn],hp_t); // successful reducer to the front
                        common_subfacet=false;
                        break;
                    }
                }
                
            } // else

            /* if(verbose && time_measured){
                cl=clock()-cl;
                verboseOutput() << cl << endl;
             }
             
             if(common_subfacet && verbose)
                 verboseOutput() << "Subfacet" << endl;*/
        } // !simplicial
        
        /* if(negative_predicted && common_subfacet){
            #pragma omp atomic
            wrong_negative++;
        }
        if(positive_predicted && !common_subfacet){
            #pragma omp atomic
            wrong_positive++;
        }*/      
        
        
        if(common_subfacet){
            add_hyperplane(new_generator,*hp_j,hyp,NewHyps,false);  // simplicial set in add_hyperplane
        }
    } // for           

    if(multithreaded_pyramid)
        #pragma omp critical(GIVEBACKHYPS)
        Facets.splice(Facets.end(),NewHyps);
    else
        Facets.splice(Facets.end(),NewHyps);

}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::collect_pos_supphyps(list<FACETDATA<Integer>*>& PosHyps, boost::dynamic_bitset<>& Zero_P, size_t& nr_pos){
           
    // positive facets are collected in a list
    
    auto ii = Facets.begin();
    nr_pos=0;
    
    for (size_t ij=0; ij< old_nr_supp_hyps; ++ij, ++ii)
        if (ii->ValNewGen>0) {
            Zero_P |= ii->GenInHyp;
            PosHyps.push_back(&(*ii));
            nr_pos++;
        }
}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::evaluate_large_rec_pyramids(size_t new_generator){
    
    size_t nrLargeRecPyrs=LargeRecPyrs.size();
    if(nrLargeRecPyrs==0)
        return;
    
    vector<list<boost::dynamic_bitset<> > > Facets_0_1(omp_get_max_threads());
    
    size_t nr_non_simplicial=0;

    auto Fac=Facets.begin();
    for(size_t i=0;i<old_nr_supp_hyps;++i,++Fac){
        if(Fac->simplicial)
            continue;            
        Facets_0_1[0].push_back(Fac->GenInHyp);
        nr_non_simplicial++;
    }
    for(size_t j=1;j<omp_get_max_threads();++j)
         Facets_0_1[j]= Facets_0_1[0];        
        
    if(verbose)
        verboseOutput() << "large pyramids " << nrLargeRecPyrs << endl;
    
    list<FACETDATA<Integer>*> PosHyps;
    boost::dynamic_bitset<> Zero_P(nr_gen);
    size_t nr_pos;
    collect_pos_supphyps(PosHyps,Zero_P,nr_pos);
    
    nrTotalComparisons+=nr_pos*nrLargeRecPyrs;
    std::exception_ptr tmp_exception;
    
    const long VERBOSE_STEPS = 50;
    long step_x_size = nrLargeRecPyrs-VERBOSE_STEPS;
    const size_t RepBound=100;
    
    // clock_t cl=clock();
    /* auto LP=LargeRecPyrs.begin();
    for(size_t i=0; i<nrLargeRecPyrs; i++,++LP){
        if(i%100 ==1)
            match_neg_hyp_with_pos_hyps(*LP,new_generator,PosHyps,Zero_P,true,nr_cand,Facets_0_1);
    }
    cl=clock()-cl;
    ticks_per_cand=(double) cl/(double) nr_cand;
    if(verbose) 
        verboseOutput() << "Ticks per cand " << ticks_per_cand << endl;*/
    
    ticks_per_cand=ticks_comp_per_supphyp*nr_non_simplicial; // estimated time for testing an irreducible by comparison
    
    bool skip_remaining=false;
    
    #pragma omp parallel if(!take_time_of_large_pyr)
    {
    size_t ppos=0;
    auto p=LargeRecPyrs.begin(); 
    
    #pragma omp for schedule(dynamic)
    for(size_t i=0; i<nrLargeRecPyrs; i++){
        
        if(skip_remaining)
            continue;
        
        for(; i > ppos; ++ppos, ++p) ;
        for(; i < ppos; --ppos, --p) ;

        if(verbose && nrLargeRecPyrs>=RepBound){
            #pragma omp critical(VERBOSE)
            while ((long)(i*VERBOSE_STEPS) >= step_x_size) {
                step_x_size += nrLargeRecPyrs;
                verboseOutput() << "." <<flush;
            }
        }
        
        // cout << "=================================" << endl;
        
        // cout << "Neg Hyp " << i << endl;
        
        // clock_t cl;
        // cl= clock();
        
        
        
        try {

            INTERRUPT_COMPUTATION_BY_EXCEPTION
            
            clock_t cl_large = 0;
            if(take_time_of_large_pyr){
                cl_large=clock();
            }
            
            match_neg_hyp_with_pos_hyps(*p,new_generator,PosHyps,Zero_P,Facets_0_1);
            
            if(take_time_of_large_pyr){
                cl_large=clock()-cl_large;
                size_t nr_pyr_gens=0;
                for(size_t i=0;i<nr_gen;++i)
                    if(p->GenInHyp[i])
                        nr_pyr_gens++;
                nr_pyr_gens++; // for the apex of the pyramid
                time_of_large_pyr[nr_pyr_gens]+=cl_large;                
            }
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
        }
        
        // cl=clock()-cl;
        // cout << "Neg Hyp " << i << " ticks " << cl << endl; 
    }
    } // parallel
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
    
    if(verbose && nrLargeRecPyrs>=RepBound)
        verboseOutput() << endl;

    LargeRecPyrs.clear();
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::check_pyr_buffer(const size_t level){
    if(level==0)
        return(nrPyramids[0] > EvalBoundLevel0Pyr);
    else
        return(nrPyramids[level] > EvalBoundPyr);
}

//---------------------------------------------------------------------------

#ifdef NMZ_MIC_OFFLOAD
template<>
void Full_Cone<long long>::try_offload(size_t max_level) {

    if (!is_pyramid && _Offload_get_device_number() < 0) // dynamic check for being on CPU (-1)
    {
        
        if (max_level >= nrPyramids.size()) max_level = nrPyramids.size()-1;
        for (size_t level = 0; level <= max_level; ++level) {                  
            if (nrPyramids[level] >= 100) {
                // cout << "XXX: Try offload of level " << level << " pyramids ..." << endl;
                mic_offloader.offload_pyramids(*this, level);
                break;
            }
        }
    }
}

template<typename Integer>
void Full_Cone<Integer>::try_offload(size_t max_level) {
}
//else it is implemented in the header


template<typename Integer>
void Full_Cone<Integer>::try_offload_loc(long place,size_t max_level){
    verboseOutput() << "From place " << place << " " << "level " << max_level << endl;
    try_offload(max_level);
}

#endif // NMZ_MIC_OFFLOAD

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_stored_pyramids(const size_t level){
// evaluates the stored non-recursive pyramids
    
#ifdef NMZ_MIC_OFFLOAD
    Pyramids_scrambled[level]=false;
    
    if(level==0 && _Offload_get_device_number() >=  0){
        verboseOutput() << "Start evaluation of " << nrPyramids[level] << " pyrs on level " << level << endl;
        // verboseOutput() << "In parallel " << omp_in_parallel() << endl;
}
#endif // NMZ_MIC_OFFLOAD

    if(Pyramids[level].empty())
        return;

    assert(omp_get_level()==omp_start_level); // assert(!omp_in_parallel());
    assert(!is_pyramid);
   
    if (Pyramids.size() < level+2) {
        Pyramids.resize(level+2);      // provide space for a new generation
        nrPyramids.resize(level+2, 0);
        Pyramids_scrambled.resize(level+2, false);
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
    size_t ppos;
    bool skip_remaining;
    std::exception_ptr tmp_exception;

    while (nrPyramids[level] > eval_down_to) {

       auto p = Pyramids[level].begin();
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

           try {
               INTERRUPT_COMPUTATION_BY_EXCEPTION
               
               Full_Cone<Integer> Pyramid(*this,*p);
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
           } catch(const std::exception& ) {
               tmp_exception = std::current_exception();
               skip_remaining = true;
               #pragma omp flush(skip_remaining)
           }
        } //end parallel for
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);

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

        try_offload(level+1);

        if (check_evaluation_buffer_size()) {
            if (verbose)
                verboseOutput() << nrPyramids[level] <<
                    " pyramids remaining on level " << level << ", ";
            Top_Cone->evaluate_triangulation();
            try_offload(level+1);
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
template<typename Integer>
void Full_Cone<Integer>::build_cone() {
    // if(dim>0){            //correction needed to include the 0 cone;
    
    // cout << "Pyr " << pyr_level << endl;
    
    if(start_from==0)    
        in_triang = vector<bool> (nr_gen,false);

    size_t RecBoundSuppHyp;
    RecBoundSuppHyp = dim*SuppHypRecursionFactor;
    if(using_GMP<Integer>())
        RecBoundSuppHyp*=GMP_time_factor; // pyramid building is more difficult for complicated arithmetic
    if(using_renf<Integer>())
        RecBoundSuppHyp*=renf_time_factor_pyr*renf_degree*renf_degree;
        
    size_t RecBoundTriang=1000000;   //  if number(supphyps)*size(triang) > RecBoundTriang pass to pyramids
    if(using_GMP<Integer>())
        RecBoundTriang*=GMP_time_factor;
    if(using_renf<Integer>())
        RecBoundTriang*=renf_time_factor;
    
    tri_recursion=false; 
    
    multithreaded_pyramid=(omp_get_level()==omp_start_level);
    
    size_t nr_original_gen=0;
    size_t steps_in_approximation = 0;
    if (!is_pyramid && is_approximation)
    {
        nr_original_gen = OriginalGenerators.nr_of_rows();    
        vector<size_t> nr_approx_points; // how many points are in the approximation
        for (size_t j=0;j<nr_original_gen;++j) {
            nr_approx_points.push_back(approx_points_keys[j].size());
        }
        // for every vertex sort the approximation points via: number of positive halfspaces / index
        vector<key_t> overall_perm;
        // stores the perm of every list 
        vector<vector<key_t> > local_perms(nr_original_gen);
        
        for (size_t current_gen = 0 ; current_gen<nr_original_gen;++current_gen){
            vector<key_t> local_perm;
            if (approx_points_keys[current_gen].size()>0){
                auto jt=approx_points_keys[current_gen].begin();
                list<pair<size_t,key_t> > max_halfspace_index_list;
                size_t tmp_hyp=0;
                // TODO: collect only those which belong to the current generator?
                for (;jt!=approx_points_keys[current_gen].end();++jt){
                    // cout << dim << " " << Support_Hyperplanes.nr_of_columns()<< " " << Generators[*jt].size() << endl;
                    tmp_hyp = v_nr_negative(Support_Hyperplanes.MxV(Generators[*jt])); // nr of negative halfspaces
                    max_halfspace_index_list.insert(max_halfspace_index_list.end(),make_pair(tmp_hyp,*jt));
                }
                max_halfspace_index_list.sort([](const pair<size_t,key_t> &left, const pair<size_t,key_t> &right) {
                    return right.first < left.first;
                });
                for(const auto& list_it: max_halfspace_index_list){
                    local_perm.push_back(list_it.second);
                }
            }
            local_perms[current_gen]=local_perm;
        }
        // concatenate the permutations
        size_t local_perm_counter=0;
        bool not_done = true;
        while (not_done){
            not_done=false;
            for (size_t current_gen=0;current_gen<nr_original_gen;++current_gen){
                if (local_perm_counter<nr_approx_points[current_gen]){
                    not_done=true;
                    overall_perm.push_back(local_perms[current_gen][local_perm_counter]);
                }
            }    
            ++local_perm_counter;
        }
        assert(overall_perm.size()==nr_gen);
        // sort the generators according to the permutations
        Generators.order_rows_by_perm(overall_perm);
    }
    
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
    
    bool check_original_gens=true;

    for (size_t i=start_from;i<nr_gen;++i) {
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION 
        
        time_t start,end;
        time (&start);
    
        start_from=i;
    
        if (in_triang[i])
            continue;
        
        if (!is_pyramid && is_approximation) steps_in_approximation++;    
        // we check whether all original generators are contained in the current cone
        if (!is_pyramid && is_approximation && check_original_gens){
            if (verbose)
                verboseOutput() << "Check...";
            size_t current_gen=0;
            auto l=Facets.begin();
            for (; l!=Facets.end(); ++l){
                if (l->is_positive_on_all_original_gens) continue;
                for (current_gen=0;current_gen<nr_original_gen;++current_gen){
                    if (!(v_scalar_product(l->Hyp,OriginalGenerators[current_gen])>=0)) {
                        l->is_negative_on_some_original_gen=true;
                        check_original_gens=false;
                        break;
                    }
                }
                if (current_gen==nr_original_gen){
                    l->is_positive_on_all_original_gens=true;    
                } else {
                    break;
                }
            } 
            if(verbose)   
                verboseOutput() << " done." << endl; 
            // now we need to stop
            if (l==Facets.end()){
                if(verbose)
                    verboseOutput() << "The original cone is now contained." << endl;
                break;
            }
        }
            
        if(do_triangulation && TriangulationBufferSize > 2*RecBoundTriang) // emermergency brake
            tri_recursion=true;               // to switch off production of simplices in favor
                                              // of non-recursive pyramids
        Integer scalar_product;                                              
        is_new_generator=false;
        auto l=Facets.begin();
        old_nr_supp_hyps=Facets.size(); // Facets will be xtended in the loop 

        long long nr_pos=0, nr_neg=0;
        long long nr_neg_simp=0, nr_pos_simp=0;
        vector<Integer> L;           
        std::exception_ptr tmp_exception;
        
        size_t lpos=0;
        #pragma omp parallel for private(L,scalar_product) firstprivate(lpos,l) reduction(+: nr_pos, nr_neg)
        for (size_t k=0; k<old_nr_supp_hyps; k++) {
            try {
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
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
            }
        }  //end parallel for
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);

        if(!is_new_generator)
            continue;

        // the i-th generator is used in the triangulation
        // in_triang[i]=true; // now at end of loop
        if (deg1_triangulation && isComputed(ConeProperty::Grading))
            deg1_triangulation = (gen_degrees[i] == 1);
        
        if (!omp_in_parallel())
            try_offload(0);
        /* if(!is_pyramid && verbose ) 
            verboseOutput() << "Neg " << nr_neg << " Pos " << nr_pos << " NegSimp " <<nr_neg_simp << " PosSimp " <<nr_pos_simp << endl; */
        // First we test whether to go to recursive pyramids because of too many supphyps
        if (recursion_allowed && nr_neg*nr_pos-(nr_neg_simp*nr_pos_simp) > RecBoundSuppHyp) {  // use pyramids because of supphyps
            if(!is_pyramid && verbose )
                verboseOutput() << "Building pyramids" << endl;
            if (do_triangulation)
                tri_recursion = true; // We can not go back to classical triangulation
            if(check_evaluation_buffer()){
                Top_Cone->evaluate_triangulation();
            }

            process_pyramids(i,true); //recursive
            // lastGen=i;
            // nextGen=i+1; 
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
        size_t nr_new_facets = Facets.size() - old_nr_supp_hyps;
        time (&end);        
        /* double dif = difftime (end,start);

        if (verbose) {
            verboseOutput() << "Generator took " << dif << " sec " <<endl;
        }*/

        if(using_renf<Integer>()){
            // v_standardize new facets in Qnormaliz for renf_elem_class        
            l=Facets.begin();
            for (size_t j=0; j<old_nr_supp_hyps;j++)
                l++;        
            for(;l!=Facets.end();++l)
                v_standardize(l->Hyp,  Norm);
        }
        
        // removing the negative hyperplanes if necessary
        if(do_all_hyperplanes || i!=last_to_be_inserted){
            l=Facets.begin();
            for (size_t j=0; j<old_nr_supp_hyps;j++){
                if (l->ValNewGen<0) {
                    if (is_approximation && l->is_negative_on_some_original_gen){
                        check_original_gens = true;
                    }
                    l=Facets.erase(l);
                }
                else
                    ++l;
            }
        }
        
        GensInCone.push_back(i);
        nrGensInCone++;
        
        Comparisons.push_back(nrTotalComparisons);
        
        in_triang[i]=true;
        
        if(verbose) {
            verboseOutput() << "gen="<< i+1 <<", ";
            if (do_all_hyperplanes || i!=last_to_be_inserted) {
                verboseOutput() << Facets.size()<<" hyp, " << nr_new_facets << " new";
            } else {
                verboseOutput() << Support_Hyperplanes.nr_of_rows()<<" hyp";
            }
            if(nrPyramids[0]>0)
                verboseOutput() << ", " << nrPyramids[0] << " pyr"; 
            if(do_triangulation||do_partial_triangulation)
                verboseOutput() << ", " << TriangulationBufferSize << " simpl";
            verboseOutput()<< endl;
        }
        
    }  // loop over i
    
    start_from=0; // in order that we can restart the primal algorithm again
    
    if (is_pyramid && do_all_hyperplanes)  // must give supphyps back to mother
        Mother->select_supphyps_from(Facets, apex, Mother_Key, in_triang);
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    // transfer Facets --> SupportHyperplanes
    if (do_all_hyperplanes) {
        nrSupport_Hyperplanes = Facets.size();
        Support_Hyperplanes = Matrix<Integer>(nrSupport_Hyperplanes,0);
        auto IHV=Facets.begin();
        for (size_t i=0; i<nrSupport_Hyperplanes; ++i, ++IHV) {
            if(keep_convex_hull_data)
                Support_Hyperplanes[i]=IHV->Hyp;
            else
                swap(Support_Hyperplanes[i],IHV->Hyp);
        }
        is_Computed.set(ConeProperty::SupportHyperplanes);
    } 
    Support_Hyperplanes.set_nr_of_columns(dim);
   
    
    if(do_extreme_rays && do_all_hyperplanes)
        compute_extreme_rays(true);
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    transfer_triangulation_to_top(); // transfer remaining simplices to top
    if(check_evaluation_buffer()){
        Top_Cone->evaluate_triangulation();
    }  
    
    if (!is_pyramid && is_approximation && verbose){
        verboseOutput() << "Performed " << steps_in_approximation << "/" << nr_gen << " steps." << endl;
    }
    // } // end if (dim>0)
    
    if(!keep_convex_hull_data)
        Facets.clear(); 

}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_bottom_facets() {

    if(verbose)
        verboseOutput() << "Computing bottom decomposition" << endl;

    vector<key_t> start_simpl=Generators.max_rank_submatrix_lex();
    Order_Vector = vector<Integer>(dim,0);
    for(size_t i=0;i<dim;++i)
        for(size_t j=0;j<dim;++j)
            Order_Vector[j]+=((unsigned long) (1+i%10))*Generators[start_simpl[i]][j];

    // First the generators for the recession cone = our cone
    Matrix<Integer> BottomGen(0,dim+1);
    vector<Integer> help(dim+1);
    for(size_t i=0;i<nr_gen;++i){
        for(size_t j=0;j<dim; ++j)
            help[j]=Generators[i][j];
        help[dim]=0;
        BottomGen.append(help);
    }
    // then the same vectors as generators of the bottom polyhedron
    for(size_t i=0;i<nr_gen;++i){
        for(size_t j=0;j<dim; ++j)
            help[j]=Generators[i][j];
        help[dim]=1;
        BottomGen.append(help);
    }
    
    Full_Cone BottomPolyhedron(BottomGen);
    BottomPolyhedron.verbose=verbose;
    BottomPolyhedron.do_extreme_rays=true;
    BottomPolyhedron.keep_order = true;
    try {     
        BottomPolyhedron.dualize_cone();  // includes finding extreme rays
    } catch(const NonpointedException& ){};

    // transfer pointedness
    assert( BottomPolyhedron.isComputed(ConeProperty::IsPointed) );
    pointed = BottomPolyhedron.pointed;
    is_Computed.set(ConeProperty::IsPointed);

    // BottomPolyhedron.Support_Hyperplanes.pretty_print(cout);

    help.resize(dim);

    // find extreme rays of Bottom among the generators
    vector<key_t> BottomExtRays;
    for(size_t i=0;i<nr_gen;++i)
        if(BottomPolyhedron.Extreme_Rays_Ind[i+nr_gen])
            BottomExtRays.push_back(i);
    /* vector<key_t> BottomExtRays; // can be used if the bool vector should not exist anymore
    size_t start_search=0;
    for(size_t i=0;i<ExtStrahl.nr_of_rows();++i){
        if(BottomPolyhedron.ExtStrahl[i][dim]==1){
            BottomPolyhedron.ExtStrahl[i].resize(dim);
            for(size_t j=0;j<nr_gen;++j){
                size_t k=(j+start_search) % nr_gen;
                if(BottomPolyhedron.ExtStrahl[i]==Generators[k]){
                    BottomExtRays.push_back(k);
                    start_search++;
                }
            }
        }
    }*/

    if(verbose)
        verboseOutput() << "Bottom has " << BottomExtRays.size() << " extreme rays" << endl;
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
 
    Matrix<Integer> BottomFacets(0,dim);
    vector<Integer> BottomDegs(0,dim);
    if (!isComputed(ConeProperty::SupportHyperplanes)) {
        Support_Hyperplanes = Matrix<Integer>(0,dim);
        nrSupport_Hyperplanes=0;
    }
    for(size_t i=0;i<BottomPolyhedron.nrSupport_Hyperplanes;++i){
        Integer test=BottomPolyhedron.Support_Hyperplanes[i][dim];
        for(size_t j=0;j<dim;++j)
            help[j]=BottomPolyhedron.Support_Hyperplanes[i][j];     
        if(test==0 && !isComputed(ConeProperty::SupportHyperplanes)){
            Support_Hyperplanes.append(help);
            nrSupport_Hyperplanes++;
        }
        if (test < 0){ 
            BottomFacets.append(help);
            BottomDegs.push_back(-test);
        }
    }
    
    is_Computed.set(ConeProperty::SupportHyperplanes);
    
    if (!pointed)
        throw NonpointedException();
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION

    vector<key_t> facet;
    for(size_t i=0;i<BottomFacets.nr_of_rows();++i){
        facet.clear();
        for(unsigned int & BottomExtRay : BottomExtRays)
            if(v_scalar_product(Generators[BottomExtRay],BottomFacets[i])==BottomDegs[i])
                facet.push_back(BottomExtRay);
        Pyramids[0].push_back(facet);
        nrPyramids[0]++;
    }
    if(verbose)
        verboseOutput() << "Bottom decomposition computed, " << nrPyramids[0] << " subcones" << endl;
}

template<typename Integer>
void Full_Cone<Integer>::start_message() {

       if (verbose) {
        verboseOutput()<<"************************************************************"<<endl;
        verboseOutput()<<"starting full cone computation" << endl;
    }   
}

template<typename Integer>
void Full_Cone<Integer>::end_message() {
    
       if (verbose) {
        verboseOutput() << "------------------------------------------------------------"<<endl;
    }   
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::build_top_cone() {
    
    primal_algorithm_initialize();

    OldCandidates.verbose=verbose;
    NewCandidates.verbose=verbose;
    
    if(dim==0)
        return;

    if( ( !do_bottom_dec || deg1_generated || dim==1 || (!do_triangulation && !do_partial_triangulation))) {        
        build_cone();
    }
    else{
        find_bottom_facets();
        start_from=nr_gen;
        deg1_triangulation=false;

        vector<list<vector<key_t> >::iterator > level0_order;
        level0_order.reserve(nrPyramids[0]);
        auto p=Pyramids[0].begin();
        for(size_t k=0;k<nrPyramids[0];++k,++p){
            level0_order.push_back(p);
        }
        for(size_t k=0;k<5*nrPyramids[0];++k){
            swap(level0_order[rand()%nrPyramids[0]],level0_order[rand()%nrPyramids[0]]);
        }
        list<vector<key_t> > new_order;
        for(size_t k=0;k<nrPyramids[0];++k){
            new_order.push_back(*level0_order[k]);
        }
        Pyramids[0].clear();
        Pyramids[0].splice(Pyramids[0].begin(),new_order);

    }   

    // try_offload(0); // superfluous since tried immediately in evaluate_stored_pyramids(0)
    
    evaluate_stored_pyramids(0);  // force evaluation of remaining pyramids

#ifdef NMZ_MIC_OFFLOAD
    if (_Offload_get_device_number() < 0) // dynamic check for being on CPU (-1)
    {
        evaluate_stored_pyramids(0);  // previous run could have left over pyramids
        mic_offloader.evaluate_triangulation();
    }
#endif // NMZ_MIC_OFFLOAD   

}


//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::check_evaluation_buffer(){

    return(omp_get_level()==omp_start_level && check_evaluation_buffer_size());
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::check_evaluation_buffer_size(){

    return(!Top_Cone->keep_triangulation && 
               Top_Cone->TriangulationBufferSize > EvalBoundTriang);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::transfer_triangulation_to_top(){

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
        tn = omp_get_ancestor_thread_num(omp_start_level+1);
  
    auto pyr_simp=TriangulationBuffer.begin();
    while (pyr_simp!=TriangulationBuffer.end()) {
        if (pyr_simp->height == 0) { // it was marked to be skipped
            Top_Cone->FS[tn].splice(Top_Cone->FS[tn].end(), TriangulationBuffer, pyr_simp++);
            --TriangulationBufferSize;
        } else {
            for (i=0; i<dim; i++)  // adjust key to topcone generators
                pyr_simp->key[i]=Top_Key[pyr_simp->key[i]];
            sort(pyr_simp->key.begin(),pyr_simp->key.end());
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
template<typename Integer>
void Full_Cone<Integer>::get_supphyps_from_copy(bool from_scratch, bool with_extreme_rays){

    if(isComputed(ConeProperty::SupportHyperplanes)){ // we have them already
        if(with_extreme_rays)
            compute_extreme_rays();        
        return;
    }
    
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
        if(isComputed(ConeProperty::ExtremeRays)){
            copy.is_Computed.set(ConeProperty::ExtremeRays);
            with_extreme_rays=false;
        }
        copy.GensInCone=GensInCone;
        copy.nrGensInCone=nrGensInCone;
        copy.Comparisons=Comparisons;
        if(!Comparisons.empty())
            copy.nrTotalComparisons=Comparisons[Comparisons.size()-1];
        
        typename list< FACETDATA<Integer> >::const_iterator l=Facets.begin();
        
        for(size_t i=0;i<old_nr_supp_hyps;++i){
            copy.Facets.push_back(*l);
            ++l;
        }
    }
    
    copy.dualize_cone();
    if(with_extreme_rays){
        copy.do_extreme_rays=true;
        copy.compute();
        Extreme_Rays_Ind=copy.Extreme_Rays_Ind;
        is_Computed.set(ConeProperty::ExtremeRays);
    }
    
    std::swap(Support_Hyperplanes,copy.Support_Hyperplanes);
    nrSupport_Hyperplanes = copy.nrSupport_Hyperplanes;
    is_Computed.set(ConeProperty::SupportHyperplanes);
    do_all_hyperplanes = false;
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::update_reducers(bool forced){
    
    if((!do_Hilbert_basis || do_module_gens_intcl) && !forced)
        return;

    if(NewCandidates.Candidates.empty())
        return;
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    if(hilbert_basis_rec_cone_known){
        NewCandidates.sort_by_deg();
        NewCandidates.reduce_by(HBRC);
        ModuleGensDepot.merge(NewCandidates);
        return;
    }

    if(nr_gen==dim)  // no global reduction in the simplicial case
        NewCandidates.sort_by_deg(); 
    if(nr_gen!=dim || forced){  // global reduction in the nonsimplicial case (or forced)
        NewCandidates.auto_reduce();
        if(verbose){
            verboseOutput() << "reducing " << OldCandidates.Candidates.size() << " candidates by " << NewCandidates.Candidates.size() << " reducers" << endl;
        }
        OldCandidates.reduce_by(NewCandidates);
    }
    OldCandidates.merge(NewCandidates);
    CandidatesSize=OldCandidates.Candidates.size();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::prepare_old_candidates_and_support_hyperplanes(){

    if(!isComputed(ConeProperty::SupportHyperplanes)){
        if (verbose) {
            verboseOutput() << "**** Computing support hyperplanes for reduction:" << endl;
        }
        get_supphyps_from_copy(false);
    }
    
    check_pointed();
    if(!pointed){
        throw NonpointedException();
    }

    int max_threads = omp_get_max_threads();
    size_t Memory_per_gen=8*nrSupport_Hyperplanes;
    size_t max_nr_gen=RAM_Size/(Memory_per_gen*max_threads);
    // cout << "max_nr_gen " << max_nr_gen << endl;
    AdjustedReductionBound=max_nr_gen;
    if(AdjustedReductionBound < 2000)
        AdjustedReductionBound=2000;


    Sorting=compute_degree_function();
    if (!is_approximation) {
        bool save_do_module_gens_intcl=do_module_gens_intcl;
        do_module_gens_intcl=false; // to avoid multiplying sort_deg by 2 for the original generators
                                    // sort_deg of new candiadtes will be multiplied by 2
                                    // so that all old candidates are tested for reducibility
        for (size_t i = 0; i <nr_gen; i++) {               
            // cout << gen_levels[i] << " ** " << Generators[i];
            if(!inhomogeneous || gen_levels[i]==0 || (!save_do_module_gens_intcl && gen_levels[i]<=1)) {
                OldCandidates.Candidates.push_back(Candidate<Integer>(Generators[i],*this));
                OldCandidates.Candidates.back().original_generator=true;
            }
        }
        for(size_t i=0;i<HilbertBasisRecCone.nr_of_rows();++i){
            HBRC.Candidates.push_back(Candidate<Integer>(HilbertBasisRecCone[i],*this));            
        }
        do_module_gens_intcl=save_do_module_gens_intcl; // restore
        if(HilbertBasisRecCone.nr_of_rows()>0){ // early enough to avoid multiplictaion of sort_deg by 2 for the elements 
                                               // in HilbertBasisRecCone
            hilbert_basis_rec_cone_known=true;
            HBRC.sort_by_deg();
        }
        if(!do_module_gens_intcl) // if do_module_gens_intcl we don't want to change the original monoid
            OldCandidates.auto_reduce();
        else
            OldCandidates.sort_by_deg();
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_triangulation(){

    // prepare reduction 
    if (do_Hilbert_basis && OldCandidates.Candidates.empty()) {
        prepare_old_candidates_and_support_hyperplanes();
    }
    
    if (TriangulationBufferSize == 0)
        return;
    
    assert(omp_get_level()==omp_start_level);

    const long VERBOSE_STEPS = 50;
    long step_x_size = TriangulationBufferSize-VERBOSE_STEPS;
    if (verbose) {
        verboseOutput() << "evaluating "<<TriangulationBufferSize<<" simplices" <<endl;
        /* verboseOutput() << "---------+---------+---------+---------+---------+"
                        << " (one | per 2%)" << endl;*/
    }
    
    totalNrSimplices += TriangulationBufferSize;
    
    if(do_Stanley_dec || keep_triangulation){ // in these cases sorting is necessary
        for(auto& simp : TriangulationBuffer)
            sort(simp.key.begin(),simp.key.end());                
    }    

    if(do_evaluation && !do_only_multiplicity) {
    
    deque<bool> done(TriangulationBufferSize,false);
    bool skip_remaining;
    std::exception_ptr tmp_exception;

    do{ // allows multiple run of loop below in case of interruption for the update of reducers
    
    skip_remaining=false;
    step_x_size = TriangulationBufferSize-VERBOSE_STEPS;

    #pragma omp parallel
    {
        auto s = TriangulationBuffer.begin();
        size_t spos=0;
        int tn = omp_get_thread_num();
        #pragma omp for schedule(dynamic) nowait
        for(size_t i=0; i<TriangulationBufferSize; i++){
            try {
                if(skip_remaining)
                    continue;
                
                for(; i > spos; ++spos, ++s) ;
                for(; i < spos; --spos, --s) ;

           INTERRUPT_COMPUTATION_BY_EXCEPTION
                
                if(done[spos])
                    continue;

                done[spos]=true;

                /* if(keep_triangulation || do_Stanley_dec)  -- now done above
                    sort(s->key.begin(),s->key.end()); */
                if(!SimplexEval[tn].evaluate(*s)){
                    #pragma omp critical(LARGESIMPLEX)
                    LargeSimplices.push_back(SimplexEval[tn]);
                }
                if (verbose) {
                    #pragma omp critical(VERBOSE)
                    while ((long)(i*VERBOSE_STEPS) >= step_x_size) {
                        step_x_size += TriangulationBufferSize;
                        verboseOutput() << "|" <<flush;
                    }
                }

                if(do_Hilbert_basis && Results[tn].get_collected_elements_size() > AdjustedReductionBound)
                    skip_remaining=true;
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
                #pragma omp flush(skip_remaining)
            }
        }
        Results[tn].transfer_candidates();
    } // end parallel
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);

    if (verbose)
        verboseOutput()  << endl;
        
    update_reducers();
        
    } while(skip_remaining);
        
    } // do_evaluation
            
    if (verbose)
    {
        verboseOutput() << totalNrSimplices << " simplices";
        if(do_Hilbert_basis)
            verboseOutput() << ", " << CandidatesSize << " HB candidates";
        if(do_deg1_elements)
            verboseOutput() << ", " << CandidatesSize << " deg1 vectors";
        verboseOutput() << " accumulated." << endl;
    }
    
    if (keep_triangulation) {
        Triangulation.splice(Triangulation.end(),TriangulationBuffer);
    } else {
        // #pragma omp critical(FREESIMPL)
        FreeSimpl.splice(FreeSimpl.begin(),TriangulationBuffer);
    }
    TriangulationBufferSize=0;

    if (verbose && use_bottom_points) {
        size_t lss=LargeSimplices.size();
        if(lss>0)
            verboseOutput() << lss << " large simplices stored" << endl;
    }

    for(size_t i=0;i<Results.size();++i)
        Results[i].transfer_candidates(); // any remaining ones
    
    update_reducers();
}

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::evaluate_triangulation(){

    assert(omp_get_level()==0);
        
    if (TriangulationBufferSize == 0)
        return;
    
    totalNrSimplices += TriangulationBufferSize;
    
    if(do_determinants){
 
        bool dummy;
        bool skip_remaining=false;
    std::exception_ptr tmp_exception;

        long nr_simplices_done=0;
        #pragma omp parallel
        {
        Matrix<renf_elem_class> work;
        auto t=TriangulationBuffer.begin();
        size_t spos=0;
        #pragma omp for
        for(size_t i=0; i<TriangulationBufferSize; i++){
            try {
            if(skip_remaining)
                continue;
            
            for(; i > spos; ++spos, ++t) ;
            for(; i < spos; --spos, --t) ;
            
            INTERRUPT_COMPUTATION_BY_EXCEPTION

            work=Generators.submatrix(t->key);
            work.row_echelon_inner_elem(dummy);
            t->vol=1;
            for(size_t i=0;i<dim;++i)
                t->vol*=work[i][i];
            
            t->vol_for_detsum=t->vol;
            
            if(do_multiplicity){
                renf_elem_class deg_prod=1;
                for(size_t j=0;j<dim;++j)
                    deg_prod*=gen_degrees[t->key[j]];
                t->vol/=deg_prod;                
            }
            
            #pragma omp atomic
            nr_simplices_done++;
            
            if(verbose && nr_simplices_done%1000==0){
                #pragma omp critical(PROGRESS)
                verboseOutput() << nr_simplices_done << " simplices done" << endl;
            }
            
            } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
            }

        } // for
        
        } // parallel
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
        
        auto t=TriangulationBuffer.begin();
        for(;t!=TriangulationBuffer.end();++t){ 
            
            INTERRUPT_COMPUTATION_BY_EXCEPTION
    
            t->vol=Iabs(t->vol);    
            // t->vol=Generators.submatrix(t->key).vol();
            detSum+=Iabs(t->vol_for_detsum);
            if(do_multiplicity){
                renf_multiplicity+=t->vol;                    
            }            
        }
    }
    
    if (keep_triangulation) {
        Triangulation.splice(Triangulation.end(),TriangulationBuffer);
    } else {
        // #pragma omp critical(FREESIMPL)
        FreeSimpl.splice(FreeSimpl.begin(),TriangulationBuffer);
    }
    TriangulationBufferSize=0;

}
#endif

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_large_simplices(){

    size_t lss = LargeSimplices.size();
    if (lss == 0)
        return;
    
    assert(omp_get_level()==omp_start_level);

    if (verbose) {
        verboseOutput() << "Evaluating " << lss << " large simplices" << endl;
    }
    size_t j;
    for (j = 0; j < lss; ++j) {
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        evaluate_large_simplex(j, lss);
    }

    // decomposition might have created new simplices  -- NO LONGER, now in Pyramids[0]
    // evaluate_triangulation();

    // also new large simplices are possible
    /* if (!LargeSimplices.empty()) {
        use_bottom_points = false;
        lss += LargeSimplices.size();
        if (verbose) {
            verboseOutput() << "Continue evaluation of " << lss << " large simplices without new decompositions of simplicial cones." << endl;
        }
        for (; j < lss; ++j) {
        
            INTERRUPT_COMPUTATION_BY_EXCEPTION
        
            evaluate_large_simplex(j, lss);
        }
    }*/ 
    assert(LargeSimplices.empty());

    for(size_t i=0;i<Results.size();++i)
        Results[i].transfer_candidates(); // any remaining ones

    update_reducers();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_large_simplex(size_t j, size_t lss) {
    if (verbose) {
        verboseOutput() << "Large simplex " << j+1 << " / " << lss << endl;
    }

    if (do_deg1_elements && !do_h_vector && !do_Stanley_dec && !deg1_triangulation) {
        compute_deg1_elements_via_projection_simplicial(LargeSimplices.front().get_key());
    }
    else {
        LargeSimplices.front().Simplex_parallel_evaluation();
        if (do_Hilbert_basis && Results[0].get_collected_elements_size() > AdjustedReductionBound) {
            Results[0].transfer_candidates();
            update_reducers();
        }
    }
    LargeSimplices.pop_front();
}

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::evaluate_large_simplex(size_t j, size_t lss) {
    assert(false);
}
#endif

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_deg1_elements_via_projection_simplicial(const vector<key_t>& key){
    
    Matrix<Integer> Gens=Generators.submatrix(key);
    Sublattice_Representation<Integer>  NewCoordinates=LLL_coordinates<Integer,Integer>(Gens);
    Matrix<Integer> Gred=NewCoordinates.to_sublattice(Gens);
    vector<Integer> GradT=NewCoordinates.to_sublattice_dual(Grading);
    
    Matrix<Integer> GradMat(0,dim);
    GradMat.append(GradT);
    Cone<Integer> ProjCone(Type::cone,Gred,Type::grading, GradMat);
    ConeProperties ForDeg1;
    ForDeg1.set(ConeProperty::Projection);
    ForDeg1.set(ConeProperty::NoLLL);
    if(using_GMP<Integer>())
            ForDeg1.set(ConeProperty::BigInt);
    ForDeg1.set(ConeProperty::Deg1Elements);
    ProjCone.compute(ForDeg1);
        
    /*if(using_GMP<Integer>())
        ProjCone.compute(ConeProperty::Projection,ConeProperty::NoLLL,ConeProperty::BigInt,);
    else
        ProjCone.compute(ConeProperty::Projection,ConeProperty::NoLLL);*/
    Matrix<Integer> Deg1=ProjCone.getDeg1ElementsMatrix();
    Deg1=NewCoordinates.from_sublattice(Deg1);   
    
    Matrix<Integer> Supp=ProjCone.getSupportHyperplanesMatrix();
    Supp=NewCoordinates.from_sublattice_dual(Supp);
    
    /*for(size_t i=0;i<dim;++i)
        for(size_t j=0;j<dim;++j)
            assert(v_scalar_product(Supp[i],Gens[j])>=0); */         
    
    vector<bool> Excluded(dim,false); // we want to discard duplicates
    for(size_t i=0;i<dim;++i){
        Integer test=v_scalar_product(Supp[i],Order_Vector);
        if(test>0)
            continue;
        if(test<0){
            Excluded[i]=true;
            continue;
        }
        size_t j;
        for(j=0;j<dim;++j){
            if(Supp[i][j]!=0)
                break;        
        }
        if(Supp[i][j]<0)
            Excluded[i]=true;        
    }
    
    for(const auto& E : Deg1.get_elements()){
        size_t i;
        for(i=0;i<dim;++i)
            if(v_scalar_product(E,Supp[i])==0 && Excluded[i])
                break;
        if(i<dim)
            continue;
            
        for(i=0;i<dim;++i)  // exclude original generators
            if(E==Gens[i])
                 break;
        if(i==dim){    
            Results[0].Deg1_Elements.push_back(E);
            Results[0].collected_elements_size++;
        }        
    }
    Results[0].transfer_candidates();
}

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::compute_deg1_elements_via_projection_simplicial(const vector<key_t>& key){
    assert(false);
}
#endif   

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::remove_duplicate_ori_gens_from_HB(){
    
return; //TODO reactivate!

    set<vector<Integer> > OriGens;
    for(auto c=OldCandidates.Candidates.begin();c!=OldCandidates.Candidates.end();){
        if(!c->original_generator){
            ++c;
            continue;
        }
        auto found=OriGens.find(c->cand);
        if(found!=OriGens.end()){
            c=OldCandidates.Candidates.erase(c);
        }
        else{
            if(c->original_generator)
                OriGens.insert(c->cand);
            ++c;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm(){
    
    if( !(do_deg1_elements || do_Hilbert_basis || do_h_vector || do_multiplicity || do_determinants || do_triangulation_size) )
        return;

    // primal_algorithm_initialize();

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

template<typename Integer>
void Full_Cone<Integer>::set_primal_algorithm_control_variables(){

    do_triangulation = false;
    do_partial_triangulation = false;
    stop_after_cone_dec=false;
    do_evaluation = false;
    do_only_multiplicity=false;
    use_bottom_points = true;
    triangulation_is_nested = false;
    triangulation_is_partial = false;

    if(do_multiplicity)     do_determinants=true;
    if (do_determinants)    do_triangulation = true;
    if(do_triangulation_size) do_triangulation = true;
    if (do_h_vector)        do_triangulation = true;
    if (do_deg1_elements)   do_partial_triangulation = true;
    if (do_Hilbert_basis)   do_partial_triangulation = true;
    // activate 
    do_only_multiplicity = do_determinants || do_multiplicity;

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
    
    // cout << "DOM " << do_only_multiplicity << " Tri " << do_triangulation << " Wit " << do_integrally_closed << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm_initialize() {
    
    set_primal_algorithm_control_variables();
    
    if(verbose){
        verboseOutput() << "Starting primal algorithm ";    
        if (do_partial_triangulation) 
            verboseOutput()<<"with partial triangulation ";
         if (do_triangulation)
                verboseOutput()<<"with full triangulation ";
         if (!do_triangulation && !do_partial_triangulation) verboseOutput()<<"(only support hyperplanes) ";
         verboseOutput()<<"..."<<endl;
    }
    
    prepare_inclusion_exclusion();

    SimplexEval = vector< SimplexEvaluator<Integer> >(omp_get_max_threads(),SimplexEvaluator<Integer>(*this));
    for(size_t i=0;i<SimplexEval.size();++i)
        SimplexEval[i].set_evaluator_tn(i);
    Results = vector< Collector<Integer> >(omp_get_max_threads(),Collector<Integer>(*this));
    Hilbert_Series.setVerbose(verbose);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm_finalize() {

    if (isComputed(ConeProperty::Grading) && !deg1_generated) {
        deg1_triangulation = false;
    }
    if (keep_triangulation) {
        is_Computed.set(ConeProperty::Triangulation);
    }
    if (do_cone_dec) {
        is_Computed.set(ConeProperty::ConeDecomposition);
    }

    evaluate_triangulation();
    assert(nrPyramids[0]==0);
    evaluate_large_simplices(); // can produce level 0 pyramids
    use_bottom_points=false; // block new attempts for subdivision
    evaluate_stored_pyramids(0); // in case subdivision took place
    evaluate_triangulation();
    FreeSimpl.clear();
    
    // collect accumulated data from the SimplexEvaluators
    for (int zi=0; zi<omp_get_max_threads(); zi++) {
        detSum += Results[zi].getDetSum();
        multiplicity += Results[zi].getMultiplicitySum();
        if (do_h_vector) {
            Hilbert_Series += Results[zi].getHilbertSeriesSum();
        }
    }
#ifdef NMZ_MIC_OFFLOAD
    // collect accumulated data from mics
    if (_Offload_get_device_number() < 0) // dynamic check for being on CPU (-1)
    {
        mic_offloader.finalize();
    }
#endif // NMZ_MIC_OFFLOAD
    if (do_h_vector) {
        Hilbert_Series.collectData();
    }
    
    if(verbose) {
        verboseOutput() << "Total number of pyramids = "<< totalNrPyr << ", among them simplicial " << nrSimplicialPyr << endl;
        // cout << "Uni "<< Unimod << " Ht1NonUni " << Ht1NonUni << " NonDecided " << NonDecided << " TotNonDec " << NonDecidedHyp<< endl;
        if(do_only_multiplicity)
            verboseOutput() << "Determinants computed = " << TotDet << endl;
        /* if(NrCompVect>0)
            cout << "Vector comparisons " << NrCompVect << " Value comparisons " << NrCompVal 
                    << " Average " << NrCompVal/NrCompVect+1 << endl; */
    }
    
    if(verbose && GMP_hyp+GMP_scal_prod+GMP_mat>0)
        verboseOutput() << "GMP transitions: matrices " << GMP_mat << " hyperplanes " << GMP_hyp << " vector operations " << GMP_scal_prod << endl; 

}
    
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::make_module_gens(){

    if(!inhomogeneous){
        NewCandidates.extract(ModuleGeneratorsOverOriginalMonoid);
        vector<Integer> Zero(dim,0);
        ModuleGeneratorsOverOriginalMonoid.push_front(Zero);
        // cout << "Mod " << endl;
        // Matrix<Integer>(ModuleGeneratorsOverOriginalMonoid).pretty_print(cout);
        // cout << "--------" << endl;
        is_Computed.set(ConeProperty::ModuleGeneratorsOverOriginalMonoid,true);
        return;
    }
    
    CandidateList<Integer> Level1OriGens;
    for(size_t i=0;i<nr_gen;++i){
            if(gen_levels[i]==1){
                Level1OriGens.push_back(Candidate<Integer>(Generators[i],*this));    
            }
    }
    CandidateList<Integer> Level1Generators=Level1OriGens;
    Candidate<Integer> new_cand(dim,Support_Hyperplanes.nr_of_rows());
    for(const auto& lnew : NewCandidates.Candidates){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        Integer level=v_scalar_product(lnew.cand,Truncation);
        if(level==1){
            new_cand=lnew;
            Level1Generators.reduce_by_and_insert(new_cand,OldCandidates);
        }
        else{
            for(const auto& l1 : Level1OriGens.Candidates){
                new_cand=sum(l1,lnew);
                Level1Generators.reduce_by_and_insert(new_cand,OldCandidates);
            }
        }        
    }
    Level1Generators.extract(ModuleGeneratorsOverOriginalMonoid);
    ModuleGeneratorsOverOriginalMonoid.sort();
    ModuleGeneratorsOverOriginalMonoid.unique();
    is_Computed.set(ConeProperty::ModuleGeneratorsOverOriginalMonoid,true);
    
    for (size_t i = 0; i <nr_gen; i++) { // the level 1 input generators have not yet ben inserted into OldCandidates              
        if(gen_levels[i]==1) {          // but they are needed for the truncated Hilbert basis com�putation
            NewCandidates.Candidates.push_back(Candidate<Integer>(Generators[i],*this));
            NewCandidates.Candidates.back().original_generator=true;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::make_module_gens_and_extract_HB(){
       
    make_module_gens();
    
    NewCandidates.divide_sortdeg_by2(); // was previously multplied by 2    
    NewCandidates.sort_by_deg();
    
    OldCandidates.merge(NewCandidates);
    OldCandidates.auto_reduce(); 
}


//---------------------------------------------------------------------------

 template<typename Integer>
 void Full_Cone<Integer>::finish_Hilbert_series(){
    
    if (do_h_vector) {
        Hilbert_Series.setShift(convertTo<long>(shift));
        Hilbert_Series.adjustShift();
        // now the shift in the HilbertSeries may change and we would have to adjust
        // the shift, the grading and more in the Full_Cone to continue to add data!
            // COMPUTE HSOP here
        if (do_hsop){
            compute_hsop();
            is_Computed.set(ConeProperty::HSOP);
        }
        Hilbert_Series.simplify();
        is_Computed.set(ConeProperty::HilbertSeries);
    }
 }
 
template<>
 void Full_Cone<renf_elem_class>::finish_Hilbert_series(){
        assert(false);
 }

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm_set_computed() {

    extreme_rays_and_deg1_check();
    if(!pointed){
        throw NonpointedException();
    }

    if (do_triangulation || do_partial_triangulation) {
        is_Computed.set(ConeProperty::TriangulationSize,true);
        if (do_evaluation) {
            is_Computed.set(ConeProperty::TriangulationDetSum,true);
        }
    }
    if ( (do_triangulation && do_evaluation && isComputed(ConeProperty::Grading))
        || (do_multiplicity && using_renf<Integer>())    )
        is_Computed.set(ConeProperty::Multiplicity,true);
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
                
    if (do_Hilbert_basis) {
        if(hilbert_basis_rec_cone_known){
            // OldCandidates.Candidates.clear();
            OldCandidates.merge(HBRC);
            OldCandidates.merge(ModuleGensDepot);            
        }
        if(do_module_gens_intcl){
                make_module_gens_and_extract_HB();
        }
        OldCandidates.sort_by_val();
        OldCandidates.extract(Hilbert_Basis);
        OldCandidates.Candidates.clear();
        Hilbert_Basis.unique();
        is_Computed.set(ConeProperty::HilbertBasis,true);
    }
    
    if (isComputed(ConeProperty::Grading) && isComputed(ConeProperty::HilbertBasis)) {
        select_deg1_elements();
        check_deg1_hilbert_basis();
    }
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    if (do_deg1_elements) {
        for(size_t i=0;i<nr_gen;i++)
            if(v_scalar_product(Grading,Generators[i])==1 && (!(is_approximation || is_global_approximation) 
                            || subcone_contains(Generators[i])))
                Deg1_Elements.push_front(Generators[i]);
        is_Computed.set(ConeProperty::Deg1Elements,true);
        Deg1_Elements.sort();
        Deg1_Elements.unique();
    }
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    if(do_h_vector)
        finish_Hilbert_series();
    

    if(do_Stanley_dec){
        is_Computed.set(ConeProperty::StanleyDec);
    }
    
    // If the grading has gcd > 1 on the recession monoid,
    // we must multiply the multiplicity by it.
    // Without this correction, the multiplicity (relative to deg/g)
    // is divided by g^r, but it must be g^{r-1}.
    // We determine g and multiply by it.
    //
    // The reason behind this correction is that the determinants
    // are computed with respect to a basis in which the
    // basic simplex has volume 1/g instead of 1.
    // The correction above takes care of this "mistake"
    // that we are forced to make in order to keep data integral.

    if(isComputed(ConeProperty::Multiplicity)){        
        Integer corr_factor;       
        if(!inhomogeneous)
            corr_factor=v_gcd(Grading);
        if(inhomogeneous && level0_dim==0)
            corr_factor=1;
        if(inhomogeneous && level0_dim>0){
            Matrix<Integer> Level0Space=ProjToLevel0Quot.kernel();
            corr_factor=0;
            for(size_t i=0;i<Level0Space.nr_of_rows();++i) 
                corr_factor=libnormaliz::gcd(corr_factor,v_scalar_product(Grading,Level0Space[i]));           
        }
        multiplicity*=convertTo<mpz_class>(corr_factor);        
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::set_degrees() {

    // Generators.pretty_print(cout);
    // cout << "Grading " << Grading;
    if (gen_degrees.size() != nr_gen && isComputed(ConeProperty::Grading)) // now we set the degrees
    {
        gen_degrees.resize(nr_gen);
        if(do_h_vector || (!using_GMP<Integer>() && !using_renf<Integer>()) )
                gen_degrees_long.resize(nr_gen);
        gen_degrees=Generators.MxV(Grading);
        for (size_t i=0; i<nr_gen; i++) {
            if (gen_degrees[i] <= 0) {
                throw BadInputException("Grading gives non-positive value "
                        + toString(gen_degrees[i])
                        + " for generator " + toString(i+1) + ".");
            }
            if(do_h_vector || (!using_GMP<Integer>() && !using_renf<Integer>()))
                convert(gen_degrees_long[i], gen_degrees[i]);
                
        }
    }
}

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::set_degrees() {
    
    if(!isComputed(ConeProperty::Grading) && !inhomogeneous)
        return;
    
    vector<renf_elem_class> GradHelp=Grading;
    if(inhomogeneous)
        GradHelp=Truncation;

    gen_degrees=Generators.MxV(GradHelp);
    for(size_t i=0;i<Generators.nr_of_rows();++i)
        if(gen_degrees[i]<=0 && (do_multiplicity || do_automorphisms))
            throw BadInputException("Volume or automorphism group not computable for unbounded nalgebraic polyhedra");
    
}
#endif
   
//---------------------------------------------------------------------------
// Normaliz modes (public)
//---------------------------------------------------------------------------

// check the do_* bools, they must be set in advance
// this method (de)activates them according to dependencies between them
template<typename Integer>
void Full_Cone<Integer>::set_implications() {

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
    if (do_module_gens_intcl) use_bottom_points= false; // extra bottom points change the originalmonoid
    if (do_Stanley_dec)       keep_triangulation = true;
    if (do_cone_dec)          keep_triangulation = true;
    if (keep_triangulation)   do_determinants = true;
    // if (do_multiplicity)    do_determinants = true; // removed because of automorphisms
    if ((do_multiplicity || do_h_vector) && inhomogeneous)    do_module_rank = true;
    if (do_Hilbert_basis)     do_deg1_elements = false; // after the Hilbert basis computation, deg 1 elements will be extracted
    
    no_descent_to_facets=do_h_vector || do_module_gens_intcl || keep_triangulation  // we must use the primal algorithm directly
                || do_triangulation_size || do_Stanley_dec || do_cone_dec || do_determinants || do_excluded_faces || do_bottom_dec;
                           
    do_only_supp_hyps_and_aux= !no_descent_to_facets && !do_multiplicity && !do_deg1_elements &&
                !do_Hilbert_basis &&!do_subdivision_points;
}

// We set the do* variables to false if the corresponding task has been done
template<typename Integer>
void Full_Cone<Integer>::deactivate_completed_tasks() {

    if(isComputed(ConeProperty::IsPointed))                    do_pointed=false;    
    if(isComputed(ConeProperty::ExtremeRays))                  do_extreme_rays=false;
    if(isComputed(ConeProperty::HilbertBasis)){                
        do_Hilbert_basis=false;
        do_integrally_closed=false;
    }
    if(isComputed(ConeProperty::Deg1Elements))                 do_deg1_elements=false;
    if(isComputed(ConeProperty::ModuleGeneratorsOverOriginalMonoid))      do_module_gens_intcl=false;
    
    if(isComputed(ConeProperty::HilbertSeries))                 do_h_vector=false;
    if(isComputed(ConeProperty::Multiplicity))                  do_multiplicity=false;
    
    if(isComputed(ConeProperty::StanleyDec))                   do_Stanley_dec=false;
    if(isComputed(ConeProperty::ConeDecomposition))             do_cone_dec=false;
    if(isComputed(ConeProperty::Triangulation))                 keep_triangulation=false;
    if(isComputed(ConeProperty::TriangulationDetSum))           do_determinants=false;   
       
    if(isComputed(ConeProperty::ModuleRank))                  do_module_rank=false;
    
    if(isComputed(ConeProperty::ClassGroup))                    do_class_group=false;
}

//---------------------------------------------------------------------------

// do computations using automorphisms
template<typename Integer>
void Full_Cone<Integer>::compute_by_automorphisms() {
    
    if((!exploit_automs_mult && !exploit_automs_vectors) || no_descent_to_facets)
        return;

    if(descent_level==0){
        if(do_Hilbert_basis){
            for(size_t i=0;i<nr_gen;++i)
                Generator_Set.insert(Generators[i]);
        }
        
        if(autom_codim_vectors<0) //set default values if not set by Cone
            autom_codim_vectors=1;
        if(autom_codim_mult<0)
            autom_codim_mult=min((int) dim/4,6);
    }
    
    if(exploit_automs_mult && do_multiplicity){
        if(descent_level< autom_codim_mult && nr_gen>= dim+4){ // otherwise direct computation
            if(inhomogeneous)
                compute_multiplicity_via_recession_cone();
            else    
                compute_multiplicity_via_automs();
        }
        is_Computed.set(ConeProperty::ExploitAutomsMult);
    }
    deactivate_completed_tasks();
    
    if(exploit_automs_vectors && do_Hilbert_basis){
        if(descent_level< autom_codim_vectors && nr_gen>= dim+4){ // otherwise direct computation
            compute_HB_via_automs();
        }
        is_Computed.set(ConeProperty::ExploitAutomsVectors);
    }
    deactivate_completed_tasks();
    
    if(exploit_automs_vectors && do_deg1_elements){
        if(descent_level< God_Father->autom_codim_mult && nr_gen>= dim+4){ // otherwise direct computation
            compute_Deg1_via_automs();
        }
        is_Computed.set(ConeProperty::ExploitAutomsVectors);
    }
    deactivate_completed_tasks();
}

size_t nr_revlex_simpl=0;

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::recursive_revlex_triangulation(vector<key_t> simplex_so_far, 
                const vector<key_t>& face_key, const vector<typename list<FACETDATA<Integer>>::const_iterator>& mother_facets,
                size_t dim ){
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    /* cout << "FACE KEY "<< face_key;
    cout << "SIMPLex " << simplex_so_far; */
    
    // handle simplex case first since no further descent is necessary
    
    if(face_key.size()==dim){
        simplex_so_far.insert(simplex_so_far.end(),face_key.begin(),face_key.end());
        nr_revlex_simpl++;
        if(nr_revlex_simpl % 10000==0){
            cout << "NR REVLEX SIMPL " << nr_revlex_simpl << endl;
        }
        return;
    }
    
    // We first find the support hyperplanes of our top cone that cut out the 
    // facets of our face
    
    vector<vector<bool> > facet_candidates;
    
    vector<typename list<FACETDATA<Integer>>::const_iterator> candidates_iterators;
    
        
    for(size_t i=0;i<mother_facets.size();++i){
        
        auto F=mother_facets[i];
        
        vector<bool> intersection(nr_gen);
        size_t nr_intersection=0;
        for(unsigned int j : face_key){
            if(F->GenInHyp[j]){
                intersection[j]=true;
                nr_intersection++;
            }
        }
        // cout << "NR " << nr_intersection << endl;
        if(nr_intersection<dim-1 || nr_intersection==face_key.size()) // too small or everything
            continue;        
        facet_candidates.push_back(intersection);
        candidates_iterators.push_back(F);
    }    

    vector<bool> the_facets(facet_candidates.size(),true);
    maximal_subsets(facet_candidates,the_facets);
    
    vector<typename list<FACETDATA<Integer>>::const_iterator> facets_of_this_face;
    for(size_t i=0;i<the_facets.size();++i)
        if(the_facets[i])
            facets_of_this_face.push_back(candidates_iterators[i]);
        
    // now we have the facets of our face via support hyperplanes of the top cone
        
    // Next we go over those facets that are opposite to next_vert 
    
    key_t next_vert=face_key[0];
    simplex_so_far.push_back(next_vert);

    for(size_t i=0;i<facets_of_this_face.size();++i){
        
        auto F=facets_of_this_face[i];
        
        if(F->GenInHyp[next_vert]) // want only those opposite to next_vert
            continue;
        vector<key_t> intersection;
        for(unsigned int j : face_key){
            if(F->GenInHyp[j])
                intersection.push_back(j);                
        }
        
        recursive_revlex_triangulation(simplex_so_far,intersection,facets_of_this_face,dim-1);
    }
}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::make_facets(){
    
    if(!isComputed(ConeProperty::SupportHyperplanes))
        support_hyperplanes();
    assert(Facets.empty());
    for(size_t i=0;i<Support_Hyperplanes.nr_of_rows();++i){
        FACETDATA<Integer> NewFacet; NewFacet.Hyp.resize(dim); NewFacet.GenInHyp.resize(nr_gen);
        for(size_t j=0;j<nr_gen;++j)
            if(v_scalar_product(Support_Hyperplanes[i],Generators[j])==0)
                NewFacet.GenInHyp[j]=true;
        NewFacet.Hyp=Support_Hyperplanes[i];
        Facets.push_back(NewFacet);     
    }
}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::revlex_triangulation(){
    
    make_facets();
    compute_extreme_rays(true);
    vector<key_t> simplex_so_far;
    simplex_so_far.reserve(dim);
    boost::dynamic_bitset<> face_ind(nr_gen);
    vector<key_t> Extreme_Rays_Key;
    for(size_t i=0;i< nr_gen;++i)
        if(Extreme_Rays_Ind[i])
            Extreme_Rays_Key.push_back(i);
        
    vector<typename list<FACETDATA<Integer>>::const_iterator> mother_facets;
    
    typename list<FACETDATA<Integer>>::const_iterator F;
    for(F=Facets.begin();F!=Facets.end();++F)
        mother_facets.push_back(F);
    
    recursive_revlex_triangulation(simplex_so_far,Extreme_Rays_Key,mother_facets,dim); 
    
            cout << "FINAL NR REVLEX SIMPL " << nr_revlex_simpl << endl;
            
    exit(0);
}


//---------------------------------------------------------------------------
// general purpose compute method
// do_* bools must be set in advance, this method does sanity checks for it
// if no bool is set it does support hyperplanes and extreme rays
template<typename Integer>
void Full_Cone<Integer>::compute() {
    
    omp_start_level=omp_get_level();
    
    if(dim==0){
        set_zero_cone();
        return;
    }


    set_implications();
    start_message();
    
    if(!do_Hilbert_basis && !do_h_vector && !do_multiplicity && !do_deg1_elements
        && !do_Stanley_dec && !keep_triangulation && !do_determinants)
        assert(Generators.max_rank_submatrix_lex().size() == dim);

    minimize_support_hyperplanes(); // if they are given
    if (inhomogeneous)
        set_levels();
    
    check_given_grading();
    // look for a grading if it is needed
    find_grading();        
    if(isComputed(ConeProperty::IsPointed) && !pointed){
        end_message();
        return;
    }
    if (!isComputed(ConeProperty::Grading))
        disable_grading_dep_comp();
    
    // revlex_triangulation(); was here for test 
    
    if (do_only_supp_hyps_and_aux || (Grading.size()>0 && !isComputed(ConeProperty::Grading))){
            // in the last case there are only two possibilities:
            // either nonpointed or bad grading

        // primal_algorithm_initialize();
        support_hyperplanes();
        compute_class_group();
        compute_automorphisms();
        deactivate_completed_tasks();
        end_message(); 
        return;
    }
       
    set_degrees();
    sort_gens_by_degree(true);
    

    bool polyhedron_is_polytope=inhomogeneous;
    if(inhomogeneous){
        find_level0_dim();
        for(size_t i=0;i<nr_gen;++i)
            if(gen_levels[i]==0){
                polyhedron_is_polytope=false;
                break;
            }                
    }
    if(polyhedron_is_polytope && (do_Hilbert_basis || do_h_vector)){ // inthis situation we must just find the 
            convert_polyhedron_to_polytope();                        // degree 1 points
            deactivate_completed_tasks();
    }

    if(do_approximation && !deg1_generated){        
        if(!isComputed(ConeProperty::ExtremeRays) || !isComputed(ConeProperty::SupportHyperplanes)){
            do_extreme_rays=true;
            dualize_cone(false);// no start or end message
        }
        assert(!(do_deg1_elements && do_subdivision_points));
        if(do_deg1_elements){
            if(!isComputed(ConeProperty::Deg1Elements)){
                if(verbose)
                    verboseOutput() << "Approximating rational by lattice polytope" << endl;
                compute_deg1_elements_via_approx_global();
                is_Computed.set(ConeProperty::Deg1Elements,true);
                deactivate_completed_tasks();
            }
        } 
        
        if(do_subdivision_points){ // now we want subdividing elements for a simplicial cone

            do_Hilbert_basis=true;
            compute_elements_via_approx(Hilbert_Basis);
            return; // in this case really done
        }        
    }
    
    compute_by_automorphisms();
    deactivate_completed_tasks();

    primal_algorithm();
    deactivate_completed_tasks();
        
    if(inhomogeneous && descent_level==0){
        find_module_rank();
    }
    
    compute_class_group();
    compute_automorphisms();
    deactivate_completed_tasks();
    
    end_message();  
    
}

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::compute() {
    
    if(dim==0){
        set_zero_cone();
        return;
    }
    
    assert(Truncation.size()==0 || Grading.size()==0);
    
    Norm=Truncation;
    if(Grading.size()>0)
        Norm=Grading;

    set_implications();
    set_degrees();

    start_message();
    
    if(!do_Hilbert_basis && !do_h_vector && !do_multiplicity && !do_deg1_elements
        && !do_Stanley_dec && !keep_triangulation && !do_determinants)
        assert(Generators.max_rank_submatrix_lex().size() == dim);

    minimize_support_hyperplanes(); // if they are given
    if (inhomogeneous)
        set_levels();
    
    check_given_grading();
    
    compute_by_automorphisms();

    if (do_only_supp_hyps_and_aux){
        support_hyperplanes();
        compute_automorphisms();
        deactivate_completed_tasks();
        end_message();
        return;
    }

    if(isComputed(ConeProperty::IsPointed) && !pointed){
        end_message();
        return;
    }

    sort_gens_by_degree(true);

    primal_algorithm(); 
    compute_automorphisms();
    deactivate_completed_tasks();

    end_message();
}
#endif

// compute the degree vector of a hsop
template<typename Integer>
vector<Integer> degrees_hsop(const vector<Integer> gen_degrees,const vector<size_t> heights){
    vector<Integer> hsop(heights.back());
    hsop[0]=gen_degrees[0];
    size_t k=1;
    while (k<heights.size() && heights[k]>heights[k-1]){
        hsop[k]=gen_degrees[k];
        k++;
    }
    size_t j=k;
    for (size_t i=k;i<heights.size();i++){
            if (heights[i]>heights[i-1]){
                hsop[j]=v_lcm_to(gen_degrees,k,i);
                j++;
            }
    }
    return hsop;
}

/*
//---------------------------------------------------------------------------
template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::copy_basic_data_from(const Full_Cone<Integer>& C){
    if(C.isComputed(ConeProperty::SupportHyperplanes)){
        Support_Hyperplanes=C.Support_Hyperplanes;
        nrSupport_Hyperplanes=C.nrSupport_Hyperplanes;
        is_Computed.set(ConeProperty::SupportHyperplanes);        
    }
    if(C.isComputed(ConeProperty::ExtremeRays)){
        Extreme_Rays_Ind=C.Extreme_Rays_Ind;
        is_Computed.set(ConeProperty::ExtremeRays);
    }
    if(C.isComputed(ConeProperty::Automorphisms)){
        Automs=C.Automs;
        is_Computed.set(ConeProperty::Automorphisms);
    }
    exploit_automorphisms=C.exploit_automorphisms;
    keep_order=true;
    verbose=C.verbose;
    descent_level=C.descent_level;


        Facet_2.Grading=Facet_Sub.to_sublattice_dual_no_div(Grading);
        Facet_2.is_Computed.set(ConeProperty::Grading);
        Facet_2.Mother=&(*this);
        Facet_2.God_Father=God_Father;
        Facet_2.do_multiplicity=true;

}
*/

//---------------------------------------------------------------------------
template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::push_supphyps_to_cone_over_facet(const vector<Integer>& fixed_point, const key_t facet_nr){
  
  Matrix<Integer> SuppHyps(0,dim);
  vector<Integer> Facet=Support_Hyperplanes[facet_nr];
  SuppHyps.append(Facet);
  Integer h=v_scalar_product(fixed_point,Facet);
  vector<Integer> NewFacet(dim);
  for(key_t i=0;i<nrSupport_Hyperplanes;++i){
    if(i==facet_nr)
        continue;
    Integer hN=v_scalar_product(fixed_point,Support_Hyperplanes[i]);
    NewFacet=FM_comb(Facet,h,Support_Hyperplanes[i],hN);
    SuppHyps.append(NewFacet);
  }
  return SuppHyps;
  
}


//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::copy_autom_params(const Full_Cone<Integer>& C){
    exploit_automs_mult=C.exploit_automs_mult;
    exploit_automs_vectors=C.exploit_automs_vectors;
    quality_of_automorphisms=C.quality_of_automorphisms;
    do_automorphisms=C.do_automorphisms;
    keep_order=true;
}
//---------------------------------------------------------------------------
// We want to replace the fixed point by a generator of the cone that has smaller height 
// over the base facet of the pyramid such that the fixed point is contained in the_facets
// pyramid with base the facet and apex the generator
/*
template<typename Integer>
vector<Integer> Full_Cone<Integer>::replace_fixed_point_by_generator(const vector<Integer>& fixed_point,
            const key_t facet_nr, const vector<Integer>& help_grading){

    Integer height_fixed_pt=v_scalar_product(Support_Hyperplanes[facet_nr],fixed_point);
    if(height_fixed_pt<=1)
        return fixed_point;

    Integer deg_fp=v_scalar_product(fixed_point,help_grading);
    Integer height_fp=v_scalar_product(fixed_point,Support_Hyperplanes[facet_nr]);
    
    bool first=true;
    Integer min_height;
    vector<Integer> min_ht_gen;
    
    for(size_t i=0;i<nr_gen;++i){
        Integer height_gen=v_scalar_product(Support_Hyperplanes[facet_nr],Generators[i]);
        Integer deg_gen=v_scalar_product(Generators[i],help_grading);
        if(deg_fp*height_gen<=deg_gen*height_fp)
            continue;
        vector<Integer> test=FM_comb(fixed_point,height_fp,Generators[i],height_gen,false);
        bool in_cone=true;
        for(size_t j=0;j<Support_Hyperplanes.nr_of_rows();++j){
            if(v_scalar_product(test,Support_Hyperplanes[j])<0){
                in_cone=false;
                break;
            }
        }
        if(!in_cone)
            continue;
        if(first || height_gen<min_height){
            first=false;
            min_ht_gen=Generators[i];
            min_height=height_gen;            
        }        
    }
    
    if(!first && min_height<height_fp)
        return min_ht_gen;
    else{
        cout << "No generator found" << endl;
        return fixed_point;
    }
}*/
//---------------------------------------------------------------------------
// version without iso classes
template<typename Integer>
void Full_Cone<Integer>::get_cone_over_facet_vectors(const vector<Integer>& fixed_point, const vector<key_t>& facet_key, 
                                      const key_t facet_nr, list<vector<Integer> >& Facet_vectors){
  
    vector<Integer> help_grading=compute_degree_function();
    
    Matrix<Integer> Facet_Gens(0,dim);
    // vector<Integer> selected_gen=replace_fixed_point_by_generator(fixed_point,facet_nr,help_grading);
         // cpuld be the fixed point
    Facet_Gens.append(fixed_point);
    Facet_Gens.append(Generators.submatrix(facet_key)); 
    
    if(verbose){
        verboseOutput() << "Finding Hilbert basis/deg 1 elements for cone over codim " << descent_level+1 <<" face" << endl;
         verboseOutput()  << "Height of pyramid apex  over face " << v_scalar_product(fixed_point,Support_Hyperplanes[facet_nr]) << endl;
    }
    
    Full_Cone ConeOverFacet(Facet_Gens);
    ConeOverFacet.verbose=verbose;

    if(isComputed(ConeProperty::Grading)){
      ConeOverFacet.Grading=Grading;
      ConeOverFacet.is_Computed.set(ConeProperty::Grading);
    }
    ConeOverFacet.descent_level=descent_level+1;
    ConeOverFacet.Mother=&(*this);
    ConeOverFacet.God_Father=God_Father;
    if(ConeOverFacet.descent_level < God_Father->autom_codim_vectors){ // otherwise dirct computation of HB
        ConeOverFacet.copy_autom_params(*this);
        ConeOverFacet.Embedding=Embedding;
    }
    ConeOverFacet.Support_Hyperplanes=push_supphyps_to_cone_over_facet(fixed_point,facet_nr);
    ConeOverFacet.do_Hilbert_basis=do_Hilbert_basis;
    ConeOverFacet.do_deg1_elements=do_deg1_elements;
    ConeOverFacet.inhomogeneous=inhomogeneous;
    if(inhomogeneous){
        ConeOverFacet.Truncation=Truncation;
    }
    ConeOverFacet.autom_codim_vectors=autom_codim_vectors;
    ConeOverFacet.compute();
    Facet_vectors.clear();
    if(do_Hilbert_basis)
        Facet_vectors.splice(Facet_vectors.begin(),ConeOverFacet.Hilbert_Basis);
    else
        Facet_vectors.splice(Facet_vectors.begin(),ConeOverFacet.Deg1_Elements); 
}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::compute_Deg1_via_automs(){
    
    compute_automorphisms(descent_level);
    
    if(!do_deg1_elements || isComputed(ConeProperty::Deg1Elements)
               || !isComputed(ConeProperty::Automorphisms) || Automs.getOrder() ==1)
        return;

    list<vector<Integer> > union_of_facets; // collects all candidates from the orbits of the HBs of the facets
    vector<Integer> fixed_point=get_fixed_point(descent_level);

    if(verbose){
        verboseOutput()  << "Computing deg1 elements via automorphisms in codim " << descent_level << endl;
        verboseOutput() << "Fixed point " << fixed_point;
    }
    
    vector<vector<key_t> > facet_keys = get_facet_keys_for_orbits(fixed_point,false);

    for(auto & facet_key : facet_keys){
        list<vector<Integer> > facet_Deg1;
        key_t facet_nr=facet_key.back();
        facet_key.pop_back();
        get_cone_over_facet_vectors(fixed_point,facet_key,facet_nr, facet_Deg1);
            
        list<vector<Integer> > union_of_orbits; // we must spread the deg 1 elements over their orbit
        for(const auto& c : facet_Deg1){
            list<vector<Integer> > orbit_of_deg1=Automs.orbit_primal(c);
            union_of_orbits.splice(union_of_orbits.end(),orbit_of_deg1);
        }
        union_of_orbits.sort();
        union_of_facets.merge(union_of_orbits);
    }
    union_of_facets.unique(); // necesary since dupocates cannot be avoided
    Deg1_Elements.splice(Deg1_Elements.begin(),union_of_facets);

    is_Computed.set(ConeProperty::Deg1Elements);     
}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::compute_HB_via_automs(){
    
    compute_automorphisms(descent_level);
    
    if(!do_Hilbert_basis || isComputed(ConeProperty::HilbertBasis)
               || !isComputed(ConeProperty::Automorphisms)  || Automs.getOrder() ==1)
        return;
    
    prepare_old_candidates_and_support_hyperplanes();

    set<vector<Integer> > union_of_facets; // collects all candidates from the orbits of the HBs of the facets
    vector<Integer> fixed_point=get_fixed_point(descent_level); // this is the number of cone points so far

    if(verbose){
        verboseOutput()  << "Computing Hilbert basis via automorphisms in codim " << descent_level << endl;
        verboseOutput() << "Fixed point " << fixed_point;
    }
    
    vector<vector<key_t> > facet_keys = get_facet_keys_for_orbits(fixed_point,false);

    for(auto & facet_key : facet_keys){
        list<vector<Integer> > facet_HB;
        key_t facet_nr=facet_key.back();
        facet_key.pop_back();
        get_cone_over_facet_vectors(fixed_point,facet_key,facet_nr, facet_HB);
        
        CandidateList<Integer> Cands_from_facet; // first we sort out the reducible elements
        for(const auto& jj : facet_HB)
            Cands_from_facet.reduce_by_and_insert(jj,*this,OldCandidates);
            
        // set<vector<Integer> > union_of_orbits; // we must spread the irreducibles over their orbit
        for(const auto& c : Cands_from_facet.Candidates){
            auto fc=union_of_facets.find(c.cand);
            if(fc!=union_of_facets.end())
                continue;
            list<vector<Integer> > orbit_of_cand=Automs.orbit_primal(c.cand);
            for(const auto& cc : orbit_of_cand)
                union_of_facets.insert(cc);
        }
    }
    cout << "Union unique size " << union_of_facets.size() << endl;
    for(const auto& v : union_of_facets)
        NewCandidates.push_back(Candidate<Integer>(v,*this));
    update_reducers(true); // we always want reduction
    OldCandidates.extract(Hilbert_Basis);
    Hilbert_Basis.sort();
    Hilbert_Basis.unique();

    is_Computed.set(ConeProperty::HilbertBasis);
    
    if (isComputed(ConeProperty::Grading)) {
        select_deg1_elements();
        check_deg1_hilbert_basis();
    }   
}

//---------------------------------------------------------------------------
template<typename Integer>
vector<Integer> Full_Cone<Integer>::get_fixed_point(size_t nr_cone_points){
// find fixed ppoint of low degree
        
    size_t mini=0;
    key_t min_orbit=0;
    for(size_t i=0;i<Automs.GenOrbits.size();++i)
        if((mini==0 || Automs.GenOrbits[i].size()<mini)
                && Automs.GenOrbits[i][0]>=nr_cone_points){
            mini=Automs.GenOrbits[i].size();
            min_orbit=i;
        }
    vector<Integer> fixed_point(dim);
    Matrix<Integer> Extreme_Rays=Generators.submatrix(Extreme_Rays_Ind);
        // Extreme_Rays.pretty_print(cout);
    for(size_t i=0;i<Automs.GenOrbits[min_orbit].size();++i){
        fixed_point=v_add(fixed_point,Extreme_Rays[Automs.GenOrbits[min_orbit][i]]);
    }
    v_make_prime(fixed_point);
    return fixed_point;
}

//---------------------------------------------------------------------------
template<typename Integer>
vector<vector<key_t> > Full_Cone<Integer>::get_facet_keys_for_orbits(const vector<Integer>& fixed_point,bool with_orbit_sizes){
// We collect only the facets that do not contain the fixed point.
// The last one (or two) entries of each key vector are abused for 
// (the orbit size and )  the number of the suport hyperplane.
// Everything for the first hyperplane in the orbit.

    vector<vector<key_t> > facet_keys;
    for(size_t k=0;k<Automs.LinFormOrbits.size();++k){
        key_t facet_nr=Automs.LinFormOrbits[k][0];
        assert(facet_nr<nrSupport_Hyperplanes); // for safety
        Integer ht=v_scalar_product(fixed_point,Support_Hyperplanes[facet_nr]);
        if(ht==0)   // fixed point in facet, does not contribute to multiplicity
            continue;
        vector<key_t> facet_gens;
        for(size_t i=0; i<nr_gen;++i){
            if(Extreme_Rays_Ind[i] && v_scalar_product(Generators[i],Support_Hyperplanes[facet_nr])==0)
                facet_gens.push_back(i);        
        }
        facet_keys.push_back(facet_gens);
        if(with_orbit_sizes)
            facet_keys.back().push_back(Automs.LinFormOrbits[k].size());
        facet_keys.back().push_back(facet_nr);
    }
    return facet_keys;
}
//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::compute_multiplicity_via_automs(){
    
    compute_automorphisms(0);
    
    if(!do_multiplicity || isComputed(ConeProperty::Multiplicity) || 
        !isComputed(ConeProperty::Grading) || !isComputed(ConeProperty::Automorphisms)  || Automs.getOrder() ==1)
        return;
    
    vector<Integer> fixed_point=get_fixed_point(0); // no cone points in this case
    Integer deg_fixed_point=v_scalar_product(fixed_point,Grading);

    vector<vector<key_t> > facet_keys = get_facet_keys_for_orbits(fixed_point,true);
    
    if(verbose){
        verboseOutput() << "Computing multiplicity via automorphisms in codim " << descent_level << endl;
        verboseOutput() << "Fixed point " << fixed_point;
    }
    
    for(auto & facet_key : facet_keys){
        key_t facet_nr=facet_key.back();
        facet_key.pop_back();
        Integer ht=v_scalar_product(fixed_point,Support_Hyperplanes[facet_nr]);
        long long orbit_size=facet_key.back();
        facet_key.pop_back();
        multiplicity+=convertTo<mpz_class>(orbit_size)*convertTo<mpz_class>(ht)
                        *facet_multiplicity(facet_key)/convertTo<mpz_class>(deg_fixed_point);
    }
    is_Computed.set(ConeProperty::Multiplicity);    
}

//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::compute_multiplicity_via_recession_cone(){
    Matrix<Integer> Level0Gens(0,dim);
    for(size_t i=0;i<nr_gen;++i){
        if(gen_levels[i]==0)
            Level0Gens.append(Generators[i]);
    }
    Sublattice_Representation<Integer> Level0Sub(Level0Gens,true);
    Matrix<Integer> RecGens=Level0Sub.to_sublattice(Level0Gens);
    Full_Cone<Integer> RecCone(RecGens);
    RecCone.Grading=Level0Sub.to_sublattice_dual_no_div(Grading);
    RecCone.is_Computed.set(ConeProperty::Grading);
    RecCone.do_multiplicity=true;
    RecCone.verbose=verbose;
    RecCone.copy_autom_params(*this);
    if(quality_of_automorphisms==AutomParam::ambient){
        RecCone.Embedding=Level0Sub.getEmbeddingMatrix().multiplication(Embedding);
    }
    RecCone.compute();
    multiplicity=RecCone.multiplicity;
    is_Computed.set(ConeProperty::Multiplicity);
}

//---------------------------------------------------------------------------
template<typename Integer>
mpq_class Full_Cone<Integer>::facet_multiplicity(const vector<key_t>& facet_key){
    
    Matrix<Integer> Facet_Gens=Generators.submatrix(facet_key);
    
        if(verbose){
        verboseOutput() << "Finding multiplicity for face of codim " << descent_level+1 <<  endl;
    }

    
    Sublattice_Representation<Integer> Facet_Sub(Facet_Gens,false);
    // By this choice we guarantee that the extreme Rays that generate the facet also
    // generate the lattice in which the multiplicity is computed.
    // This allows for more efficient isomorphism check.
    // The lattice can be smaller than the intersection of the facet with the full lattice.
    // We take care of this by mutiplying the computed multiplicity with the external index (see below).
    
    Matrix<Integer> Transformed_Facet_Gens=Facet_Sub.to_sublattice(Facet_Gens);
    Full_Cone Facet(Transformed_Facet_Gens);
    Facet.verbose=verbose;
    
    Facet.Grading=Facet_Sub.to_sublattice_dual_no_div(Grading);
    Facet.is_Computed.set(ConeProperty::Grading);
    Facet.Mother=&(*this);
    Facet.God_Father=God_Father;
    Facet.copy_autom_params(*this);
    if(quality_of_automorphisms==AutomParam::ambient){
        Facet.Embedding=Facet_Sub.getEmbeddingMatrix().multiplication(Embedding);
    }
    Facet.inhomogeneous=inhomogeneous;
    if(inhomogeneous){
        Facet.Truncation=Facet_Sub.to_sublattice_dual_no_div(Truncation);
    }
    Facet.descent_level=descent_level+1;
    Facet.keep_order=true;
    Facet.Support_Hyperplanes=Facet_Sub.to_sublattice_dual(Support_Hyperplanes);
    Facet.compute();
    bool found;
    const IsoType<Integer>& face_class=God_Father->FaceClasses.find_type(Facet,found);
    if(found){
        if(verbose){
            verboseOutput() << "Found isomorphism class" << endl;
        }
        mpq_class mmm=face_class.getMultiplicity();       
        return mmm*Facet_Sub.getExternalIndex();        
    } else{
        if(verbose){
            verboseOutput() << "New isomorphism class" << endl;
        }
        Full_Cone Facet_2(Transformed_Facet_Gens);
        Facet_2.Automs=Facet.Automs;
        Facet_2.is_Computed.set(ConeProperty::Automorphisms);
        Facet_2.Extreme_Rays_Ind=Facet.Extreme_Rays_Ind;
        Facet_2.is_Computed.set(ConeProperty::ExtremeRays);
        Facet_2.Support_Hyperplanes=Facet.Support_Hyperplanes;
        Facet_2.nrSupport_Hyperplanes=Facet.nrSupport_Hyperplanes;
        Facet_2.is_Computed.set(ConeProperty::SupportHyperplanes);
        Facet_2.copy_autom_params(*this);
        Facet_2.inhomogeneous=inhomogeneous;
        Facet_2.Truncation=Facet.Truncation;
        if(quality_of_automorphisms==AutomParam::ambient){
            Facet_2.Embedding=Facet_Sub.getEmbeddingMatrix().multiplication(Embedding);
        }
        Facet_2.verbose=verbose;
        Facet_2.descent_level=descent_level+1;
        Facet_2.Grading=Facet_Sub.to_sublattice_dual_no_div(Grading);
        Facet_2.is_Computed.set(ConeProperty::Grading);
        Facet_2.Mother=&(*this);
        Facet_2.God_Father=God_Father;
        Facet_2.do_multiplicity=true;
        Facet_2.verbose=true;
        Facet_2. autom_codim_mult=autom_codim_mult;
        Facet_2.compute();
        mpq_class mult_before=Facet_2.multiplicity;
        bool added;
        God_Father->FaceClasses.add_type(Facet_2, added);
        assert(mult_before==Facet_2.multiplicity);
        return Facet_2.multiplicity*Facet_Sub.getExternalIndex();
    }    
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_hsop(){
        vector<long> hsop_deg(dim,1);
        // if all extreme rays are in degree one, there is nothing to compute
        if (!isDeg1ExtremeRays()){
            if(verbose){
            verboseOutput() << "Computing heights ... " << flush;
            }
            
            vector<bool> choice = Extreme_Rays_Ind;
            if (inhomogeneous){
                for (size_t i=0; i<Generators.nr_of_rows(); i++) {
                    if (Extreme_Rays_Ind[i] && v_scalar_product(Generators[i],Truncation) != 0) {
                        choice[i]=false;
                    }
                }
            }
            Matrix<Integer> ER = Generators.submatrix(choice);
            Matrix<Integer> SH = getSupportHyperplanes();
            if (inhomogeneous){
                    Sublattice_Representation<Integer> recession_lattice(ER,true);
                    Matrix<Integer> SH_raw = recession_lattice.to_sublattice_dual(SH);
                    Matrix<Integer> ER_embedded = recession_lattice.to_sublattice(ER);
                    Full_Cone<Integer> recession_cone(ER_embedded);
                    recession_cone.Support_Hyperplanes = SH_raw;
                    recession_cone.dualize_cone();
                    SH = recession_lattice.from_sublattice_dual(recession_cone.getSupportHyperplanes());
            }
            vector<size_t> ideal_heights(ER.nr_of_rows(),1);
            // the heights vector is clear in the simplicial case
            if (is_simplicial){
                    for (size_t j=0;j<ideal_heights.size();j++) ideal_heights[j]=j+1;
            } else {
                list<pair<boost::dynamic_bitset<> , size_t> > facet_list;
                list<vector<key_t> > facet_keys;
                vector<key_t> key;
                size_t d = dim;
                if (inhomogeneous) d = level0_dim;
                assert(d>0); // we want to use d-1
                for (size_t i=SH.nr_of_rows();i-->0;){
                    boost::dynamic_bitset<> new_facet(ER.nr_of_rows());
                    key.clear();
                    for (size_t j=0;j<ER.nr_of_rows();j++){
                        if (v_scalar_product(SH[i],ER[j])==0){
                            new_facet[new_facet.size()-1-j]=1;
                        } else {
                            key.push_back(j);
                        }
                    }
                    facet_list.push_back(make_pair(new_facet,d-1));
                    facet_keys.push_back(key);
                }
                facet_list.sort(); // should be sorted lex
                //~ cout << "FACETS:" << endl;
                //~ //cout << "size: " << facet_list.size() << " | " << facet_list << endl;
                //~ for (auto jt=facet_list.begin();jt!=facet_list.end();++jt){
                        //~ cout << jt->first << " | " << jt->second << endl;
                //~ }
                //cout << "facet non_keys: " << facet_keys << endl;
                heights(facet_keys,facet_list,ER.nr_of_rows()-1,ideal_heights,d-1);
            }
        if(verbose){
            verboseOutput() << "done." << endl;
            if(!inhomogeneous)
                assert(ideal_heights[ER.nr_of_rows()-1]==dim);
            else
                assert(ideal_heights[ER.nr_of_rows()-1]==level0_dim);
            verboseOutput() << "Heights vector: " << ideal_heights;   
        }
        vector<Integer> er_deg = ER.MxV(Grading);
        hsop_deg = convertTo<vector<long> >(degrees_hsop(er_deg,ideal_heights));
        } 
        if(verbose){
            verboseOutput() << "Degrees of HSOP: " << hsop_deg;   
        }
        Hilbert_Series.setHSOPDenom(hsop_deg);
}

template<>
void Full_Cone<renf_elem_class>::compute_hsop(){
    assert(false);
}



// recursive method to compute the heights
// TODO: at the moment: facets are a parameter. global would be better
template<typename Integer>
void Full_Cone<Integer>::heights(list<vector<key_t> >& facet_keys,list<pair<boost::dynamic_bitset<>,size_t> > faces, size_t index,vector<size_t>& ideal_heights,size_t max_dim){
    // since we count the index backwards, this is the actual nr of the extreme ray
    
    size_t ER_nr = ideal_heights.size()-index-1;
    //~ cout << "starting calculation for extreme ray nr " << ER_nr << endl;
    list<pair<boost::dynamic_bitset<>,size_t> > not_faces;
    auto face_it=faces.begin();
    for (;face_it!=faces.end();++face_it){
        if (face_it->first.test(index)){ // check whether index is set
            break;
        }
        // resize not_faces
        face_it->first.resize(index);
    }
    not_faces.splice(not_faces.begin(),faces,faces.begin(),face_it);
    //~ cout << "faces not containing it:" << endl;
    //~ for (auto jt=not_faces.begin();jt!=not_faces.end();++jt){
                    //~ cout << jt->first << " | " << jt->second << endl;
    //~ }
    //~ cout << "faces containing it:" << endl;
    //~ for (auto jt=faces.begin();jt!=faces.end();++jt){
                    //~ cout << jt->first << " | " << jt->second << endl;
    //~ }
    auto not_faces_it=not_faces.begin();
    // update the heights
    if (ER_nr>0){
        if (!not_faces.empty()){
            ideal_heights[ER_nr] = ideal_heights[ER_nr-1];
            // compute the dimensions of not_faces
            vector<bool> choice = Extreme_Rays_Ind;
            if (inhomogeneous){
                for (size_t i=0; i<Generators.nr_of_rows(); i++) {
                    if (Extreme_Rays_Ind[i] && v_scalar_product(Generators[i],Truncation) != 0) {
                        choice[i]=false;
                    }
                }
            }
            Matrix<Integer> ER = Generators.submatrix(choice);
            int tn;
            if(omp_get_level()==omp_start_level)
                tn=0;
            else  tn = omp_get_ancestor_thread_num(omp_start_level+1);
            Matrix<Integer>& Test = Top_Cone->RankTest[tn];
            vector<key_t> face_key;
            for (;not_faces_it!=not_faces.end();++not_faces_it){
                if (not_faces_it->second==0){ // dimension has not yet been computed
                    // generate the key vector
                    face_key.resize(0);
                    for (size_t i=0;i<not_faces_it->first.size();++i){
                        if (not_faces_it->first.test(i)){
                            face_key.push_back(ER.nr_of_rows()-1-i);
                        }
                    }
                    not_faces_it->second = Test.rank_submatrix(ER,face_key);
                }
                if (not_faces_it->second==max_dim) break;
            }
            if (not_faces_it==not_faces.end()) {
                --max_dim;
                ideal_heights[ER_nr] = ideal_heights[ER_nr-1]+1;
            }
        } else {
            ideal_heights[ER_nr] = ideal_heights[ER_nr-1]+1;
            --max_dim;
        }
    }
    // we computed all the heights
    if (index==0) return;
    // if inner point, we can skip now
    
    // take the union of all faces not containing the current extreme ray
    boost::dynamic_bitset<> union_faces(index);
    not_faces_it = not_faces.begin();
    for (;not_faces_it!=not_faces.end();++not_faces_it){
        union_faces |= not_faces_it->first; // take the union
    }
    //cout << "Their union: " << union_faces << endl;
    // the not_faces now already have a size one smaller
    union_faces.resize(index+1);
    list<pair<boost::dynamic_bitset<>,size_t> > new_faces;
    // delete all facets which only consist of the previous extreme rays
    auto facet_it=facet_keys.begin();
    size_t counter=0;
    while(facet_it!=facet_keys.end()){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        counter=0;
        for (size_t i=0;i<facet_it->size();i++){
            if (facet_it->at(i)<=ER_nr) continue;
            // now we only have new extreme rays
            counter = i;
            break;
        }
        size_t j=ER_nr+1;
        for (;j<ideal_heights.size();j++){
            if (facet_it->at(counter)!=j){ // facet contains the element j
                    break;
            } else if (counter < facet_it->size()-1) counter++;
        }
        if (j==ideal_heights.size()){
            facet_it = facet_keys.erase(facet_it);
        } else ++facet_it;
    }
    facet_it=facet_keys.begin();
    
    // main loop
    for (;facet_it!=facet_keys.end();++facet_it){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        // check whether the facet is contained in the faces not containing the generator
        // and the previous generators
        // and check whether the generator is in the facet    
        // check whether intersection with facet contributes
        bool not_containing_el =false;
        // bool whether the facet contains an element which is NOT in the faces not containing the current extreme ray
        bool containing_critical_el=false; 
        counter=0;
        //cout << "check critical for facet " << *it << endl;
        for (size_t i=0;i<facet_it->size();i++){
            if (facet_it->at(i)==ER_nr){
                not_containing_el = true;
            }
            if (facet_it->at(i)<=ER_nr && i<facet_it->size()-1) continue;
            counter=i; // now we have elements which are bigger than the current extreme ray
            if (not_containing_el){
                for (size_t j=ER_nr+1;j<ideal_heights.size();j++){
                    if (facet_it->at(counter)!=j){ // i.e. j is in the facet
                        if (!union_faces.test(ideal_heights.size()-1-j)){
                            containing_critical_el = true;
                            break;
                        }
                    } else if (counter<facet_it->size()-1) counter++;
                }
            }
            break;
        }
        if(not_containing_el && containing_critical_el){ //facet contributes
            //cout << "Taking intersections with the facet " << *facet_it << endl;
            face_it =faces.begin();
            for (;face_it!=faces.end();++face_it){
                boost::dynamic_bitset<> intersection(face_it->first);
                for (unsigned int i : *facet_it){
                    if (i>ER_nr) intersection.set(ideal_heights.size()-1-i,false);
                }
                intersection.resize(index);
                if (intersection.any()){
                    // check whether the intersection lies in any of the not_faces
                    not_faces_it = not_faces.begin();
                    for (;not_faces_it!=not_faces.end();++not_faces_it){
                            if (intersection.is_subset_of(not_faces_it->first)) break;
                    }
                    if (not_faces_it== not_faces.end()) new_faces.push_back(make_pair(intersection,0)); 
                }
            }
       }
    }
    // the new faces need to be sort in lex order anyway. this can be used to reduce operations
    // for subset checking
    new_faces.sort();
    auto outer_it = new_faces.begin();
    auto inner_it = new_faces.begin();
    for (;outer_it!=new_faces.end();++outer_it){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        // work with a not-key vector
        vector<key_t> face_not_key;
        for (size_t i=0;i<outer_it->first.size();i++){
            if (!outer_it->first.test(i)){
                face_not_key.push_back(i);
            }
        }
        inner_it=new_faces.begin();
        size_t i=0;
        while (inner_it!=outer_it){
            i=0;
            for (;i<face_not_key.size();++i){
                if (inner_it->first.test(face_not_key[i])) break; //inner_it has an element which is not in outer_it
            }
            if (i==face_not_key.size()){
                inner_it = new_faces.erase(inner_it); //inner_it is a subface of outer_it
            } else ++inner_it;
        }
    }
    new_faces.merge(not_faces);
    /*cout << "The new faces: " << endl;
    for (const auto& jt : new_faces){
        cout << jt.first << " | " << jt.second << endl;
    }*/
    
    heights(facet_keys,new_faces,index-1,ideal_heights,max_dim);
}



template<typename Integer>
void Full_Cone<Integer>::convert_polyhedron_to_polytope() {
    
    if(verbose){
        verboseOutput() << "Converting polyhedron to polytope" << endl;
        verboseOutput() << "Pointed since cone over polytope" << endl;      
    }
    pointed=true;
    is_Computed.set(ConeProperty::IsPointed);    
    Full_Cone Polytope(Generators);
    Polytope.pointed=true;
    Polytope.is_Computed.set(ConeProperty::IsPointed);    
    Polytope.keep_order=true;
    Polytope.Grading=Truncation;
    Polytope.is_Computed.set(ConeProperty::Grading);
    if(isComputed(ConeProperty::SupportHyperplanes)){
        Polytope.Support_Hyperplanes=Support_Hyperplanes;
        Polytope.nrSupport_Hyperplanes=nrSupport_Hyperplanes;
        Polytope.is_Computed.set(ConeProperty::SupportHyperplanes);     
    }
    if(isComputed(ConeProperty::ExtremeRays)){
        Polytope.Extreme_Rays_Ind=Extreme_Rays_Ind;
        Polytope.is_Computed.set(ConeProperty::ExtremeRays);        
    }
    Polytope.do_deg1_elements=true;
    Polytope.verbose=verbose;
    Polytope.compute();

    if(Polytope.isComputed(ConeProperty::SupportHyperplanes) && 
        !isComputed(ConeProperty::SupportHyperplanes)){
        Support_Hyperplanes=Polytope.Support_Hyperplanes;
        nrSupport_Hyperplanes=Polytope.nrSupport_Hyperplanes;
        is_Computed.set(ConeProperty::SupportHyperplanes);  
    }
    if(Polytope.isComputed(ConeProperty::ExtremeRays) &&
        !isComputed(ConeProperty::ExtremeRays)){
        Extreme_Rays_Ind=Polytope.Extreme_Rays_Ind;
        is_Computed.set(ConeProperty::ExtremeRays);     
    }
    if(Polytope.isComputed(ConeProperty::Deg1Elements)){
        module_rank=Polytope.Deg1_Elements.size();
        if(do_Hilbert_basis){
            Hilbert_Basis=Polytope.Deg1_Elements;
            is_Computed.set(ConeProperty::HilbertBasis);
        }
        is_Computed.set(ConeProperty::ModuleRank);
        if(isComputed(ConeProperty::Grading)){
            multiplicity=1; // of the recession cone;
            is_Computed.set(ConeProperty::Multiplicity);
            if(do_h_vector){
                vector<num_t> hv(1);
                typename list<vector<Integer> >::const_iterator hb=Polytope.Deg1_Elements.begin();
                for(;hb!=Polytope.Deg1_Elements.end();++hb){
                    size_t deg = convertTo<long>(v_scalar_product(Grading,*hb));
                    if(deg+1>hv.size())
                        hv.resize(deg+1);
                    hv[deg]++;                        
                }    
                Hilbert_Series.add(hv,vector<denom_t>());
                Hilbert_Series.setShift(convertTo<long>(shift));
                Hilbert_Series.adjustShift();
                Hilbert_Series.simplify();
                is_Computed.set(ConeProperty::HilbertSeries);
            }
        }  
    }
}

template<>
void Full_Cone<renf_elem_class>::convert_polyhedron_to_polytope() {
    assert(false);    
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_deg1_elements_via_approx_global() {
    
    compute_elements_via_approx(Deg1_Elements);
    
    /*// now already done in simplex.cpp and directly for generators
    for(auto e=Deg1_Elements.begin(); e!=Deg1_Elements.end();)
        if(!contains(*e))
            e=Deg1_Elements.erase(e);
        else
            ++e; */
    if(verbose)
            verboseOutput() << Deg1_Elements.size() << " deg 1 elements found" << endl; 
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_elements_via_approx(list<vector<Integer> >& elements_from_approx) {

    if (!isComputed(ConeProperty::Grading)){
        support_hyperplanes(); // the only thing we can do now
        return;
    }
    assert(elements_from_approx.empty());
    vector<list<vector<Integer> > > approx_points = latt_approx();
    vector<vector<key_t> > approx_points_indices;
    key_t current_key =0;
    //cout << "Approximation points: " << endl;
    //for (size_t j=0;j<dim;++j){
        ////cout << "Original generator " << j << ": " << Generators[j] << endl;
        ////cout << approx_points[j] << endl;
    //}
    //cout << "Nr of ER: " << getExtremeRays().size() << endl;
    //Matrix<Integer> all_approx_points(Generators);
    Matrix<Integer> all_approx_points(0,dim);
    for (size_t i=0;i<nr_gen;i++){
        vector<key_t> indices(approx_points[i].size());
        if (!approx_points[i].empty()){

             all_approx_points.append(approx_points[i]);
             for (size_t j=0;j<approx_points[i].size();++j){
                 indices[j]=current_key;
                 current_key++;
             }
             
         }
         approx_points_indices.push_back(indices);
    }
    if(verbose){
        verboseOutput() << "Nr original generators: " << nr_gen << endl;
        verboseOutput() << "Nr approximation points: " << all_approx_points.nr_of_rows() << endl;
    }
    Full_Cone C_approx(all_approx_points);
    C_approx.OriginalGenerators = Generators;
    C_approx.approx_points_keys = approx_points_indices;
    C_approx.verbose = verbose;
    //C_temp.build_cone_approx(*this,approx_points_indices);
    

    //Full_Cone C_approx(C_temp.getGenerators()); 
    //Full_Cone C_approx(all_approx_points); // latt_approx computes a matrix of generators
   
    C_approx.is_approximation=true;
    // C_approx.Generators.pretty_print(cout);
    C_approx.do_deg1_elements=do_deg1_elements;  // for supercone C_approx that is generated in degree 1
    C_approx.do_Hilbert_basis=do_Hilbert_basis;
    C_approx.do_all_hyperplanes=false;    // not needed
    C_approx.Subcone_Support_Hyperplanes=Support_Hyperplanes; // *this is a subcone of C_approx, used to discard elements
    C_approx.Support_Hyperplanes=Support_Hyperplanes; // UNFORTUNATELY NEEDED IN REDUCTION FOR SUBDIVIUSION BY APPROX
    C_approx.is_Computed.set(ConeProperty::SupportHyperplanes);
    C_approx.nrSupport_Hyperplanes = nrSupport_Hyperplanes;
    C_approx.Grading = Grading;
    C_approx.is_Computed.set(ConeProperty::Grading);
    C_approx.Truncation=Truncation;
    C_approx.TruncLevel=TruncLevel;

    //if(verbose)
        //verboseOutput() << "Computing elements in approximating cone with "
                        //<< C_approx.Generators.nr_of_rows() << " generators." << endl;
    if(verbose){
        verboseOutput() << "Computing elements in approximating cone." << endl;
    }

    bool verbose_tmp = verbose;
    verbose =false;
    C_approx.compute();
    verbose = verbose_tmp;
    
    //vector<key_t> used_gens;
    //for (size_t j=0;j<nr_gen;++j){
        //if (C_approx.in_triang[j]) used_gens.push_back(j);
    //}
    //C_approx.Generators = C_approx.Generators.submatrix(used_gens);
    //if (verbose){
        //verboseOutput() << "Used "<< C_approx.Generators.nr_of_rows() << " / " << C_approx.nr_gen << " generators." << endl;
    //}
    //C_approx.nr_gen=C_approx.Generators.nr_of_rows();
    
    
    // TODO: with the current implementation, this is always the case!
    if(!C_approx.contains(*this) || Grading!=C_approx.Grading){
        throw FatalException("Wrong approximating cone.");
    }

    if(verbose)
        verboseOutput() << "Sum of dets of simplicial cones evaluated in approximation = " << C_approx.detSum << endl;   
    
    if(verbose)
        verboseOutput() << "Returning to original cone" << endl;
    if(do_deg1_elements)
        elements_from_approx.splice(elements_from_approx.begin(),C_approx.Deg1_Elements);
    if(do_Hilbert_basis)
        elements_from_approx.splice(elements_from_approx.begin(),C_approx.Hilbert_Basis);
}


// -s
template<typename Integer>
void Full_Cone<Integer>::support_hyperplanes() {

    if(!isComputed(ConeProperty::SupportHyperplanes)){
        sort_gens_by_degree(false); // we do not want to triangulate here
        build_top_cone();           
    }
    extreme_rays_and_deg1_check();
    if(inhomogeneous){
        find_level0_dim();
        if(do_module_rank) 
            find_module_rank();
    }
}

//---------------------------------------------------------------------------
// Checks and auxiliary algorithms
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::extreme_rays_and_deg1_check() {
    check_pointed();
    if(!pointed){
        throw NonpointedException();
    }
    //cout << "Generators" << endl;
    //Generators.pretty_print(cout);
    //cout << "SupportHyperplanes" << endl;
    //Support_Hyperplanes.pretty_print(cout);
    compute_extreme_rays();
    deg1_check();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_given_grading(){
    
    if(Grading.size()==0)
        return;

    bool positively_graded=true;
    
    if(!isComputed(ConeProperty::Grading)){
        size_t neg_index=0;
        Integer neg_value;
        bool nonnegative=true;
        vector<Integer> degrees = Generators.MxV(Grading);
        for (size_t i=0; i<degrees.size(); ++i) {
            if (degrees[i]<=0 && (!inhomogeneous || gen_levels[i]==0)) { 
                // in the inhomogeneous case: test only generators of tail cone
                positively_graded=false;;
                if(degrees[i]<0){
                    nonnegative=false;
                    neg_index=i;
                    neg_value=degrees[i];
                }
            }
        }

        if(!nonnegative){
            throw BadInputException("Grading gives negative value "
                    + toString(neg_value) + " for generator "
                    + toString(neg_index+1) + "!");
        }
    }
    
    if(positively_graded){
        is_Computed.set(ConeProperty::Grading);    
        if(inhomogeneous)
            find_grading_inhom();
        set_degrees();
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_grading(){
    
    if(inhomogeneous) // in the inhomogeneous case we do not allow implicit grading
        return;

    deg1_check(); // trying to find grading under which all generators have the same degree
    if (!isComputed(ConeProperty::Grading) && (do_multiplicity || do_deg1_elements || do_h_vector)) {
        if (!isComputed(ConeProperty::ExtremeRays)) {
            if (verbose) {
                verboseOutput() << "Cannot find grading s.t. all generators have the degree 1! Computing Extreme rays first:" << endl;
            }
            get_supphyps_from_copy(true);
            extreme_rays_and_deg1_check();
            if(!pointed){
                throw NonpointedException();
            };

            // We keep the SupportHyperplanes, so we do not need to recompute them
            // for the last generator, and use them to make a global reduction earlier
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_level0_dim(){
    
    if(isComputed(ConeProperty::RecessionRank))
        return;

    if(!isComputed(ConeProperty::Generators)){
        throw FatalException("Missing Generators.");
    }
    
    Matrix<Integer> Help(nr_gen,dim);
    for(size_t i=0; i<nr_gen;++i)
        if(gen_levels[i]==0)
            Help[i]=Generators[i];
        
    ProjToLevel0Quot=Help.kernel(false); // Necessary for the module rank
                                    // For level0_dim the rank of Help would be enough
    
    level0_dim=dim-ProjToLevel0Quot.nr_of_rows();
    // level0_dim=Help.rank();
    
    is_Computed.set(ConeProperty::RecessionRank);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_level0_dim_from_HB(){
// we use the Hilbert basis if we don't have the extreme reys.
// This is possible if the HB was computed by the dual algorithm.
// Would be enough if we would take the extreme reys of the recession cone,
// but they have not been extracted from the HB
    
    if(isComputed(ConeProperty::RecessionRank))
        return;

    assert(isComputed(ConeProperty::HilbertBasis));
    
    Matrix<Integer> Help(0,dim);
    for(const auto& H : Hilbert_Basis)
        if(v_scalar_product(H,Truncation)==0)
            Help.append(H);
        
    ProjToLevel0Quot=Help.kernel(); // Necessary for the module rank
                                    // For level0_dim the rank of Help would be enough
    
    level0_dim=dim-ProjToLevel0Quot.nr_of_rows();
    
    is_Computed.set(ConeProperty::RecessionRank);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_module_rank(){
    
    if(isComputed(ConeProperty::ModuleRank))
        return;   

    if(level0_dim==dim){
        module_rank=0;
        is_Computed.set(ConeProperty::ModuleRank);
        return;
    }     
    if(isComputed(ConeProperty::HilbertBasis)){
        find_module_rank_from_HB();
        return;
    }

    // size_t HBrank = module_rank;

    if(do_module_rank)
        find_module_rank_from_proj();
    
    /* if(isComputed(ConeProperty::HilbertBasis))
        assert(HBrank==module_rank);
    */
    
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_module_rank_from_proj(){
    
    if(verbose){
        verboseOutput() << "Computing projection to quotient mod level 0" << endl;
    } 

    Matrix<Integer> ProjGen(nr_gen,dim-level0_dim);
    for(size_t i=0;i<nr_gen;++i){
        ProjGen[i]=ProjToLevel0Quot.MxV(Generators[i]);
    }
    
    vector<Integer> GradingProj=ProjToLevel0Quot.transpose().solve_ZZ(Truncation);
    
    Full_Cone<Integer> Cproj(ProjGen);
    Cproj.verbose=false;
    Cproj.Grading=GradingProj;
    Cproj.is_Computed.set(ConeProperty::Grading);
    Cproj.do_deg1_elements=true;
    Cproj.compute();
    
    module_rank=Cproj.Deg1_Elements.size();
    is_Computed.set(ConeProperty::ModuleRank);
    return;
}
        
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_module_rank_from_HB(){
    
    if(level0_dim==0){
        module_rank=Hilbert_Basis.size();
        is_Computed.set(ConeProperty::ModuleRank);
        return;        
    }

    
    set<vector<Integer> > Quotient;
    vector<Integer> v;
    
    // cout << "=======================" << endl;
    // ProjToLevel0Quot.print(cout);
    // cout << "=======================" << endl;
    
    for(const auto& h : Hilbert_Basis){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        v=ProjToLevel0Quot.MxV(h);
        bool zero=true;
        for(size_t j=0;j<v.size();++j)
            if(v[j]!=0){
                zero=false;
                break;
            }
        if(!zero)
            Quotient.insert(v);
    }
    
    module_rank=Quotient.size();
    is_Computed.set(ConeProperty::ModuleRank);

}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_grading_inhom(){

    if(Grading.size()==0 || Truncation.size()==0){
         throw FatalException("Cannot find grading in the inhomogeneous case!");
    }
    
    if(shift!=0)  // to avoid double computation
        return;

    bool first=true;
    Integer level,degree,quot=0,min_quot=0;
    for(size_t i=0;i<nr_gen;++i){
        level=v_scalar_product(Truncation,Generators[i]);
        if(level==0)
            continue;
        degree=v_scalar_product(Grading,Generators[i]);
        quot=degree/level;
        // cout << Generators[i];
        // cout << "*** " << degree << " " << level << " " << quot << endl;
        if(level*quot>=degree)
            quot--;
        if(first){
            min_quot=quot;
            first=false;
        }
        if(quot<min_quot)
            min_quot=quot;
        // cout << "+++ " << min_quot << endl;
    }
    shift = min_quot;
    for(size_t i=0;i<dim;++i) // under this grading all generators have positive degree
        Grading[i] = Grading[i] - shift * Truncation[i];
        
    // shift--;  // NO LONGER correction for the Hilbert series computation to have it start in degree 0
}

#ifdef ENFNORMALIZ
template<>
void Full_Cone<renf_elem_class>::find_grading_inhom(){
    assert(false);
}
#endif

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::set_levels() {
    if(inhomogeneous && Truncation.size()!=dim){
        throw FatalException("Truncation not defined in inhomogeneous case.");
    }    

    if(gen_levels.size()!=nr_gen) // now we compute the levels
    {
        gen_levels.resize(nr_gen);
        vector<Integer> gen_levels_Integer=Generators.MxV(Truncation);
        for (size_t i=0; i<nr_gen; i++) {
            if (gen_levels_Integer[i] < 0) {
                throw FatalException("Truncation gives non-positive value "
                        + toString(gen_levels_Integer[i]) + " for generator "
                        + toString(i+1) + ".");
            }
            convert(gen_levels[i], gen_levels_Integer[i]);
            // cout << "Gen " << Generators[i];
            // cout << "level " << gen_levels[i] << endl << "----------------------" << endl;
        }
    }    
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::sort_gens_by_degree(bool triangulate) {
    // if(deg1_extreme_rays)  // gen_degrees.size()==0 || 
    // return;
    
    if(keep_order)
        return;

    // we first order the generaors by "support hyperplanes" for computations using automorphisms
    // in order to have an intrinsic useful sorting
    if(isComputed(ConeProperty::SupportHyperplanes) && descent_level>0){
        Matrix<Integer> TranspType(Support_Hyperplanes.nr_of_rows(),Generators.nr_of_rows());
 
        #pragma omp parallel for
        for(size_t i=0;i<Support_Hyperplanes.nr_of_rows();++i)
            for(size_t j=0;j<Generators.nr_of_rows();++j)
                TranspType[i][j]=v_scalar_product(Support_Hyperplanes[i],Generators[j]);

        Matrix<Integer> OrderSupps=TranspType.sort_by_nr_of_zeroes();
        Matrix<Integer> Type=TranspType.transpose();
        vector<key_t> new_perm=Type.perm_by_lex();
        Generators.order_rows_by_perm(new_perm);
        compose_perm_gens(new_perm);
        if(verbose)        
            verboseOutput() << "Generators sorted lexicographically by scalar products with support hyperplanes" << endl;
    }
    
    Matrix<Integer> Weights(0,dim);
    vector<bool> absolute;
    if(triangulate){
        if(isComputed(ConeProperty::Grading)){
            Weights.append(Grading);
            absolute.push_back(false);
        }
        else{
            Weights.append(vector<Integer>(dim,1));
            absolute.push_back(true);
        }
    }
    
    vector<key_t> perm=Generators.perm_by_weights(Weights,absolute);

    Generators.order_rows_by_perm(perm);
    order_by_perm_bool(Extreme_Rays_Ind,perm);

    if(isComputed(ConeProperty::Grading)){
        order_by_perm(gen_degrees,perm);
        if(do_h_vector || (!using_GMP<Integer>() && !using_renf<Integer>()))
            order_by_perm(gen_degrees_long,perm);
    }

    if(inhomogeneous && gen_levels.size()==nr_gen)
        order_by_perm(gen_levels,perm);

    if(triangulate){
        Integer roughness;
        if(isComputed(ConeProperty::Grading)){
                roughness=gen_degrees[nr_gen-1]/gen_degrees[0];
        }
        else{
            Integer max_norm=0, min_norm=0;
            for(size_t i=0;i<dim;++i){
            max_norm+=Iabs(Generators[nr_gen-1][i]);
            min_norm+=Iabs(Generators[0][i]); 
            }
            roughness=max_norm/min_norm;
        }
        if(verbose){
            verboseOutput() << "Roughness " << roughness <<  endl;
        }
        
        if(roughness >= 10 && !suppress_bottom_dec){
            do_bottom_dec=true;
            if(verbose){
                    verboseOutput() << "Bottom decomposition activated" << endl;
            }
        }
    }
    
    if(exploit_automs_vectors && descent_level==0 && isComputed(ConeProperty::Grading)){
        vector<key_t> inverse_order(nr_gen); 
        for(size_t i=0;i<nr_gen;++i)
            inverse_order[i]=nr_gen-1-i;
        vector<key_t> largest_simplex=Generators.max_rank_submatrix_lex(inverse_order);
        HB_bound=-1;
        for(size_t i=0;i<dim;++i)
            HB_bound+=convertTo<Integer>(gen_degrees[largest_simplex[i]]);
    }
    
    if (verbose) {
        if(triangulate){ 
            if(isComputed(ConeProperty::Grading)){
                verboseOutput() <<"Generators sorted by degree and lexicographically" << endl;
                verboseOutput() << "Generators per degree:" << endl;
                verboseOutput() << count_in_map<Integer,long>(gen_degrees);
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

template<typename Integer>
void Full_Cone<Integer>::compose_perm_gens(const vector<key_t>& perm) {
    order_by_perm(PermGens,perm);
}

//---------------------------------------------------------------------------

// an alternative to compute() for the basic tasks that need no triangulation
template<typename Integer>
void Full_Cone<Integer>::dualize_cone(bool print_message){
    
    omp_start_level=omp_get_level();
    
    if(dim==0){
        set_zero_cone();
        return;
    }

    // DO NOT CALL do_vars_check!!

    bool save_tri      = do_triangulation;
    bool save_part_tri = do_partial_triangulation;
    /* do_triangulation         = false;
    do_partial_triangulation = false; */
    
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

template<typename Integer>
vector<key_t> Full_Cone<Integer>::find_start_simplex() const {
        return Generators.max_rank_submatrix_lex();
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
    for(auto ll=S.begin();ll!=S.end()&&i<k;++ll){
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

void Full_Cone<Integer>::minimize_support_hyperplanes(){
    if(Support_Hyperplanes.nr_of_rows() == 0){
        return;
    }
    if(isComputed(ConeProperty::SupportHyperplanes)){
        nrSupport_Hyperplanes=Support_Hyperplanes.nr_of_rows();
        return;
    }
    if (verbose) {
        verboseOutput() << "Minimize the given set of support hyperplanes by "
                        << "computing the extreme rays of the dual cone" << endl;
    }
    Full_Cone<Integer> Dual(Support_Hyperplanes);
    Dual.verbose=false; // verbose;
    Dual.Support_Hyperplanes = Generators;
    Dual.is_Computed.set(ConeProperty::SupportHyperplanes);
    Dual.compute_extreme_rays();
    Support_Hyperplanes = Dual.Generators.submatrix(Dual.Extreme_Rays_Ind); //only essential hyperplanes
    is_Computed.set(ConeProperty::SupportHyperplanes);
    nrSupport_Hyperplanes=Support_Hyperplanes.nr_of_rows();
    do_all_hyperplanes=false;
}
    

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_extreme_rays(bool use_facets){

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

template<typename Integer>
void Full_Cone<Integer>::compute_extreme_rays_rank(bool use_facets){

    if (verbose) verboseOutput() << "Select extreme rays via rank ... " << flush;

    size_t i;
    vector<key_t> gen_in_hyperplanes;
    gen_in_hyperplanes.reserve(Support_Hyperplanes.nr_of_rows());
    Matrix<Integer> M(Support_Hyperplanes.nr_of_rows(),dim);

    deque<bool> Ext(nr_gen,false);
    #pragma omp parallel for firstprivate(gen_in_hyperplanes,M)
    for(i=0;i<nr_gen;++i){
//        if (isComputed(ConeProperty::Triangulation) && !in_triang[i])
//            continue;
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        gen_in_hyperplanes.clear();
        if(use_facets){
            typename list<FACETDATA<Integer>>::const_iterator IHV=Facets.begin();            
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

template<typename Integer>
void Full_Cone<Integer>::compute_extreme_rays_compare(bool use_facets){

    if (verbose) verboseOutput() << "Select extreme rays via comparison ... " << flush;

    size_t i,j,k;
    // Matrix<Integer> SH=getSupportHyperplanes().transpose();
    // Matrix<Integer> Val=Generators.multiplication(SH);
    size_t nc=Support_Hyperplanes.nr_of_rows();
    
    vector<vector<bool> > Val(nr_gen);
    for (i=0;i<nr_gen;++i)
       Val[i].resize(nc);
        
    // In this routine Val[i][j]==1, i.e. true, indicates that
    // the i-th generator is contained in the j-th support hyperplane
    
    vector<key_t> Zero(nc);
    vector<key_t> nr_ones(nr_gen);

    for (i = 0; i <nr_gen; i++) {
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        k=0;
        Extreme_Rays_Ind[i]=true;
        if(use_facets){
            typename list<FACETDATA<Integer>>::const_iterator IHV=Facets.begin();            
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

template<typename Integer>
void Full_Cone<Integer>::compute_class_group() { // from the support hyperplanes
    if(!do_class_group || !isComputed(ConeProperty::SupportHyperplanes) || isComputed(ConeProperty::ClassGroup)
        || descent_level!=0)
        return;
    Matrix<Integer> Trans=Support_Hyperplanes; // .transpose();
    size_t rk;
    Trans.SmithNormalForm(rk);
    ClassGroup.push_back(Support_Hyperplanes.nr_of_rows()-rk);
    for(size_t i=0;i<rk;++i)
        if(Trans[i][i]!=1)
            ClassGroup.push_back(Trans[i][i]);
    is_Computed.set(ConeProperty::ClassGroup);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::select_deg1_elements() { // from the Hilbert basis

    if(inhomogeneous || descent_level>0)
        return;
    for(const auto & h : Hilbert_Basis){
        if(v_scalar_product(Grading,h)==1)
            Deg1_Elements.push_back(h);
    }
    is_Computed.set(ConeProperty::Deg1Elements,true);
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::subcone_contains(const vector<Integer>& v) {
    for (size_t i=0; i<Subcone_Support_Hyperplanes.nr_of_rows(); ++i)
        if (v_scalar_product(Subcone_Support_Hyperplanes[i],v) < 0)
            return false;
    for (size_t i=0; i<Subcone_Equations.nr_of_rows(); ++i)
        if (v_scalar_product(Subcone_Equations[i],v) != 0)
            return false;
    if(is_global_approximation)
        if(v_scalar_product(Subcone_Grading,v)!=1)
            return false;

    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::contains(const vector<Integer>& v) {
    for (size_t i=0; i<Support_Hyperplanes.nr_of_rows(); ++i)
        if (v_scalar_product(Support_Hyperplanes[i],v) < 0)
            return false;
    return true;
}
//---------------------------------------------------------------------------

template<typename Integer>
bool Full_Cone<Integer>::contains(const Full_Cone& C) {
    for(size_t i=0;i<C.nr_gen;++i)
        if(!contains(C.Generators[i])){
            cerr << "Missing generator " << C.Generators[i] << endl;
            return(false);
    }
    return(true);
}
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::select_deg1_elements(const Full_Cone& C) {  // from vectors computed in 
                                                              // the auxiliary cone C
    assert(isComputed(ConeProperty::SupportHyperplanes));
    assert(C.isComputed(ConeProperty::Deg1Elements));
    for(const auto& h : C.Deg1_Elements){
        if(contains(h))
            Deg1_Elements.push_back(h);
    }
    is_Computed.set(ConeProperty::Deg1Elements,true);
}

//---------------------------------------------------------------------------


// so far only for experiments
/*
template<typename Integer>
void Full_Cone<Integer>::select_Hilbert_Basis(const Full_Cone& C) {  // from vectors computed in 
                                                              // the auxiliary cone C
    assert(isComputed(ConeProperty::SupportHyperplanes));
    assert(C.isComputed(ConeProperty::Deg1Elements));
    typename list<vector<Integer> >::const_iterator h = C.Hilbert_Basis.begin();
    for(;h!=C.Hilbert_Basis.end();++h){
        if(contains(*h))
            // Deg1_Elements.push_back(*h);
            cout << *h;
    }
    exit(0);
    is_Computed.set(ConeProperty::Deg1Elements,true);
}
*/

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_pointed() {
    if (isComputed(ConeProperty::IsPointed))
        return;
    assert(isComputed(ConeProperty::SupportHyperplanes));
    if (isComputed(ConeProperty::Grading)){
        pointed=true;
        if (verbose) verboseOutput() << "Pointed since graded" << endl;
        is_Computed.set(ConeProperty::IsPointed);
        return;
    }
    if (verbose) verboseOutput() << "Checking pointedness ... " << flush;
    if(Support_Hyperplanes.nr_of_rows() <= dim*dim/2){
        pointed = (Support_Hyperplanes.rank() == dim);
    }
    else
        pointed = (Support_Hyperplanes.max_rank_submatrix_lex().size() == dim);
    is_Computed.set(ConeProperty::IsPointed);
    if(pointed && Grading.size()>0){
        throw BadInputException("Grading not positive on pointed cone.");
    }
    if (verbose) verboseOutput() << "done." << endl;
}


//---------------------------------------------------------------------------
template<typename Integer>
void Full_Cone<Integer>::disable_grading_dep_comp() {

    if (do_multiplicity || do_deg1_elements || do_h_vector) {
        if (do_default_mode) {
            do_deg1_elements = false;
            do_h_vector = false;            
            if(!explicit_full_triang){
                do_triangulation=false;
                if(do_Hilbert_basis)
                    do_partial_triangulation=true;
            }

        } else {
            throw NotComputableException("No grading specified and cannot find one. Cannot compute some requested properties!");
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::deg1_check() {

    if(inhomogeneous)  // deg 1 check disabled since it makes no sense in this case
        return;
        
    if (!isComputed(ConeProperty::Grading) && Grading.size()==0          // we still need it and
     && !isComputed(ConeProperty::IsDeg1ExtremeRays)) { // we have not tried it
        if (isComputed(ConeProperty::ExtremeRays)) {
            Matrix<Integer> Extreme=Generators.submatrix(Extreme_Rays_Ind);
            if (has_generator_with_common_divisor) 
                Extreme.make_prime();
            try{
                Grading = Extreme.find_linear_form();
            } catch(const ArithmeticException& e) { // if the exception has been thrown, the grading has 
                Grading.resize(0);   // we consider the grafing as non existing -- though this may not be true
                if(verbose)
                    verboseOutput() << "Giving up the check for a grading" << endl;
            }

            if (Grading.size() == dim && v_scalar_product(Grading,Extreme[0])==1) {
                is_Computed.set(ConeProperty::Grading);
            } else {
                deg1_extreme_rays = false;
                Grading.clear();
                is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
            }
        } else // extreme rays not known
        if (!deg1_generated_computed) {
            Matrix<Integer> GenCopy = Generators;
            if (has_generator_with_common_divisor)
                GenCopy.make_prime();
            try{
            Grading = GenCopy.find_linear_form();
            } catch(const ArithmeticException& e) { // if the exception has been thrown,
                Grading.resize(0);   // we consider the grafing as non existing-- though this may not be true
                if(verbose)
                    verboseOutput() << "Giving up the check for a grading" << endl;
            }
            if (Grading.size() == dim && v_scalar_product(Grading,GenCopy[0])==1) {
                is_Computed.set(ConeProperty::Grading);
            } else {
                deg1_generated = false;
                deg1_generated_computed = true;
                Grading.clear();
            }
        }
    }

    //now we hopefully have a grading

    if (!isComputed(ConeProperty::Grading)) {
        if (isComputed(ConeProperty::ExtremeRays)) {
            // there is no hope to find a grading later
            deg1_generated = false;
            deg1_generated_computed = true;
            deg1_extreme_rays = false;
            is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
            disable_grading_dep_comp();
        }
        return; // we are done
    }
    
    set_degrees();

    vector<Integer> divided_gen_degrees = gen_degrees;
    if (has_generator_with_common_divisor) {
        Matrix<Integer> GenCopy = Generators;
        GenCopy.make_prime();
        convert(divided_gen_degrees, GenCopy.MxV(Grading));
    }

    if (!deg1_generated_computed) {
        deg1_generated = true;
        for (size_t i = 0; i < nr_gen; i++) {
            if (divided_gen_degrees[i] != 1) {
                deg1_generated = false;
                break;
            }
        }
        deg1_generated_computed = true;
        if (deg1_generated) {
            deg1_extreme_rays = true;
            is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
        }
    }
    if (!isComputed(ConeProperty::IsDeg1ExtremeRays)
      && isComputed(ConeProperty::ExtremeRays)) {
        deg1_extreme_rays = true;
        for (size_t i = 0; i < nr_gen; i++) {
            if (Extreme_Rays_Ind[i] && divided_gen_degrees[i] != 1) {
                deg1_extreme_rays = false;
                break;
            }
        }
        is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_deg1_hilbert_basis() {
    if (isComputed(ConeProperty::IsDeg1HilbertBasis) || inhomogeneous || descent_level>0)
        return;

    if ( !isComputed(ConeProperty::Grading) || !isComputed(ConeProperty::HilbertBasis)) {
        if (verbose) {
            errorOutput() << "WARNING: unsatisfied preconditions in check_deg1_hilbert_basis()!" <<endl;
        }
        return;
    }

    if (isComputed(ConeProperty::Deg1Elements)) {
        deg1_hilbert_basis = (Deg1_Elements.size() == Hilbert_Basis.size());
    } else {
        deg1_hilbert_basis = true;
        for (const auto& h : Hilbert_Basis) {
            if (v_scalar_product(h,Grading)!=1) {
                deg1_hilbert_basis = false;
                break;
            }
        }
    }
    is_Computed.set(ConeProperty::IsDeg1HilbertBasis);
}

//---------------------------------------------------------------------------

// Computes the generators of a supercone approximating "this" by a cone over a lattice polytope
// for every vertex of the simplex, we get a matrix with the integer points of the respective Weyl chamber
template<typename Integer>
vector<list<vector<Integer> > > Full_Cone<Integer>::latt_approx() {
    
    assert(isComputed(ConeProperty::Grading));
    assert(isComputed(ConeProperty::ExtremeRays));
    Matrix<Integer> G(1,dim);
    G[0]=Grading;
    Matrix<Integer> G_copy=G;
    
    // Lineare_Transformation<Integer> NewBasis(G); // gives a new basis in which the grading is a coordinate
    size_t dummy;
    Matrix<Integer> U=G_copy.SmithNormalForm(dummy);   // the basis elements are the columns of U

    Integer dummy_denom;                             
    // vector<Integer> dummy_diag(dim); 
    Matrix<Integer> T=U.invert(dummy_denom);       // T is the coordinate transformation
                                                            // to the new basis: v --> Tv (in this case)
                                                    // for which the grading is the FIRST coordinate

    assert(dummy_denom==1);  // for safety 

    // It can happen that -Grading has become the first row of T, but we want Grading. If necessary we replace the
    // first row by its negative, and correspondingly the first column of U by its negative

    if(T[0]!=Grading){
        for(size_t i=0;i<dim;++i){
            U[i][0]*=-1;
            T[0][i]*=-1;
        }
    }
    assert(T[0] == Grading);
    
    Matrix<Integer> M;
    vector<list<vector<Integer> > > approx_points;
    size_t nr_approx=0;
    for(size_t i=0;i<nr_gen;++i){
        list<vector<Integer> > approx;
        if(Extreme_Rays_Ind[i]){
            //cout << "point before transformation: " << Generators[i];
            approx_simplex(T.MxV(Generators[i]),approx,1);
            // TODO: NECESSARY?
            //approx.unique()
            //cout << "Approximation points for generator " << i << ": " << Generators[i] << endl;
            nr_approx = 0;
            for(auto& jt : approx){  // reverse transformation
                jt=U.MxV(jt);
                v_make_prime(jt);
                ++nr_approx;
                //cout << jt << endl;
            }
            //~ M=Matrix<Integer>(approx);
            //~ 
            //~ // remove duplicates
            //~ for (size_t j=0;j<approx_points.size();j++){
                //~ M.remove_duplicate(approx_points[j]);
            //~ }
            // +++ TODO: REMOVE DUPLICATES MORE EFFICIENT +++
            if (nr_approx>1){
                for (size_t j=0;j<approx_points.size();j++){
                    for (const auto& jt : approx_points[j]){
                        approx.remove(jt);
                    }
                }
            }
        }
        approx_points.push_back(approx);
    }
    
    
    //cout << "-------" << endl;
     //M.print(cout);
     //cout << "-------" << endl;
    
    return(approx_points);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::prepare_inclusion_exclusion() {

    if (ExcludedFaces.nr_of_rows() == 0)
        return;

    do_excluded_faces = do_h_vector || do_Stanley_dec;

    if (verbose && !do_excluded_faces) {
        errorOutput() << endl << "WARNING: excluded faces, but no h-vector computation or Stanley decomposition"
                      << endl << "Therefore excluded faces will be ignored" << endl;
    }

    if (isComputed(ConeProperty::ExcludedFaces) &&
            (isComputed(ConeProperty::InclusionExclusionData) || !do_excluded_faces) ) {
        return;
    }

    // indicates which generators lie in the excluded faces
    vector<boost::dynamic_bitset<> > GensInExcl(ExcludedFaces.nr_of_rows());

    for(size_t j=0;j<ExcludedFaces.nr_of_rows();++j){ // now we produce these indicators
        bool first_neq_0=true;           // and check whether the linear forms in ExcludedFaces
        bool non_zero=false;             // have the cone on one side
        GensInExcl[j].resize(nr_gen,false);
        for(size_t i=0; i< nr_gen;++i){
            Integer test=v_scalar_product(ExcludedFaces[j],Generators[i]);
            if(test==0){
                GensInExcl[j].set(i);
                continue;
            }
            non_zero=true;
            if(first_neq_0){
                first_neq_0=false;
                if(test<0){
                    for(size_t k=0;k<dim;++k)     // replace linear form by its negative
                        ExcludedFaces[j][k]*=-1;  // to get cone in positive halfspace 
                    test*=-1;                     // (only for error check)
                }    
            }
            if(test<0){
                throw FatalException("Excluded hyperplane does not define a face.");
            }
                
        }
        if(!non_zero){  // not impossible if the hyperplane contains the vector space spanned by the cone
            throw FatalException("Excluded face contains the full cone.");
        }       
    }
    
    vector<bool> essential(ExcludedFaces.nr_of_rows(),true);
    bool remove_one=false;
    for(size_t i=0;i<essential.size();++i)
        for(size_t j=i+1;j<essential.size();++j){
            if(GensInExcl[j].is_subset_of(GensInExcl[i])){
                essential[j]=false;
                remove_one=true;
                continue;
            }
            if(GensInExcl[i].is_subset_of(GensInExcl[j])){
                essential[i]=false;
                remove_one=true;
            }
        }
    if(remove_one){
        Matrix<Integer> Help(0,dim);
        vector<boost::dynamic_bitset<> > HelpGensInExcl;
        for(size_t i=0;i<essential.size();++i)
            if(essential[i]){
                Help.append(ExcludedFaces[i]);
                HelpGensInExcl.push_back(GensInExcl[i]);
            }
        ExcludedFaces=Help;
        GensInExcl=HelpGensInExcl;
    }
    is_Computed.set(ConeProperty::ExcludedFaces);



    if (isComputed(ConeProperty::InclusionExclusionData) || !do_excluded_faces) {
        return;
    }

    vector< pair<boost::dynamic_bitset<> , long> > InExScheme;  // now we produce the formal 
    boost::dynamic_bitset<> all_gens(nr_gen);             // inclusion-exclusion scheme
    all_gens.set();                         // by forming all intersections of
                                           // excluded faces
    InExScheme.push_back(pair<boost::dynamic_bitset<> , long> (all_gens, 1));
    size_t old_size=1;
    
    for(size_t i=0;i<ExcludedFaces.nr_of_rows();++i){
        for(size_t j=0;j<old_size;++j)
            InExScheme.push_back(pair<boost::dynamic_bitset<> , long>
                   (InExScheme[j].first & GensInExcl[i], -InExScheme[j].second));
        old_size*=2;
    }
    
    InExScheme.erase(InExScheme.begin()); // remove full cone
    
    // map<boost::dynamic_bitset<>, long> InExCollect;
    InExCollect.clear();
    
    for(size_t i=0;i<old_size-1;++i){               // we compactify the list of faces
        auto F=InExCollect.find(InExScheme[i].first);    // obtained as intersections
        if(F!=InExCollect.end())                    // by listing each face only once
            F->second+=InExScheme[i].second;        // but with the right multiplicity
        else
            InExCollect.insert(InExScheme[i]);
    }
     
    for(auto F=InExCollect.begin();F!=InExCollect.end();){   // faces with multiplicity 0
        if(F->second==0)                                 // can be erased
            InExCollect.erase(F++);
        else{
            ++F;
        }    
    }
     
    if(verbose){
        verboseOutput() << endl;
        verboseOutput() << "InEx complete, " << InExCollect.size() << " faces involved" << endl;
    }
     
    is_Computed.set(ConeProperty::InclusionExclusionData);
}

//---------------------------------------------------------------------------

template<typename Integer>

bool Full_Cone<Integer>::check_extension_to_god_father(){

    assert(dim==God_Father->dim);
    vector<Integer> test(dim);
    for(size_t k=0;k<Automs.LinMaps.size();++k){
        for(size_t i=0;i<God_Father->nr_gen;++i){
            test=Automs.LinMaps[k].MxV(God_Father->Generators[i]);
            if(God_Father->Generator_Set.find(test)==God_Father->Generator_Set.end())
                return false;
        }        
    }
    return true;
}

//---------------------------------------------------------------------------

/* computes a degree function, s.t. every generator has value >0 */
template<typename Integer>
vector<Integer> Full_Cone<Integer>::compute_degree_function() const {
    size_t i;  
    vector<Integer> degree_function(dim,0);
    if (isComputed(ConeProperty::Grading)) { //use the grading if we have one
        for (i=0; i<dim; i++) {
            degree_function[i] = Grading[i];
        }
    } else { // add hyperplanes to get a degree function
        if(verbose) {
            verboseOutput()<<"computing degree function... "<<flush;
        }
        size_t h;
        for (h=0; h<Support_Hyperplanes.nr_of_rows(); ++h) {
            for (i=0; i<dim; i++) {
                degree_function[i] += Support_Hyperplanes.get_elem(h,i);
            }
        }
        if(using_renf<Integer>())
            v_standardize(degree_function,Norm);
        else
            v_make_prime(degree_function);
        if(verbose) {
            verboseOutput()<<"done."<<endl;
        }
    }
    return degree_function;
}

//---------------------------------------------------------------------------

/* adds generators, they have to lie inside the existing cone */
template<typename Integer>
void Full_Cone<Integer>::add_generators(const Matrix<Integer>& new_points) {
    is_simplicial = false;
    int nr_new_points = new_points.nr_of_rows();
    int nr_old_gen = nr_gen;
    Generators.append(new_points);
    nr_gen += nr_new_points;
    set_degrees();
    Top_Key.resize(nr_gen);
    Extreme_Rays_Ind.resize(nr_gen);
    for (size_t i=nr_old_gen; i<nr_gen; ++i) {
        Top_Key[i] = i;
        Extreme_Rays_Ind[i] = false;
    }
    // inhom cones
    if (inhomogeneous) {
        set_levels();
    }
    // excluded faces have to be reinitialized
    is_Computed.set(ConeProperty::ExcludedFaces, false);
    is_Computed.set(ConeProperty::InclusionExclusionData, false);
    prepare_inclusion_exclusion();

    if (do_Hilbert_basis) {
        // add new points to HilbertBasis
        for (size_t i = nr_old_gen; i < nr_gen; ++i) {
            if(!inhomogeneous || gen_levels[i]<=1) {
                NewCandidates.reduce_by_and_insert(Generators[i],*this,OldCandidates);
                NewCandidates.Candidates.back().original_generator=true;
            }
        }
        // OldCandidates.auto_reduce();
    }
}

//---------------------------------------------------------------------------
// Constructors
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::reset_tasks(){

    do_triangulation_size = false;
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
    do_subdivision_points=false;    
    do_extreme_rays=false;
    do_pointed=false;
    do_all_hyperplanes=true;
        
    do_bottom_dec=false;
    keep_order=false;

    nrSimplicialPyr=0;
    totalNrPyr=0;
    is_pyramid = false;


    exploit_automs_vectors=false;
    exploit_automs_mult=false;
    do_automorphisms=false;
    autom_codim_vectors=-1;
    autom_codim_mult=-1;
    
    use_existing_facets=false;

    triangulation_is_nested = false;
    triangulation_is_partial = false;
    
    hilbert_basis_rec_cone_known=false;
    time_measured=false;
    
    keep_convex_hull_data=false;
}


//---------------------------------------------------------------------------

template<typename Integer>
Full_Cone<Integer>::Full_Cone(){
    reset_tasks();
}

template<typename Integer>
Full_Cone<Integer>::Full_Cone(const Matrix<Integer>& M, bool do_make_prime){ // constructor of the top cone
    
    omp_start_level=omp_get_level();
        
    dim=M.nr_of_columns();
    if(dim>0)
        Generators=M;
    
    /*
    cout << "------------------" << endl;
    cout << "dim " << dim << endl;
    M.pretty_print(cout);
    cout << "------------------" << endl;
    M.transpose().pretty_print(cout);
    cout << "==================" << endl;
    */
    
    // assert(M.row_echelon()== dim); rank check now done later 
    
    /*index=1;                      // not used at present
    for(size_t i=0;i<dim;++i)
        index*=M[i][i];
    index=Iabs(index); */

    //make the generators coprime, remove 0 rows and duplicates
    has_generator_with_common_divisor = false;
    if (do_make_prime) {
        Generators.make_prime();
    } else {
        nr_gen = Generators.nr_of_rows();
        for (size_t i = 0; i < nr_gen; ++i) {
            if (v_gcd(Generators[i]) != 1) {
                has_generator_with_common_divisor = true;
                break;
            }
        }
    }
    Generators.remove_duplicate_and_zero_rows();
    
    nr_gen = Generators.nr_of_rows();

    if (nr_gen != static_cast<size_t>(static_cast<key_t>(nr_gen))) {
        throw FatalException("Too many generators to fit in range of key_t!");
    }
    
    multiplicity = 0;
#ifdef ENFNORMALIZ
    renf_multiplicity=0;
#endif
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
    // in_triang = vector<bool> (nr_gen,false); // now in build_cone
    deg1_triangulation = true;
/*
    if(dim==0){            //correction needed to include the 0 cone;
        multiplicity = 1;
#ifdef ENFNORMALIZ
    renf_multiplicity=1;
#endif
        Hilbert_Series.add(vector<num_t>(1,1),vector<denom_t>());
        is_Computed.set(ConeProperty::HilbertSeries);
        is_Computed.set(ConeProperty::Triangulation);
    }
*/
    pyr_level=-1;
    descent_level=0;
    Top_Cone=this;
    God_Father=this;
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=i;
    totalNrSimplices=0;
    TriangulationBufferSize=0;
    CandidatesSize=0;
    detSum = 0;
    shift = 0;
    
    FS.resize(omp_get_max_threads());
    
    Pyramids.resize(20);  // prepare storage for pyramids
    nrPyramids.resize(20,0);
    Pyramids_scrambled.resize(20,false);
      
    recursion_allowed=true;
    
    // nextGen=0;
    store_level=0;
    
    Comparisons.reserve(nr_gen);
    nrTotalComparisons=0;

    inhomogeneous=false;
    
    level0_dim=dim; // must always be defined

    start_from=0;
    old_nr_supp_hyps=0;
    
    verbose=false;
    OldCandidates.dual=false;
    OldCandidates.verbose=verbose;
    NewCandidates.dual=false;
    NewCandidates.verbose=verbose;
    
    RankTest = vector< Matrix<Integer> >(omp_get_max_threads(), Matrix<Integer>(0,dim));
    RankTest_float = vector< Matrix<nmz_float> >(omp_get_max_threads(), Matrix<nmz_float>(0,dim));
    
    do_bottom_dec=false;
    suppress_bottom_dec=false;
    keep_order=false;

    is_approximation=false;
    is_global_approximation=false;
    
    PermGens.resize(nr_gen);
    for(size_t i=0;i<nr_gen;++i)
        PermGens[i]=i;
    
    Mother=&(*this);

    don_t_add_hyperplanes=false;
    take_time_of_large_pyr=false;
    
    renf_degree=2; // default value to prevent desasters
}

//---------------------------------------------------------------------------
// converts a Cone_Dual_Mode into a Full_Cone
template<typename Integer>
Full_Cone<Integer>::Full_Cone(Cone_Dual_Mode<Integer> &C) {
    
    omp_start_level=omp_get_level();

    is_Computed = bitset<ConeProperty::EnumSize>();  //initialized to false

    dim = C.dim;
    Generators.swap(C.Generators);
    nr_gen = Generators.nr_of_rows();
    if (Generators.nr_of_rows() > 0) 
        is_Computed.set(ConeProperty::Generators);
    has_generator_with_common_divisor = false;
    Extreme_Rays_Ind.swap(C.ExtremeRaysInd);
    if (!Extreme_Rays_Ind.empty()) is_Computed.set(ConeProperty::ExtremeRays);

    multiplicity = 0;

#ifdef ENFNORMALIZ
    renf_multiplicity=0;
#endif
    in_triang = vector<bool>(nr_gen,false);
 
    Basis_Max_Subspace=C.BasisMaxSubspace;
    is_Computed.set(ConeProperty::MaximalSubspace);    
    pointed = (Basis_Max_Subspace.nr_of_rows()==0);
    is_Computed.set(ConeProperty::IsPointed);
    is_simplicial = nr_gen == dim;
    deg1_extreme_rays = false;
    deg1_generated = false;
    deg1_generated_computed = false;
    deg1_triangulation = false;
    deg1_hilbert_basis = false;
    
    reset_tasks();
    
    if (!Extreme_Rays_Ind.empty()) { // only then we can assume that all entries on C.Supp.. are relevant
        Support_Hyperplanes.swap(C.SupportHyperplanes);
        // there may be duplicates in the coordinates of the Full_Cone
        Support_Hyperplanes.remove_duplicate_and_zero_rows();
        is_Computed.set(ConeProperty::SupportHyperplanes);
    }
    if(!C.do_only_Deg1_Elements){
        Hilbert_Basis.swap(C.Hilbert_Basis);
        is_Computed.set(ConeProperty::HilbertBasis);
    }
    else {
        Deg1_Elements.swap(C.Hilbert_Basis);
        is_Computed.set(ConeProperty::Deg1Elements);
    }
    if(dim==0){            //correction needed to include the 0 cone;
        multiplicity = 1;
#ifdef ENFNORMALIZ
        renf_multiplicity=1;
#endif
        Hilbert_Series.add(vector<num_t>(1,1),vector<denom_t>());
        is_Computed.set(ConeProperty::HilbertSeries);
    }
    pyr_level=-1;
    Top_Cone=this;
    God_Father=this;
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=i;
    totalNrSimplices=0;
    TriangulationBufferSize=0;
    CandidatesSize=0;
    detSum = 0;
    shift = 0;
    
    tri_recursion=false;
    
    // nextGen=0;
    
    inhomogeneous=C.inhomogeneous;
    
    level0_dim=dim; // must always be defined
    
    use_existing_facets=false;
    start_from=0;
    old_nr_supp_hyps=0;
    verbose = C.verbose;
    OldCandidates.dual=false;
    OldCandidates.verbose=verbose;
    NewCandidates.dual=false;
    NewCandidates.verbose=verbose;
    
    descent_level=0;
    // approx_level = 1; ???? Noch gebraucht ???

    is_approximation=false;
    
    don_t_add_hyperplanes=false;
    take_time_of_large_pyr=false;
    
    verbose=C.verbose;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_grading_after_dual_mode(){

    if(dim>0 && Grading.size()>0 && !isComputed(ConeProperty::Grading)) {
        if(isComputed(ConeProperty::Generators)){
            vector<Integer> degrees=Generators.MxV(Grading);
            vector<Integer> levels;
            if(inhomogeneous)
                levels=Generators.MxV(Truncation);
            size_t i=0;
            for(;i<degrees.size();++i){
                if(degrees[i]<=0 &&(!inhomogeneous || levels[i]==0))
                    break;
            }
            if(i==degrees.size())
                is_Computed.set(ConeProperty::Grading);
        }
        else if(isComputed(ConeProperty::HilbertBasis)){
            auto hb=Hilbert_Basis.begin();
            for(;hb!=Hilbert_Basis.end();++hb){
                if(v_scalar_product(*hb,Grading)<=0 && (!inhomogeneous || v_scalar_product(*hb,Truncation)==0))
                    break;
            }
            if(hb==Hilbert_Basis.end())
                is_Computed.set(ConeProperty::Grading);
        }   
    }
    if(isComputed(ConeProperty::Deg1Elements)){
        auto hb=Deg1_Elements.begin();
        for(;hb!=Deg1_Elements.end();++hb){
            if(v_scalar_product(*hb,Grading)<=0)
                break;
        }
        if(hb==Deg1_Elements.end())
            is_Computed.set(ConeProperty::Grading);
    }

    if(Grading.size()>0 && !isComputed(ConeProperty::Grading)){
        throw BadInputException("Grading not positive on pointed cone.");
    }
}

template<typename Integer>
void Full_Cone<Integer>::dual_mode() {
    
    omp_start_level=omp_get_level();
    
    if(dim==0){
        set_zero_cone();
        return;
    }

    use_existing_facets=false; // completely irrelevant here
    start_from=0;
    old_nr_supp_hyps=0;
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    compute_class_group();
    
    check_grading_after_dual_mode();
    
    compute_automorphisms();
        
    if(dim>0 && !inhomogeneous){
        deg1_check();
        if (isComputed(ConeProperty::Grading) && !isComputed(ConeProperty::Deg1Elements)) {
            if (verbose) { 
                verboseOutput() << "Find degree 1 elements" << endl;
            }
            select_deg1_elements();
        }
    }
    
    if(dim==0){
        deg1_extreme_rays = deg1_generated = true;
        Grading=vector<Integer>(dim);
        is_Computed.set(ConeProperty::IsDeg1ExtremeRays);
        deg1_generated_computed = true;
        is_Computed.set(ConeProperty::Grading);
    }
    if(!inhomogeneous && isComputed(ConeProperty::HilbertBasis)){
        if (isComputed(ConeProperty::Grading)) 
            check_deg1_hilbert_basis();
    }

    if (inhomogeneous && isComputed(ConeProperty::Generators)) {
       set_levels();
       find_level0_dim();
       find_module_rank();
    }
    
    if (inhomogeneous && !isComputed(ConeProperty::Generators) && isComputed(ConeProperty::HilbertBasis)) {
       find_level0_dim_from_HB();
       find_module_rank();
    }
    
    use_existing_facets=false;
    start_from=0;
}

//---------------------------------------------------------------------------

/* constructor for pyramids */
template<typename Integer>
Full_Cone<Integer>::Full_Cone(Full_Cone<Integer>& C, const vector<key_t>& Key) {
    
    omp_start_level=C.omp_start_level;

    Generators = C.Generators.submatrix(Key);
    if(using_renf<Integer>() && C.Generators_float.nr_of_rows()>0)
        Generators_float = C.Generators_float.submatrix(Key);        
    dim = Generators.nr_of_columns();
    nr_gen = Generators.nr_of_rows();
    has_generator_with_common_divisor = C.has_generator_with_common_divisor;
    is_simplicial = nr_gen == dim;
    
    Top_Cone=C.Top_Cone; // relate to top cone
    C.God_Father=C.God_Father;
    Top_Key.resize(nr_gen);
    for(size_t i=0;i<nr_gen;i++)
        Top_Key[i]=C.Top_Key[Key[i]];

  
    multiplicity = 0;
#ifdef ENFNORMALIZ
    renf_multiplicity=0;
#endif
    
    Extreme_Rays_Ind = vector<bool>(nr_gen,false);
    is_Computed.set(ConeProperty::ExtremeRays, C.isComputed(ConeProperty::ExtremeRays));
    if(isComputed(ConeProperty::ExtremeRays))
        for(size_t i=0;i<nr_gen;i++)
            Extreme_Rays_Ind[i]=C.Extreme_Rays_Ind[Key[i]];
    // in_triang = vector<bool> (nr_gen,false); // now in build_cone
    deg1_triangulation = true;
    
    Grading=C.Grading;
    is_Computed.set(ConeProperty::Grading, C.isComputed(ConeProperty::Grading));
    Order_Vector=C.Order_Vector;

    // Note: For the computation of pyramids we do not call primal_algorithm.
    // Therefore it is necessary to set do_triangulation etc. here (and not only
    // the computation goals).
    do_triangulation=C.do_triangulation;
    do_partial_triangulation=C.do_partial_triangulation;
    do_only_multiplicity=C.do_only_multiplicity;
    stop_after_cone_dec=C.stop_after_cone_dec;
    
    // now the computation goals
    do_extreme_rays=false;
    do_triangulation_size=C.do_triangulation;
    do_determinants=C.do_determinants;
    do_multiplicity=C.do_multiplicity;
    do_deg1_elements=C.do_deg1_elements;
    do_h_vector=C.do_h_vector;
    do_Hilbert_basis=C.do_Hilbert_basis;
    keep_triangulation=C.keep_triangulation;
    do_evaluation=C.do_evaluation;
    do_Stanley_dec=C.do_Stanley_dec;
    do_bottom_dec=false;
    keep_order=true;
    do_all_hyperplanes=true; //  must be reset for non-recursive pyramids
    use_existing_facets=false;

    // not used in a pyramid, but set for precaution
    deg1_extreme_rays = false;
    deg1_generated = false;
    deg1_generated_computed = false;
    deg1_hilbert_basis = false;    
    inhomogeneous=C.inhomogeneous;   // at present not used in proper pyramids
    
    is_pyramid=true;
    
    pyr_level=C.pyr_level+1;
    descent_level=0; // should never be used in pyramids
    store_level = C.store_level;

    totalNrSimplices=0;
    detSum = 0;
    multiplicity = 0;
    shift = C.shift;
    level0_dim = C.level0_dim; // must always be defined
    
    if(C.gen_degrees.size()>0){ // now we copy the degrees
        gen_degrees.resize(nr_gen);
        if(C.do_h_vector || (!using_GMP<Integer>() && !using_renf<Integer>()) )
            gen_degrees_long.resize(nr_gen);
        for (size_t i=0; i<nr_gen; i++) {
            gen_degrees[i] = C.gen_degrees[Key[i]];
            if(C.do_h_vector || (!using_GMP<Integer>() && !using_renf<Integer>()) )
                gen_degrees_long[i] = C.gen_degrees_long[Key[i]];
        }
    }
    if(C.gen_levels.size()>0){ // now we copy the levels
        gen_levels.resize(nr_gen);
        for (size_t i=0; i<nr_gen; i++) {
            gen_levels[i] = C.gen_levels[Key[i]];
        }
    }
    
    TriangulationBufferSize=0; // not used in pyramids
    CandidatesSize=0; // ditto
    
    recursion_allowed=C.recursion_allowed; // must be reset if necessary 
    // multithreaded_pyramid=false; // SEE ABOVE
    
    Comparisons.reserve(nr_gen);
    nrTotalComparisons=0;
    
    start_from=0;
    old_nr_supp_hyps=0;
    verbose=false;
    OldCandidates.dual=false;
    OldCandidates.verbose=verbose;
    NewCandidates.dual=false;
    NewCandidates.verbose=verbose;

    is_approximation = C.is_approximation;
    
    is_global_approximation = C.is_global_approximation;
    
    do_bottom_dec=false;
    suppress_bottom_dec=false;
    keep_order=true;
    
    time_measured=C.time_measured;
    ticks_comp_per_supphyp=C.ticks_comp_per_supphyp;
    ticks_rank_per_row=C.ticks_rank_per_row;
    
    don_t_add_hyperplanes=false;
    take_time_of_large_pyr=false;
    renf_degree=C.renf_degree;
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
bool Full_Cone<Integer>::isDeg1ExtremeRays() const{
    return deg1_extreme_rays;
}

template<typename Integer>
bool Full_Cone<Integer>::isDeg1HilbertBasis() const{
    return deg1_hilbert_basis;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Full_Cone<Integer>::getGrading() const{
    return Grading;
}

//---------------------------------------------------------------------------

template<typename Integer>
mpq_class Full_Cone<Integer>::getMultiplicity()const{
    return multiplicity;
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Full_Cone<Integer>::getShift()const{
    return shift;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Full_Cone<Integer>::getModuleRank()const{
    return module_rank;
}


//---------------------------------------------------------------------------

template<typename Integer>
const Matrix<Integer>& Full_Cone<Integer>::getGenerators()const{
    return Generators;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<bool> Full_Cone<Integer>::getExtremeRays()const{
    return Extreme_Rays_Ind;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Full_Cone<Integer>::getNrExtremeRays() const{
    size_t n=0;
    for(size_t j=0;j<nr_gen;++j)
        if(Extreme_Rays_Ind[j])
            ++n;
    return n;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getSupportHyperplanes()const{
    return Support_Hyperplanes;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getHilbertBasis()const{
    if(Hilbert_Basis.empty())
        return Matrix<Integer>(0,dim);
    else
        return Matrix<Integer>(Hilbert_Basis);
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getModuleGeneratorsOverOriginalMonoid()const{
    if(ModuleGeneratorsOverOriginalMonoid.empty())
        return Matrix<Integer>(0,dim);
    else
        return Matrix<Integer>(ModuleGeneratorsOverOriginalMonoid);
}


//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getDeg1Elements()const{
    if(Deg1_Elements.empty())
        return Matrix<Integer>(0,dim);
    else
        return Matrix<Integer>(Deg1_Elements);
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getExcludedFaces()const{
    return(ExcludedFaces);
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
    verboseOutput()<<"\nHilbert Series  is:\n";
    verboseOutput()<<Hilbert_Series;
}

#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template class Full_Cone<long>;
#endif
template class Full_Cone<long long>;
template class Full_Cone<mpz_class>;

#ifdef ENFNORMALIZ
template class Full_Cone<renf_elem_class>;
#endif

} //end namespace


/*

//---------------------------------------------------------------------------
// version with isomorphism classes -- has no real effect
template<typename Integer>
void Full_Cone<Integer>::get_cone_over_facet_HB(const vector<Integer>& fixed_point, const vector<key_t>& facet_key, 
                                      const key_t facet_nr, list<vector<Integer> >& Facet_HB){
  
  
    Matrix<Integer> Facet_Gens(0,dim);
    Facet_Gens.append(fixed_point);
    Facet_Gens.append(Generators.submatrix(facet_key));   
    
    for(long i=0;i<descent_level+1;++i)
        cout << "$$$$$$  ";
    cout << " " << Facet_Gens.nr_of_rows() << endl;
    cout << "Height FP over facet " << v_scalar_product(fixed_point,Support_Hyperplanes[facet_nr]) << endl;
    
    Full_Cone ConeOverFacet(Facet_Gens);
    ConeOverFacet.verbose=verbose;

    if(isComputed(ConeProperty::Grading)){
      ConeOverFacet.Grading=Grading;
      ConeOverFacet.is_Computed.set(ConeProperty::Grading);
    }
     ConeOverFacet.descent_level=descent_level+1;
    ConeOverFacet.Mother=&(*this);
    ConeOverFacet.God_Father=God_Father;
    ConeOverFacet.exploit_automorphisms=true;
    ConeOverFacet.full_automorphisms=full_automorphisms;
    ConeOverFacet.ambient_automorphisms=ambient_automorphisms;
    ConeOverFacet.input_automorphisms=input_automorphisms;
    ConeOverFacet.Embedding=Embedding;
    ConeOverFacet.keep_order=true;
    ConeOverFacet.Support_Hyperplanes=push_supphyps_to_cone_over_facet(fixed_point,facet_nr);
    // ConeOverFacet.do_Hilbert_basis=true;
    ConeOverFacet.compute();
    if(ConeOverFacet.isComputed(ConeProperty::HilbertBasis)){
        Facet_HB.splice(Facet_HB.begin(),ConeOverFacet.Hilbert_Basis);
        return;        
    }
    bool found;
    const IsoType<Integer>& face_class=God_Father->FaceClasses.find_type(ConeOverFacet,found);
    if(found){
        ConeOverFacet.import_HB_from(face_class);
        Facet_HB.clear();
        Facet_HB.splice(Facet_HB.begin(),ConeOverFacet.Hilbert_Basis);
        if(ConeOverFacet.isComputed(ConeProperty::HilbertBasis))
            return;
    } 

    Full_Cone Facet_2(Facet_Gens);
    Facet_2.Automs=ConeOverFacet.Automs;
    Facet_2.is_Computed.set(ConeProperty::Automorphisms);
    Facet_2.Embedding=Embedding;
    Facet_2.full_automorphisms=full_automorphisms;
    Facet_2.ambient_automorphisms=ambient_automorphisms;
    Facet_2.input_automorphisms=input_automorphisms;
    Facet_2.exploit_automorphisms=true;
    Facet_2.keep_order=true;
    Facet_2.Extreme_Rays_Ind=ConeOverFacet.Extreme_Rays_Ind;
    Facet_2.is_Computed.set(ConeProperty::ExtremeRays);
    Facet_2.Support_Hyperplanes=ConeOverFacet.Support_Hyperplanes;
    Facet_2.nrSupport_Hyperplanes=ConeOverFacet.nrSupport_Hyperplanes;
    Facet_2.is_Computed.set(ConeProperty::SupportHyperplanes);
    Facet_2.verbose=verbose;
    Facet_2.descent_level=descent_level+1;
    Facet_2.full_automorphisms=full_automorphisms;
    Facet_2.ambient_automorphisms=ambient_automorphisms;
    Facet_2.input_automorphisms=input_automorphisms;
    if(isComputed(ConeProperty::Grading)){
        Facet_2.Grading=Grading;
        Facet_2.is_Computed.set(ConeProperty::Grading);
    }
    Facet_2.Mother=&(*this);
    Facet_2.God_Father=God_Father;
    Facet_2.do_Hilbert_basis=true;
    Facet_2.compute();
    bool added;
    God_Father->FaceClasses.add_type(Facet_2, added);
    Facet_HB.clear();
    Facet_HB.splice(Facet_HB.begin(),Facet_2.Hilbert_Basis);
    return;
}
*/

//---------------------------------------------------------------------------

/*
// not used at present
template<typename Integer>
void Full_Cone<Integer>::import_HB_from(const IsoType<Integer>& copy){

    assert(isComputed(ConeProperty::Automorphisms));
    
    size_t N=copy.getHilbertBasis().nr_of_rows();
    if(N==0){
        is_Computed.set(ConeProperty::HilbertBasis);
        return;    
    }
        
    assert(Hilbert_Basis.empty());
    for(size_t i=0;i<nr_gen;++i)
        Hilbert_Basis.push_back(Generators[i]);
    
    vector<key_t> CanBasisKey=Generators.max_rank_submatrix_lex(Automs.CanLabellingGens);    
    Matrix<Integer> Transform=copy.getCanTransform().multiplication(Generators.submatrix(CanBasisKey));
    Integer D=Transform.matrix_gcd();
    if(D!=copy.getCanDenom()) // not liftable
        return;
    Transform.scalar_division(D);
    for(size_t i=0;i<N;++i){
        Hilbert_Basis.push_back(Transform.VxM(copy.getHilbertBasis()[i]));        
    }
    
    is_Computed.set(ConeProperty::HilbertBasis);
    return;     
}
*/
