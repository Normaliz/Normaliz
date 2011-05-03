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

#include "full_cone.h"
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
void Full_Cone<Integer>::add_hyperplane(const size_t& ind_gen, const FMDATA & positive,const FMDATA & negative){
	size_t k;
	
	// NEW: indgen is the index of the generator being inserted    
	
	// vector<Integer> hyperplane(hyp_size,0); // initialized with 0
	
	FMDATA NewHypIndVal; NewHypIndVal.Hyp.resize(dim); NewHypIndVal.GenInHyp.resize(nr_gen);
	
	Integer used_for_tests;
	if (test_arithmetic_overflow) {  // does arithmetic tests
		for (k = 0; k <dim; k++) {
			NewHypIndVal.Hyp[k]=positive.ValNewGen*negative.Hyp[k]-negative.ValNewGen*positive.Hyp[k];
			used_for_tests =(positive.ValNewGen%overflow_test_modulus)*(negative.Hyp[k]%overflow_test_modulus)-(negative.ValNewGen%overflow_test_modulus)*(positive.Hyp[k]%overflow_test_modulus);
			if (((NewHypIndVal.Hyp[k]-used_for_tests) % overflow_test_modulus)!=0) {
				errorOutput()<<"Arithmetic failure in Full_cone::add_hyperplane. Possible arithmetic overflow.\n";
				throw ArithmeticException();
			}
		}
	}
	else  {                      // no arithmetic tests
		for (k = 0; k <dim; k++) {
			NewHypIndVal.Hyp[k]=positive.ValNewGen*negative.Hyp[k]-negative.ValNewGen*positive.Hyp[k];
		}
	}
	NewHypIndVal.Hyp=v_make_prime(NewHypIndVal.Hyp);
	NewHypIndVal.ValNewGen=0; // not really needed, only for completeness
	
	NewHypIndVal.GenInHyp=positive.GenInHyp & negative.GenInHyp; // new hyperplane contains old gen iff both pos and neg do
	NewHypIndVal.GenInHyp.set(ind_gen);  // new hyperplane contains new generator
	
	#pragma omp critical(HYPERPLANE)
	HypIndVal.push_back(NewHypIndVal);
}


//---------------------------------------------------------------------------


template<typename Integer>
void Full_Cone<Integer>::transform_values(const size_t& ind_gen){

	//to see if possible to replace the function .end with constant iterator since push-back is performed.

	// NEW: ind_gen is the index of the generator being inserted

	register size_t i,k,nr_zero_i;
	register size_t subfacet_dim=dim-2; // NEW dimension of subfacet
	register size_t facet_dim=dim-1; // NEW dimension of facet
	
	const bool tv_verbose = false; //verbose && Support_Hyperplanes.size()>10000; //verbose in this method call
	
		
	// preparing the computations
	list <FMDATA*> l_Pos_Simp,l_Pos_Non_Simp;
	list <FMDATA*> l_Neg_Simp,l_Neg_Non_Simp;
	list <FMDATA*> l_Neutral_Simp, l_Neutral_Non_Simp;
	
	boost::dynamic_bitset<> Zero_Positive(nr_gen),Zero_Negative(nr_gen);

	bool simplex;
	bool ranktest;
	
	if (tv_verbose) verboseOutput()<<"transform_values: create SZ,Z,PZ,P,NS,N"<<endl<<flush;
	size_t ipos=0;
	
	typename list<FMDATA>::iterator ii = HypIndVal.begin();
	
	size_t listsize=HypIndVal.size();

	for (size_t kk=0; kk<listsize; ++kk) {
		for(;kk > ipos; ++ipos, ++ii) ;
		for(;kk < ipos; --ipos, --ii) ;
		simplex=false;
		
		nr_zero_i=ii->GenInHyp.count();
		if(ii->ValNewGen>0)
			Zero_Positive|=ii->GenInHyp;
		else if(ii->ValNewGen<0)
			Zero_Negative|=ii->GenInHyp;        
		if (nr_zero_i==dim-1)
			simplex=true;
			
		if (ii->ValNewGen==0) {
			ii->GenInHyp.set(ind_gen);  // Must be set explicitly !!
			if (simplex) {
				l_Neutral_Simp.push_back(&(*ii));
			}   else {
				l_Neutral_Non_Simp.push_back(&(*ii));
			}
		}
		else if (ii->ValNewGen>0) {
			if (simplex) {
				l_Pos_Simp.push_back(&(*ii));
			} else {
				l_Pos_Non_Simp.push_back(&(*ii));
			}
		} 
		else if (ii->ValNewGen<0) {
			if (simplex) {
				l_Neg_Simp.push_back(&(*ii));
			} else {
				l_Neg_Non_Simp.push_back(&(*ii));
			}
		}
	}
	
	boost::dynamic_bitset<> Zero_PN(nr_gen);
	Zero_PN=Zero_Positive & Zero_Negative;
	
	if (tv_verbose) verboseOutput()<<"transform_values: copy to vector"<<endl;

	size_t nr_PosSimp  = l_Pos_Simp.size();
	size_t nr_PosNSimp = l_Pos_Non_Simp.size();
	size_t nr_NegSimp  = l_Neg_Simp.size();
	size_t nr_NegNSimp = l_Neg_Non_Simp.size();
	size_t nr_NeuSimp  = l_Neutral_Simp.size();
	size_t nr_NeuNSimp = l_Neutral_Non_Simp.size();

	ranktest=false;
	if (nr_PosNSimp+nr_NegNSimp+nr_NeuNSimp > dim*dim*dim/6)
		ranktest=true;

	vector <FMDATA*> Pos_Simp(nr_PosSimp);
	vector <FMDATA*> Pos_Non_Simp(nr_PosNSimp);
	vector <FMDATA*> Neg_Simp(nr_NegSimp);
	vector <FMDATA*> Neg_Non_Simp(nr_NegNSimp);
	vector <FMDATA*> Neutral_Simp(nr_NeuSimp);
	vector <FMDATA*> Neutral_Non_Simp(nr_NeuNSimp);

	for (k = 0; k < Pos_Simp.size(); k++) {
		Pos_Simp[k]=l_Pos_Simp.front();
		l_Pos_Simp.pop_front();
	}

	for (k = 0; k < Pos_Non_Simp.size(); k++) {
		Pos_Non_Simp[k]=l_Pos_Non_Simp.front();
		l_Pos_Non_Simp.pop_front();
	}

	for (k = 0; k < Neg_Simp.size(); k++) {
		Neg_Simp[k]=l_Neg_Simp.front();
		l_Neg_Simp.pop_front();
	}

	for (k = 0; k < Neg_Non_Simp.size(); k++) {
		Neg_Non_Simp[k]=l_Neg_Non_Simp.front();
		l_Neg_Non_Simp.pop_front();
	}

	for (k = 0; k < Neutral_Simp.size(); k++) {
		Neutral_Simp[k]=l_Neutral_Simp.front();
		l_Neutral_Simp.pop_front();
	}

	for (k = 0; k < Neutral_Non_Simp.size(); k++) {
		Neutral_Non_Simp[k]=l_Neutral_Non_Simp.front();
		l_Neutral_Non_Simp.pop_front();
	}

	if (tv_verbose) verboseOutput()<<"PS "<<nr_PosSimp<<" P "<<nr_PosNSimp<<" NS "<<nr_NegSimp<<" N "<<nr_NegNSimp<<" ZS "<<nr_NeuSimp<<" Z "<<nr_NeuNSimp<<endl<<flush;

	if (tv_verbose) verboseOutput()<<"transform_values: fill multimap with subfacets of NS"<<endl<<flush;
	
	multimap < boost::dynamic_bitset<>, int> Neg_Subfacet_Multi;

	boost::dynamic_bitset<> zero_i(nr_gen);
	boost::dynamic_bitset<> subfacet(nr_gen);

	for (i=0; i<nr_NegSimp;i++){
		zero_i=Zero_PN & Neg_Simp[i]->GenInHyp;
		nr_zero_i=zero_i.count();
		
		if(nr_zero_i==subfacet_dim) // NEW This case treated separately
			Neg_Subfacet_Multi.insert(pair <boost::dynamic_bitset<>, int> (zero_i,i));
			
		else{       
			for (k =0; k<nr_gen; k++) {  //TODO use BOOST ROUTINE
				if(zero_i.test(k)) {              
					subfacet=zero_i;
					subfacet.reset(k);  // remove k-th element from facet to obtain subfacet
					Neg_Subfacet_Multi.insert(pair <boost::dynamic_bitset<>, int> (subfacet,i));
				}
			}
		}
	}


	if (tv_verbose) verboseOutput()<<"transform_values: go over multimap of size "<< Neg_Subfacet_Multi.size() <<endl<<flush;

	multimap < boost::dynamic_bitset<>, int > ::iterator jj;
	multimap < boost::dynamic_bitset<>, int > ::iterator del;
	jj =Neg_Subfacet_Multi.begin();                               // remove negative subfecets shared
	while (jj!= Neg_Subfacet_Multi.end()) {                       // by two neg simpl facets
		del=jj++;
		if (jj!=Neg_Subfacet_Multi.end() && (*jj).first==(*del).first) {   //delete since is the intersection of two negative simplicies
			Neg_Subfacet_Multi.erase(del);
			del=jj++;
			Neg_Subfacet_Multi.erase(del);
		}
	}

	if (tv_verbose) verboseOutput()<<"transform_values: singlemap size "<<Neg_Subfacet_Multi.size()<<endl<<flush;
	
	size_t nr_NegSubfMult = Neg_Subfacet_Multi.size();
	size_t nr_NegSubf;
	map < boost::dynamic_bitset<>, int > Neg_Subfacet;
	
	#pragma omp parallel private(jj)
	{
	size_t i,j,k,t,nr_zero_i;
	boost::dynamic_bitset<> subfacet(dim-2);
	jj = Neg_Subfacet_Multi.begin();
	size_t jjpos=0;

	map < boost::dynamic_bitset<>, int > ::iterator last_inserted=Neg_Subfacet.begin(); // used to speedup insertion into the new map
	bool found;
	#pragma omp for schedule(dynamic)
	for (size_t j=0; j<nr_NegSubfMult; ++j) {             // remove negative subfacets shared
		for(;j > jjpos; ++jjpos, ++jj) ;                // by non-simpl neg or neutral facets 
		for(;j < jjpos; --jjpos, --jj) ;

		subfacet=(*jj).first;
		found=false; 
		for (i = 0; i <Neutral_Simp.size(); i++) {
			found=subfacet.is_subset_of(Neutral_Simp[i]->GenInHyp);
			if(found)
				break;
		}
		if (!found) {
			for (i = 0; i <Neutral_Non_Simp.size(); i++) {
				found=subfacet.is_subset_of(Neutral_Non_Simp[i]->GenInHyp);
				if(found)
					break;                    
				}
			if(!found) {
				for (i = 0; i <Neg_Non_Simp.size(); i++) {
					found=subfacet.is_subset_of(Neg_Non_Simp[i]->GenInHyp);
					if(found)
						break; 
				}
			}
		}
		if (!found) {
			#pragma omp critical(NEGATIVE_SUBFACET)
			{last_inserted=Neg_Subfacet.insert(last_inserted,*jj);}
		}
	}
	
	#pragma omp single
	{nr_NegSubf = Neg_Subfacet.size();}

	#pragma omp single nowait
	if (tv_verbose) {
		verboseOutput()<<"transform_values: reduced map size "<<nr_NegSubf<<endl<<flush;
	} 
	#pragma omp single nowait
	{Neg_Subfacet_Multi.clear();}

	#pragma omp single nowait
	if (tv_verbose) {
		verboseOutput()<<"transform_values: PS vs NS"<<endl<<flush;
	}
	
	boost::dynamic_bitset<> zero_i(nr_gen);
	map <boost::dynamic_bitset<>, int> ::iterator jj_map;

	#pragma omp for schedule(dynamic) nowait           // Now matching positive and negative (sub)facets
	for (i =0; i<nr_PosSimp; i++){ //Positive Simp vs.Negative Simp
		zero_i=Pos_Simp[i]->GenInHyp & Zero_PN;
		nr_zero_i=zero_i.count();
		
		if (nr_zero_i==subfacet_dim) {                 // NEW slight change in logic. Positive simpl facet shared at most
			jj_map=Neg_Subfacet.find(zero_i);           // one subfacet with negative simpl facet
			if (jj_map!=Neg_Subfacet.end()) {
				add_hyperplane(ind_gen,*Pos_Simp[i],*Neg_Simp[(*jj_map).second]);
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
						add_hyperplane(ind_gen,*Pos_Simp[i],*Neg_Simp[(*jj_map).second]);
						(*jj_map).second = -1;
					}
				}
			}
		}
	}

	#pragma omp single nowait
	if (tv_verbose) {
		verboseOutput()<<"transform_values: NS vs P"<<endl<<flush;
	}

//  for (jj_map = Neg_Subfacet.begin(); jj_map != Neg_Subfacet.end(); ++jj_map)  //Neg_simplex vs. Pos_Non_Simp
	jj_map = Neg_Subfacet.begin();
	jjpos=0;
	#pragma omp for schedule(dynamic) nowait
	for (size_t j=0; j<nr_NegSubf; ++j) {
		for( ; j > jjpos; ++jjpos, ++jj_map) ;
		for( ; j < jjpos; --jjpos, --jj_map) ;

		if ( (*jj_map).second != -1 ) {  // skip used subfacets
			for (i = 0; i <Pos_Non_Simp.size(); i++) {
				if(jj_map->first.is_subset_of(Pos_Non_Simp[i]->GenInHyp)){
					add_hyperplane(ind_gen,*Pos_Non_Simp[i],*Neg_Simp[(*jj_map).second]);
					break;
				}
			}
		}
	}
	
	#pragma omp single nowait
	if (tv_verbose) {
		verboseOutput()<<"transform_values: PS vs N"<<endl<<flush;
	}

	vector<size_t> key(nr_gen);
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
			for (j=0; j<Neg_Non_Simp.size(); j++){ // search negative facet with common subfacet
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
					add_hyperplane(ind_gen,*Pos_Simp[i],*Neg_Non_Simp[j]);
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
	
	bool exactly_two;
	FMDATA *hp_i, *hp_j, *hp_t; // pointers to current hyperplanes
	
	size_t missing_bound, nr_common_zero;
	boost::dynamic_bitset<> common_zero(nr_gen);
	vector<size_t> common_key(nr_gen);
	
	#pragma omp for schedule(dynamic) nowait
	for (size_t i =0; i<nr_PosNSimp; i++){ //Positive Non Simp vs.Negative Non Simp

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
				for (j=0; j<nr_NegNSimp; j++){
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
							Test.write(k+1,Generators.read(common_key[k]+1)); 

						if (Test.rank_destructiv()<subfacet_dim) {
							exactly_two=false;
						}
					} // ranktest
					else{                 // now the comparison test
						for (t=0;t<nr_PosNSimp;t++){
							hp_t=Pos_Non_Simp[t];
							if (t!=i && common_zero.is_subset_of(hp_t->GenInHyp)) {                                
								exactly_two=false;
								break;
							}
						}
						if (exactly_two) {
							for (t=0;t<Neg_Non_Simp.size();t++){
								hp_t=Neg_Non_Simp[t];
								if (t!=j && common_zero.is_subset_of(hp_t->GenInHyp)) {  
									exactly_two=false;
									break;
								}
							}
						}
						if (exactly_two) {
							for (t=0;t<nr_NeuNSimp;t++){
								hp_t=Neutral_Non_Simp[t];
								if(common_zero.is_subset_of(hp_t->GenInHyp)) {  
									exactly_two=false;
									break;
								}
							}
							
						}                        
					} // else
					if (exactly_two) {  //intersection of i and j is a subfacet
						add_hyperplane(ind_gen,*hp_i,*hp_j);
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
void Full_Cone<Integer>::add_simplex(const size_t& new_generator){
	typename list<FMDATA>::const_iterator i=HypIndVal.begin();
	typename list< pair< vector<size_t>, Integer> >::const_iterator j;
	size_t nr_zero_i, nr_nonzero_i, not_in_i=0, l, k; 
	size_t s, Triangulation_size=Triangulation.size(); // the size we start with
	vector<size_t> key(dim);

	size_t ipos=0;
	size_t listsize=HypIndVal.size();

	#pragma omp parallel for private(j,nr_zero_i,nr_nonzero_i,l,k,s) firstprivate(ipos, i, key, not_in_i) schedule(dynamic)
	for (size_t kk=0; kk<listsize; ++kk) {
		for(;kk > ipos; ++ipos, ++i) ;
		for(;kk < ipos; --ipos, --i) ;

		if ((*i).ValNewGen<0) {
			nr_zero_i=(*i).GenInHyp.count();

			if (nr_zero_i==dim-1) { //simplicial
				l=0;
				for (k = 0; k <nr_gen; k++) {
					if ((*i).GenInHyp[k]==1) {
						key[l]=k+1;
						l++;
					}
				}
				key[dim-1]=new_generator+1;

				store_key(key,i->ValNewGen);
				if(!keep_triangulation){
					Simplex<Integer> simp(key);        
					simp.evaluate(*this,i->ValNewGen);
				}
			}
			else {
				j =Triangulation.begin();
				for (s=0; s<Triangulation_size; s++){
					key=j->first;
					nr_nonzero_i=0;
					k=0;
					do{
						if ( !(*i).GenInHyp.test(key[k]-1)) {
							nr_nonzero_i++;
							not_in_i=k;
						}
						k++;
					} while((k<dim)&&(nr_nonzero_i<2));
					
					if (nr_nonzero_i<=1){
						key[not_in_i]=new_generator+1;

						store_key(key,i->ValNewGen);
						if(!keep_triangulation){
							Simplex<Integer> simp(key);        
							simp.evaluate(*this,i->ValNewGen);
						}
					}

					j++; 
									   
				}  // s
				
				
			} // else
			
		} // if < 0
		
	} // for kk
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::store_key(const vector<size_t>& key, const Integer& height){
	pair<vector<size_t>,Integer> newsimplex;
	newsimplex.first=key;
	newsimplex.second=height;
	#pragma omp critical(TRIANG)
	Triangulation.push_back(newsimplex);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::process_pyramids(const size_t ind_gen,const bool recursive){

	typename list< FMDATA >::iterator l=HypIndVal.begin();

#ifdef _WIN32 //for 32 and 64 bit windows

	size_t lpos=0, listsize=HypIndVal.size();

   #pragma omp parallel for firstprivate(lpos,l) schedule(dynamic) 
	for (size_t k=0; k<listsize; k++) {
		for(;k > lpos; lpos++, l++) ;
		for(;k < lpos; lpos--, l--) ;

		// only triangulation of Pyramids of height >=2 needeed
		if(l->ValNewGen>=0 ||(!recursive && l->ValNewGen>=-1))
			continue;
	
		process_pyramid((*l), ind_gen, recursive);
	} //end for

#else         // all other systems

	if (is_pyramid) { //do not create a new parallel region
		for (; l != HypIndVal.end(); ++l) {
			// only triangulation of Pyramids of height >=2 needeed
			if(l->ValNewGen>=0 ||(!recursive && l->ValNewGen>=-1))
				continue;
		
			#pragma omp task firstprivate(l)
			{
				process_pyramid((*l), ind_gen, recursive);
			}
		} //end for
		#pragma omp taskwait
	} else {
		#pragma omp parallel if(!is_pyramid)
		{
			#pragma omp single
			{
				for (; l != HypIndVal.end(); ++l) {
	
					// only triangulation of Pyramids of height >=2 needeed
					if(l->ValNewGen>=0 ||(!recursive && l->ValNewGen>=-1))
						continue;
				
					#pragma omp task firstprivate(l)
					{
						process_pyramid((*l), ind_gen, recursive);
					}
			
				} //end for
			} //end single
			#pragma omp taskwait
		} //end parallel
	}
#endif //else _WIN32
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::process_pyramid(FMDATA& l, const size_t ind_gen,const bool recursive){

	vector<size_t> Pyramid_key;
	Pyramid_key.reserve(nr_gen);
	size_t i; 
	boost::dynamic_bitset<> in_Pyramid(nr_gen,false);           
	Pyramid_key.push_back(ind_gen+1);
	in_Pyramid.set(ind_gen);
	for(i=0;i<nr_gen;i++){
		if(in_triang[i] && v_scalar_product(l.Hyp,Generators.read(i+1))==0){ // we cannot trust that HypIndVal is
			Pyramid_key.push_back(i+1);                                      // up-to-date
			in_Pyramid.set(i);
		}
	}
	Full_Cone<Integer> Pyramid(*this,Generators.submatrix(Pyramid_key));
	Pyramid.do_triangulation= !recursive || do_triangulation;
	if(Pyramid.do_triangulation)
		Pyramid.do_partial_triangulation=false;
	Pyramid.build_cone(); 
	
	if(recursive && keep_triangulation){        
		typename list<pair<vector<size_t>,Integer> >::iterator pyr_simp=Pyramid.Triangulation.begin();
		pair<vector<size_t>,Integer> newsimplex;
		newsimplex.first=vector<size_t> (dim);
		for(;pyr_simp!=Pyramid.Triangulation.end();pyr_simp++){
			for(i=0;i<dim;i++)
				newsimplex.first[i]=Pyramid_key[pyr_simp->first[i]-1];
			newsimplex.second=pyr_simp->second;
			#pragma omp critical(TRIANG)
			Triangulation.push_back(newsimplex);          
		}
	}
	Pyramid.Triangulation.clear();        
	 
	if(recursive){         
		typename list<vector<Integer> >::iterator pyr_hyp = Pyramid.Support_Hyperplanes.begin();
		bool new_global_hyp;
		FMDATA NewHypIndVal;
		Integer test;   
		for(;pyr_hyp!=Pyramid.Support_Hyperplanes.end();pyr_hyp++){
			if(v_scalar_product(Generators.read(ind_gen+1),*pyr_hyp)>0)
				continue;
			new_global_hyp=true;
			for(i=0;i<nr_gen;i++){
				if(in_Pyramid.test(i) || !in_triang[i])
					continue;
				test=v_scalar_product(Generators.read(i+1),*pyr_hyp);
				if(test<=0){
					new_global_hyp=false;
					break;
				}
			
			}
			if(new_global_hyp){
				NewHypIndVal.Hyp=*pyr_hyp;                
				#pragma omp critical(HYPERPLANE)
				HypIndVal.push_back(NewHypIndVal);
			}
		}
	}
	Pyramid.Support_Hyperplanes.clear();
	
	if(do_h_vector) {
		#pragma omp critical(HVECTOR)
		H_Vector=v_add(H_Vector,Pyramid.H_Vector);
	}
	
	#pragma omp critical(MULTIPLICITY)
	multiplicity += Pyramid.multiplicity;

	#pragma omp critical(CANDIDATES)
	Candidates.splice(Candidates.begin(),Pyramid.Candidates);
	#pragma omp critical(HT1ELEMENTS)
	Ht1_Elements.splice(Ht1_Elements.begin(),Pyramid.Ht1_Elements);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::find_and_evaluate_start_simplex(){

	size_t i,j;
	Integer factor;

	
	Simplex<Integer> S = find_start_simplex();
	vector<size_t> key=S.read_key();   // generators indexed from 1
		
	// vector<bool> in_triang(nr_gen,false);
	for (i = 0; i < dim; i++) {
		in_triang[key[i]-1]=true;
	}
	   
	Matrix<Integer> H=S.read_support_hyperplanes();
	for (i = 0; i <dim; i++) {
		FMDATA NewHypIndVal; NewHypIndVal.Hyp.resize(dim); NewHypIndVal.GenInHyp.resize(nr_gen);
		NewHypIndVal.Hyp=H.read(i+1);
		for(j=0;j < dim;j++)
			if(j!=i)
				NewHypIndVal.GenInHyp.set(key[j]-1);
		HypIndVal.push_back(NewHypIndVal);
	}
	
	if(!is_pyramid)
		Order_Vector = vector<Integer>(dim,0);
	
	if(!is_pyramid && do_h_vector){
		Matrix<Integer> G=S.read_generators();
		//srand(12345);
		for(i=0;i<dim;i++){
			factor=(unsigned long)(2*(rand()%(2*dim))+3);
			for(j=0;j<dim;j++)
				Order_Vector[j]+=factor*G.read(i+1,j+1);        
		} 
	}

	//the volume is an upper bound for the height
	if(do_triangulation)
		store_key(key,-S.read_volume());
	if(!keep_triangulation) {
		Simplex<Integer> simp(key);
		simp.evaluate(*this, -S.read_volume());
	}
}

//int RekTiefe=-1;

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
		if (do_triangulation) verboseOutput()<<"with full triangulation ";
		if (!do_triangulation && !do_partial_triangulation) verboseOutput()<<"(only support hyperplanes) ";
		verboseOutput()<<"..."<<endl;
	}
	size_t i;
	
	bool pyramid_recursion=false;

//#pragma omp critical(REKTIEFE)
	//RekTiefe++;
	//cout << RekTiefe;

	// DECIDE WHETHER TO USE RECURSION
	size_t RecBoundSuppHyp = dim*dim*dim;
	RecBoundSuppHyp *= RecBoundSuppHyp * 10; //dim^6 * 10
	size_t bound_div = nr_gen-dim+1;
	if(bound_div > 3*dim) bound_div = 3*dim;
	RecBoundSuppHyp /= bound_div;

	size_t RecBoundTriang=RecBoundFactor; //dim*dim*dim*dim*RecBoundFactor;

//if(!is_pyramid) cout << "RecBoundSuppHyp = "<<RecBoundSuppHyp<<endl;

	find_and_evaluate_start_simplex();
	
	Integer scalar_product;
	bool new_generator;
	short round;

	for (round = 0; round <= 1; round++) {  //two times, first only KNOWN extreme rays are considered
		for (i = 0; i < nr_gen; i++) {
			if ((in_triang[i]==false)&&((round==1)||(Extreme_Rays[i]==true))) {
				new_generator=false;

				typename list< FMDATA >::iterator l=HypIndVal.begin();

				size_t nr_pos=0; size_t nr_neg=0;
				vector<Integer> L;
				size_t old_nr_supp_hyps=HypIndVal.size();                
				
				size_t lpos=0;
				size_t listsize=HypIndVal.size();
				#pragma omp parallel for private(L,scalar_product) firstprivate(lpos,l) schedule(dynamic)
				for (size_t k=0; k<listsize; k++) {
					for(;k > lpos; lpos++, l++) ;
					for(;k < lpos; lpos--, l--) ;

					L=Generators.read(i+1);
					scalar_product=v_scalar_product(L,(*l).Hyp);
					if (test_arithmetic_overflow && !v_test_scalar_product(L,(*l).Hyp,scalar_product,overflow_test_modulus)) {
						error_msg("error: Arithmetic failure in Full_cone::support_hyperplanes. Possible arithmetic overflow.\n");
						throw ArithmeticException();
					}

					l->ValNewGen=scalar_product;
					if (scalar_product<0) {
						new_generator=true;
						nr_neg++;
					}
					if (scalar_product>0) 
						nr_pos++;
				}  //end parallel for
				
				if(!new_generator)
					continue;
					
//				if(pyramid_recursion || nr_neg*nr_pos>RecBoundSuppHyp){
				if(pyramid_recursion || nr_neg*nr_pos>RecBoundSuppHyp || nr_neg*Triangulation.size()>RecBoundTriang){
					if(!pyramid_recursion && !keep_triangulation)
						Triangulation.clear();
					pyramid_recursion=true;
					process_pyramids(i,true); //recursive
				}
				else{
					if(do_triangulation)
						add_simplex(i); 
					if(do_partial_triangulation)
						process_pyramids(i,false); // non-recursive
					transform_values(i);
				}
				
				//removing the negative hyperplanes
				l=HypIndVal.begin();
				for (size_t j=0; j<old_nr_supp_hyps;j++){
					if (l->ValNewGen<0) 
						l=HypIndVal.erase(l);
					else 
						l++;
				}
				
				in_triang[i]=true;

				if (verbose && !is_pyramid) {
					verboseOutput() << "generator="<< i+1 <<" and "<<HypIndVal.size()<<" hyperplanes... ";
					if(keep_triangulation)
						verboseOutput() << Triangulation.size() << " simplices ";
					verboseOutput()<<endl;
				}
			}
		}
	}

	typename list<FMDATA>::const_iterator IHV=HypIndVal.begin();
	for(;IHV!=HypIndVal.end();IHV++){
		Support_Hyperplanes.push_back(IHV->Hyp);
	}

	} // end if (dim>0)
	
	HypIndVal.clear();
	
//#pragma omp critical(REKTIEFE)
//	RekTiefe--;

	is_Computed.set(ConeProperty::SupportHyperplanes);
	if(keep_triangulation) {
		if(!is_pyramid) {
			//sort the keys
			typename list< pair<vector<size_t>, Integer> >::iterator it = Triangulation.begin();
			while (it!=Triangulation.end()) {
				sort(it->first.begin(),it->first.end());
				++it;
			}
		}
		is_Computed.set(ConeProperty::Triangulation);
	}
	
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::extreme_rays_and_ht1_check() {
	check_pointed();
	if(!pointed) return;
	compute_extreme_rays();
	check_ht1_extreme_rays();
	check_ht1_generated();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_support_hyperplanes(){
	if(is_Computed.test(ConeProperty::SupportHyperplanes))
		return;
	bool save_tri=do_triangulation;
	bool save_part_tri=do_partial_triangulation;
	do_triangulation=false;
	do_partial_triangulation=false;
	build_cone();
	do_triangulation=save_tri;
	do_partial_triangulation=save_part_tri;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_support_hyperplanes_triangulation(){
	keep_triangulation=true;
	do_triangulation=true;   
	build_cone();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::evaluate_triangulation(){

	size_t listsize = Triangulation.size();

	const long VERBOSE_STEPS = 50;
	long step_x_size = listsize-VERBOSE_STEPS;
	if (verbose) {
		verboseOutput() << "evaluating "<<listsize<<" simplices" <<endl;
		verboseOutput() << "---------+---------+---------+---------+---------+"
		                << " (one | per 2%)" << endl;
	}
	

	#pragma omp parallel 
	{
		typename list<pair<vector<size_t>,Integer> >::iterator s = Triangulation.begin();
		size_t spos=0;
		#pragma omp for schedule(dynamic) 
		for(size_t i=0; i<listsize; i++){
			for(; i > spos; ++spos, ++s) ;
			for(; i < spos; --spos, --s) ;

			Simplex<Integer> simp(s->first);
			simp.evaluate(*this,s->second);
			s->second=simp.read_volume();
			if (verbose) {
				#pragma omp critical(VERBOSE)
				while ((long)(i*VERBOSE_STEPS) >= step_x_size) {
					step_x_size += listsize;
					verboseOutput() << "|" <<flush;
				}
			}
		}
	}

	if (verbose) {
		verboseOutput() << endl << listsize<<" simplices evaluated." <<endl;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm_main(){

	if (is_Computed.test(ConeProperty::IsHt1ExtremeRays) && !ht1_extreme_rays) {
		if (do_ht1_elements)
			return;
		if (do_h_vector)
			do_h_vector=false;
	}
		
	if (keep_triangulation) {
		evaluate_triangulation();
	} else {
		Support_Hyperplanes.clear();
		is_Computed.reset(ConeProperty::SupportHyperplanes);
		for(size_t i=0;i<nr_gen;i++)
			in_triang[i]=false;
//		cout << "New build " << endl;
		build_cone();
		extreme_rays_and_ht1_check();
		if(!pointed) return;
	}
//	cout << "Nr Invert " << NrInvert << endl;
	
	if (ht1_extreme_rays && do_triangulation)
		is_Computed.set(ConeProperty::Multiplicity,true);
		
	if (do_Hilbert_basis) {
		global_reduction();
		is_Computed.set(ConeProperty::HilbertBasis,true);
		check_integrally_closed();
	}
	
	if (ht1_extreme_rays && do_Hilbert_basis) {
		select_ht1_elements();
		check_ht1_hilbert_basis();
	}
	if (do_ht1_elements) {
		for(size_t i=0;i<nr_gen;i++)
			if(v_scalar_product(Linear_Form,Generators.read(i+1))==1)
				Ht1_Elements.push_front(Generators.read(i+1));
		Ht1_Elements.sort();
		Ht1_Elements.unique();
		is_Computed.set(ConeProperty::Ht1Elements,true);
	}
	if (do_h_vector) {
		is_Computed.set(ConeProperty::HVector);
		compute_polynomial();
	}
}
//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm_keep_triang() {
	compute_support_hyperplanes_triangulation();
	extreme_rays_and_ht1_check();    
	if(!pointed) return;
	if (ht1_extreme_rays && !ht1_generated) { //TODO ht1_triangulated einbauen und nutzen
		if (verbose) {
			cout << "not all generators have height 1, but extreme rays have"<<endl
			     << "making a new triangulation with only extreme rays" <<endl;
		}
		Support_Hyperplanes.clear();
		is_Computed.set(ConeProperty::SupportHyperplanes,false);
		Triangulation.clear();
		in_triang = vector<bool>(nr_gen,false);
		is_Computed.set(ConeProperty::Triangulation,false);
		compute_support_hyperplanes_triangulation();
	}
		
	primal_algorithm_main();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::primal_algorithm_immediate_evaluation(){
	if (do_triangulation || do_ht1_elements || do_h_vector) {
		check_ht1_generated();
		if(!ht1_generated) {
			compute_support_hyperplanes();
			extreme_rays_and_ht1_check();
			if(!pointed) return;
		}
	}

	primal_algorithm_main();
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
	compute_support_hyperplanes();
	extreme_rays_and_ht1_check();
	reset_tasks();
}

// -v
template<typename Integer>
void Full_Cone<Integer>::support_hyperplanes_triangulation() { 
	primal_algorithm_keep_triang();
	reset_tasks();
}


// -V
template<typename Integer>
void Full_Cone<Integer>::support_hyperplanes_triangulation_pyramid() {   
	do_triangulation=true; 
	primal_algorithm_immediate_evaluation();
	reset_tasks();
}

//-n
template<typename Integer>
void Full_Cone<Integer>::triangulation_hilbert_basis() {
	do_Hilbert_basis=true;
	primal_algorithm_keep_triang();
	reset_tasks();
}

// -N
template<typename Integer>
void Full_Cone<Integer>::hilbert_basis() {
	do_Hilbert_basis=true;
	do_partial_triangulation=true;
	primal_algorithm_immediate_evaluation();
	reset_tasks();
}

// -h
template<typename Integer>
void Full_Cone<Integer>::hilbert_basis_polynomial() {
	do_Hilbert_basis=true;
	do_h_vector=true;
	primal_algorithm_keep_triang();   
	reset_tasks();    
}

// -H
template<typename Integer>
void Full_Cone<Integer>::hilbert_basis_polynomial_pyramid() {
	do_Hilbert_basis=true;
	do_h_vector=true;
	do_triangulation=true;
	primal_algorithm_immediate_evaluation();
	reset_tasks();    
}

// -p
template<typename Integer>
void Full_Cone<Integer>::hilbert_polynomial() {
	do_ht1_elements=true;
	do_h_vector=true;
	primal_algorithm_keep_triang();
	reset_tasks();
}

// -P
template<typename Integer>
void Full_Cone<Integer>::hilbert_polynomial_pyramid() {
	do_ht1_elements=true;
	do_h_vector=true;
	do_triangulation=true;
	primal_algorithm_immediate_evaluation();
	reset_tasks();
}

// -1
template<typename Integer>
void Full_Cone<Integer>::ht1_elements() {
	do_ht1_elements=true;
	do_partial_triangulation=true;
	primal_algorithm_immediate_evaluation();
	reset_tasks();
}

template<typename Integer>
void Full_Cone<Integer>::dual_mode() {
	Support_Hyperplanes.sort();
	Support_Hyperplanes.unique();
	Support_Hyperplanes.remove(vector<Integer>(dim,0));

	if(dim>0) {            //correction needed to include the 0 cone;
		Linear_Form = Generators.homogeneous(ht1_extreme_rays);
		ht1_generated = ht1_extreme_rays;
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
		is_Computed.set(ConeProperty::IsHt1Generated);

		if (ht1_extreme_rays) {
			is_Computed.set(ConeProperty::LinearForm);
			if (verbose) { 
				cout << "Find height 1 elements" << endl;
			}
			typename list < vector <Integer> >::const_iterator h;
			for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
				if (v_scalar_product((*h),Linear_Form)==1) {
					Ht1_Elements.push_back((*h));
				}
			}
			is_Computed.set(ConeProperty::Ht1Elements);
		}
	} else {
		ht1_extreme_rays = ht1_generated = true;
		Linear_Form=vector<Integer>(dim);
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
		is_Computed.set(ConeProperty::LinearForm);
	}
	if (ht1_extreme_rays) check_ht1_hilbert_basis();
	check_integrally_closed();
}

//---------------------------------------------------------------------------
// Checks and auxiliary algorithms
//---------------------------------------------------------------------------


template<typename Integer>
Simplex<Integer> Full_Cone<Integer>::find_start_simplex() const {
	if (is_Computed.test(ConeProperty::ExtremeRays)) {
		vector<size_t> marked_extreme_rays(0);
		for (size_t i=0; i<nr_gen; i++) {
			if (Extreme_Rays[i])
				marked_extreme_rays.push_back(i+1);
		}
		vector<size_t> key_extreme = Generators.submatrix(Extreme_Rays).max_rank_submatrix_lex();
		assert(key_extreme.size() == dim);
		vector<size_t> key(dim);
		for (size_t i=0; i<dim; i++) {
			key[i] = marked_extreme_rays[key_extreme[i]-1];
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
void Full_Cone<Integer>::compute_extreme_rays(){
	size_t i,j,k,l,t;
	Matrix<Integer> SH=getSupportHyperplanes().transpose();
	Matrix<Integer> Val=Generators.multiplication(SH);
	size_t nc=Val.nr_of_columns();
	vector<size_t> Zero(nc);
	vector<size_t> nr_zeroes(nr_gen);

	for (i = 0; i <nr_gen; i++) {
		k=0;
		Extreme_Rays[i]=true;
		for (j = 0; j <nc; j++) {
			if (Val.get_elem(i+1,j+1)==0) {
				k++;
			}
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
			if (Val.get_elem(i+1,j+1)==0) {
				Zero[k]=j;
				k++;
			}
		}

		for (j = 0; j <nr_gen; j++) {
			if (i!=j && Extreme_Rays[j]                // not compare with itself or a known nonextreme ray
					 && nr_zeroes[i]<nr_zeroes[j]) {   // or something whose zeroes cannot be a superset
				l=0;
				for (t = 0; t < nr_zeroes[i]; t++) {
					if (Val.get_elem(j+1,Zero[t]+1)==0)
						l++;
					if (l>=nr_zeroes[i]) {
						Extreme_Rays[i]=false;
						break;
					}
				}
			}
		}
	}

	// cout << "Extr durch" << endl;
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
	if (is_Computed.test(ConeProperty::IsPointed))
		return;
	Matrix<Integer> SH = getSupportHyperplanes();
	pointed = (SH.rank_destructiv() == dim);
	is_Computed.set(ConeProperty::IsPointed);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::check_ht1_generated() {
	if (is_Computed.test(ConeProperty::IsHt1Generated))
		return;

	if (is_Computed.test(ConeProperty::ExtremeRays)) {
		check_ht1_extreme_rays();
		if (ht1_extreme_rays) {
			ht1_generated = true;
			for (size_t i = 0; i < nr_gen; i++) {
				if (!Extreme_Rays[i] && v_scalar_product(Generators[i], Linear_Form) != 1) {
					ht1_generated = false;
					break;
				}
			}
		}
	} else {
		Linear_Form = Generators.homogeneous(ht1_generated);
		if (ht1_generated) {
			ht1_extreme_rays=true;
			is_Computed.set(ConeProperty::IsHt1ExtremeRays);
			is_Computed.set(ConeProperty::LinearForm);
		}
	}
	is_Computed.set(ConeProperty::IsHt1Generated);
}

template<typename Integer>
void Full_Cone<Integer>::check_ht1_extreme_rays() {
	if (is_Computed.test(ConeProperty::IsHt1ExtremeRays))
		return;

	if (ht1_generated) {
		ht1_extreme_rays=true;
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
		return;
	}
	assert(is_Computed.test(ConeProperty::ExtremeRays));
	vector<size_t> key;
	for (size_t i = 0; i < nr_gen; i++) {
		if (Extreme_Rays[i])
			key.push_back(i+1);
	}
	Matrix<Integer> Extreme=Generators.submatrix(key);
	Linear_Form = Extreme.homogeneous(ht1_extreme_rays);
	is_Computed.set(ConeProperty::IsHt1ExtremeRays);
	if (ht1_extreme_rays) {
		is_Computed.set(ConeProperty::LinearForm);
	}
}

template<typename Integer>
void Full_Cone<Integer>::check_ht1_hilbert_basis() {
	if (is_Computed.test(ConeProperty::IsHt1HilbertBasis))
		return;

	if ( !is_Computed.test(ConeProperty::IsHt1ExtremeRays) || !is_Computed.test(ConeProperty::HilbertBasis)) {
		errorOutput() << "Warning: unsatisfied preconditions in check_ht1_hilbert_basis()!" <<endl;
		return;
	}
	
	if (is_Computed.test(ConeProperty::Ht1Elements)) {
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
	if (is_Computed.test(ConeProperty::IsIntegrallyClosed))
		return;

	if ( !is_Computed.test(ConeProperty::HilbertBasis)) {
		errorOutput() << "Warning: unsatisfied preconditions in check_integrally_closed()!" <<endl;
		return;
	}
	integrally_closed = false;
	if (Hilbert_Basis.size() <= nr_gen) {
		integrally_closed = true;
		typename list< vector<Integer> >::iterator h;
		for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
			integrally_closed = false;
			for (size_t i=1; i<= nr_gen; i++) {
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
	register size_t i;
	size_t s=Support_Hyperplanes.size();
	vector <Integer> candidate=v_cut_front(new_element,dim);
	vector <Integer> scalar_product=l_multiplication(Support_Hyperplanes,candidate);
	typename list< vector<Integer>* >::iterator j;
	vector<Integer> *reducer;
	for (j =Irred.begin(); j != Irred.end(); j++) {
		reducer=(*j);
		for (i = 1; i <= s; i++) {
			if ((*reducer)[i]>scalar_product[i-1]){
				break;
			}
		}
		if (i==s+1) {
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
	
	list <vector<Integer> > Candidates_with_Scalar_Product;
	list <vector<Integer> > HB;
	typename list <vector<Integer> >::iterator c;
	typename list <vector<Integer> >::const_iterator h;
	typename list <vector<Integer> >::iterator cit;
	
	for (size_t i = 0; i <nr_gen; i++) {
		Candidates.push_front(Generators.read(i+1));
	}
	Candidates.sort();
	Candidates.unique();

	if (nr_gen == dim) { // cone is simplicial, therefore no global reduction is necessary
		if (verbose) {
			verboseOutput()<<"Cone is simplicial, no global reduction necessary."<<endl;
		}
		for (cit = Candidates.begin(); cit != Candidates.end(); ++cit) {
			Hilbert_Basis.push_back(v_cut_front(*cit,dim));
		}
		Candidates.clear();
		return;
	}
	

	vector<Integer> degree_function=compute_degree_function();

	cit = Candidates.begin();
	size_t cpos = 0;
	size_t listsize=Candidates.size();
	
	if(verbose) {
		verboseOutput()<<"computing the degrees of the candidates... "<<flush;
	}
	//go over candidates: do single scalar product
	//for (c = Candidates.begin(); c != Candidates.end(); c++) { 
	vector<Integer> scalar_product;
	for (size_t j=0; j<listsize; ++j) {
		for(;j > cpos; ++cpos, ++cit) ;
		for(;j < cpos; --cpos, --cit) ;

		norm=v_scalar_product(degree_function,(*cit));

		vector <Integer> new_element(1);
		new_element[0]=norm;
		new_element=v_merge(new_element,(*cit));
		Candidates_with_Scalar_Product.push_back(new_element);
	}
	Candidates.clear();         //delete old set
	if(verbose) {
		verboseOutput()<<"sorting the list... "<<endl<<flush;
	}
	Candidates_with_Scalar_Product.sort();
	if (verbose) {
		verboseOutput()<< Candidates_with_Scalar_Product.size() <<" candidate vectors sorted."<<endl;
	}
	
	// do global reduction
	list< vector<Integer> > HBtmp(0);
	Integer norm_crit;
	while ( !Candidates_with_Scalar_Product.empty() ) {
		//use norm criterion to find irreducible elements
		c=Candidates_with_Scalar_Product.begin();
		norm_crit=(*c)[0]*2;  //candidates with smaller norm are irreducible
		if ( Candidates_with_Scalar_Product.back()[0] < norm_crit) { //all candidates are irreducible
			if (verbose) {
				verboseOutput()<<Hilbert_Basis.size()+Candidates_with_Scalar_Product.size();
				verboseOutput()<<" Hilbert Basis elements of degree <= "<<norm_crit-1<<"; done"<<endl;
			}
			while ( !Candidates_with_Scalar_Product.empty()) {
				Hilbert_Basis.push_back(v_cut_front(*c,dim)); // already of the final type 
				c=Candidates_with_Scalar_Product.erase(c);
			}
			break;
		}
		while ( (*c)[0] < norm_crit ) { //can't go over the end because of the previous if
			// push to HBtmp with scalar products
			vector <Integer> candidate=v_cut_front(*c,dim);
			vector <Integer> scalar_products=l_multiplication(Support_Hyperplanes,candidate);
			vector <Integer> new_HB_element(1);
			new_HB_element[0]=(*c)[0];
			new_HB_element=v_merge(new_HB_element,scalar_products);
			new_HB_element=v_merge(new_HB_element,candidate);
			HBtmp.push_back(new_HB_element);
			Hilbert_Basis.push_back(candidate); // already of the final type 
			c=Candidates_with_Scalar_Product.erase(c);
		}
		size_t csize=Candidates_with_Scalar_Product.size();
		if (verbose) {
			verboseOutput()<<Hilbert_Basis.size()<< " Hilbert Basis elements of degree <= "<<norm_crit-1<<"; "<<csize<<" candidates left"<<endl;
		}

		// reduce candidates against HBtmp
		// fill pointer list
		list < vector <Integer>* >  HBpointers;  // used to put "reducer" to the front
		c=HBtmp.begin();
		while (c!=HBtmp.end()) {
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
		
	//  list< vector<Integer>* > HBcopy(HBpointers); //one copy for each thread

		c=Candidates_with_Scalar_Product.begin();
		cpos=0;
		#pragma omp for schedule(dynamic)
		for (size_t k=0; k<csize; ++k) {
			for(;k > cpos; ++cpos, ++c) ;
			for(;k < cpos; --cpos, --c) ;
			
			if ( is_reducible(HBpointers, *c) ) {
				(*c)[0]=-1; //mark as reducible
			}

			if (verbose) {
				#pragma omp critical(VERBOSE)
				{
				counter++;

				while (counter*VERBOSE_STEPS >= step_x_size) {
					steps_done++;
					step_x_size += csize;
					verboseOutput() << "|" <<flush;
//					cout<<counter<<" ";
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
		c=Candidates_with_Scalar_Product.begin();
		while(c != Candidates_with_Scalar_Product.end()) {
			if((*c)[0]==-1) {
				c=Candidates_with_Scalar_Product.erase(c);
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
		verboseOutput()<<"computing degree function... ";
	}
	size_t i;  
	vector<Integer> degree_function(dim,0);
	if(is_Computed.test(ConeProperty::LinearForm)){ //use Linear_From in homogeneous case
		for (i=0; i<dim; i++) {
			degree_function[i] = Linear_Form[i];
		}
		if(verbose) {
			verboseOutput()<<"using homogenous linear form."<<endl<<flush;
		}
	} else { // add hyperplanes to get a degree function
		typename list< vector<Integer> >::const_iterator h;
		for (h=Support_Hyperplanes.begin(); h!=Support_Hyperplanes.end(); ++h) {
			for (i=0; i<dim; i++) {
				degree_function[i]+=(*h)[i];
			}
		} //TODO parallel addition in each thread and final addition at the end
		//TODO make_prime()?
		if(verbose) {
			verboseOutput()<<"done."<<endl<<flush;
		}
	}
	return degree_function;
}

//--------------------------------------------------------------------------
// Hilbert polynomial
//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Full_Cone<Integer>::compute_e_vector(){
	size_t i,j;
	vector <Integer> E_Vector(dim,0);
	vector <Integer> Q=H_Vector;
	Q.push_back(0);
	for (i = 0; i <dim; i++) {
		for (j = 0; j <dim; j++) {
			E_Vector[i]+=Q[j];
		}
		E_Vector[i]/=permutations<Integer>(1,i);
		for (j = 1; j <=dim; j++) {
			Q[j-1]=(unsigned long)j*Q[j];
		}
	}
	return E_Vector;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::compute_polynomial(){
	size_t i,j;
	Integer factorial=permutations<Integer>(1,dim);
	if ((factorial-permutations_modulo<Integer>(1,dim,overflow_test_modulus))%overflow_test_modulus != 0) {
		errorOutput() << "Hilbert polynom has too big coefficients. Its computation is omitted." <<endl;
		return;
	}
	Integer mult_factor = factorial;
	vector <Integer> E_Vector=compute_e_vector();
	vector <Integer> C(dim,0);
	C[0]=1;
	for (i = 0; i <dim; i++) {
		mult_factor=permutations<Integer>(i,dim);
		if (((dim-1-i)%2)==0) {
			for (j = 0; j <dim; j++) {
				Hilbert_Polynomial[2*j]+=mult_factor*E_Vector[dim-1-i]*C[j];
			}
		}
		else {
			for (j = 0; j <dim; j++) {
				Hilbert_Polynomial[2*j]-=mult_factor*E_Vector[dim-1-i]*C[j];
			}
		}
		for (j = dim-1; 0 <j; j--) {
			C[j]=(unsigned long)(i+1)*C[j]+C[j-1];
		}
		C[0]=permutations<Integer>(1,i+1);
	}
	for (i = 0; i <dim; i++) {
		mult_factor=gcd<Integer>(Hilbert_Polynomial[2*i],factorial);
		Hilbert_Polynomial[2*i]/= mult_factor;
		Hilbert_Polynomial[2*i+1]= factorial/mult_factor;
	}
	is_Computed.set(ConeProperty::HilbertPolynomial);
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Full_Cone<Integer>::primary_multiplicity() const{
	size_t i,j,k;
	Integer primary_multiplicity=0;
	vector <size_t> key,new_key(dim-1);
	Matrix<Integer> Projection(nr_gen,dim-1);
	for (i = 1; i <= nr_gen; i++) {
		for (j = 1; j <= dim-1; j++) {
			Projection.write(i,j,Generators.read(i,j));
		}
	}
	typename list< vector<Integer> >::const_iterator h;
	typename list<pair<vector<size_t>,Integer> >::const_iterator t;
	for (h =Support_Hyperplanes.begin(); h != Support_Hyperplanes.end(); ++h){
		if ((*h)[dim-1]!=0) {
			for (t =Triangulation.begin(); t!=Triangulation.end(); ++t){
				key=t->first;
				for (i = 0; i <dim; i++) {
					k=0;
					for (j = 0; j < dim; j++) {
						if (j!=i && Generators.read(key[j],dim)==1) {
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
							Simplex<Integer> S(new_key,Projection);
							primary_multiplicity+=S.read_volume();
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
Full_Cone<Integer>::Full_Cone(){
	dim=0;
	nr_gen=0;
	hyp_size=0;
}

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
}

//---------------------------------------------------------------------------

template<typename Integer>
Full_Cone<Integer>::Full_Cone(Matrix<Integer> M){
	dim=M.nr_of_columns();
	if (dim!=M.rank()) {
		error_msg("error: Matrix with rank = number of columns needed in the constructor of the object Full_Cone<Integer>.\nProbable reason: Cone not full dimensional (<=> dual cone not pointed)!");
		throw NormalizException();
	}
	Generators = M;
	nr_gen=Generators.nr_of_rows();
	//make the generators coprime and remove 0 rows
	vector<Integer> gcds = Generators.make_prime();
	vector<size_t> key=v_non_zero_pos(gcds);
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
	Support_Hyperplanes = list< vector<Integer> >();
	Triangulation = list<pair<vector<size_t>,Integer> >();
	in_triang = vector<bool> (nr_gen,false);
	HypIndVal = list<FMDATA>();
	Hilbert_Basis = list< vector<Integer> >();
	Candidates = list< vector<Integer> >();
	Ht1_Elements = list< vector<Integer> >();  
	if(dim>0){            //correction needed to include the 0 cone;
		H_Vector = vector<Integer>(dim);
		Hilbert_Polynomial = vector<Integer>(2*dim);
	} else {
		multiplicity = 1;
		H_Vector = vector<Integer>(1,1);
		Hilbert_Polynomial = vector<Integer>(2,1);
		Hilbert_Polynomial[0] = 0;
		is_Computed.set(ConeProperty::HilbertPolynomial);
		is_Computed.set(ConeProperty::HVector);
		is_Computed.set(ConeProperty::Triangulation);
	}
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
	ht1_hilbert_basis = false;
	integrally_closed = false;
	
	reset_tasks();
	
	Extreme_Rays = vector<bool>(nr_gen,true); //all generators are extreme rays
	is_Computed.set(ConeProperty::ExtremeRays);
	Matrix<Integer> SH = C.SupportHyperplanes;
	Support_Hyperplanes = list< vector<Integer> >();
	for (size_t i=1; i<= SH.nr_of_rows(); i++) {
		Support_Hyperplanes.push_back(SH.read(i));
	}
	is_Computed.set(ConeProperty::SupportHyperplanes);
	Triangulation = list<pair<vector<size_t>,Integer> >();
	in_triang = vector<bool>(nr_gen,false);
	HypIndVal = list<FMDATA>();
	Hilbert_Basis = C.Hilbert_Basis;
	is_Computed.set(ConeProperty::HilbertBasis);
	Ht1_Elements = list< vector<Integer> >();
	if(dim>0){            //correction needed to include the 0 cone;
		H_Vector = vector<Integer>(dim);
		Hilbert_Polynomial = vector<Integer>(2*dim);
	} else {
		multiplicity = 1;
		H_Vector = vector<Integer>(1,1);
		Hilbert_Polynomial = vector<Integer>(2,1);
		Hilbert_Polynomial[0] = 0;
		is_Computed.set(ConeProperty::HVector);
		is_Computed.set(ConeProperty::HilbertPolynomial);
	}
}
//---------------------------------------------------------------------------

/* constructor for pyramids */
template<typename Integer>
Full_Cone<Integer>::Full_Cone(const Full_Cone<Integer>& C, Matrix<Integer> M) {
	dim = M.nr_of_columns();
	Generators = M;
	nr_gen = Generators.nr_of_rows();
	multiplicity = 0;
	is_Computed =  bitset<ConeProperty::EnumSize>();
	Extreme_Rays = vector<bool>(nr_gen,false);
	Support_Hyperplanes = list< vector<Integer> >();
	Triangulation = list< pair<vector<size_t>,Integer> >();
	Hilbert_Basis = list< vector<Integer> >();
	Ht1_Elements = list< vector<Integer> >();
	Candidates = list< vector<Integer> >();
	H_Vector = vector<Integer>(dim);
	in_triang = vector<bool> (nr_gen,false);
	HypIndVal = list<FMDATA>();
	
	Linear_Form=C.Linear_Form;
	Order_Vector=C.Order_Vector;
	
	do_triangulation=C.do_triangulation;
	do_partial_triangulation=C.do_partial_triangulation;
	do_ht1_elements=C.do_ht1_elements;
	do_h_vector=C.do_h_vector;
	do_Hilbert_basis=C.do_Hilbert_basis;
	keep_triangulation=C.keep_triangulation;
	is_pyramid=true;
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
Integer Full_Cone<Integer>::getMultiplicity()const{
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
	size_t i=1;
	typename list< vector<Integer> >::const_iterator l;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++) {
		M.write(i,(*l));
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Full_Cone<Integer>::getTriangulation(list< vector<size_t> >& Triang, list<Integer>& TriangVol) const {
	Triang.clear();
	TriangVol.clear();
	vector<size_t> key(dim);
	typename list<pair<vector<size_t>,Integer> >::const_iterator l;
	for (l =Triangulation.begin(); l != Triangulation.end(); l++) {
		key=l->first;
		sort(key.begin(),key.end());
		Triang.push_back(key);
		TriangVol.push_back(l->second);
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Full_Cone<Integer>::getHilbertBasis()const{
	size_t s= Hilbert_Basis.size();
	Matrix<Integer> M(s,dim);
	size_t i=1;
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
	size_t i=1;
	typename list< vector<Integer> >::const_iterator l;
	for (l =Ht1_Elements.begin(); l != Ht1_Elements.end(); l++) {
		M.write(i,(*l));
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Full_Cone<Integer>::getHVector() const{
	return H_Vector;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Full_Cone<Integer>::getHilbertPolynomial() const{
	return Hilbert_Polynomial;
}


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
	verboseOutput()<<"\nh-vector is:\n";
	v_read(H_Vector);
	verboseOutput()<<"\nHilbert polvnomial is:\n";
	v_read(Hilbert_Polynomial);
}

} //end namespace

