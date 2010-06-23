/*
 * Normaliz 2.2
 * Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim
 * With contributions by Christof Soeger
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
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>
#include <time.h>
using namespace std;

//---------------------------------------------------------------------------

#include "full_cone.h"
#include "vector_operations.h"
#include "lineare_transformation.h"
#include "list_operations.h"

//---------------------------------------------------------------------------
extern bool test_arithmetic_overflow;
extern int overflow_test_modulus;
extern bool verbose;
extern bool optimize_speed;
extern void global_error_handling();
struct v_compare_shelling {
	bool operator () (const vector<Integer>& u,const vector<Integer>& v) const 	{
		int dim=u.size()-1;
		Integer a,b;
		a=u[dim]*v[dim-1];
		b=u[dim-1]*v[dim];
		if (a<b)
			return true;
		else if (a>b)
			return false;
		else
			return (v_scalar_multiplication_two(u,v[dim-1])<v_scalar_multiplication_two(v,u[dim-1]));
	}
};


//---------------------------------------------------------------------------
//private
//---------------------------------------------------------------------------

void Full_Cone::add_hyperplane(const int& size, const vector<Integer>& positive,const vector<Integer>& negative){
	int k;
	vector<Integer> hyperplane(hyp_size,0); // initialized with 0
	Integer used_for_tests;
	if (test_arithmetic_overflow==true) {  // does arithmetic tests
		for (k = 0; k <size; k++) {
			hyperplane[k]=positive[size]*negative[k]-negative[size]*positive[k];
			used_for_tests =(positive[size]%overflow_test_modulus)*(negative[k]%overflow_test_modulus)-(negative[size]%overflow_test_modulus)*(positive[k]%overflow_test_modulus);
			if (((hyperplane[k]-used_for_tests) % overflow_test_modulus)!=0) {
				error("error: Arithmetic failure in Full_cone::add_hyperplane. Possible arithmetic overflow.\n");
			}
		}
	}
	else  {                      // no arithmetic tests
		for (k = 0; k <size; k++) {
			hyperplane[k]=positive[size]*negative[k]-negative[size]*positive[k];
		}
	}
	hyperplane=v_make_prime(hyperplane);
	
	#pragma omp critical(HYPERPLANE)
	Support_Hyperplanes.push_back(hyperplane);
}


//---------------------------------------------------------------------------


void Full_Cone::transform_values(const int& size, const vector <int> & test_key){

	//to see if possible to replace the function .end with constant iterator since push-back is performed.

	vector<Integer> hyperplane(hyp_size,0); // initialized with 0
	register int i,j,k,t,nr_zero_i,nr_zero_i_and_j,sub=dim-3;
	
	const bool tv_verbose = false; //verbose && Support_Hyperplanes.size()>10000; //verbose in this method call
		
	// preparing the computations
	list < vector<Integer>* > l_Positive_Simplex,l_Positive_Non_Simplex;
	list < vector<Integer>* > l_Negative_Simplex,l_Negative_Non_Simplex;
	list < vector<Integer>* > l_Neutral_Simplex, l_Neutral_Non_Simplex;
	vector <bool> Zero_Positive(hyp_size,false),Zero_Negative(hyp_size,false);
	list < vector<Integer> > Non_Simplex;
	bool simplex;
	bool rangtest;
	
	if (tv_verbose) cout<<"transform_values: create SZ,Z,PZ,P,NS,N"<<endl<<flush;
	int ipos=0;
	list< vector<Integer> >::iterator ii = Support_Hyperplanes.begin();
	int listsize=Support_Hyperplanes.size();
	//for (ii =Support_Hyperplanes.begin();ii!= Support_Hyperplanes.end();ii++){
//	#pragma omp parallel for private(simplex, nr_zero_i, k) firstprivate(ipos, ii) schedule(dynamic)
	for (int kk=0; kk<listsize; ++kk) {
		for(;kk > ipos; ++ipos, ++ii) ;
		for(;kk < ipos; --ipos, --ii) ;
		simplex=false;
		nr_zero_i=0;
		for (k = dim; k < size; k++) {
			if ((*ii)[k]==0) {
				nr_zero_i++;
				if ((*ii)[size]>0) {
//					#pragma omp atomic
					Zero_Positive[k]=true;
				} else if ((*ii)[size]<0) {
//					#pragma omp atomic
					Zero_Negative[k]=true;
				}
			}
		}
		if (nr_zero_i==dim-1) {
			simplex=true;
		}
		if ((*ii)[size]==0) {
			if (simplex) {
//				#pragma omp critical(NeutS)
				l_Neutral_Simplex.push_back(&(*ii));
			}	else {
//				#pragma omp critical(NeutNS)
				l_Neutral_Non_Simplex.push_back(&(*ii));
			}
		} else 
		if ((*ii)[size]>0) {
			if (simplex) {
//				#pragma omp critical(PosS)
				l_Positive_Simplex.push_back(&(*ii));
			} else {
//				#pragma omp critical(PosNS)
				l_Positive_Non_Simplex.push_back(&(*ii));
			}
		} else 
		if ((*ii)[size]<0) {
			if (simplex) {
//				#pragma omp critical(NegS)
				l_Negative_Simplex.push_back(&(*ii));
			} else {
//				#pragma omp critical(NegNS)
				l_Negative_Non_Simplex.push_back(&(*ii));
			}
		}
	}

	vector <bool> Zero_PN(hyp_size,false);
	for (k = dim; k < size; k++)
		if (Zero_Positive[k]&&Zero_Negative[k])
			Zero_PN[k]=true;

	int PosNSsize = l_Positive_Non_Simplex.size();

	if (tv_verbose) cout<<"transform_values: copy to vector"<<endl;
	vector < vector<Integer>* > Positive_Simplex(l_Positive_Simplex.size());
	vector < vector<Integer>* > Positive_Non_Simplex(PosNSsize);
	vector < vector<Integer>* > Negative_Simplex(l_Negative_Simplex.size());
	vector < vector<Integer>* > Negative_Non_Simplex(l_Negative_Non_Simplex.size());
	vector < vector<Integer>* > Neutral_Simplex(l_Neutral_Simplex.size());
	vector < vector<Integer>* > Neutral_Non_Simplex(l_Neutral_Non_Simplex.size());

	for (k = 0; k < Positive_Simplex.size(); k++) {
		Positive_Simplex[k]=l_Positive_Simplex.front();
		l_Positive_Simplex.pop_front();
	}

	for (k = 0; k < Positive_Non_Simplex.size(); k++) {
		Positive_Non_Simplex[k]=l_Positive_Non_Simplex.front();
		l_Positive_Non_Simplex.pop_front();
	}

	for (k = 0; k < Negative_Simplex.size(); k++) {
		Negative_Simplex[k]=l_Negative_Simplex.front();
		l_Negative_Simplex.pop_front();
	}

	for (k = 0; k < Negative_Non_Simplex.size(); k++) {
		Negative_Non_Simplex[k]=l_Negative_Non_Simplex.front();
		l_Negative_Non_Simplex.pop_front();
	}

	for (k = 0; k < Neutral_Simplex.size(); k++) {
		Neutral_Simplex[k]=l_Neutral_Simplex.front();
		l_Neutral_Simplex.pop_front();
	}

	for (k = 0; k < Neutral_Non_Simplex.size(); k++) {
		Neutral_Non_Simplex[k]=l_Neutral_Non_Simplex.front();
		l_Neutral_Non_Simplex.pop_front();
	}
	if (tv_verbose) cout<<"PS "<<Positive_Simplex.size()<<" P "<<Positive_Non_Simplex.size()<<" NS "<<Negative_Simplex.size()<<" N "<<Negative_Non_Simplex.size()<<" ZS "<<Neutral_Simplex.size()<<" Z "<<Neutral_Non_Simplex.size()<<endl<<flush;
	 
	/*
	   possible improvement using the fact that in the lifted version all
	   hyperplanes hyp[dim-1]!=0 are simplicies???
	 */
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: fill multimap with subfacets of NS"<<endl<<flush;
	multimap < vector< int >, int > Negative_Subfacet_Multi;

	
//	#pragma omp parallel private(i,k,nr_zero_i)
	{
	vector< int > zero_i(nr_gen);
	vector< int > subfacet(dim-2);
//	#pragma omp for schedule(dynamic)
	for (i=0; i<Negative_Simplex.size();i++){
		nr_zero_i=0;
		for (k = dim; k < size; k++) {
			if (Zero_PN[k]&& ((*Negative_Simplex[i])[k]==0)){
				zero_i[nr_zero_i]=k;
				nr_zero_i++;
			}
		}
		if(nr_zero_i>sub){
			for (k = 0; k <dim-2; k++) {
				subfacet[k]=zero_i[k];
			}
//			#pragma omp critical(MULTISET)
			Negative_Subfacet_Multi.insert(pair<vector< int >, int>(subfacet,i));
			if (nr_zero_i==dim-1){
				for (k = dim-2; k >0; k--) {
					subfacet[k-1]=zero_i[k];
//					#pragma omp critical(MULTISET)
					Negative_Subfacet_Multi.insert(pair<vector< int >, int>(subfacet,i));
				}
			}
		}
	}
	}

	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: go over multimap of size "<< Negative_Subfacet_Multi.size() <<endl<<flush;
	multimap < vector< int >, int > ::iterator jj;
	multimap < vector< int >, int > ::iterator del;
	jj =Negative_Subfacet_Multi.begin();
	while (jj!= Negative_Subfacet_Multi.end()) {
		del=jj++;
		if (jj!=Negative_Subfacet_Multi.end() && (*jj).first==(*del).first) {   //delete since is the intersection of two negative simplicies
			Negative_Subfacet_Multi.erase(del);
			del=jj++;
			Negative_Subfacet_Multi.erase(del);
		}
	}
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: singlemap size "<<Negative_Subfacet_Multi.size()<<endl<<flush;
	
	listsize = Negative_Subfacet_Multi.size();
	map < vector< int >, int > Negative_Subfacet;
	#pragma omp parallel private(i, j, k, jj, nr_zero_i)
	{
	vector< int > subfacet(dim-2);
	jj = Negative_Subfacet_Multi.begin();
	int jjpos=0;
	map < vector< int >, int > ::iterator last_inserted=Negative_Subfacet.begin(); // used to speedup insertion into the new map
	bool found;
	#pragma omp for schedule(dynamic)
	for (int j=0; j<listsize; ++j) {
		for(;j > jjpos; ++jjpos, ++jj) ;
		for(;j < jjpos; --jjpos, --jj) ;

		subfacet=(*jj).first;
		found=false; 
		for (i = 0; i <Neutral_Simplex.size(); i++) {
			for (k = 0; k < dim-2; k++)
				if((*Neutral_Simplex[i])[subfacet[k]]!=0)
					break;
			if (k==dim-2) {
				found=true;
				break;
			}
		}
		if (!found) {
			for (i = 0; i <Neutral_Non_Simplex.size(); i++) {
				for (k = 0; k < dim-2; k++)
					if((*Neutral_Non_Simplex[i])[subfacet[k]]!=0)
						break;
				if (k==dim-2) {
					found=true;
					break;
				}
			}
			if(!found) {
				for (i = 0; i <Negative_Non_Simplex.size(); i++) {
					for (k = 0; k < dim-2; k++)
						if((*Negative_Non_Simplex[i])[subfacet[k]]!=0)
								break;
					if (k==dim-2) {
						found=true;
						break;
					}
				}
			}
		}
		if (!found) {
			#pragma omp critical(NEGATIVE_SUBFACET)
			{last_inserted=Negative_Subfacet.insert(last_inserted,*jj);}
		}
	}
	
	#pragma omp single
	if (tv_verbose) cout<<"transform_values: reduced map size "<<Negative_Subfacet.size()<<endl<<flush;
	#pragma omp single nowait
	Negative_Subfacet_Multi.clear();


	#pragma omp single
	if (tv_verbose) cout<<"transform_values: PS vs NS"<<endl<<flush;
	
	vector< int > zero_i(nr_gen);
	map < vector< int >, int > ::iterator jj_map;
	#pragma omp for schedule(dynamic)
	for (i =0; i<Positive_Simplex.size(); i++){ //Positive Simplex vs.Negative Simplex
		nr_zero_i=0;
		for (k = dim; k < size; k++) {
			if (Zero_PN[k] && (*Positive_Simplex[i])[k]==0) {
				zero_i[nr_zero_i]=k; //contains the indices where *i is 0
				nr_zero_i++;
			}
		}
		if (nr_zero_i>sub) {     //sub=dim-3, else can not contain an effective subfacet
			for (k = 0; k <dim-2; k++) {
				subfacet[k]=zero_i[k];
			}
			jj_map=Negative_Subfacet.find(subfacet);
			if (jj_map!=Negative_Subfacet.end()) {
				add_hyperplane(size,*Positive_Simplex[i],*Negative_Simplex[(*jj_map).second]);
				//Negative_Subfacet.erase(jj_map);
				(*jj_map).second = -1;
			}
			if (nr_zero_i==dim-1){
				for (k = dim-2; k >0; k--) {
					subfacet[k-1]=zero_i[k];
					jj_map=Negative_Subfacet.find(subfacet);
					if (jj_map!=Negative_Subfacet.end()) {
						add_hyperplane(size,*Positive_Simplex[i],*Negative_Simplex[(*jj_map).second]);
						//Negative_Subfacet.erase(jj_map);
						(*jj_map).second = -1;
					}
				}
			}
		}
	}

	#pragma omp single
	if (tv_verbose) cout<<"transform_values: NS vs P"<<endl<<flush;
	#pragma omp single
	listsize = Negative_Subfacet.size();

//	for (jj_map = Negative_Subfacet.begin(); jj_map != Negative_Subfacet.end(); ++jj_map) { //Negative_simplex vs. Positive_Non_Simplex
	jj_map = Negative_Subfacet.begin();
	jjpos=0;
	#pragma omp for schedule(dynamic)
	for (int j=0; j<listsize; ++j) {
		for( ; j > jjpos; ++jjpos, ++jj_map) ;
		for( ; j < jjpos; --jjpos, --jj_map) ;

		if ( (*jj_map).second != -1 ) {
			subfacet=(*jj_map).first;
			for (i = 0; i <Positive_Non_Simplex.size(); i++) {
				for (k = 0; k <dim-2; k++)
					if ((*Positive_Non_Simplex[i])[subfacet[k]]!=0)
						break;
				if (k==dim-2) {
					add_hyperplane(size,*Positive_Non_Simplex[i],*Negative_Simplex[(*jj_map).second]);
					break;
				}
			}
		}
	}
	} //END parallel

	
	
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: PS vs N"<<endl<<flush;
	listsize = Positive_Simplex.size();
	#pragma omp parallel private(k,j,t,nr_zero_i,nr_zero_i_and_j)
	{
	vector< int > zero_i(nr_gen);
	#pragma omp for schedule(dynamic)
	for (int i =0; i<listsize; i++){ //Positive Simplex vs.Negative Non Simplex
		nr_zero_i=0;
		for (k = dim; k < size; k++) {
			if (Zero_PN[k] && (*Positive_Simplex[i])[k]==0) {
				zero_i[nr_zero_i]=k; //contains the indices where i is 0
				nr_zero_i++;
			}
		}
		if (nr_zero_i>sub) {
			for (j=0; j<Negative_Non_Simplex.size(); j++){
				nr_zero_i_and_j=0;
				for (k = 0; k < nr_zero_i; k++)
					if ((*Negative_Non_Simplex[j])[zero_i[k]]==0)
						nr_zero_i_and_j++;
				if(nr_zero_i_and_j==dim-2){
					add_hyperplane(size,*Positive_Simplex[i],*Negative_Non_Simplex[j]);
					if (nr_zero_i==dim-2) {
						break;
					}
				}
			}
		}
	}

	#pragma omp single
	{
		rangtest=false;
		if (Positive_Non_Simplex.size()+Negative_Non_Simplex.size()+Neutral_Non_Simplex.size()>dim*dim*dim/6) {
			rangtest=true;
		}
	}

	#pragma omp single
	if (tv_verbose) cout<<"transform_values: P vs N"<<endl<<flush;
	#pragma omp single
	PosNSsize = Positive_Non_Simplex.size();
	
	bool exactly_two;
	vector< int > zero_i_and_j(nr_gen);
	vector<Integer> hp_i, *hp_j, *hp_t; // pointers to current hyperplanes
	
	#pragma omp for schedule(dynamic)
	for (int i =0; i<PosNSsize; i++){ //Positive Non Simplex vs.Negative Non Simplex
		nr_zero_i=0;
		hp_i=(*Positive_Non_Simplex[i]);
		for (k = dim; k < size; k++) {
			if (Zero_PN[k] && (hp_i)[k]==0) {
				zero_i[nr_zero_i]=k; //contains the indices where i is 0
				nr_zero_i++;
			}
		}
		if (nr_zero_i>sub) {
			for (j=0; j<Negative_Non_Simplex.size(); j++){
				hp_j=(Negative_Non_Simplex[j]);
				nr_zero_i_and_j=0;
				for (k = 0; k < nr_zero_i; k++){
					if ((*hp_j)[zero_i[k]]==0){
						zero_i_and_j[nr_zero_i_and_j]=zero_i[k]; //contains the indices where both *i and *j are 0
						nr_zero_i_and_j++;
					}
				}
				if(nr_zero_i_and_j>sub){//intersection of *i and *j may be a subfacet
					exactly_two=true;
					if (rangtest) {
						Matrix Test(nr_zero_i_and_j,dim);
						for (k = 0; k < nr_zero_i_and_j; k++) {
							Test.write(k+1,Generators.read(test_key[zero_i_and_j[k]]));
						}
						if (Test.rank_destructiv()<dim-2) {
							exactly_two=false;
						}
					}
					else{
						for (t=0;t<PosNSsize;t++){
							hp_t=(Positive_Non_Simplex[t]);
							if (t!=i) {
								k=0;
								while((k<nr_zero_i_and_j)&&((*hp_t)[zero_i_and_j[k]]==0))
									k++;
								if (k==nr_zero_i_and_j) {
									exactly_two=false;
									break;
								}
							}
						}
						if (exactly_two) {
							for (t=0;t<Negative_Non_Simplex.size();t++){
								hp_t=(Negative_Non_Simplex[t]);
								if (t!=j) {
									k=0;
									while((k<nr_zero_i_and_j)&&((*hp_t)[zero_i_and_j[k]]==0))
										k++;
									if (k==nr_zero_i_and_j) {
										exactly_two=false;
										break;
									}
								}
							}
							if (exactly_two) {
								for (t=0;t<Neutral_Non_Simplex.size();t++){
									hp_t=(Neutral_Non_Simplex[t]);
									k=0;
									while((k<nr_zero_i_and_j)&&((*hp_t)[zero_i_and_j[k]]==0))
										k++;
									if (k==nr_zero_i_and_j) {
										exactly_two=false;
										break;
									}
								}
							}
						}
					}
					if (exactly_two) {  //intersection of i and j is a subfacet
						add_hyperplane(size,hp_i,*hp_j);
					}
				}
			}
		}
	}
	} //END parallel
	#pragma omp critical(VERBOSE)



	//removing the negative hyperplanes
	if (tv_verbose) cout<<"transform_values: remove negative hyperplanes"<<endl<<flush;
	list< vector<Integer> >::iterator l;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); ){
		if ((*l)[size]<0) {
			l=Support_Hyperplanes.erase(l);
		} else {
			++l;
		}
	}
	if (tv_verbose) cout<<"transform_values: done"<<endl;
}

//---------------------------------------------------------------------------

void Full_Cone::add_simplex(const int& new_generator,const int& size,const vector<int>& col, const vector<int>& col_inv){
	list< vector<Integer> >::const_iterator i=Support_Hyperplanes.begin();
	list< Simplex >::const_iterator j;
	int nr_zero_i, nr_nonzero_i, not_in_i=0, l, k, s, Triangulation_size=Triangulation.size();
	vector<int> key(dim);

	int ipos=0;
	int listsize=Support_Hyperplanes.size();
	//for (i = Support_Hyperplanes.begin(); i != Support_Hyperplanes.end(); i++){
	#pragma omp parallel for private(j,nr_zero_i,nr_nonzero_i,l,k,s) firstprivate(ipos, i, key, not_in_i) schedule(dynamic)
	for (int kk=0; kk<listsize; ++kk) {
		for(;kk > ipos; ++ipos, ++i) ;
		for(;kk < ipos; --ipos, --i) ;

		if ((*i)[size]<0) {
			nr_zero_i=0;
			for (k = dim; k <size; k++) {
				if ((*i)[k]==0) {
					nr_zero_i++;
				}
			}
			if (nr_zero_i==dim-1) { //simplicial
				l=0;
				for (k = dim; k <size; k++) {
					if ((*i)[k]==0) {
						key[l]=col_inv[k-dim]+1;
						l++;
					}
				}
				key[dim-1]=new_generator+1;
				Simplex simp(key);
				#pragma omp critical(TRIANG)
				Triangulation.push_back(simp);
			}
			else {
				j =Triangulation.begin();
				for (s=0; s<Triangulation_size; s++){
					key=(*j).read_key();
					nr_nonzero_i=0;
					k=0;
					do{
						if ((*i)[col[key[k]-1]] !=0) {
							nr_nonzero_i++;
							not_in_i=k;
						}
						k++;
					} while((k<dim)&&(nr_nonzero_i<2));
					if (nr_nonzero_i<=1){
						key[not_in_i]=new_generator+1;
						Simplex simp(key);
						#pragma omp critical(TRIANG)
						Triangulation.push_back(simp);
					}
					j++;
				}
			}
		}
	}
}

//---------------------------------------------------------------------------

bool Full_Cone::is_reducible(list< vector<Integer> >& Ired, const vector< Integer >& new_element){
	register int i;
	int s=Support_Hyperplanes.size();
	vector <Integer> candidate=v_cut_front(new_element,dim);
	vector <Integer> scalar_product=l_multiplication(Support_Hyperplanes,candidate);
	list< vector<Integer> >::iterator j;
	for (j =Ired.begin(); j != Ired.end(); j++) {
		for (i = 1; i <= s; i++) {
			if ((*j)[i]>scalar_product[i-1]){
				break;
			}
		}
		if (i==s+1) {
			//found a "reducer" and move it to the front
			Ired.push_front(*j);
			Ired.erase(j);
			return true;
		}
	}
	return false;
}

//---------------------------------------------------------------------------

bool Full_Cone::is_reducible(list< vector<Integer>* >& Ired, const vector< Integer >& new_element){
	register int i;
	int s=Support_Hyperplanes.size();
	vector <Integer> candidate=v_cut_front(new_element,dim);
	vector <Integer> scalar_product=l_multiplication(Support_Hyperplanes,candidate);
	list< vector<Integer>* >::iterator j;
	vector<Integer> *reducer;
	for (j =Ired.begin(); j != Ired.end(); j++) {
		reducer=(*j);
		for (i = 1; i <= s; i++) {
			if ((*reducer)[i]>scalar_product[i-1]){
				break;
			}
		}
		if (i==s+1) {
			//found a "reducer" and move it to the front
			Ired.push_front(*j);
			Ired.erase(j);
			return true;
		}
	}
	return false;
}

//---------------------------------------------------------------------------

void Full_Cone::reduce_and_insert(const vector< Integer >& new_element){
	if (new_element[0]==0) {
		return; // new_element=0
	}
	else {
		register int i,c=1,s=Support_Hyperplanes.size()+1;
		vector <Integer> candidate(dim);
		for (i = 0; i <dim; i++) {
			candidate[i]=new_element[i+1];
		}
		vector <Integer> scalar_product=l_multiplication(Support_Hyperplanes,candidate);
		list< vector<Integer> >::iterator j;
		for (j =Hilbert_Basis.begin(); j != Hilbert_Basis.end(); j++) {
			if ( new_element[0]<2*(*j)[0] ) {
				break; //new_element is not reducible since the norm is too small, goes into HB
			} else {
				if ((*j)[c]<=scalar_product[c-1]){
					for (i = 1; i < s; i++) {
						if ((*j)[i]>scalar_product[i-1]){
							c=i;
							break;
						}
					}
					if (i==s) {  //found a "reducer" and move it to the front
						Hilbert_Basis.push_front(*j);
						Hilbert_Basis.erase(j);
						return;
					}
					//new_element is not in the Hilbert Basis
				}
			}
		}
		vector <Integer> new_HB_element(1);
		new_HB_element[0]=new_element[0];
		new_HB_element=v_merge(new_HB_element,scalar_product);
		new_HB_element=v_merge(new_HB_element,candidate);
		Hilbert_Basis.push_back(new_HB_element);
	}
}

//---------------------------------------------------------------------------

void Full_Cone::reduce_and_insert_speed(const vector< Integer >& new_element){
	if (new_element[0]==0) {
		return; // new_element=0
	}
	else {
		register int i,c=1,s=Support_Hyperplanes.size()+1;
		list< vector<Integer> >::iterator j;
		for (j =Hilbert_Basis.begin(); j != Hilbert_Basis.end(); j++) {
			if (new_element[0]<2*(*j)[0]) {
				break; //new_element is not reducible (not the sum of 2 or more basis elements)
			}
			else  {
				if ((*j)[c]<=new_element[c]){
					for (i = 1; i < s; i++) {
						if ((*j)[i]>new_element[i]){
							c=i;
							break;
						}
					}
					if (i==s) { 				//new_element is not in the Hilbert Basis
						Hilbert_Basis.push_front(*j); //put j to the front, it is a promising candidate to reduce with
						Hilbert_Basis.erase(j);
						return;
					}
				}
			}
		}
		Hilbert_Basis.push_back(new_element);
	}
}


//---------------------------------------------------------------------------

void Full_Cone::reduce( list< vector< Integer > >& Ired, list< vector< Integer > >& Red, const int& size){
	Ired.sort();
	register int i,c=1;
	vector<Integer> dummy(size+3,0);
	Red.push_front(dummy);
	Red.push_back(dummy);
	list< vector<Integer> >::iterator j;
	list< vector<Integer> >::iterator s;
	for (s = Red.begin(); s != Red.end(); s++) {
		for (j =Ired.begin(); j != Ired.end(); j++) {
			if ((*s)[0]<2*(*j)[0]) {
				break; //element is not reducible;
			}
			else  {

				if ((*j)[c]<=(*s)[c]){
					for (i = 1; i <= size; i++) {
						if ((*j)[i]>(*s)[i]){
							c=i;
							break;
						}
					}
					if (i==size+1) {
						Ired.push_front(*j);
						Ired.erase(j);
						s=Red.erase(s);
						//	if(s!=Red.begin())
						s--;
						break;
					}
				}
				//new_element is not in the Hilbert Basis
			}
		}
	}
	Red.pop_front();
	Red.pop_back();
}

//---------------------------------------------------------------------------

void Full_Cone::reduce_and_insert(const Matrix& New_Elements){
	int i, j=New_Elements.nr_of_rows();
	for (i = 1; i <= j; i++) {
		reduce_and_insert(New_Elements.read(i));
	}
}

//---------------------------------------------------------------------------

void Full_Cone::reduce_and_insert(const list< vector<Integer> >& New_Elements){
	list< vector<Integer> >::const_iterator l;
	for (l =New_Elements.begin(); l != New_Elements.end(); l++) {
		reduce_and_insert(*l);
	}
}

//---------------------------------------------------------------------------

void Full_Cone::find_new_face(){
	if (verbose==true) {
		cout<<"combinatorial evaluation of the shelling ..."<<endl<<flush;
	}
	int i;
	vector<int> facet(dim-1),key;
	list< Simplex >::iterator l;
	list< int > help;
	list< int >::const_iterator m;
	set< vector<int> > Facets;
	pair<set< vector<int> >::iterator, bool> ret;
	for (l =Triangulation.begin(); l!=Triangulation.end(); l++){
		help.clear();
		key=(*l).read_key();
		for (i = 0; i <dim-1; i++) {
			facet[i]=key[i];
		}
		ret = Facets.insert(facet);
		if (!ret.second) { // facet was in the set before
			help.push_back(key[dim-1]);
			Facets.erase(ret.first);
		}
		for (i = dim-1; i >0; i--) {
			facet[i-1]=key[i];
			ret=Facets.insert(facet);
			if (!ret.second) {  // facet was in the set before
				help.push_back(key[i-1]);
				Facets.erase(ret.first);
			}
		}
		i=0;
		vector< int > new_face(help.size());
		for (m = help.begin(); m != help.end(); m++) {
			new_face[i]=(*m);
			i++;
		}
		(*l).write_new_face(new_face);
	}
	if (verbose==true) {
		cout<<"done."<<endl<<flush; 
	}
}

//---------------------------------------------------------------------------
//public
//---------------------------------------------------------------------------

Full_Cone::Full_Cone(){
	dim=0;
	nr_gen=0;
	hyp_size=0;
}

//---------------------------------------------------------------------------

Full_Cone::Full_Cone(Matrix M){
	dim=M.nr_of_columns();
	if (dim!=M.rank()) {
		error("error: Matrix with rank = number of columns needed in the constructor of the object Full_Cone.\nProbable reason: Cone not full dimensional!");	
	}
	Generators = M;
	nr_gen=Generators.nr_of_rows();
	//make the generators coprime and remove 0 rows
	vector<Integer> gcds = Generators.make_prime();
	vector<int> key=v_non_zero_pos(gcds);
	if (key.size() < nr_gen) {
		Generators=Generators.submatrix(key);
		nr_gen=Generators.nr_of_rows();
	}
	hyp_size=dim+nr_gen;
	multiplicity = 0;
	is_Computed =  bitset<ConeProperty::EnumSize>();
	is_pointed = false;
	is_ht1_extreme_rays = false;
	is_ht1_generated = false;
	is_ht1_hilbert_basis = false;
	is_integrally_closed = false;
	Extreme_Rays = vector<bool>(nr_gen,false);
	Support_Hyperplanes = list< vector<Integer> >();
	Triangulation = list< Simplex >();
	Hilbert_Basis = list< vector<Integer> >();
	Homogeneous_Elements = list< vector<Integer> >();
	if(dim>0){            //correction needed to include the 0 cone;
		H_Vector = vector<Integer>(dim);
		Hilbert_Polynomial = vector<Integer>(2*dim);
	} else {
		multiplicity = 1;
		H_Vector = vector<Integer>(1,1);
		Hilbert_Polynomial = vector<Integer>(2,1);
		Hilbert_Polynomial[0] = 0;
	}
}

//---------------------------------------------------------------------------

Full_Cone::Full_Cone(const Cone_Dual_Mode &C) {

	dim = C.dim;
	Generators = C.get_generators();
	nr_gen = Generators.nr_of_rows();

	hyp_size = dim+nr_gen;
	multiplicity = 0;
	is_Computed =  bitset<ConeProperty::EnumSize>();
	is_pointed = true;
	is_Computed.set(ConeProperty::IsPointed);
	is_ht1_extreme_rays = false;
	is_ht1_generated = false;
	is_ht1_hilbert_basis = false;
	is_integrally_closed = false;
	Extreme_Rays = vector<bool>(nr_gen,true); //all generators are extreme rays
	is_Computed.set(ConeProperty::ExtremeRays);
	Matrix SH = C.SupportHyperplanes;
	Support_Hyperplanes = list< vector<Integer> >();
	for (int i=1; i<= SH.nr_of_rows(); i++) {
		Support_Hyperplanes.push_back(SH.read(i));
	}
	is_Computed.set(ConeProperty::SupportHyperplanes);
	Triangulation = list< Simplex >();
	Hilbert_Basis = C.Hilbert_Basis;
	is_Computed.set(ConeProperty::HilbertBasis);
	Homogeneous_Elements = list< vector<Integer> >();
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

Full_Cone::Full_Cone(const Full_Cone& C){
	dim=C.dim;
	nr_gen=C.nr_gen;
	hyp_size=C.hyp_size;
	is_Computed = C.is_Computed;
	is_pointed=C.is_pointed;
	is_ht1_generated=C.is_ht1_generated;
	is_ht1_extreme_rays=C.is_ht1_extreme_rays;
	is_ht1_hilbert_basis = C.is_ht1_hilbert_basis;
	is_integrally_closed = C.is_integrally_closed;
	Linear_Form=C.Linear_Form;
	multiplicity=C.multiplicity;
	Generators=C.Generators;
	Extreme_Rays=C.Extreme_Rays;
	Support_Hyperplanes=C.Support_Hyperplanes;
	Triangulation=C.Triangulation;
	Hilbert_Basis=C.Hilbert_Basis;
	Homogeneous_Elements=C.Homogeneous_Elements;
	H_Vector=C.H_Vector;
	Hilbert_Polynomial=C.Hilbert_Polynomial;
}

//---------------------------------------------------------------------------

/* constructor for recursively generated subcones */
Full_Cone::Full_Cone(Matrix M, int i) {  
	dim = M.nr_of_columns();
	Generators = M;
	nr_gen = Generators.nr_of_rows();
	hyp_size = dim+nr_gen;
	multiplicity = 0;
	is_Computed =  bitset<ConeProperty::EnumSize>();
	Extreme_Rays = vector<bool>(nr_gen,false);
	Support_Hyperplanes = list< vector<Integer> >();
	Triangulation = list< Simplex >();
	Hilbert_Basis = list< vector<Integer> >();
	Homogeneous_Elements = list< vector<Integer> >();
	H_Vector = vector<Integer>(dim);
}

//---------------------------------------------------------------------------

Full_Cone::~Full_Cone(){
	//automatic destructor
}

//---------------------------------------------------------------------------

void Full_Cone::print()const{
	cout<<"\ndim="<<dim<<".\n";
	cout<<"\nnr_gen="<<nr_gen<<".\n";
	cout<<"\nhyp_size="<<hyp_size<<".\n";
	cout<<"\nHomogeneous is "<<is_ht1_generated<<".\n";
	cout<<"\nLinear_Form is:\n";
	v_read(Linear_Form);
	cout<<"\nMultiplicity is "<<multiplicity<<".\n";
	cout<<"\nGenerators are:\n";
	Generators.read();
	cout<<"\nExtreme_rays are:\n";
	v_read(Extreme_Rays);
	cout<<"\nSupport Hyperplanes are:\n";
	l_read(Support_Hyperplanes);
	cout<<"\nTriangulation is:\n";
	l_read(Triangulation);
	cout<<"\nHilbert basis is:\n";
	l_read(Hilbert_Basis);
	cout<<"\nHomogeneous elements are:\n";
	l_read(Homogeneous_Elements);
	cout<<"\nh-vector is:\n";
	v_read(H_Vector);
	cout<<"\nHilbert polvnomial is:\n";
	v_read(Hilbert_Polynomial);
}

//---------------------------------------------------------------------------

bool Full_Cone::isComputed(ConeProperty::Enum prop) const{
	return is_Computed[prop];
}

//---------------------------------------------------------------------------

int Full_Cone::read_dimension()const{
	return dim;
}

//---------------------------------------------------------------------------

int Full_Cone::read_nr_generators()const{
	return nr_gen;
}

//---------------------------------------------------------------------------

bool Full_Cone::read_homogeneous() const{
	return is_ht1_extreme_rays;
}

bool Full_Cone::isHt1HilbertBasis() const{
	return is_ht1_hilbert_basis;
}

bool Full_Cone::isIntegrallyClosed() const{
	return is_integrally_closed;
}

//---------------------------------------------------------------------------

vector<Integer> Full_Cone::read_linear_form() const{
	return Linear_Form;
}

//---------------------------------------------------------------------------

Integer Full_Cone::read_multiplicity()const{
	return multiplicity;
}

//---------------------------------------------------------------------------

Matrix Full_Cone::read_generators()const{
	return Generators;
}

//---------------------------------------------------------------------------

vector<bool> Full_Cone::read_extreme_rays()const{
	return Extreme_Rays;
}

//---------------------------------------------------------------------------

Matrix Full_Cone::read_support_hyperplanes()const{
	int s= Support_Hyperplanes.size();
	Matrix M(s,dim);
	int i=1;
	list< vector<Integer> >::const_iterator l;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++) {
		M.write(i,(*l));
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

Matrix Full_Cone::read_triangulation()const{
	int s= Triangulation.size();
	Matrix M(s,dim);
	vector<int> key(dim);
	int i=1;
	list< Simplex >::const_iterator l;
	for (l =Triangulation.begin(); l != Triangulation.end(); l++) {
		key=(*l).read_key();
		sort(key.begin(),key.end());
		M.write(i,key);
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

Matrix Full_Cone::read_triangulation_volume()const{
	int s= Triangulation.size();
	Matrix M(s,dim+1);
	vector<int> key(dim);
	int k,i=1;
	list< Simplex >::const_iterator l;
	for (l =Triangulation.begin(); l != Triangulation.end(); l++) {
		key=(*l).read_key();
		sort(key.begin(),key.end());
		for (k = 0; k <dim; k++) {
			M.write(i,k+1,key[k]);
		}
		M.write(i,dim+1,(*l).read_volume());
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

Matrix Full_Cone::read_hilbert_basis()const{
	int s= Hilbert_Basis.size();
	Matrix M(s,dim);
	int i=1;
	list< vector<Integer> >::const_iterator l;
	for (l =Hilbert_Basis.begin(); l != Hilbert_Basis.end(); l++) {
		M.write(i,(*l));
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

Matrix Full_Cone::read_homogeneous_elements()const{
	int s= Homogeneous_Elements.size();
	Matrix M(s,dim);
	int i=1;
	list< vector<Integer> >::const_iterator l;
	for (l =Homogeneous_Elements.begin(); l != Homogeneous_Elements.end(); l++) {
		M.write(i,(*l));
		i++;
	}
	return M;
}

//---------------------------------------------------------------------------

vector<Integer> Full_Cone::read_h_vector() const{
	return H_Vector;
}

//---------------------------------------------------------------------------

vector<Integer> Full_Cone::read_hilbert_polynomial() const{
	return Hilbert_Polynomial;
}


//---------------------------------------------------------------------------
// control methods (public)
//---------------------------------------------------------------------------

void Full_Cone::support_hyperplane_common() {
	 check_pointed();
    if(!is_pointed) return;
    compute_extreme_rays();
    check_ht1_extreme_rays();
    check_ht1_generated();
}

void Full_Cone::support_hyperplanes() {
	compute_support_hyperplanes();
	support_hyperplane_common();
}

void Full_Cone::support_hyperplanes_pyramid() {
	compute_support_hyperplanes_pyramid();
	support_hyperplane_common();
}

void Full_Cone::support_hyperplanes_partial_triangulation() {
	compute_support_hyperplanes(true);
	support_hyperplane_common();
	//process_non_compressed(non_compressed); //TODO in the method again (nicht gut)
}

void Full_Cone::support_hyperplanes_triangulation() {
	compute_support_hyperplanes_triangulation();
	support_hyperplane_common();
	if(!is_pointed) return;
	if (is_ht1_extreme_rays && !is_ht1_generated) {
		if (verbose) {
			cout << "not all generators have height 1, but extreme rays have"<<endl
			     << "making a new triangulation with only extreme rays" <<endl;
			v_read(Extreme_Rays);
		}
		Support_Hyperplanes.clear();
		is_Computed.set(ConeProperty::SupportHyperplanes,false);
		Triangulation.clear();
		is_Computed.set(ConeProperty::Triangulation,false);
		compute_support_hyperplanes_triangulation();
	}
	compute_multiplicity();
}

void Full_Cone::support_hyperplanes_triangulation_pyramid() {
	compute_support_hyperplanes_pyramid(true);
	support_hyperplane_common();
   if(!is_pointed) return;
	if (is_ht1_extreme_rays && !is_ht1_generated) {
		if (verbose) {
			cout << "not all generators have height 1, but extreme rays have"<<endl
			     << "making a new triangulation with only extreme rays" <<endl;
			v_read(Extreme_Rays);
		}
		Support_Hyperplanes.clear();
		is_Computed.set(ConeProperty::SupportHyperplanes,false);
		Triangulation.clear();
		is_Computed.set(ConeProperty::Triangulation,false);
		compute_support_hyperplanes_pyramid(true);
	}
	compute_multiplicity();
}

void Full_Cone::triangulation_hilbert_basis(){
	compute_support_hyperplanes_triangulation();
	support_hyperplane_common();
   if(!is_pointed) return;
	if (is_ht1_extreme_rays && !is_ht1_generated) {
		if (verbose) {
			cout << "not all generators have height 1, but extreme rays have"<<endl
			     << "making a new triangulation with only extreme rays" <<endl;
			v_read(Extreme_Rays);
		}
		Support_Hyperplanes.clear();
		is_Computed.set(ConeProperty::SupportHyperplanes,false);
		Triangulation.clear();
		is_Computed.set(ConeProperty::Triangulation,false);
		compute_support_hyperplanes_triangulation();
	}
	compute_hilbert_basis();
	if (is_ht1_extreme_rays) check_ht1_hilbert_basis();
	check_integrally_closed();
}

void Full_Cone::hilbert_basis(){
	support_hyperplanes_partial_triangulation();
	if(!is_pointed) return;
	compute_hilbert_basis();
	if (is_ht1_extreme_rays) check_ht1_hilbert_basis();
	check_integrally_closed();
}

void Full_Cone::ht1_elements(){
	support_hyperplanes_partial_triangulation();
	if(!is_pointed) return;
	compute_ht1_elements();
}

void Full_Cone::hilbert_basis_polynomial(){
	compute_support_hyperplanes();
	check_pointed();
	if(!is_pointed) return;
	compute_extreme_rays();
	
	check_ht1_extreme_rays();
	if ( !is_ht1_extreme_rays ) {
		if (verbose) {
			cout << "************************************************************" << endl;
			cout << "extreme rays not in height 1, using computation type hilbert_basis" << endl;
		}
		Support_Hyperplanes.clear();
		is_Computed.set(ConeProperty::SupportHyperplanes,false);
		hilbert_basis();
	} else {
		check_ht1_generated();
		if(dim>0) {            //correction needed to include the 0 cone;
			triangulation_lift();
			find_new_face();
			compute_hilbert_basis_polynomial();
			compute_polynomial();
		}
	}
	if (is_ht1_extreme_rays) check_ht1_hilbert_basis();
	check_integrally_closed();
}

void Full_Cone::hilbert_polynomial(){
	compute_support_hyperplanes();
	check_pointed();
	if(!is_pointed) return;
	compute_extreme_rays();
	
	check_ht1_extreme_rays();
	if ( !is_ht1_extreme_rays ) {
		if (verbose) {
			cout << "************************************************************" << endl;
			cout << "extreme rays not in height 1, using computation type hilbert_basis" << endl;
		}
		Support_Hyperplanes.clear();
		is_Computed.set(ConeProperty::SupportHyperplanes,false);
		hilbert_basis();
	} else {
		check_ht1_generated();
		if(dim>0) {            //correction needed to include the 0 cone;
			triangulation_lift();
			find_new_face();
			compute_hilbert_polynomial();
			compute_polynomial();
		}
	}
}

void Full_Cone::dual_mode() {
		Support_Hyperplanes.sort();
		Support_Hyperplanes.unique();
		Support_Hyperplanes.remove(vector<Integer>(dim,0));

	if(dim>0) {            //correction needed to include the 0 cone;
		Linear_Form = Generators.homogeneous(is_ht1_extreme_rays);
		is_ht1_generated = is_ht1_extreme_rays;
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
		is_Computed.set(ConeProperty::IsHt1Generated);

		if (is_ht1_extreme_rays) {
			list < vector <Integer> >::const_iterator h;
			for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
				if (v_scalar_product((*h),Linear_Form)==1) {
					Homogeneous_Elements.push_back((*h));
				}
			}
			is_Computed.set(ConeProperty::Ht1Elements);
		}
	} else {
		is_ht1_extreme_rays = is_ht1_generated = true;
		Linear_Form=vector<Integer>(dim);
		is_Computed.set(ConeProperty::IsHt1ExtremeRays);
		is_Computed.set(ConeProperty::IsHt1Generated);
	}
	if (is_ht1_extreme_rays) check_ht1_hilbert_basis();
	check_integrally_closed();
}

//---------------------------------------------------------------------------
// compute methods (private)
//---------------------------------------------------------------------------

void Full_Cone::compute_extreme_rays(){
	int i,j,k,l,t;
	Matrix SH=read_support_hyperplanes();
	SH=SH.transpose();
	Matrix Val=Generators.multiplication(SH);
	int nc=Val.nr_of_columns();
	vector<int> Zero(nc);
	for (i = 0; i <nr_gen; i++) {
		Extreme_Rays[i]=true;
	}
	for (i = 1; i <=nr_gen; i++) {
		k=0;
		for (j = 1; j <=nc; j++) {
			if (Val.get_elem(i,j)==0) {
				Zero[k]=j;
				k++;
			}
		}
		if (k<dim-1||k==nc) {     // not contained in enough facets  or in all (0 as generator)
			Extreme_Rays[i-1]=false;
		}
		else {
			for (j = 1; j <=nr_gen; j++) {
				if (i!=j && Extreme_Rays[j-1]!=false) { // not compare with itself or a known nonextreme ray
					l=0;
					for (t = 0; t < k; t++) {
						if (Val.get_elem(j,Zero[t])==0) {
							l++;
						}
						if (l>=k) {
							Extreme_Rays[i-1]=false;
							break;
						}
					}
				}
			}

		}
	}
	is_Computed.set(ConeProperty::ExtremeRays);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_support_hyperplanes(const bool do_partial_triangulation) {
	if(dim>0){            //correction needed to include the 0 cone;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing support hyperplanes ..."<<endl;
	}
	int i,j;
	//Initialization of the list of support hyperplanes
	vector<Integer> hyperplane(hyp_size,0),L,R; // initialized with 0
	Simplex S;
	if (is_Computed.test(ConeProperty::ExtremeRays)) {
		S = Simplex(Generators.submatrix(Extreme_Rays));
	} else {
		S = Simplex(Generators);
	}
	vector<int> key=S.read_key();
	vector<bool> in_triang(nr_gen,false);
	vector<int> test_key(hyp_size);
	for (i = 0; i < dim; i++) {
		in_triang[key[i]-1]=true;
	}
	bool first_simplex_compressed = do_partial_triangulation;

	Matrix G=S.read_generators();
	G=G.transpose();
	Matrix H=S.read_support_hyperplanes();
	Matrix P=H.multiplication(G);
	for (i = 1; i <=dim; i++) {
		L=H.read(i);
		R=P.read(i);
		for (j = 0; j < dim; j++) {
			hyperplane[j]=L[j];
			hyperplane[j+dim]=R[j];
			test_key[j+dim]=key[j];
			if (first_simplex_compressed && hyperplane[i]>1)
				first_simplex_compressed = false;
		}
		Support_Hyperplanes.push_back(hyperplane);
	}

	
	int size=2*dim;
	//test if the first simplex is compressed
	if(do_partial_triangulation && !first_simplex_compressed){
//		for (i=dim; i<size; i++) {
//			if (hyperplane[i]>1) {
				Triangulation.push_back(key);
//				break;
//			}
//		}
	}

	//computation of support hyperplanes
	Integer scalar_product;
	bool new_generator;
	list< vector<int> > non_compressed;
	for (j = 0; j <= 1; j++) {  //two times, first only extreme rays are considered
		for (i = 0; i < nr_gen; i++) {
			if ((in_triang[i]==false)&&((j==1)||(Extreme_Rays[i]==true))) {
				new_generator=false;
				list< vector<Integer> >::iterator l=Support_Hyperplanes.begin();
				int lpos=0;
				int listsize=Support_Hyperplanes.size();
				//	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++){
				#pragma omp parallel for private(L,scalar_product) firstprivate(lpos,l) schedule(dynamic)
				for (int k=0; k<listsize; k++) {
					for(;k > lpos; lpos++, l++) ;
					for(;k < lpos; lpos--, l--) ;

					L=Generators.read(i+1);
					scalar_product=v_scalar_product(L,(*l));
					if (test_arithmetic_overflow && v_test_scalar_product(L,(*l),scalar_product,overflow_test_modulus)==false) {
						error("error: Arithmetic failure in Full_cone::support_hyperplanes. Possible arithmetic overflow.\n");
					}

					(*l)[size]=scalar_product;
					if (scalar_product<0) {
						new_generator=true;
						if (do_partial_triangulation && scalar_product<-1) { //found non-compressed pyramid
							// make new subcone (gens in hyperplane + new gen)
							vector<int> pyramid;
							pyramid.reserve(size-dim+1);
							pyramid.push_back(i+1);
							for (int g=dim; g<size; g++) {
								if ((*l)[g]==0) {
									pyramid.push_back(test_key[g]);
								}
							}
							if (pyramid.size()==dim) { //simplicial case
								Simplex simp(pyramid);
								#pragma omp critical(TRIANG)
								Triangulation.push_back(simp);
							}
							else {
								#pragma omp critical(NON_COMP)
								non_compressed.push_back(pyramid);
							}
						}
					}
				}
				if (new_generator) {
					in_triang[i]=true;
					test_key[size]=i+1;
					transform_values(size,test_key);
					size++;
				}
				if (verbose==true) {
					cout<<"generator="<< i+1 <<" and "<<Support_Hyperplanes.size()<<" hyperplanes... ";
					cout<<endl;
				}
			}
		}
	}	
	
	l_cut(Support_Hyperplanes,dim);
	if(do_partial_triangulation && non_compressed.size()>0) process_non_compressed(non_compressed);
	} // end if (dim>0)
	is_Computed.set(ConeProperty::SupportHyperplanes);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_support_hyperplanes_pyramid(const bool do_triang) {
	if (verbose) {
		cout << "\n************************************************************\n";
		cout << "computing support hyperplanes ";
		if (do_triang) cout << "and triangulation ";
		cout << "(pyramid)..." << endl;
	}

	int i, j;
	vector<Integer> hyperplane(hyp_size, 0), L, R;
	Simplex S;
	if (is_Computed.test(ConeProperty::ExtremeRays)) {
		S = Simplex(Generators.submatrix(Extreme_Rays));
	} else {
		S = Simplex(Generators);
	}
	vector<int> key = S.read_key();
	if (do_triang) {
		Triangulation.push_back(key);
	}

	vector<bool> in_triang(nr_gen, false);
	vector<int> test_key(hyp_size);
	for (i = 0; i < dim; i++) {
		in_triang[key[i]-1] = true;
	}
	Matrix G = S.read_generators();
	G = G.transpose();
	Matrix H = S.read_support_hyperplanes();
	Matrix P = H.multiplication(G);
	for (i = 1; i <= dim; i++) {
		L = H.read(i);
		R = P.read(i);
		for (j = 0; j < dim; j++) {
			hyperplane[j] = L[j];
			hyperplane[j + dim] = R[j];
			test_key[j + dim] = key[j];
		}
		Support_Hyperplanes.push_back(hyperplane);
	}

	int size = 2 * dim;
	int g;
	bool verbose_bak = verbose;
	verbose = false;
	Integer scalar_product;
	bool new_generator;
	for (j = 0; j <= 1; j++) {  //two times, first only extreme rays are considered
		for (i = 0; i < nr_gen; i++) {
			if ((in_triang[i] == false) && ((j == 1) || (Extreme_Rays[i] == true))) {
				new_generator = false;
				list<vector<Integer> >::iterator l = Support_Hyperplanes.begin();
				int lpos = 0;
				int listsize = Support_Hyperplanes.size();
				//      for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++){
				#pragma omp parallel for private(L,scalar_product,g) firstprivate(lpos,l) schedule(dynamic)
				for (int k = 0; k < listsize; ++k) {
					for (; k > lpos; ++lpos, ++l) ;
					for (; k < lpos; --lpos, --l) ;

					L = Generators.read(i + 1);
					scalar_product = v_scalar_product(L, (*l));
					if (test_arithmetic_overflow && v_test_scalar_product(L,(*l),scalar_product,overflow_test_modulus)==false) {
						error("error: Arithmetic failure in Full_cone::support_hyperplanes. Possible arithmetic overflow.\n");
					}
					(*l)[size] = scalar_product;
					if (scalar_product < 0) {
						new_generator = true;
						vector<int> pyramid;
						pyramid.reserve(size - dim + 1);
						pyramid.push_back(i + 1);
						for (g = dim; g < size; g++) {
							if ((*l)[g] == 0) {
								pyramid.push_back(test_key[g]);
							}
						}

						list<vector<Integer> > subconeSH;
						if (pyramid.size() == dim) {   //do simplical subcones directly
							vector<Integer> hyperplane(hyp_size, 0);
							Simplex S(pyramid, Generators);
							Matrix H = S.read_support_hyperplanes();
							for (int k = 1; k <= dim; k++) {
								subconeSH.push_back(H.read(k));
							}
							if (do_triang) {
								#pragma omp critical(TRIANG_P)
								Triangulation.push_back(S);
							}
						} else {   //compute support hyperplanes of the subcone
							Full_Cone subcone(Generators.submatrix(pyramid), 1);
							if (do_triang) {
								//TODO use test to decide which function is called in recursion
								subcone.compute_support_hyperplanes_triangulation();
								list<Simplex>::const_iterator sub_it = subcone.Triangulation.begin();
								list<Simplex>::const_iterator sub_end = subcone.Triangulation.end();
								vector<int> subkey;
								 // adjust indices and add to Triangulation
								while (sub_it != sub_end) {
									vector<int> simplex(dim, 0);
									subkey = (*sub_it).read_key();
									for (int j = 0; j < dim; j++) {
										simplex[j] = pyramid[subkey[j] - 1];
									}
									#pragma omp critical(TRIANG_P)
									Triangulation.push_back(simplex);
									++sub_it;
								}

							} else {
								subcone.compute_support_hyperplanes();
							}
							subconeSH = subcone.Support_Hyperplanes;
						}

						list<vector<Integer> >::const_iterator sub_it = subconeSH.begin();
						list<vector<Integer> >::const_iterator sub_end = subconeSH.end();
						vector<Integer> scalar_prods = vector<Integer> (nr_gen);
						// add support hyperplanes when needed
						while (sub_it != sub_end) {
							//check if all old generators are in the key or >0
							for (g = dim; g < size; g++) {
								if (in_triang[test_key[g] - 1] && (*l)[g] != 0) {  //old gen && not in key
									scalar_product = v_scalar_product(Generators.read(test_key[g]), (*sub_it));
									if (scalar_product <= 0) {	 //hyperplane is no new support hyperplane
										break;
									}
									scalar_prods[g-dim] = scalar_product;
								}
							}

							if (g == size) {  //compute remaining scalar products
								for (g = dim; g < size; g++) {
									if (in_triang[test_key[g]-1] && (*l)[g]==0) {  //old gen && in key
										scalar_product = v_scalar_product(Generators.read(test_key[g]),(*sub_it));
										scalar_prods[g-dim] = scalar_product;
									}
								}

								scalar_product = v_scalar_product(L, (*sub_it));
								scalar_prods[size - dim] = scalar_product;
								vector<Integer> hyperplane = v_merge((*sub_it), scalar_prods);
								#pragma omp critical(HYPERPLANE_P)
								Support_Hyperplanes.push_back(hyperplane);
							}
							++sub_it;
						}
					}
				}

				if (new_generator) {
					in_triang[i] = true;
					list<vector<Integer> >::iterator l;
					//remove negative hyperplanes
					for (l = Support_Hyperplanes.begin(); l != Support_Hyperplanes.end();) {
						if ((*l)[size] < 0) {
							l = Support_Hyperplanes.erase(l);
						} else
							l++;
					}
					test_key[size] = i + 1;
					size++;
				}

				if (verbose_bak == true) {
					cout << "generator=" << i + 1 << " and " << Support_Hyperplanes.size() << " hyperplanes... "<< endl;
				}
			}

		}

	}

	verbose = verbose_bak;
	l_cut(Support_Hyperplanes, dim);
	is_Computed.set(ConeProperty::SupportHyperplanes);
	is_Computed.set(ConeProperty::Triangulation, do_triang);
}



//---------------------------------------------------------------------------

void Full_Cone::support_hyperplanes_dynamic(){
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing support hyperplanes ..."<<endl;
	}
	int i,j,k;
	bool simplicial;
	set<Integer> test_simplicial;
	Integer scalar_product, scalar_product_small;
	//intialization of the list of support hyperplanes
	vector<Integer> hyperplane(hyp_size,0),L,R; // initialized with 0
	Simplex S(Generators);
	vector<int> key=S.read_key();
	vector<bool> in_triang(nr_gen,false);
	vector <int> test_key(hyp_size);
	for (i = 0; i < dim; i++) {
		 in_triang[key[i]-1]=true;
	}
	Matrix G=S.read_generators();
	G=G.transpose();
	Matrix H=S.read_support_hyperplanes();
	Matrix P=H.multiplication(G);
	for (i = 1; i <=dim; i++) {
		 L=H.read(i);
		 R=P.read(i);
		 for (j = 0; j < dim; j++) {
			 hyperplane[j]=L[j];
			 hyperplane[j+dim]=R[j];
			 test_key[j+dim]=key[j];
		}
		Support_Hyperplanes.push_back(hyperplane);
	}
	//computation of support hyperplanes
	int size=2*dim;
	bool new_generator;
	for (j = 0; j <= 1; j++) {  //two times, first only extreme rays are considered
		for (i = 0; i < nr_gen; i++) {
			if ((in_triang[i]==false)&&((j==1)||(Extreme_Rays[i]==true))) {
				new_generator=false;
				simplicial=true;
				list< vector<Integer> >::iterator l;
				for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++){
					L=Generators.read(i+1);
					scalar_product=v_scalar_product(L,(*l));
					scalar_product_small=scalar_product-L[dim-1]*(*l)[dim-1];
					if ((*l)[dim-1]!=0) {
						if (scalar_product_small % (*l)[dim-1]==0) { 
							test_simplicial.insert(scalar_product_small / (*l)[dim-1]);
						}
					}
					if (test_arithmetic_overflow && v_test_scalar_product(L,(*l),scalar_product,overflow_test_modulus)==false) {
							error("error: Arithmetic failure in Full_cone::support_hyperplanes_dynamic. Possible arithmetic overflow.\n");
					}
					
					(*l)[size]=scalar_product;
					if (scalar_product<0) {
					   new_generator=true;
					}
					if (scalar_product==0) {
					   simplicial=false;
					}
				}
				if (simplicial==false) {
					k=0;
					while(test_simplicial.find(-L[dim-1])!=test_simplicial.end()){
						 L[dim-1]++;
						 k++;
					}
					Generators.write(i+1,dim,L[dim-1]);
					for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++){
						(*l)[size]=(*l)[size]+k* (*l)[dim-1];
					}
				    
				}
				if (new_generator) {
					in_triang[i]=true;
					test_key[size]=i+1;
					transform_values(size,test_key);
					size++;
				}
				if (verbose==true) {
					cout<<"generator="<< i+1 <<" and "<<Support_Hyperplanes.size()<<" hyperplanes..."<<endl;
				}
			}
		}
	}
	
	l_cut(Support_Hyperplanes,dim);
	is_Computed.set(ConeProperty::SupportHyperplanes);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_support_hyperplanes_triangulation(){
	if(dim>0){            //correction needed to include the 0 cone;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing support hyperplanes and triangulation ..."<<endl;
	}
	int i,j;
	//intialization of the lists of support hyperplanes and triangulation
	vector<Integer> hyperplane(hyp_size,0),L,R; // initialized with 0
	Simplex S;
	if (is_Computed.test(ConeProperty::ExtremeRays)) {
		S = Simplex(Generators.submatrix(Extreme_Rays));
	} else {
		S = Simplex(Generators);
	}
	Triangulation.push_back(S);
	vector<int> key=S.read_key();
	vector<bool> in_triang(nr_gen,false);
	vector<int> test_key(hyp_size);
	vector<int> col(nr_gen,0),col_inv(nr_gen,0); //col[i] contains the position in the hyperplane
	//of the scalar product the hyperplane
	//and the i-th generator
	//col_inv is the inverse of col
	for (i = 0; i < dim; i++) {
		in_triang[key[i]-1]=true;
		col[key[i]-1]=dim+i;
		col_inv[i]=key[i]-1;
	}
	Matrix G=S.read_generators();
	G=G.transpose();
	Matrix H=S.read_support_hyperplanes();
	Matrix P=H.multiplication(G);
	for (i = 1; i <=dim; i++) {
		L=H.read(i);
		R=P.read(i);
		for (j = 0; j < dim; j++) {
			hyperplane[j]=L[j];
			hyperplane[j+dim]=R[j];
			test_key[j+dim]=key[j];
		}
		Support_Hyperplanes.push_back(hyperplane);
	}
	//computation of support hyperplanes and triangulation
	int size=2*dim;
	Integer scalar_product;
	bool new_generator;
	for (j = 0; j <= 1; j++) {  //two times, first only extreme rays are considered
		for (i = 0; i < nr_gen; i++) {
			if ((in_triang[i]==false)&&((j==1)||(Extreme_Rays[i]==true))) {
				new_generator=false;
				list< vector<Integer> >::iterator l=Support_Hyperplanes.begin();
				int lpos=0;
				int listsize=Support_Hyperplanes.size();
				//	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); l++){
				#pragma omp parallel for private(L,scalar_product) firstprivate(lpos,l) schedule(dynamic)
				for (int k=0; k<listsize; ++k) {
					for(;k > lpos; ++lpos, ++l) ;
					for(;k < lpos; --lpos, --l) ;

					L=Generators.read(i+1);
					scalar_product=v_scalar_product(L,(*l));
					if (test_arithmetic_overflow && v_test_scalar_product(L,(*l),scalar_product,overflow_test_modulus)==false) {
						error("error: Arithmetic failure in Full_cone::support_hyperplanes_triangulation. Possible arithmetic overflow.\n");
					}
					(*l)[size]=scalar_product;
					if (scalar_product<0) {
						new_generator=true;
					}
				}
				if (new_generator) {
					in_triang[i]=true;
					test_key[size]=i+1;
					add_simplex(i,size,col,col_inv);
					transform_values(size,test_key);
					col[i]=size;
					col_inv[size-dim]=i;
					size++;
				}
				if (verbose==true) {
					cout<<"generator="<< i+1 <<", "<<Support_Hyperplanes.size()<<" hyperplanes..."<<" and "<<Triangulation.size()<<" simplicies"<<endl;
				}
			}
		}
	}
	
	if (verbose) {
		cout<<Support_Hyperplanes.size()<<" hyperplanes and " <<Triangulation.size()<<" simplicies"<<endl;
	}
	l_cut(Support_Hyperplanes,dim);
	} // end if (dim>0)
	is_Computed.set(ConeProperty::SupportHyperplanes);
	is_Computed.set(ConeProperty::Triangulation);
}

//---------------------------------------------------------------------------

void Full_Cone::check_pointed() {
	Matrix SH = read_support_hyperplanes();
	is_pointed = (SH.rank() == dim);
	is_Computed.set(ConeProperty::IsPointed);
}

void Full_Cone::check_ht1_generated() {
	if (is_ht1_extreme_rays) {
		is_ht1_generated = true;
	 	for (int i = 0; i < nr_gen; i++) {
			if (v_scalar_product(Generators.read(i+1), Linear_Form) != 1) {
				is_ht1_generated = false;
				return ;
			}
		}
	} else {
		Linear_Form = Generators.homogeneous(is_ht1_generated);
	}
	is_Computed.set(ConeProperty::IsHt1Generated);

}

void Full_Cone::check_ht1_extreme_rays() {
	if (is_ht1_generated) {
		is_ht1_extreme_rays=true;
		return ;
	}
	vector<int> key;
	for (int i = 0; i < nr_gen; i++) {
		if (Extreme_Rays[i])
			key.push_back(i+1);
	}
	Matrix Extreme=Generators.submatrix(key);
	Linear_Form = Extreme.homogeneous(is_ht1_extreme_rays);
	is_Computed.set(ConeProperty::IsHt1ExtremeRays);
}

void Full_Cone::check_ht1_hilbert_basis() {
	if ( !is_Computed.test(ConeProperty::IsHt1ExtremeRays) || !is_Computed.test(ConeProperty::HilbertBasis)) {
		cerr << "Warning: unsatisfied preconditions in check_ht1_hilbert_basis()!" <<endl;
		return;
	}
	
	if (is_Computed.test(ConeProperty::Ht1Elements)) {
		is_ht1_hilbert_basis = (Homogeneous_Elements.size() == Hilbert_Basis.size());
	} else {
		is_ht1_hilbert_basis = true;
		list< vector<Integer> >::iterator h;
		for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
			if (v_scalar_product((*h),Linear_Form)!=1) {
				is_ht1_hilbert_basis = false;
				break;
			}
		}
	}
	is_Computed.set(ConeProperty::IsHt1HilbertBasis);
}

void Full_Cone::check_integrally_closed() {
	if ( !is_Computed.test(ConeProperty::HilbertBasis)) {
		cerr << "Warning: unsatisfied preconditions in check_integrally_closed()!" <<endl;
		return;
	}
	is_integrally_closed = false;
	if (Hilbert_Basis.size() <= nr_gen) {
		is_integrally_closed = true;
		list< vector<Integer> >::iterator h;
		for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
			is_integrally_closed = false;
			for (int i=1; i<= nr_gen; i++) {
				if ((*h) == Generators.read(i)) {
					is_integrally_closed = true;
					break;
				}
			}
			if (!is_integrally_closed) {
				break;
			}
		}
	}
	is_Computed.set(ConeProperty::IsIntegrallyClosed);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_multiplicity(){
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing multiplicity ..."<<endl;
	}
	multiplicity=0;
	int listsize=Triangulation.size();
	//for (l =Triangulation.begin(); l!=Triangulation.end(); l++) {
	#pragma omp parallel 
	{	
	Integer volume;
	Integer mult=0;
	list< Simplex >::iterator l=Triangulation.begin();
	int lpos=0;
	#pragma omp for schedule(dynamic)
	for (int k=0; k<listsize; ++k) {
		for(;k > lpos; ++lpos, ++l) ;
		for(;k < lpos; --lpos, --l) ;

		Simplex S=(*l);
		S.initialize(Generators);
		volume=S.read_volume();
		(*l).write_volume(volume);
		mult+=volume;
		if (verbose==true && (k+1)%10000==0) {
			cout<<"simplex="<<k+1<<endl;
		}
	}
	#pragma omp critical(MULT)
	multiplicity+=mult;
	} //END parallel
}

//---------------------------------------------------------------------------

void Full_Cone::compute_ht1_elements() {
	if(!is_ht1_extreme_rays)
		return;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing height 1 elements..."<<endl;
	}

	vector<Integer> generator;
	for (int i = 0; i <nr_gen; i++) {
		generator = Generators.read(i+1);
		if (is_ht1_generated || Extreme_Rays[i] || v_scalar_product(generator,Linear_Form)==1)
			Homogeneous_Elements.push_back(generator);
	}

	multiplicity=0;
	int listsize=Triangulation.size();
	//for (l =Triangulation.begin(); l!=Triangulation.end(); l++) {
	#pragma omp parallel
	{
	Integer volume;
	Integer mult=0;
	list <vector<Integer> > HE;
	list <vector<Integer> >::const_iterator h;

	list< Simplex >::iterator l=Triangulation.begin();
	int lpos=0;
	#pragma omp for schedule(dynamic)
	for (int k=0; k<listsize; ++k) {
		for(;k > lpos; ++lpos, ++l) ;
		for(;k < lpos; --lpos, --l) ;

		Simplex S=(*l);
		S.initialize(Generators);
		volume=S.read_volume();
		(*l).write_volume(volume);
		mult+=volume;

		S.ht1_elements(Linear_Form);
		HE=S.read_homogeneous_elements();
		for (h = HE.begin(); h != HE.end(); ++h) {
			#pragma omp critical(HT1ELEMENTS)
			Homogeneous_Elements.push_back((*h));
		}

		if (verbose==true && (k+1)%5000==0) {
			cout<<"simplex="<<k+1<<endl;
		}
	}
	#pragma omp critical(MULT)
	multiplicity+=mult;
	} //END parallel
	Homogeneous_Elements.sort();
	Homogeneous_Elements.unique();
	is_Computed.set(ConeProperty::Ht1Elements);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_hilbert_basis(){
	if(dim>0){            //correction needed to include the 0 cone;
	int counter=0,i;
	Integer volume;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing Hilbert basis ..."<<endl;
	}

	// local hilbert basis
	set <vector<Integer> > Candidates;
	set <vector<Integer> >::iterator cit;
	list <vector<Integer> > HB;
	list <vector<Integer> >::const_iterator h;
	for (i = 1; i <=nr_gen; i++) {
		Candidates.insert(Generators.read(i));
	}
	multiplicity=0;
	list< Simplex >::iterator l=Triangulation.begin();
	int lpos=0;
	int listsize=Triangulation.size();
	//for (l =Triangulation.begin(); l!=Triangulation.end(); l++) {
	#pragma omp parallel for private(volume,HB,h) firstprivate(lpos,l) schedule(dynamic)
	for (int k=0; k<listsize; ++k) {
		for(;k > lpos; ++lpos, ++l) ;
		for(;k < lpos; --lpos, --l) ;

		Simplex S=(*l);
		S.hilbert_basis_interior(Generators);
		volume=S.read_volume();
		(*l).write_volume(volume);
		#pragma omp critical(MULT)
		multiplicity += volume;
		HB=S.acces_hilbert_basis();
		for (h = HB.begin(); h != HB.end(); ++h) {
			#pragma omp critical(CANDI)
			Candidates.insert((*h));
		}
		if (verbose==true) {
			#pragma omp critical(VERBOSE)
			{
				counter++;
				if (counter%1000==0) {
					cout<<"simplex="<<counter<<" and "<< Candidates.size() <<" candidate vectors to be globally reduced."<<endl;
				}
			}
		}
	}
	if (verbose) {
		cout<<"simplex="<<counter<<" and "<< Candidates.size() <<" candidate vectors to be globally reduced."<<endl;
	}

	global_reduction(Candidates);

	if(is_ht1_generated==true){
		for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); ++h) {
			if (v_scalar_product((*h),Linear_Form)==1) {
				Homogeneous_Elements.push_back((*h));
			}
		}
		is_Computed.set(ConeProperty::Ht1Elements);
	}
	} // end if (dim>0)
	is_Computed.set(ConeProperty::HilbertBasis);
}

//---------------------------------------------------------------------------

void Full_Cone::global_reduction(set < vector<Integer> >& Candidates) {
	Integer norm;
	
	list <vector<Integer> > Candidates_with_Scalar_Product;
	list <vector<Integer> > HB;
	list <vector<Integer> >::iterator c;
	list <vector<Integer> >::const_iterator h;
	set <vector<Integer> >::iterator cit;

	if (nr_gen == dim) { // cone is simplicial, therefore no global reduction is necessary
		if (verbose) {
			cout<<"Cone is simplicial, no global reduction necessary."<<endl;
		}
		for (cit = Candidates.begin(); cit != Candidates.end(); ++cit) {
			Hilbert_Basis.push_back(v_cut_front(*cit,dim));
		}
		Candidates.clear();
		return;
	}

	vector<Integer> degree_function=compute_degree_function();

	cit = Candidates.begin();
	int cpos = 0;
	int listsize=Candidates.size();
	
	if(verbose) {
		cout<<"Computing the degrees of the candidates, "<<flush;
	}
	//go over candidates: do single scalar product
	//for (c = Candidates.begin(); c != Candidates.end(); c++) { 
	vector<Integer> scalar_product;
	for (int j=0; j<listsize; ++j) {
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
		cout<<"sorting the list, "<<flush;
	}
	Candidates_with_Scalar_Product.sort();
	if (verbose) {
		cout<< Candidates_with_Scalar_Product.size() <<" candidate vectors sorted."<<endl;
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
				cout<<Hilbert_Basis.size()+Candidates_with_Scalar_Product.size();
				cout<<" Hilbert Basis elements of degree <= "<<norm_crit-1<<", done"<<endl;
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
		int csize=Candidates_with_Scalar_Product.size();
		if (verbose) {
			cout<<Hilbert_Basis.size()<< " Hilbert Basis elements of degree <= "<<norm_crit-1<<"; "<<csize<<" candidates left"<<endl;
		}

		// reduce candidates against HBtmp
		// fill pointer list
		list < vector <Integer>* >  HBpointers;  // used to put "reducer" to the front
		c=HBtmp.begin();
		while (c!=HBtmp.end()) {
			HBpointers.push_back(&(*(c++)));
		}

		#pragma omp parallel private(c,cpos) firstprivate(HBpointers)
		{
		
	//	list< vector<Integer>* > HBcopy(HBpointers); //one copy for each thread

		c=Candidates_with_Scalar_Product.begin();
		cpos=0;
		#pragma omp for schedule(dynamic)
		for (int k=0; k<csize; ++k) {
			for(;k > cpos; ++cpos, ++c) ;
			for(;k < cpos; --cpos, --c) ;
			if (verbose && k%10000==0 && k!=0) {
				cout<<k<<" / "<<csize<<endl<<flush;
			}
			
			if ( is_reducible(HBpointers, *c) ) {
				(*c)[0]=-1;	//mark as reducible
			}
		}
		} //end parallel

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
		cout<<Hilbert_Basis.size()<< " Hilbert Basis elements"<<endl;
	}
}


//---------------------------------------------------------------------------

vector<Integer> Full_Cone::compute_degree_function() const {
	if(verbose) {
		cout<<"computing degree function: ";
	}
	int i;	
	vector<Integer> degree_function(dim,0);
	if(is_ht1_generated==true){ //use Linear_From in homogeneous case
		for (i=0; i<dim; i++) {
			degree_function[i] = Linear_Form[i];
		}
		if(verbose) {
			cout<<" using homogenous linear form."<<endl<<flush;
		}
	} else { // add hyperplanes to get a degree function
		list< vector<Integer> >::const_iterator h;
		for (h=Support_Hyperplanes.begin(); h!=Support_Hyperplanes.end(); ++h) {
			for (i=0; i<dim; i++) {
				degree_function[i]+=(*h)[i];
			}
		} //TODO parallel addition in each thread and final addition at the end
		if(verbose) {
			cout<<"done."<<endl<<flush;
		}
	}
	return degree_function;
}

//---------------------------------------------------------------------------

bool Full_Cone::low_part_simplicial(){
	support_hyperplanes_dynamic();		//change needed for dynamic lifting
	vector<Integer> val;
	//int i,counter;
	list< vector<Integer> >::iterator l;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end();){
		if ((*l)[dim-1]>0) {         // consider just the lower facets
/*			val=Generators.MxV((*l));
			counter=0;
			for (i = 0; i < nr_gen; i++) {
				if (val[i]==0) {
					counter++;
				}
			}
			if (counter!=dim-1) {   // more then dim vertices in one lower facet, the facet is not simplicial
				return false;
			}
*/			++l;
		}
		else {                     //delete the upper facets
			l=Support_Hyperplanes.erase(l);  //only this should remain, other test not needed anymore
		}
	}
	return true;
}

//---------------------------------------------------------------------------

void Full_Cone::line_shelling(){  //try shelling with a line of direction (0,0 ... 0,1)
	if (verbose==true) {
		 cout<<"computing a line shelling ..."<<endl;
	}
  int i,j;
	vector<Integer> Support_Hyperplane_With_Intersection(dim+1);
	list< vector<Integer> >::const_iterator l;
	set < vector<Integer>, v_compare_shelling>::const_iterator k;
	vector<Integer> line(dim,0);
	for (i = 0; i <dim; i++){
	  for (j = 1; j <=nr_gen; j++) {
		 line[i]+=Generators.read(j,i+1);
	  }
	}
	line[dim-1]=0;
	set < vector < Integer >, v_compare_shelling > Shelling ;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); ++l){
	  for (i = 0; i <dim; i++){
		 Support_Hyperplane_With_Intersection[i]=(*l)[i];
	  }
	  Support_Hyperplane_With_Intersection[dim]=v_scalar_product((*l),line);
	  Shelling.insert(Support_Hyperplane_With_Intersection);
	}
	Support_Hyperplanes.clear();
	for (k = Shelling.begin(); k != Shelling.end(); ++k) {
	  Support_Hyperplanes.push_back((*k));
	}
	l_cut(Support_Hyperplanes,dim);
	if (verbose==true) {
	  cout<<"done."<<endl;
	}
}

//---------------------------------------------------------------------------

void Full_Cone::triangulation_lift(){
	if ( ! is_Computed.test(ConeProperty::SupportHyperplanes) ) {
		error("error: Status support hyperplanes needed before using Full_Cone::triangulation_lift.");
		return;
	}
	int i,j,counter=0,nr_extreme_rays=0;
	for (i = 0; i < nr_gen; i++) {
		if (Extreme_Rays[i]==true) {
			nr_extreme_rays++;
		}
	}
	vector<int> Extreme(nr_extreme_rays);
	for (i = 0; i < nr_gen; i++) {
		if (Extreme_Rays[i]==true) {
			Extreme[counter]=i+1;
			counter++;
		}
	}
	Matrix Extreme_Generators=Generators.submatrix(Extreme);
	if (Extreme_Generators.nr_of_columns()==Extreme_Generators.nr_of_rows()) {
		Simplex S(Extreme);
		Triangulation.push_back(S); 
	}
	else {
		Full_Cone Lifted;
		lift(Lifted,Extreme_Generators);
		if ( !Lifted.isComputed(ConeProperty::SupportHyperplanes) ) {
			error("error: Status support hyperplanes for the lifted cone needed when using Full_Cone::triangulation_lift.");
			return;
		}
		Lifted.line_shelling();
		if (verbose==true) {
			cout<<"computing triangulation ..."<<endl;
		}
		vector<int> key(dim);
		vector<Integer> v;
		list < vector<Integer> >::const_iterator l;
		for (l = Lifted.Support_Hyperplanes.begin(); l != Lifted.Support_Hyperplanes.end(); ++l) {
			v=Lifted.Generators.MxV((*l));
			counter=0;
			for (j = 0; j < nr_extreme_rays; j++) {
				if (v[j+1]==0) {  //j+1 because of the extra unit vector
					key[counter]=Extreme[j];
					counter++;
				}
			}
			Simplex S(key);
			Triangulation.push_back(S);
		}
	}
	if (verbose==true) {
		cout<<"computed triangulation has "<<Triangulation.size()<<" simplices"<<endl;
	}
	is_Computed.set(ConeProperty::Triangulation);
}

//---------------------------------------------------------------------------

vector<Integer> Full_Cone::compute_e_vector(){
	int i,j;
	vector <Integer> E_Vector(dim,0);
	vector <Integer> Q=H_Vector;
	Q.push_back(0);
	for (i = 0; i <dim; i++) {
		for (j = 0; j <dim; j++) {
			E_Vector[i]+=Q[j];
		}
		E_Vector[i]/=permutations(1,i);
		for (j = 1; j <=dim; j++) {
			Q[j-1]=j*Q[j];
		}
	}
	return E_Vector;
}

//---------------------------------------------------------------------------

void Full_Cone::compute_polynomial(){
#ifndef normbig
	if (dim > 21) {
		cerr << "Hilbert polynom has too big coefficients. Its computation is omitted." <<endl;
		return;
	}
#endif
	int i,j;
	Integer mult_factor, factorial=permutations(1,dim);
	vector <Integer> E_Vector=compute_e_vector();
	vector <Integer> C(dim,0);
	C[0]=1;
	for (i = 0; i <dim; i++) {
		mult_factor=permutations(i,dim);
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
			C[j]=(i+1)*C[j]+C[j-1];
		}
		C[0]=permutations(1,i+1);
	}
	for (i = 0; i <dim; i++) {
		mult_factor=gcd(Hilbert_Polynomial[2*i],factorial);
		Hilbert_Polynomial[2*i]/= mult_factor;
		Hilbert_Polynomial[2*i+1]= factorial/mult_factor;
	}
	is_Computed.set(ConeProperty::HilbertPolynomial);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_hilbert_polynomial(){
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing Hilbert polynomial ..."<<endl;
	}
	int counter=0;
	Integer volume;
	multiplicity=0;
	for (int i = 1; i <=nr_gen; i++) {
		Homogeneous_Elements.push_back(Generators.read(i));
	}
	list< vector<Integer> > HE;
	list< Simplex >::iterator l=Triangulation.begin();
	int lpos=0;
	int listsize=Triangulation.size();
	//for (l =Triangulation.begin(); l!=Triangulation.end(); l++) {
	#pragma omp parallel for private(volume,HE) firstprivate(lpos,l) schedule(dynamic)
	for (int k=0; k<listsize; ++k) {
		for(;k > lpos; ++lpos, ++l) ;
		for(;k < lpos; --lpos, --l) ;

		Simplex S=(*l);
		#pragma omp critical(H_VECTOR)
		H_Vector[S.read_new_face_size()]++;
		S.h_vector(Generators,Linear_Form);
		volume=S.read_volume();
		(*l).write_volume(volume);
		#pragma omp critical(MULT)
		multiplicity += volume;
		HE=S.read_homogeneous_elements();
		#pragma omp critical(HOMOGENEOUS)
		Homogeneous_Elements.merge(HE);
		#pragma omp critical(H_VECTOR)
		H_Vector=v_add(H_Vector,S.read_h_vector());
		if (verbose==true) {
			#pragma omp critical(VERBOSE)
			{
				counter++;
				if (counter%1000==0) {
					cout<<"simplex="<<counter<<endl;
				}
			}
		}
	}
	Homogeneous_Elements.sort();
	Homogeneous_Elements.unique();
	is_Computed.set(ConeProperty::Ht1Elements);
	is_Computed.set(ConeProperty::HVector);
}

//---------------------------------------------------------------------------

void Full_Cone::compute_hilbert_basis_polynomial(){
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing Hilbert basis and polynomial ..."<<endl<<flush;
	}
	int counter=0,i;
	Integer volume;
	vector<Integer> scalar_product;
	set <vector<Integer> > Candidates;
	list <vector<Integer> >  HB,HE;
	set <vector<Integer> >::iterator c;
	list <vector<Integer> >::const_iterator h;
	for (i = 1; i <=nr_gen; i++) {
		vector<Integer> Generator = Generators.read(i);
		Candidates.insert(Generator);
		if (is_ht1_generated || v_scalar_product(Generator, Linear_Form)==1) { 
			Homogeneous_Elements.push_back(Generator);
		}
	}
	multiplicity=0;
	list< Simplex >::iterator l=Triangulation.begin();
	int lpos=0;
	int listsize=Triangulation.size();
	//for (l =Triangulation.begin(); l!=Triangulation.end(); l++) {
	#pragma omp parallel for private(volume,HE,HB) firstprivate(lpos,l) schedule(dynamic)
	for (int k=0; k<listsize; ++k) {
		for(;k > lpos; ++lpos, ++l) ;
		for(;k < lpos; --lpos, --l) ;
		
		Simplex S=(*l);
		#pragma omp critical(H_VECTOR)
		H_Vector[S.read_new_face_size()]++;
		S.hilbert_basis_interior_h_vector(Generators,Linear_Form);
		volume=S.read_volume();
		(*l).write_volume(volume);
		#pragma omp critical(MULT)
		multiplicity += volume;
		HE=S.read_homogeneous_elements();
		#pragma omp critical(HOMOGENEOUS)
		Homogeneous_Elements.merge(HE);
		#pragma omp critical(H_VECTOR)
		H_Vector=v_add(H_Vector,S.read_h_vector());
		HB=S.acces_hilbert_basis();
		#pragma omp critical(CANDI)
		Candidates.insert(HB.begin(),HB.end());
		
		if (verbose==true) {
			#pragma omp critical(VERBOSE)
			{
				counter++;
				if (counter%500==0) {
					cout<<"simplex="<<counter<<" and "<< Candidates.size() <<" candidate vectors to be globally reduced."<<endl;
				}
			}
		}
	}
		
	Homogeneous_Elements.sort();
	Homogeneous_Elements.unique();
	is_Computed.set(ConeProperty::Ht1Elements);
	is_Computed.set(ConeProperty::HilbertBasis);
	is_Computed.set(ConeProperty::HVector);
		
	global_reduction(Candidates);
}

//---------------------------------------------------------------------------

Integer Full_Cone::primary_multiplicity() const{
	int i,j,k;
	Integer primary_multiplicity=0;
	vector <int> key,new_key(dim-1);
	Matrix Projection(nr_gen,dim-1);
	for (i = 1; i <= nr_gen; i++) {
		for (j = 1; j <= dim-1; j++) {
			Projection.write(i,j,Generators.read(i,j));
		}
	}
	list< vector<Integer> >::const_iterator h;
	list< Simplex >::const_iterator t;
	for (h =Support_Hyperplanes.begin(); h != Support_Hyperplanes.end(); ++h){
		if ((*h)[dim-1]!=0) {
			for (t =Triangulation.begin(); t!=Triangulation.end(); ++t){
				key=(*t).read_key();
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
							Simplex S(new_key,Projection);
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

void Full_Cone::error(string s) const{
	cerr <<"\nFull Cone "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

void lift(Full_Cone& Lifted,Matrix Extreme_Generators){
	if (verbose) {
		cout<<"start lifting the cone ...";
	}
	int i,j,nr_extreme_gen=Extreme_Generators.nr_of_rows(),dim=Extreme_Generators.nr_of_columns();
	// add an extra dimension, and the (0,...,0,1) vector
	Matrix New_Generators(nr_extreme_gen+1,dim+1);
	New_Generators.write(1, dim+1, 1);  // (0,...,0,1)
	for (i = 1; i <= nr_extreme_gen; i++) {
		for (j = 1; j <= dim; j++) {
			New_Generators.write(i+1,j,Extreme_Generators.read(i,j));
		}
	}

	Lifted = Full_Cone(New_Generators,1);
	if (Lifted.low_part_simplicial()==true) {
		if (verbose) {
			cout<<"lifting done."<<endl;
		}
		return;
	}

	cerr<<"error: Dynamic Lifting has failed in Full_Cone::lift.";
	global_error_handling();
}

//---------------------------------------------------------------------------

bool Full_Cone::check_compressed() {
	//TODO check for "status"
	// test if the scalar products of the support hyperplanes with the generators are always <=1
	int i;
	bool is_compressed=true;
	list< vector<Integer> >::const_iterator ss = Support_Hyperplanes.begin();
	vector<Integer> scalarProduct;

	for (; ss!=Support_Hyperplanes.end()&&is_compressed; ++ss) {
		scalarProduct = Generators.MxV(*ss);
		for (i=0; i<nr_gen && is_compressed; i++) {
			if (scalarProduct[i]>1) {
				cout<<"check_compressed():" << scalarProduct[i] << endl;
				is_compressed=false;
			}
		}
	}
	return is_compressed;
}

//---------------------------------------------------------------------------

void Full_Cone::process_non_compressed(list< vector<int> > &non_compressed) {
	int listsize=non_compressed.size();
	if (verbose) {
		cout<<"\n************************************************************\n";
		cout<<"Progress "<<listsize<<" subcones"<<endl;
	}
	int itpos=0;
	list< vector<int> >::const_iterator it=non_compressed.begin();

	//override global verbose value for recursion
	bool verbose_bak = verbose;
	verbose=false;

	int verbose_step=10000;
	if (verbose_bak) {
		while (listsize/5 < verbose_step && verbose_step >= 100) {
			verbose_step/=10;
		}
	}

	#pragma omp parallel for firstprivate(itpos,it) schedule(dynamic)
	for (int i=0; i<listsize; ++i) {
		for(;i > itpos; ++itpos, ++it) ;
		for(;i < itpos; --itpos, --it) ;
//		v_read(*it); cout<<" subcone "<<i+1;
		//compute triangulations of subcones
		Full_Cone subcone(Generators.submatrix(*it),1);
		subcone.compute_support_hyperplanes_triangulation();	// use normal method to compute triangulation of subcone
//		subcone.support_hyperplanes(true);				// use "compressed test" method to compute non-compressed subcones of subcone
		list< Simplex >::const_iterator sub_it  = subcone.Triangulation.begin();
		list< Simplex >::const_iterator sub_end = subcone.Triangulation.end();
		vector<int> key;
		// adjust indices and add to Triangulation
		while (sub_it!=sub_end) {
			vector<int> simplex(dim,0);
			key=(*sub_it).read_key();
			for (int j=0; j<dim; j++) {
				simplex[j]=(*it)[key[j]-1];
			}
			#pragma omp critical(TRIANG)
			Triangulation.push_back(simplex);
	//		v_read(*sub_it); cout<<" simplex of subcone "<<i+1;
	//		v_read(simplex); cout<<" corresponding simplex in fullcone";
			++sub_it;
		}
		if (verbose_bak && (i+1)%verbose_step==0) {
			cout<<i+1<<" subcones done, found "<<Triangulation.size()<<" simplices"<<endl<<flush;
		}
	}
	if (verbose_bak) {
		cout<<listsize<<" subcones done, found "<<Triangulation.size()<<" simplices"<<endl<<flush;
	}

	//restore global verbose value
	verbose=verbose_bak;
}

