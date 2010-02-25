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
extern int lifting_bound;
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

	//to see if posible to replace the funtion .end with constant iterator since push-back is performed.

	vector<Integer> hyperplane(hyp_size,0); // initialized with 0
	register int i,j,k,t,nr_zero_i,nr_zero_i_and_j,sub=dim-3;
	
	bool tv_verbose = verbose && Support_Hyperplanes.size()>10000;  //verbose in this method call
		
	// preparing the computations
	list < vector<Integer>* > l_Positive_Simplex,l_Positive_Non_Simplex;
	list < vector<Integer>* > l_Negative_Simplex,l_Negative_Non_Simplex;
	list < vector<Integer>* > l_Neutral_Simplex, l_Neutral_Non_Simplex;
	vector <bool> Zero_Positive(hyp_size,false),Zero_Negative(hyp_size,false);
	list < vector<Integer> > Non_Simplex;
	bool simplex;

//	if (tv_verbose) cout<<"transform_values: create SZ,Z,PZ,P,NS,N"<<endl;
	int ipos=0;
	list< vector<Integer> >::iterator ii = Support_Hyperplanes.begin();
	int listsize=Support_Hyperplanes.size();
	//for (ii =Support_Hyperplanes.begin();ii!= Support_Hyperplanes.end();ii++){
//	#pragma omp parallel for private(simplex, nr_zero_i, k) firstprivate(ipos, ii) schedule(dynamic)
	for (int kk=0; kk<listsize; kk++) {
		for(;kk > ipos; ipos++, ii++) ;
		for(;kk < ipos; ipos--, ii--) ;
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

//	if (tv_verbose) cout<<"transform_values: copy to vector"<<endl;
	vector < vector<Integer>* > Positive_Simplex(l_Positive_Simplex.size());
	vector < vector<Integer>* > Positive_Non_Simplex(l_Positive_Non_Simplex.size());
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
	if (tv_verbose) cout<<"PS "<<Positive_Simplex.size()<<" P "<<Positive_Non_Simplex.size()<<" NS "<<Negative_Simplex.size()<<" N "<<Negative_Non_Simplex.size()<<" ZS "<<Neutral_Simplex.size()<<" Z "<<Neutral_Non_Simplex.size()<<endl;
	 
	/*
	   possible improvement using the fact that in the lifted version all
	   hyperplanes hyp[dim-1]!=0 are simplicies???
	 */
//	#pragma omp parallel sections private(i,j,k,nr_zero_i) 
	{
//	#pragma omp section
	{
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: fill multimap with subfacets of NS"<<endl<<flush;
	vector< int > zero_i(nr_gen);
	vector< int > subfacet(dim-2);
	multimap < vector< int >, int > Negative_Subfacet_Multi;

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
			Negative_Subfacet_Multi.insert(pair<vector< int >, int>(subfacet,i));
			if (nr_zero_i==dim-1){
				for (k = dim-2; k >0; k--) {
					subfacet[k-1]=zero_i[k];
					Negative_Subfacet_Multi.insert(pair<vector< int >, int>(subfacet,i));
				}
			}
		}
	}


	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: go over multimap of size "<< Negative_Subfacet_Multi.size() <<endl;
	multimap < vector< int >, int > ::iterator jj;
	multimap < vector< int >, int > ::iterator del;
	jj =Negative_Subfacet_Multi.begin();
	while (jj!= Negative_Subfacet_Multi.end()) {
		del=jj;
		jj++;
		if (jj!=Negative_Subfacet_Multi.end() && (*jj).first==(*del).first) {   //delete since is the intersection of two negative simplicies
			Negative_Subfacet_Multi.erase(del);
			del=jj;
			jj++;
			Negative_Subfacet_Multi.erase(del);
		}
	}
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: singlemap size "<<Negative_Subfacet_Multi.size()<<endl;
	
	int Negative_Subfacet_Multi_Size=Negative_Subfacet_Multi.size();
	map < vector< int >, int > Negative_Subfacet;
	#pragma omp parallel private(i,j,k,jj)
	{
	vector< int > subfacet(dim-2);
	jj = Negative_Subfacet_Multi.begin();
	int jjpos=0;
	map < vector< int >, int > ::iterator last_inserted=Negative_Subfacet.begin(); // used to speedup insertion into the new map
	bool found;
	#pragma omp for schedule(dynamic)
//	while (jj!= Negative_Subfacet_Multi.end()) {
	for (int j=0; j<Negative_Subfacet_Multi_Size; j++) {
		for(;j > jjpos; jjpos++, jj++) ;
		for(;j < jjpos; jjpos--, jj--) ;

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
//		jj++;
	}
	}//END parallel
	
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: reduced map size "<<Negative_Subfacet.size()<<endl;
	Negative_Subfacet_Multi.clear();
	//making computations

	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: PS vs NS"<<endl;
	
	map < vector< int >, int > ::iterator jj_map;
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
				Negative_Subfacet.erase(jj_map);
			}
			if (nr_zero_i==dim-1){
				for (k = dim-2; k >0; k--) {
					subfacet[k-1]=zero_i[k];
					jj_map=Negative_Subfacet.find(subfacet);
					if (jj_map!=Negative_Subfacet.end()) {
						add_hyperplane(size,*Positive_Simplex[i],*Negative_Simplex[(*jj_map).second]);
						Negative_Subfacet.erase(jj_map);
					}
				}
			}
		}
	}

	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: NS vs P"<<endl;
	for (jj_map = Negative_Subfacet.begin();jj_map != Negative_Subfacet.end() ; jj_map++) { //Negative_simplex vs. Positive_Non_Simplex
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
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: multimap and NS vs ... done"<<endl;
	} //END section

//	#pragma omp section
	{
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: PS vs N"<<endl;
	int PosSsize=Positive_Simplex.size();
	#pragma omp parallel private(k,j,nr_zero_i,nr_zero_i_and_j)
	{
	vector< int > zero_i(nr_gen);
	#pragma omp for schedule(dynamic)
	for (int i =0; i<PosSsize; i++){ //Positive Simplex vs.Negative Non Simplex
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
	} //END parallel
	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: PS vs N done"<<endl;
	} //END section

	
//	#pragma omp section
	{
	bool rangtest=false;
	if (Positive_Non_Simplex.size()+Negative_Non_Simplex.size()+Neutral_Non_Simplex.size()>dim*dim*dim/6) {
		rangtest=true;
	}

	#pragma omp critical(VERBOSE)
	if (tv_verbose) cout<<"transform_values: P vs N"<<endl;
	int PosNSsize=Positive_Non_Simplex.size();
	#pragma omp parallel private(k,j,t,nr_zero_i,nr_zero_i_and_j)
	{
	bool exactly_two;
	vector< int > zero_i(nr_gen);
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
	if (tv_verbose) cout<<"transform_values: P vs N done"<<endl;
	} //END section	
	} //END sections
	
	//removing the negative hyperplanes
	if (tv_verbose) cout<<"transform_values: remove negative hyperplanes"<<endl;
	list< vector<Integer> >::iterator l;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end(); ){
		if ((*l)[size]<0) {
			l=Support_Hyperplanes.erase(l);
		}
		else
			l++;
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
	for (int kk=0; kk<listsize; kk++) {
		for(;kk > ipos; ipos++, i++) ;
		for(;kk < ipos; ipos--, i--) ;

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
						Hilbert_Basis.push_front(*j); //put j to the front, it is a promissing candidate to reduce with
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

bool Full_Cone::reduce( list< vector< Integer >* >& Ired, const vector< Integer >& new_element, const int& size){
	register int i,c=1;
	list< vector<Integer>* >::iterator j;
	vector<Integer> *reducer;
	for (j =Ired.begin(); j != Ired.end(); j++) {
		reducer=(*j);
		if (new_element[0]<=(*reducer)[0])
			continue;
		if ((*reducer)[c]<=new_element[c]){
			for (i = 1; i <=size ; i++) {
				if ((*reducer)[i]>new_element[i]){
					c=i;
					break;
				}
			}
			if (i==size+1) {
				Ired.push_front(*j);
				Ired.erase(j);
				return true;
			}
		}
		//new_element is reducible
	}
	return false;
}

//---------------------------------------------------------------------------

bool Full_Cone::reduce( list< vector< Integer > >& Ired, const vector< Integer >& new_element, const int& size){
	register int i,c=1;
	list< vector<Integer> >::iterator j;
	for (j =Ired.begin(); j != Ired.end(); j++) {
		if (new_element[0]<=(*j)[0])
			continue;
		if ((*j)[c]<=new_element[c]){
			for (i = 1; i <=size ; i++) {
				if ((*j)[i]>new_element[i]){
					c=i;
					break;
				}
			}
			if (i==size+1) {
				Ired.push_front(*j);
				Ired.erase(j);
				return true;
			}
		}
		//new_element is reducible
	}
	return false;
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

void Full_Cone::reduce_and_insert(const vector< Integer >& new_element, const int& size){
	register int i,c=1;
	list< vector<Integer> >::iterator j;
	for (j =Hilbert_Basis.begin(); j != Hilbert_Basis.end(); j++) {
		if (new_element[0]<2*(*j)[0]) {
			break; //new_element is not reducible;
		}
		else  {

			if ((*j)[c]<=new_element[c]){
				for (i = 1; i <= size; i++) {
					if ((*j)[i]>new_element[i]){
						c=i;
						break;
					}
				}
				if (i==size+1) {
					Hilbert_Basis.push_front(*j);
					Hilbert_Basis.erase(j);
					return;
				}
			}
			//new_element is not in the Hilbert Basis
		}
	}
	Hilbert_Basis.push_back(new_element);
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


void Full_Cone::reduce_and_insert_extreme( const vector< Integer >& new_element){
	register int i,c=1;
	list< vector<Integer> >::iterator j;
	for (j =Support_Hyperplanes.begin(); j != Support_Hyperplanes.end(); j++) {
		if (new_element[0]<=(*j)[0])
			continue;
		if ((*j)[c]<=new_element[c]){
			for (i = 1; i <=nr_gen ; i++) {
				if ((*j)[i]>new_element[i]){
					c=i;
					break;
				}
			}
			if (i==nr_gen+1) {
				Support_Hyperplanes.push_front(*j);
				Support_Hyperplanes.erase(j);
				return; //new element is not an extreme ray
			}
		}
		//new_element is reducible
	}
	Support_Hyperplanes.push_back(new_element);
}

//---------------------------------------------------------------------------

void Full_Cone::find_new_face(){
	if (verbose==true) {
		cout<<"computing new faces using the shelling ..."<<endl;
	}
	int i;
	vector<int> facet(dim-1),key;
	list< Simplex >::iterator l;
	list< int > help;
	list< int >::const_iterator m;
	set< vector<int> > Facets;
	set< vector<int> >::iterator del;
	for (l =Triangulation.begin(); l!=Triangulation.end(); l++){
		help.clear();
		key=(*l).read_key();
		for (i = 0; i <dim-1; i++) {
			facet[i]=key[i];
		}
		del=Facets.find(facet);
		if (del!=Facets.end()) {
			help.push_back(key[dim-1]);
			Facets.erase(del);
		}
		else
			Facets.insert(facet);
		for (i = dim-1; i >0; i--) {
			facet[i-1]=key[i];
			del=Facets.find(facet);
			if (del!=Facets.end()) {
				help.push_back(key[i-1]);
				Facets.erase(del);
			}
			else
				Facets.insert(facet);
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
		cout<<"done."<<endl; 
	}
}

//---------------------------------------------------------------------------
//public
//---------------------------------------------------------------------------

Full_Cone::Full_Cone(){
	dim=0;
	nr_gen=0;
	hyp_size=0;
	status="non initialized";
}

//---------------------------------------------------------------------------

Full_Cone::Full_Cone(Matrix M){
	dim=M.nr_of_columns();
	if (dim!=M.rank()) {
		error("error: Matrix with rank = number of columns needed in the constructor of the object Full_Cone.");	
	}
	Generators = M;
	nr_gen=Generators.nr_of_rows();
	vector<Integer> gcds = Generators.make_prime();
	vector<int> key=v_non_zero_pos(gcds);
	if (key.size() < nr_gen) {
		Generators=Generators.submatrix(key);
		nr_gen=Generators.nr_of_rows();
	}
	hyp_size=dim+nr_gen;
	if (dim>0) {            //correction needed to include the 0 cone;
		Linear_Form = Generators.homogeneous(homogeneous);
	} else {
		homogeneous = true;
	}
	multiplicity = 0;
	Extreme_Rays = vector<bool>(nr_gen,false);
	Support_Hyperplanes = list< vector<Integer> >();
	Triangulation = list< Simplex >();
	Hilbert_Basis =  list< vector<Integer> >();
	Homogeneous_Elements = list< vector<Integer> >();
	if(dim>0){            //correction needed to include the 0 cone;
		vector<Integer> Help7(dim);
		H_Vector=Help7;
		vector<Integer> Help8(2*dim);
		Hilbert_Polynomial=Help8;
	} else {
		vector<Integer> Help7(1,1);
		H_Vector=Help7;
		vector<Integer> Help8(2,1);
		Hilbert_Polynomial=Help8;
		Hilbert_Polynomial[0]=0;
	}
	status="initialized, before computations";
	Generators.print("bla.egn");
}

//---------------------------------------------------------------------------

Full_Cone::Full_Cone(const Full_Cone& C){
	dim=C.dim;
	nr_gen=C.nr_gen;
	hyp_size=C.hyp_size;
	status=C.status;
	homogeneous=C.homogeneous;
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

Full_Cone::~Full_Cone(){
	//automatic destructor
}

//---------------------------------------------------------------------------

void Full_Cone::print()const{
	cout<<"\ndim="<<dim<<".\n";
	cout<<"\nnr_gen="<<nr_gen<<".\n";
	cout<<"\nhyp_size="<<hyp_size<<".\n";
	cout<<"\nStatus is "<<status<<".\n";
	cout<<"\nHomogeneous is "<<homogeneous<<".\n";
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

int Full_Cone::read_dimension()const{
	return dim;
}

//---------------------------------------------------------------------------

int Full_Cone::read_nr_generators()const{
	return nr_gen;
}

//---------------------------------------------------------------------------

int Full_Cone::read_hyp_size()const{
	return hyp_size;
}

//---------------------------------------------------------------------------

string Full_Cone::read_status()const{
	return status;
}

//---------------------------------------------------------------------------

bool Full_Cone::read_homogeneous()const{
	return homogeneous;
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

//list< vector<int> > Full_Cone::get_triangulation_list() const {
//	return Triangulation;
//}
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

void Full_Cone::extreme_rays(){
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
			if (Val.read(i,j)==0) {
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
						if (Val.read(j,Zero[t])==0) {
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
}

//---------------------------------------------------------------------------

void Full_Cone::support_hyperplanes(const bool compressed_test) {
	if(dim>0){            //correction needed to include the 0 cone;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing support hyperplanes ..."<<endl;
	}
	int i,j;
	//intialization of the list of support hyperplanes
	vector<Integer> hyperplane(hyp_size,0),L,R; // initialized with 0
	Simplex S(Generators);
	vector<int> key=S.read_key();
	vector<bool> in_triang(nr_gen,false);
	vector<int> test_key(hyp_size);
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

	
	int size=2*dim;
	//test if the first simplex is compressed
	if(compressed_test){
		for (i=dim; i<size; i++) {
			if (hyperplane[i]>1) {
				Triangulation.push_back(key);
				break;
			}
		}
	}

	//computation of support hyperplanes
	Integer scalar_product;
	bool new_generator;
	int nr_non_compressed=0, nr_non_compressed_simp=0, max_heigth=0;
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
						if (compressed_test && scalar_product<-1) { //found non-compressed piece
							// make new subcone (gens in hyperplane + new gen)
							vector<int> piece;
							for (int g=dim; g<size; g++) {
								if ((*l)[g]==0) {
									piece.push_back(test_key[g]);
								}
							}
							piece.push_back(i+1);
							if (piece.size()==dim) { //simplicial case
								Simplex simp(piece);
								#pragma omp critical(TRIANG)
								Triangulation.push_back(simp);
								#pragma omp critical
								nr_non_compressed_simp++;
							}
							else {
								#pragma omp critical(NON_COMP)
								non_compressed.push_back(piece);
							}
							if (verbose) {
								#pragma omp critical
								{
									nr_non_compressed++;
									max_heigth=min(max_heigth, (int)explicit_cast_to_long(scalar_product));
								}
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
					cout<<"generator="<< i+1 <<" and "<<Support_Hyperplanes.size()<<" hyperplanes... "<<nr_non_compressed<<" "<<nr_non_compressed_simp<<" "<<max_heigth<<endl;
				}
			}
		}
	}	
	
	l_cut(Support_Hyperplanes,dim);
	if(compressed_test) {
		//compute triangulations of the not compressed, not simplicial pieces and add them to this.Triangulation
		process_non_compressed(non_compressed);
	}
	} // end if (dim>0)
	status="support hyperplanes";
	extreme_rays();
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
							error("error: Arithmetic failure in Full_cone::support_hyperplanes. Possible arithmetic overflow.\n");
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
						 L[dim-1]--;
						 k--;
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
	status="support hyperplanes";
	extreme_rays();
}

//---------------------------------------------------------------------------

void Full_Cone::support_hyperplanes_triangulation(){
	if(dim>0){            //correction needed to include the 0 cone;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing support hyperplanes and triangulation ..."<<endl;
	}
	int i,j;
	//intialization of the lists of support hyperplanes and triangulation
	vector<Integer> hyperplane(hyp_size,0),L,R; // initialized with 0
	Simplex S(Generators);
	vector<int> key=S.read_key();
	Triangulation.push_back(key);
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
				for (int k=0; k<listsize; k++) {
					for(;k > lpos; lpos++, l++) ;
					for(;k < lpos; lpos--, l--) ;

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
	status="support hyperplanes";
	extreme_rays();
}

//---------------------------------------------------------------------------

void Full_Cone::support_hyperplanes_triangulation_multiplicity(){
	if(dim>0){            //correction needed to include the 0 cone;
	support_hyperplanes_triangulation();
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing multiplicity ..."<<endl;
	}
	int counter=0;
	Integer volume;
	multiplicity=0;
	list< Simplex >::iterator l=Triangulation.begin();
	int lpos=0;
	int listsize=Triangulation.size();
	//for (l =Triangulation.begin(); l!=Triangulation.end(); l++) {
	#pragma omp parallel for private(volume) firstprivate(lpos,l) schedule(dynamic)
	for (int k=0; k<listsize; k++) {
		for(;k > lpos; lpos++, l++) ;
		for(;k < lpos; lpos--, l--) ;

		Simplex S=(*l);
		S.initialize(Generators);
		volume=S.read_volume();
		(*l).write_volume(volume);
		#pragma omp critical(MULT)
		multiplicity=multiplicity+volume;
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
	} // end if (dim>0)
	status="triangulation";
}

//---------------------------------------------------------------------------
void Full_Cone::hilbert_basis(const bool compressed_test){
	if(compressed_test) {
		//Hilbert Basis aus den nicht komprimierten Kegeln der Triangulierung bestimmen, status nacher "normal" aber keine Triangulierung!
		support_hyperplanes(true);
	} else {
		support_hyperplanes_triangulation();
	}
	Matrix SH=read_support_hyperplanes();
	if (SH.rank()!=dim) {
		error("error: Not pointed cone detected. This program is limited to pointed cones only.");
	}
	only_hilbert_basis();
}
//---------------------------------------------------------------------------
void Full_Cone::only_hilbert_basis(){
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
	for (int k=0; k<listsize; k++) {
		for(;k > lpos; lpos++, l++) ;
		for(;k < lpos; lpos--, l--) ;

		Simplex S=(*l);
		S.hilbert_basis_interior(Generators);
		volume=S.read_volume();
		(*l).write_volume(volume);
		#pragma omp critical(MULT)
		multiplicity += volume;
		HB=S.acces_hilbert_basis();
		for (h = HB.begin(); h != HB.end(); h++) {
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

	if(homogeneous==true){
		for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); h++) {
			if (v_scalar_product((*h),Linear_Form)==1) {
				Homogeneous_Elements.push_back((*h));
			}
		}
	}
	} // end if (dim>0)
	status="normal";
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
		for (cit = Candidates.begin(); cit != Candidates.end(); cit++) {
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
		cout<<"Compute norm of the candidates, "<<flush;
	}
	//go over candidates: do single scalar product
	//for (c = Candidates.begin(); c != Candidates.end(); c++) { 
	vector<Integer> scalar_product;
	for (int j=0; j<listsize; j++) {
		for(;j > cpos; cpos++, cit++) ;
		for(;j < cpos; cpos--, cit--) ;

		norm=v_scalar_product(degree_function,(*cit));

		vector <Integer> new_element(1);
		new_element[0]=norm;
		new_element=v_merge(new_element,(*cit));
		Candidates_with_Scalar_Product.push_back(new_element);
	}
	Candidates.clear();         //delete old set
	if(verbose) {
		cout<<"sort the list, "<<flush;
	}
	Candidates_with_Scalar_Product.sort();
	if (verbose) {
		cout<< Candidates_with_Scalar_Product.size() <<" candidate vectors sorted."<<endl;
	}
	
	// do global reduction
	list< vector<Integer> > HBtmp(0);
	Integer norm_crit;
	while ( !Candidates_with_Scalar_Product.empty() ) {
		if (verbose) cout<<"new reduction round"<<endl<<flush;
		//use norm criterion to find irreducible elements
		c=Candidates_with_Scalar_Product.begin();
		norm_crit=(*c)[0]*2;  //candidates with smaller norm are irreducible
		if ( Candidates_with_Scalar_Product.back()[0] < norm_crit) { //all candidates are irreducible
			if (verbose) {
				cout<<Hilbert_Basis.size()+Candidates_with_Scalar_Product.size();
				cout<<" Hilbert Basis elements of norm <= "<<norm_crit-1<<", done"<<endl;
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
			cout<<Hilbert_Basis.size()<< " Hilbert Basis elements of norm <= "<<norm_crit-1<<"; "<<csize<<" candidates left"<<endl;
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
		for (int k=0; k<csize; k++) {
			for(;k > cpos; cpos++, c++) ;
			for(;k < cpos; cpos--, c--) ;
			if (verbose && k%10000==0) {
				cout<<k<<" / "<<csize<<endl<<flush;
			}
			
			if ( is_reducible(HBpointers, *c) ) {
				(*c)[0]=-1;	//mark as reducible
			}
		}
		} //end parallel
		if (verbose) cout<<csize<<" / "<<csize<<endl<<flush;

		// delete reducible candidates
		c=Candidates_with_Scalar_Product.begin();
		while(c != Candidates_with_Scalar_Product.end()) {
			if((*c)[0]==-1) {
				c=Candidates_with_Scalar_Product.erase(c);
			} else {
				c++;
			}
		}
		HBtmp.clear();
	}
}


//---------------------------------------------------------------------------

vector<Integer> Full_Cone::compute_degree_function() const {
	if(verbose) {
		cout<<"computing degree function: ";
	}
	int i;	
	vector<Integer> degree_function(dim,0);
	if(homogeneous==true){ //use Linear_From in homogeneous case
		if(verbose) {
			cout<<" using homogenous linear form, ";
		}
		for (i=0; i<dim; i++) {
			degree_function[i] = Linear_Form[i];
		}
	} else { // add hyperplanes to get a degree function
		list< vector<Integer> >::const_iterator h;
		for (h=Support_Hyperplanes.begin(); h!=Support_Hyperplanes.end(); h++) {
			for (i=0; i<dim; i++) {
				degree_function[i]+=(*h)[i];
			}
		} //TODO parallel addition in each thread and final addition at the end
	}
	if(verbose) {
		cout<<"done"<<endl<<flush;
	}
	return degree_function;
}

//---------------------------------------------------------------------------

bool Full_Cone::low_part_simplicial(){
	support_hyperplanes_dynamic();		//change needed for dynamic lifting
	vector<Integer> val;
	int i,counter;
	list< vector<Integer> >::iterator l;
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end();){
		if ((*l)[dim-1]>0) {         // consider just the lower facets
			val=Generators.MxV((*l));
			counter=0;
			for (i = 0; i < nr_gen; i++) {
				if (val[i]==0) {
					counter++;
				}
			}
			if (counter!=dim-1) {   // more then dim vertices in one lower facet, the facet is not simplicial
				return false;
			}
			l++;
			if (l == Support_Hyperplanes.end()) {
			}
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
	for (l =Support_Hyperplanes.begin(); l != Support_Hyperplanes.end();l++){
	  for (i = 0; i <dim; i++){
		 Support_Hyperplane_With_Intersection[i]=(*l)[i];
	  }
	  Support_Hyperplane_With_Intersection[dim]=v_scalar_product((*l),line);
	  Shelling.insert(Support_Hyperplane_With_Intersection);
	}
	Support_Hyperplanes.clear();
	for (k = Shelling.begin(); k != Shelling.end(); k++) {
	  Support_Hyperplanes.push_back((*k));
	}
	l_cut(Support_Hyperplanes,dim);
	if (verbose==true) {
	  cout<<"done."<<endl;
	}
}

//---------------------------------------------------------------------------

void Full_Cone::triangulation_lift(){
	if (status!="support hyperplanes") {
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
		if (Lifted.read_status()!="support hyperplanes") {
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
		for (l = Lifted.Support_Hyperplanes.begin(); l != Lifted.Support_Hyperplanes.end(); l++) {
			v=Lifted.Generators.MxV((*l));
			counter=0;
			for (j = 0; j < nr_extreme_rays; j++) {
				if (v[j]==0) {
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
}

//---------------------------------------------------------------------------

void Full_Cone::hilbert_polynomial(){
	if (homogeneous==false) {
		hilbert_basis();
	}
	else{
		if(dim>0){            //correction needed to include the 0 cone;
		support_hyperplanes();
		triangulation_lift();
		find_new_face();
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
		for (int k=0; k<listsize; k++) {
			for(;k > lpos; lpos++, l++) ;
			for(;k < lpos; lpos--, l--) ;

			Simplex S=(*l);
			#pragma omp critical(H_VECTOR)
			H_Vector[S.read_new_face_size()]++;
			S.h_vector(Generators,Linear_Form);
			volume=S.read_volume();
			(*l).write_volume(volume);
			#pragma omp critical(MULT)
			multiplicity+=volume;
			HE=S.read_homogeneous_elements();
			#pragma omp critical(HOMOGENEOUS)
			Homogeneous_Elements.merge(HE);
			#pragma omp critical(H_VECTOR)
			H_Vector=v_add(H_Vector,S.read_h_vector());
			if (verbose==true) {
				#pragma omp critical(VERBOSE)
				{
					counter++;
					if (counter%100==0) {
						cout<<"simplex="<<counter<<endl;
					}
				}
			}
		}
		compute_polynomial();
		Homogeneous_Elements.sort();
		Homogeneous_Elements.unique();
		} // end if (dim>0)
		status="hilbert polynomial";
	}
}

//---------------------------------------------------------------------------

void Full_Cone::hilbert_basis_polynomial(){
	if (homogeneous==false) {
		hilbert_basis();
	}
	else{
		if(dim>0){            //correction needed to include the 0 cone;
		support_hyperplanes();
		triangulation_lift();
		find_new_face();
		if (verbose==true) {
			cout<<"\n************************************************************\n";
			cout<<"computing Hilbert basis and polynomial ..."<<endl;
		}
		int counter=0,i;
		Integer volume;
		vector<Integer> scalar_product;
		set <vector<Integer> > Candidates;
		list <vector<Integer> >  HB,HE;
		set <vector<Integer> >::iterator c;
		list <vector<Integer> >::const_iterator h;
		for (i = 1; i <=nr_gen; i++) {
			Candidates.insert(Generators.read(i));
			Homogeneous_Elements.push_back(Generators.read(i));
		}
		multiplicity=0;
		list< Simplex >::iterator l;
		for (l =Triangulation.begin(); l!=Triangulation.end(); l++){
			Simplex S=(*l);
			H_Vector[S.read_new_face_size()]++;
			S.hilbert_basis_interior_h_vector(Generators,Linear_Form);
			volume=S.read_volume();
			(*l).write_volume(volume);
			multiplicity=multiplicity+volume;
			HE=S.read_homogeneous_elements();
			Homogeneous_Elements.merge(HE);
			H_Vector=v_add(H_Vector,S.read_h_vector());
			HB=S.acces_hilbert_basis();
			for (h = HB.begin(); h != HB.end(); h++) {
				Candidates.insert((*h));
			}
			if (verbose==true) {
				counter++;
				if (counter%100==0) {
					cout<<"simplex="<<counter<<" and "<< Candidates.size() <<" candidate vectors to be globally reduced."<<endl;
				}
			}
		}
		
		global_reduction(Candidates);
		
		compute_polynomial();
		Homogeneous_Elements.sort();
		Homogeneous_Elements.unique();
		} // end if (dim>0)
		status="hilbert basis polynomial";
	}
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
	for (h =Support_Hyperplanes.begin(); h != Support_Hyperplanes.end(); h++){
		if ((*h)[dim-1]!=0) {
			for (t =Triangulation.begin(); t!=Triangulation.end(); t++){
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

void Full_Cone::add_hyperplane(const int& hyp_counter, const bool& lifting, vector<Integer>& halfspace){
	if (verbose==true) {
		cout<<"adding hyperplane "<<hyp_counter<<" ..."<<endl;
	}
	int i,sign;
	bool  not_done;
	list < vector <Integer> > Positive_Ired,Negative_Ired,Neutral_Ired;
	Integer orientation, scalar_product,diff,factor;
	vector <Integer> hyperplane=Generators.read(hyp_counter);
	list< vector<Integer> >::iterator h;
	if (lifting==true) {
		orientation=v_scalar_product(hyperplane,halfspace);
		if(orientation<0){
			orientation=-orientation;
			v_scalar_multiplication(halfspace,-1); //transforming into the generator of the positive halfspace
		}
		for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); h++) { //reduction  modulo  the generator of the positive halfspace
			scalar_product=v_scalar_product_unequal_vectors_end(hyperplane,(*h));
			sign=1;
			if (scalar_product<0) {
				scalar_product=-scalar_product;
				sign=-1;
			}
			factor=scalar_product/orientation;
			for (i = 0; i < dim; i++) {
				(*h)[nr_gen+3+i]=(*h)[nr_gen+3+i]-sign*factor*halfspace[i];
			}
		}
		//adding the generators of the halfspace to negative and positive
		vector <Integer> hyp_element(hyp_size+3,0);
		for (i = 0; i < dim; i++) {
			hyp_element[nr_gen+3+i]= halfspace[i];
		}
		hyp_element[hyp_counter]=orientation;
		hyp_element[0]=orientation;
		if (orientation==0){ //never
			Neutral_Ired.push_back(hyp_element);
		}
		else{
			Positive_Ired.push_back(hyp_element);
			v_scalar_multiplication(hyp_element,-1);
			hyp_element[hyp_counter]=orientation;
			hyp_element[0]=orientation;
			Negative_Ired.push_back(hyp_element);
		}
	}
	for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); h++) { //dividing into negative and positive
		(*h)[hyp_counter]=v_scalar_product_unequal_vectors_end(hyperplane,(*h));
		if ((*h)[hyp_counter]>0) {
			(*h)[nr_gen+1]=1;     // generation
			(*h)[nr_gen+2]=0;     //not sum
			(*h)[0]+=(*h)[hyp_counter];
			Positive_Ired.push_back((*h));
		}
		if ((*h)[hyp_counter]<0) {
			(*h)[hyp_counter]=-(*h)[hyp_counter];
			(*h)[nr_gen+1]=1;
			(*h)[nr_gen+2]=0;
			(*h)[0]+=(*h)[hyp_counter];
			Negative_Ired.push_back((*h));
		}
		if ((*h)[hyp_counter]==0) {
			(*h)[nr_gen+1]=0;
			(*h)[nr_gen+2]=0;
			Neutral_Ired.push_back((*h));
		}
	}
	Neutral_Ired.sort();
	Positive_Ired.sort();
	Negative_Ired.sort();
	//long int counter=0;
	list < vector <Integer> > New_Positive,New_Negative,New_Neutral,Pos,Neg,Neu;
	list < vector<Integer> >::const_iterator p,n;
	list < vector <Integer> >::iterator c;
	not_done=true;
	while(not_done){
		not_done=false;
		New_Positive.clear();
		New_Negative.clear();
		New_Neutral.clear();
		//generating new elements
		cout<<"+"<<flush;
//		for(p = Positive_Ired.begin(); p != Positive_Ired.end(); p++){
		list < vector<Integer>* > Positive,Negative,Neutral; // pointer lists, used to move reducers to the front
		list < vector<Integer> >::iterator it;
		it=Positive_Ired.begin();
		while (it!=Positive_Ired.end()) {
			Positive.push_back(&(*(it++)));
		}
		it=Negative_Ired.begin();
		while (it!=Negative_Ired.end()) {
			Negative.push_back(&(*(it++)));
		}
		it=Neutral_Ired.begin();
		while (it!=Neutral_Ired.end()) {
			Neutral.push_back(&(*(it++)));
		}
		int psize=Positive_Ired.size();
		#pragma omp parallel private(p,n,diff) firstprivate(Positive,Negative,Neutral)
		{
		int ppos=0;
		p = Positive_Ired.begin();
		#pragma omp for schedule(dynamic)
		for(i = 0; i<psize; i++){
			for(;i > ppos; ppos++, p++) ;
			for(;i < ppos; ppos--, p--) ;

			for (n = Negative_Ired.begin(); n != Negative_Ired.end(); n++){
				if ((*p)[nr_gen+1]<=1&&(*n)[nr_gen+1]<=1&&((*p)[nr_gen+1]!=0||(*n)[nr_gen+1]!=0)) {
					if (((*p)[nr_gen+2]!=0&&(*p)[nr_gen+2]<=(*n)[hyp_counter])||((*n)[nr_gen+2]!=0&&(*n)[nr_gen+2]<=(*p)[hyp_counter]))
						continue;
					//	counter++;
					diff=(*p)[hyp_counter]-(*n)[hyp_counter];
					vector <Integer> new_candidate=v_add((*p),(*n));

					if (diff>0) {
						new_candidate[hyp_counter]=diff;
						new_candidate[0]-=2*(*n)[hyp_counter];
						if(reduce(Positive,new_candidate,hyp_counter)) {
							continue;
						}
						if(reduce(Neutral,new_candidate,hyp_counter-1)) {
							continue;
						}    
						new_candidate[nr_gen+1]=2;
						new_candidate[nr_gen+2]=(*p)[hyp_counter];
						#pragma omp critical(NEW_POSITIVE)
						New_Positive.push_back(new_candidate);
					}
					if (diff<0) {
						new_candidate[hyp_counter]=-diff;
						new_candidate[0]-=2*(*p)[hyp_counter];
						if(reduce(Negative,new_candidate,hyp_counter)) {
							continue;
						}
						if(reduce(Neutral,new_candidate,hyp_counter-1)) {
							continue;
						}    
						new_candidate[nr_gen+1]=2;
						new_candidate[nr_gen+2]=(*n)[hyp_counter];
						#pragma omp critical(NEW_NEGATIVE)
						New_Negative.push_back(new_candidate);
					}
					if (diff==0) {
						new_candidate[hyp_counter]=0;
						new_candidate[0]-=2*(*p)[hyp_counter];
						if(reduce(Neutral,new_candidate,hyp_counter-1)) {
							continue;
						}
						new_candidate[nr_gen+1]=0;
						new_candidate[nr_gen+2]=0;
						#pragma omp critical(NEW_NEUTRAL)
						New_Neutral.push_back(new_candidate);
					}
					/*	if (counter==10000000) {
						counter=0;
						New_Neutral.sort();
						New_Positive.sort();
						New_Negative.sort();
						reduce(Neutral,New_Neutral,hyp_counter-1);
						reduce(Neutral,New_Positive,hyp_counter-1);
						reduce(Neutral,New_Negative,hyp_counter-1);
						reduce(Positive,New_Positive,hyp_counter);
						reduce(Negative,New_Negative,hyp_counter);
						Neu.merge(New_Neutral);
						Pos.merge(New_Positive);
						Neg.merge(New_Negative);
						New_Neutral.clear();
						New_Positive.clear();
						New_Negative.clear();
						}     */
				}
			}
		}
		//end generation of new elements
		#pragma omp single nowait
		New_Neutral.sort();
		#pragma omp single nowait
		New_Positive.sort();
		#pragma omp single nowait
		New_Negative.sort();
		} //END PARALLEL
		cout<<"-"<<flush;
		//reducing the new vectors agains them self
		//Neutral_Ired=Neutral;
		//Positive_Ired=Positive;
		//Negative_Ired=Negative;
		/*reduce(Neutral,New_Neutral,hyp_counter-1);
		  reduce(Neutral,New_Positive,hyp_counter-1);
		  reduce(Neutral,New_Negative,hyp_counter-1);
		  New_Positive.sort();
		  New_Negative.sort();
		  reduce(Positive,New_Positive,hyp_counter);
		  reduce(Negative,New_Negative,hyp_counter);
		  New_Neutral.sort();
		  New_Positive.sort();
		  New_Negative.sort();
		  New_Neutral.merge(Neu);
		  New_Positive.merge(Pos);
		  New_Negative.merge(Neg);*/
		Hilbert_Basis.clear();
		for(c=New_Neutral.begin();c != New_Neutral.end();c++) {
			reduce_and_insert((*c),hyp_counter-1);
		}
		New_Neutral=Hilbert_Basis;
		Hilbert_Basis.clear();
		for(c=New_Positive.begin();c != New_Positive.end();c++) {
			reduce_and_insert((*c),hyp_counter);
		}
		New_Positive=Hilbert_Basis;
		Hilbert_Basis.clear();
		for(c=New_Negative.begin();c != New_Negative.end();c++) {
			reduce_and_insert((*c),hyp_counter);
		}
		New_Negative=Hilbert_Basis;
		if (New_Neutral.size()!=0) {
			New_Positive.sort();
			reduce(New_Neutral,New_Positive, hyp_counter-1);
			New_Negative.sort();
			reduce(New_Neutral,New_Negative, hyp_counter-1);
			reduce(New_Neutral,Neutral_Ired, hyp_counter-1);
			reduce(New_Neutral,Positive_Ired, hyp_counter-1);
			reduce(New_Neutral,Negative_Ired, hyp_counter-1);   
			Neutral_Ired.merge(New_Neutral);
		}
		if (New_Positive.size()!=0) {
			not_done=true;
			reduce(New_Positive,Positive_Ired, hyp_counter);
			Positive_Ired.merge(New_Positive);
		}
		if (New_Negative.size()!=0) {
			not_done=true;
			reduce(New_Negative,Negative_Ired, hyp_counter);
			Negative_Ired.merge(New_Negative);
		}
		for (c = Positive_Ired.begin(); c != Positive_Ired.end(); c++){
			if((*c)[nr_gen+1]>0) {
				(*c)[nr_gen+1]--;
			}
		}
		for (c = Negative_Ired.begin(); c != Negative_Ired.end(); c++){
			if((*c)[nr_gen+1]>0) {
				(*c)[nr_gen+1]--;
			}
		}
	}
	//still posible to have double elements in the Hilbert basis, coming from different generations

	set< vector<Integer> > Help;
	set< vector<Integer> >::iterator d;
	for (c = Positive_Ired.begin(); c != Positive_Ired.end(); c++) {
		(*c)[nr_gen+1]=0;
		(*c)[nr_gen+2]=0;
		Help.insert(*c);
	}
	for (c = Neutral_Ired.begin(); c != Neutral_Ired.end(); c++) {
		(*c)[nr_gen+1]=0;
		(*c)[nr_gen+2]=0;
		Help.insert(*c);
	}
	Hilbert_Basis.clear();
	d=Help.begin();
	while(d != Help.end()) {
		Hilbert_Basis.push_back(*d);
		Help.erase(d);
		d=Help.begin();
	}
	if (verbose==true) {
		cout<<"Hilbert basis size="<<Hilbert_Basis.size()<<endl;
	}
}

//---------------------------------------------------------------------------

Matrix Full_Cone::add_hyperplane(const int& hyp_counter, const Matrix& Basis_Max_Subspace){
	int i,j,rank_subspace=Basis_Max_Subspace.nr_of_rows();
	vector <Integer> scalar_product,hyperplane=Generators.read(hyp_counter),halfspace;
	bool lifting=false;
	Matrix New_Basis_Max_Subspace=Basis_Max_Subspace;
	if (rank_subspace!=0) {
		scalar_product=Basis_Max_Subspace.MxV(hyperplane);
		for (i = 0; i <rank_subspace; i++) 
			if (scalar_product[i]!=0)
				break;
		if (i!=rank_subspace) {    // the new hyperplane is not contained in the maximal subspace
			lifting=true;
			//computing new maximal subspace
			Matrix M(1,rank_subspace);
			M.write(1,scalar_product);
			Lineare_Transformation LT=Transformation(M);
			Matrix Lifted_Basis_Factor_Space_over_Ker_and_Ker=LT.get_right();
			Lifted_Basis_Factor_Space_over_Ker_and_Ker=Lifted_Basis_Factor_Space_over_Ker_and_Ker.transpose();
			Matrix  Ker(rank_subspace-1,rank_subspace);
			for (j = 1; j <= rank_subspace-1; j++) {
				Ker.write(j, Lifted_Basis_Factor_Space_over_Ker_and_Ker.read(j+1));
			}
			New_Basis_Max_Subspace=Ker.multiplication(Basis_Max_Subspace);
			halfspace=Basis_Max_Subspace.VxM(Lifted_Basis_Factor_Space_over_Ker_and_Ker.read(1));
		}
	}
	add_hyperplane(hyp_counter, lifting, halfspace);
	return New_Basis_Max_Subspace;
}

//---------------------------------------------------------------------------

void Full_Cone::extreme_rays_reduction(){
	list < vector <Integer> >::iterator c;
	int i,k;
	for (c=Hilbert_Basis.begin();c!=Hilbert_Basis.end();c++){
		k=0;
		for (i = 1; i <= nr_gen; i++) {
			if ((*c)[i]!=0) {
				(*c)[i]=1;
				k++;
			}

		}
		(*c)[0]=k;       //if (k>=dim-1) improuves speed much here
	}
	Hilbert_Basis.sort();
	for (c=Hilbert_Basis.begin();c!=Hilbert_Basis.end();c++){
		reduce_and_insert_extreme((*c));
	}
	l_cut_front(Support_Hyperplanes,dim);
}

//---------------------------------------------------------------------------

void Full_Cone::extreme_rays_rank(){
	list < vector <Integer> >::iterator c;
	list <int> zero_list;
	int i,j,k;
	for (c=Hilbert_Basis.begin();c!=Hilbert_Basis.end();c++){
		zero_list.clear();
		for (i = 1; i <= nr_gen; i++) {
			if ((*c)[i]==0) {
				zero_list.push_back(i);
			}
		}
		k=zero_list.size();
		if (k>=dim-1) {
			vector <int> zero_vector(k);
			for (j = 0; j < k; j++) {
				zero_vector[j]=zero_list.front();
				zero_list.pop_front();
			}
			Matrix Test=Generators.submatrix(zero_vector);
			if (Test.rank()>=dim-1) {
				Support_Hyperplanes.push_back((*c));
			}
		}
	}
	l_cut_front(Support_Hyperplanes,dim);
}

//---------------------------------------------------------------------------

void Full_Cone::hilbert_basis_dual(){
if(dim>0){            //correction needed to include the 0 cone;
	if (verbose==true) {
		cout<<"\n************************************************************\n";
		cout<<"computing Hilbert basis ..."<<endl;
	}
	int hyp_counter;      // current hyperplane
	Matrix Basis_Max_Subspace(dim);      //identity matrix
	for (hyp_counter = 1; hyp_counter <= nr_gen; hyp_counter++) {
		Basis_Max_Subspace=add_hyperplane(hyp_counter,Basis_Max_Subspace);
	}
	extreme_rays_rank();
	l_cut_front(Hilbert_Basis,dim);
	Matrix M=read_support_hyperplanes();
	if(dim>0){            //correction needed to include the 0 cone;
		Linear_Form=M.homogeneous_low_dim(homogeneous);
		if(homogeneous==true){
			list < vector <Integer> >::const_iterator h;
			for (h = Hilbert_Basis.begin(); h != Hilbert_Basis.end(); h++) {
				if (v_scalar_product((*h),Linear_Form)==1) {
					Homogeneous_Elements.push_back((*h));
				}
			}
		}
	} else {
		homogeneous=true;
		vector<Integer> Help(dim);
		Linear_Form=Help;
	}
	} // end if (dim>0)
	status="dual";
}

//---------------------------------------------------------------------------

void Full_Cone::error(string s) const{
	cerr <<"\nFull Cone "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

void lift(Full_Cone& Lifted,Matrix Extreme_Generators){
	if (verbose==true) {
		cout<<"start lifting the cone ...";
	}
	int i,j,counter=0,nr_extreme_gen=Extreme_Generators.nr_of_rows(),dim=Extreme_Generators.nr_of_columns();
	Matrix New_Generators(nr_extreme_gen,dim+1);
	for (i = 1; i <= nr_extreme_gen; i++) {
		for (j = 1; j <= dim; j++) {
			New_Generators.write(i,j,Extreme_Generators.read(i,j));
		}
	}
	// try 10 times a random lifting
	time_t t;
	srand((unsigned) time(&t));
	while(counter<10){
		for (i = 1; i <= nr_extreme_gen; i++){
			j=rand();
			j=j%lifting_bound;
			New_Generators.write(i,dim+1,j+1);
		}
		if (New_Generators.rank()==dim+1) {
			Full_Cone Help(New_Generators);
			Lifted=Help;
			if (Lifted.low_part_simplicial()==true) {
				if (verbose==true) {
					cout<<"lifting done."<<endl;
				}
				return;
			}
		}

		counter++;
	}
	cerr<<"error: Random lifting has failed in Full_Cone::lift.";
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

	for (; ss!=Support_Hyperplanes.end()&&is_compressed; ss++) {
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
		cout<<"verbose_step="<<verbose_step<<endl;
	}

	#pragma omp parallel for firstprivate(itpos,it) schedule(dynamic)
	for (int i=0; i<listsize; i++) {
		for(;i > itpos; itpos++, it++) ;
		for(;i < itpos; itpos--, it--) ;
//		v_read(*it); cout<<" subcone "<<i+1;
		//compute triangulations of subcones
		Full_Cone subcone(Generators.submatrix(*it));
		subcone.support_hyperplanes_triangulation();	// use normal method to compute triangulation of subcone
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
			sub_it++;
		}
		if (verbose_bak && (i+1)%verbose_step==0) {
			cout<<i+1<<" subcones done, found "<<Triangulation.size()<<" simplices"<<endl;
		}
	}
	if (verbose_bak) {
		cout<<listsize<<" subcones done, found "<<Triangulation.size()<<" simplices"<<endl;
	}

	//restore global verbose value
	verbose=verbose_bak;
}

