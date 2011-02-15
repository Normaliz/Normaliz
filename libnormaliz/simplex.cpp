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

#include <algorithm>
#include <string>
#include <iostream>
#include <set>

#include "integer.h"
#include "vector_operations.h"
#include "matrix.h"
#include "simplex.h"
#include "list_operations.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//Private
//---------------------------------------------------------------------------
template<typename Integer>
void Simplex<Integer>::reduce_and_insert_interior(const vector< Integer >& new_element){
	//implementing this function as a tree searching may speed up computations ...
	if (new_element[0]==0) {
		return; // new_element=0
	}
	else {
		register size_t i,c=1,d=dim+1;
		typename list< vector<Integer> >::iterator j;
		for (j =Hilbert_Basis.begin(); j != Hilbert_Basis.end(); j++) {
			if (new_element[0]<2*(*j)[0]) {
				break; //new_element is not reducible;
			}
			else  {
				if ((*j)[c]<=new_element[c]){
					for (i = 1; i < d; i++) {
						if ((*j)[i]>new_element[i]){
							c=i;
							break;
						}
					}
					if (i==d) {
						Hilbert_Basis.push_front(*j);
						Hilbert_Basis.erase(j);
						return;
					}
					//new_element is not in the Hilbert Basis
				}
			}
		}
		Hilbert_Basis.push_back(new_element);
	}
}

//---------------------------------------------------------------------------
//Public
//---------------------------------------------------------------------------

template<typename Integer>
Simplex<Integer>::Simplex(){
	status="non initialized";
}

template<typename Integer>
Simplex<Integer>::Simplex(const vector<size_t>& k){
	dim=k.size();
	key=k;
	volume=0;
	status="key initialized";
}

template<typename Integer>
Simplex<Integer>::Simplex(const Matrix<Integer>& Map){
	dim=Map.nr_of_columns();
	key=Map.max_rank_submatrix_lex(dim);
	Generators=Map.submatrix(key);
	vector< Integer > help(dim);
	Support_Hyperplanes=Invert(Generators, help, volume); //test for arithmetic
	//overflow performed
	diagonal=v_abs(help);
	Support_Hyperplanes=Support_Hyperplanes.transpose();
	multiplicators=Support_Hyperplanes.make_prime();
	Hilbert_Basis = list< vector<Integer> >();
	H_Vector = vector<Integer>(dim,0);
	status="initialized";
}

//---------------------------------------------------------------------------

template<typename Integer>
Simplex<Integer>::Simplex(const vector<size_t>& k, const Matrix<Integer>& Map){
	key=k;
	Generators=Map.submatrix(k);
	dim=k.size();
	vector< Integer > help(dim);
	Support_Hyperplanes=Invert(Generators, help, volume);  //test for arithmetic
	//overflow performed
	diagonal=v_abs(help);
	Support_Hyperplanes=Support_Hyperplanes.transpose();
	multiplicators=Support_Hyperplanes.make_prime();
	Hilbert_Basis = list< vector<Integer> >();
	H_Vector=vector<Integer>(dim,0);
	status="initialized";
}

//---------------------------------------------------------------------------

template<typename Integer>
Simplex<Integer>::Simplex(const Simplex<Integer>& S){
	dim=S.dim;
	status=S.status;
	volume=S.volume;
	key=S.key;
	Generators=S.Generators;
	diagonal=S.diagonal;
	multiplicators=S.multiplicators;
	New_Face=S.New_Face;
	Support_Hyperplanes=S.Support_Hyperplanes;
	Hilbert_Basis=S.Hilbert_Basis;
	Ht1_Elements=S.Ht1_Elements;
	H_Vector=S.H_Vector;
}

//---------------------------------------------------------------------------

template<typename Integer>
Simplex<Integer>::~Simplex(){
	//automatic destructor
}

//---------------------------------------------------------------------------

template<typename Integer>
void Simplex<Integer>::write_new_face(const vector<size_t>& Face){
	New_Face=Face;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Simplex<Integer>::read() const{
	cout<<"\nDimension="<<dim<<"\n";
	cout<<"\nStatus="<<status<<"\n";
	cout<<"\nVolume="<<volume<<"\n";
	cout<<"\nKey is:\n";
	v_read(key);
	cout<<"\nGenerators are:\n";
	Generators.read();
	cout<<"\nDiagonal is:\n";
	v_read(diagonal);
	cout<<"\nMultiplicators are:\n";
	v_read(multiplicators);
	cout<<"\nNew face is:\n";
	v_read(New_Face);
	cout<<"\nSupport Hyperplanes are:\n";
	Support_Hyperplanes.read();
	Matrix<Integer> M=read_hilbert_basis();
	cout<<"\nHilbert Basis is:\n";
	M.read();
	cout<<"\nh-vector is:\n";
	v_read(H_Vector);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Simplex<Integer>::read_k() const{
	v_read(key);
	v_read(New_Face);
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Simplex<Integer>::read_dimension() const{
	return dim;
}

//---------------------------------------------------------------------------

template<typename Integer>
string Simplex<Integer>::read_status() const{
	return status;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Simplex<Integer>::write_volume(const Integer& vol){
	volume=vol;
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Simplex<Integer>::read_volume() const{
	return volume;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t> Simplex<Integer>::read_key() const{
	return key;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Simplex<Integer>::read_generators() const{
	return Generators;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Simplex<Integer>::read_diagonal() const{
	return diagonal;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Simplex<Integer>::read_multiplicators() const{
	return multiplicators;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t> Simplex<Integer>::read_new_face() const{
	return New_Face;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Simplex<Integer>::read_new_face_size() const{
	return New_Face.size();
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Simplex<Integer>::read_support_hyperplanes() const{
	return Support_Hyperplanes;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Simplex<Integer>::read_hilbert_basis()const{
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
list< vector<Integer> > Simplex<Integer>::read_ht1_elements()const{
	list< vector<Integer> > HE=Ht1_Elements;
	return HE;
}

//---------------------------------------------------------------------------

template<typename Integer>
const list< vector<Integer> >& Simplex<Integer>::acces_hilbert_basis()const{
	const list< vector<Integer> >& HB=Hilbert_Basis;
	return HB;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Simplex<Integer>::read_h_vector() const{
	return H_Vector;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Simplex<Integer>::read_hilbert_basis_size() const{
	return Hilbert_Basis.size();
}

//---------------------------------------------------------------------------

template<typename Integer>
int Simplex<Integer>::compare(const Simplex<Integer>& S) const{
	return v_difference_ordered_fast(key,S.key);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Simplex<Integer>::initialize(const Matrix<Integer>& Map){
	assert(status != "non initialized");

	if (status=="key initialized") {
		Generators=Map.submatrix(key);
		vector< Integer > help(dim);
		Support_Hyperplanes=Invert(Generators, help, volume); //test for arithmetic
		//overflow performed
		diagonal=v_abs(help);
		Support_Hyperplanes=Support_Hyperplanes.transpose();
		multiplicators=Support_Hyperplanes.make_prime();
		Hilbert_Basis=list< vector<Integer> >();
		Ht1_Elements=list< vector<Integer> >();
		H_Vector=vector<Integer>(dim,0);
		status="initialized";
	}
}

//---------------------------------------------------------------------------

//TODO remove
size_t NrInvert=0;

/* evaluates a simplex in regard to all data, key must be initialized */
template<typename Integer>
Integer Simplex<Integer>::evaluate(Full_Cone<Integer>& C, const Integer& height) {
 
	Generators=C.Generators.submatrix(key);
	H_Vector=vector<Integer>(dim,0);
	
	bool unimodular=false;
	vector<Integer> Indicator;
	if(height >=-1 || (!C.do_h_vector && !C.do_Hilbert_basis && !C.do_ht1_elements)) {
		Matrix<Integer> RS(1,dim);  // (transpose of) right hand side
		RS.write(1,C.Order_Vector); // to hold order vector
		RS=RS.transpose();
		vector<Integer> diag(dim);
		Matrix<Integer> Sol=Generators.transpose().solve_destructiv(RS,diag,volume);
		Indicator=Sol.transpose().read(1);
		if(volume==1)
			unimodular=true;
	}
			
	// in this cases we had to add the volume and nothing else is to be done
	if ( (!C.do_h_vector && !C.do_Hilbert_basis && !C.do_ht1_elements) 
	  || (unimodular && !C.do_h_vector) ) {
		#pragma omp critical(MULTIPLICITY)
		C.multiplicity+=volume;
		return volume;
	}

	bool decided=true;
	size_t i,j;
	size_t Deg=0;        // Deg is the degree in which the 0 vector is counted
	if(unimodular){   // it remains to count the 0-vector in the h-vector 
		for(i=0;i<dim;i++){
			if(Indicator[i]<0)   // facet opposite of vertex i excluded
				Deg++;
			if(Indicator[i]==0){ // Order_Vector in facet, to be decided later
				decided=false;
				break;
			}
		}
		if(decided){
			#pragma omp critical(HVECTOR) //only change in the H-vector, so done directliy
			C.H_Vector[Deg]++;    // Done, provided decided==true
			#pragma omp critical(MULTIPLICITY)
			C.multiplicity+=volume;
			return volume;               // if not we need lex decision, see below
		}
	} // We have tried to take care of the unimodular case WITHOUT the matrix inversion

	#pragma omp atomic
	NrInvert++;
	vector< Integer > help(dim);
	Matrix<Integer> InvGen=Invert(Generators, help, volume);
	diagonal=v_abs(help);
	vector<bool> Excluded(dim,false);
	Integer Test; 
	
	if(C.do_h_vector){
		Deg=0;
		if (Indicator.size() != dim) { //it hasn't been computed yet
			Indicator = InvGen.VxM(C.Order_Vector);
		}
		for(i=0;i<dim;i++) // register excluded facets and degree shift for 0-vector
		{
			Test=Indicator[i];
			if(Test<0)
			{
				Excluded[i]=true; // the facet opposite to vertex i is excluded
				Deg++;
			}
			if(Test==0){  // Order_Vector in facet, now lexicographic decision
				for(j=0;j<dim;j++){
					if(InvGen[j][i]<0){ // COLUMNS of InvGen give supp hyps
						Excluded[i]=true;
						Deg++;
						break;
					}
					if(InvGen[j][i]>0) // facet included
						break;
				}
			}
		}
		H_Vector[Deg]++; // now the 0 vector is finally taken care of
		if(unimodular){     // and in the unimodular case nothing left to be done
			#pragma omp critical(HVECTOR)
			C.H_Vector[Deg]++; 
			#pragma omp critical(MULTIPLICITY)
			C.multiplicity+=volume;
			return volume;
		}
	}
	
	vector < Integer > norm(1);
	list < vector<Integer> > Candidates;
	typename list <vector <Integer> >::iterator c;
	size_t last;
	vector<Integer> point(dim,0);
 
	Matrix<Integer> elements(dim,dim); //all 0 matrix
	Matrix<Integer> V = InvGen; //Support_Hyperplanes.multiply_rows(multiplicators).transpose();
	V.reduction_modulo(volume); //makes reduction when adding V easier

	while (1){
		last = dim;
		for (int k = dim-1; k >= 0; k--) {
			if (point[k] < diagonal[k]-1) {
				last = k;
				break;
			}
		}
		if (last >= dim) {
			break;
		}

		point[last]++;
		elements[last] = v_add(elements[last], V[last]); //TODO angepasste operation
		v_reduction_modulo(elements[last],volume);       //für beide schritte schreiben

		for (i = last+1; i <dim; i++) {
			point[i]=0;
			elements[i] = elements[last];
		}
		
		norm[0]=0; // norm[0] is just the sum of coefficients, = volume*degree
		for (i = 0; i < dim; i++) {  // since generators have degree 1
			norm[0]+=elements[last][i];
		}

		if(C.do_h_vector){
			Deg=explicit_cast_to_long<Integer>(norm[0]/volume); // basic degree, here we use that all generators have degree 1            
			for(i=0;i<dim;i++)  // take care of excluded facets and increase degree where necessary
				if(elements[last][i]==0 && Excluded[i])
					Deg++;
			
			H_Vector[Deg]++; // count element in h-vector        
		}
		
		if(C.do_ht1_elements && norm[0]==volume) // found degree 1 element
		{        
			help=Generators.VxM(elements[last]);
			v_scalar_division(help,volume);
			Ht1_Elements.push_back(help);
			continue;
		} 
		
		// now we are left with the case of Hilbert bases

		if(C.do_Hilbert_basis){
			Candidates.push_back(v_merge(norm,elements[last]));
		}
	}
	
	if(C.do_h_vector) {
		for(size_t i=0; i<dim; i++) {
			if(H_Vector[i]!=0) {
				#pragma omp critical(HVECTOR)
				C.H_Vector[i]+=H_Vector[i];
			}
		}
	}
	
	if(C.do_ht1_elements) {
		#pragma omp critical(HT1ELEMENTS)
		C.Ht1_Elements.splice(C.Ht1_Elements.begin(),Ht1_Elements);
	}
	
	if(!C.do_Hilbert_basis) {
		#pragma omp critical(MULTIPLICITY)
		C.multiplicity+=volume;    
		return volume;
	}

	Candidates.sort();        
	typename list <vector <Integer> >::iterator cand=Candidates.begin();
	while(cand != Candidates.end()) {
		reduce_and_insert_interior((*cand));
		Candidates.pop_front();
		cand=Candidates.begin();
	}

	//inverse transformation
	//some test for arithmetic overflow may be implemented here

	l_cut_front(Hilbert_Basis,dim);
	typename list< vector<Integer> >::iterator jj;
	for (jj =Hilbert_Basis.begin(); jj != Hilbert_Basis.end(); jj++) {
		*jj=Generators.VxM(*jj);
		v_scalar_division(*jj,volume);
	} 
	
	#pragma omp critical(CANDIDATES)
	C.Candidates.splice(C.Candidates.begin(),Hilbert_Basis);
	
	#pragma omp critical(MULTIPLICITY)
	C.multiplicity+=volume;
	return volume; 
}

} /* end namespace */
