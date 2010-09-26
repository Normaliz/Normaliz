/*
 * Normaliz 2.5
 * Copyright (C) 2007-2010  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
#include <list>
#include <iostream>
#include <algorithm>
#include <map>

#include "mode.h"
#include "integer.h"
#include "vector_operations.h"
#include "simplex.h"
#include "lineare_transformation.h"
#include "sublattice_representation.h"
#include "full_cone.h"
#include "cone_dual_mode.h"

//---------------------------------------------------------------------------

namespace libnormaliz {

template<typename Integer>
void make_main_computation(const int& mode, string& computation_type,const Matrix<Integer>& Input, Output<Integer>& Out){
	if ((mode<0 || mode>6) && mode!=10) {
		cerr<<"warning: Unknown computation mode "<<mode<<". The program will run in computation mode 0."<<endl;
		run_mode_0(computation_type,Input, Out);
		return;
	}
	int dim=Input.nr_of_columns();
	switch (mode){
		case 0: run_mode_0(computation_type,Input, Out);break;
		case 1: run_mode_1(computation_type,Input, Out);break;
		case 2: run_mode_2(computation_type,Input, Out);break;
		case 3: run_mode_3(computation_type,Input, Out);break;
		case 4: run_mode_456(computation_type, Matrix<Integer>(0,dim), Matrix<Integer>(0,dim), Input, Out);break;
		case 5: run_mode_456(computation_type, Matrix<Integer>(0,dim), Input, Matrix<Integer>(0,dim), Out);break;
		case 6: run_mode_456(computation_type, Input, Matrix<Integer>(0,dim), Matrix<Integer>(0,dim), Out);break;
		case 10: run_mode_10(computation_type,Input, Out);break;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_0( string& computation_type,const Matrix<Integer>& Input, Output<Integer>& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in hilbert_basis mode."<<endl;
		computation_type="hilbert_basis";
	}
	Sublattice_Representation<Integer> Basis_Change(Input,true);
	Matrix<Integer> Full_Cone_Generators = Basis_Change.to_sublattice(Input);
	Full_Cone<Integer> Result = make_computations(computation_type, Full_Cone_Generators);

	Out.set_result(Result);
	Out.compose_basis_change(Basis_Change);
	Out.cone();
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_1( string& computation_type,const Matrix<Integer>& Input, Output<Integer>& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in computation mode hilbert_basis."<<endl;
		computation_type="hilbert_basis";
	}
	Sublattice_Representation<Integer> Basis_Change(Input,false);
	Matrix<Integer> Full_Cone_Generators = Basis_Change.to_sublattice(Input);
	Full_Cone<Integer> Result = make_computations(computation_type, Full_Cone_Generators);

	Out.set_result(Result);
	Out.compose_basis_change(Basis_Change);
	Out.cone();
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_2( string& computation_type,const Matrix<Integer>& Input, Output<Integer>& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in computation mode hilbert_basis."<<endl;
		computation_type="hilbert_basis";
	}
	int i,j,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
	Integer number;
	Matrix<Integer> Generators(nr_rows,nr_columns+1,1);
	for(i=1; i<=nr_rows; i++){
		for(j=1; j<=nr_columns; j++) {
			number=Input.read(i,j);
			Generators.write(i,j,number);
		}
	}

	Sublattice_Representation<Integer> Basis_Change(Generators,true);
	Matrix<Integer> Full_Cone_Generators = Basis_Change.to_sublattice(Generators);
	Full_Cone<Integer> Result = make_computations(computation_type, Full_Cone_Generators);

	Out.set_result(Result);
	Out.compose_basis_change(Basis_Change);
	Out.polytop();
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_3( string& computation_type,const Matrix<Integer>& Input, Output<Integer>& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in computation mode hilbert_basis."<<endl;
		computation_type="hilbert_basis";
	}
	int i,j,k,l,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
	bool primary=true;
	Integer number;
	Matrix<Integer> Full_Cone_Generators(nr_rows+nr_columns,nr_columns+1,0);
	for (i = 1; i <= nr_columns; i++) {
		Full_Cone_Generators.write(i,i,1);
	}
	for(i=1; i<=nr_rows; i++){
		Full_Cone_Generators.write(i+nr_columns,nr_columns+1,1);
		for(j=1; j<=nr_columns; j++) {
			number=Input.read(i,j);
			Full_Cone_Generators.write(i+nr_columns,j,number);
		}
	}
	Matrix<Integer> Prim_Test=Input;
	for(i=1; i<=nr_rows; i++){           //preparing the  matrix for primarity test
		k=0;
		for(j=1; j<=nr_columns; j++) {
			if (k<2) {
				if (Input.read(i,j)!=0 )
					k++;
			}
			if (k==2) {
				for (l = 1; l <= nr_columns; l++) {
					Prim_Test.write(i,l,0);
				}
				break;
			}
		}
	}
	for(j=1; j<=nr_columns; j++){         //primarity test
		for(i=1; i<=nr_rows && Prim_Test.read(i,j)==0; i++);
		if (i>nr_rows) {
			primary=false;
			break;
		}
	}
	Full_Cone<Integer> Result=make_computations(computation_type, Full_Cone_Generators);
	Out.set_result(Result);
	Out.compose_basis_change(Sublattice_Representation<Integer>(Result.read_dimension()));
	Out.rees(primary);
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_4( string& computation_type,const Matrix<Integer>& Input, const int& nr_equations, Output<Integer>& Out){
	int i;
	int dim=Input.nr_of_columns();
	Matrix<Integer> Equations(nr_equations, dim);
	int nr_inequalities=Input.nr_of_rows()-nr_equations;
	Matrix<Integer> Inequalities(nr_inequalities, dim);
	for(i=1;i<=nr_inequalities;i++){
		Inequalities.write(i,Input.read(i));
	}
	for(i=1;i<=nr_equations;i++){
		Equations.write(i,Input.read(i+nr_inequalities));
	}
	run_mode_equ_inequ(computation_type,Equations,Inequalities,Out);
	
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_5( string& computation_type,const Matrix<Integer>& Input, Output<Integer>& Out){
	int dim=Input.nr_of_columns();
	Matrix<Integer> Inequalities(dim);
	run_mode_equ_inequ(computation_type,Input,Inequalities,Out);    
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_456(string& computation_type, const Matrix<Integer>& Congruences, Matrix<Integer> Equations,  Matrix<Integer> Inequalities, Output<Integer>& Out) {
	

	int dim = 0;
	int nr_cong = Congruences.nr_of_rows();
	if (nr_cong > 0) {	
		dim = Congruences.nr_of_columns() -1;
		if (Equations.nr_of_rows() > 0 &&  Equations.nr_of_columns() != dim) {
			cerr << "Error: dimensions of input matricies do not match!";
			global_error_handling();
		}
	} else if (Equations.nr_of_rows() > 0) {
		dim = Equations.nr_of_columns();
	} else if (Inequalities.nr_of_rows() > 0) {
		dim = Inequalities.nr_of_columns();
	}
	// use positive orthant if no inequalities are given
	if (Inequalities.nr_of_rows() == 0) {
		Inequalities = Matrix<Integer>(Equations.nr_of_columns());
	}
		
	// handle Congurences
	if (nr_cong > 0) {
		int i,j;
		
		//add slack variables
		Matrix<Integer> Cong_Slack(nr_cong, dim+nr_cong);
		for (i = 1; i <= nr_cong; i++) {
			for (j = 1; j <= dim; j++) {
				Cong_Slack.write(i,j,Congruences.read(i,j));
			}
			Cong_Slack.write(i,dim+i,Congruences.read(i,dim+1));
		}

		//compute kernel
		Lineare_Transformation<Integer> Diagonalization = Transformation(Cong_Slack);
		int rank = Diagonalization.get_rank();
		Matrix<Integer> H = Diagonalization.get_right();
		Matrix<Integer> Ker_Basis_Transpose(dim, dim+nr_cong-rank);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim+nr_cong; j++) {
				Ker_Basis_Transpose.write(i,j-rank,H.read(i,j));
			}
		}

		//TODO now a new linear transformation is computed, necessary??
		Sublattice_Representation<Integer> Basis_Change(Ker_Basis_Transpose.transpose(),false);
		Out.compose_basis_change(Basis_Change);
		Equations = Basis_Change.to_sublattice_dual(Equations);
		Inequalities = Basis_Change.to_sublattice_dual(Inequalities);
	}

	run_mode_equ_inequ(computation_type, Equations, Inequalities, Out);
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_equ_inequ( string& computation_type,const Matrix<Integer>& Equations, const Matrix<Integer>& Inequalities, Output<Integer>& Out){
	int i,j,dim=Equations.nr_of_columns();
	Lineare_Transformation<Integer> Diagonalization=Transformation(Equations);
	int rank=Diagonalization.get_rank();

	Matrix<Integer> Help=Diagonalization.get_right();
	Matrix<Integer> Ker_Basis_Transpose(dim,dim-rank);
	for (i = 1; i <= dim; i++) {
		for (j = rank+1; j <= dim; j++) {
			Ker_Basis_Transpose.write(i,j-rank,Help.read(i,j));
		}
	}
	Sublattice_Representation<Integer> Basis_Change(Ker_Basis_Transpose.transpose(),true);
	Out.compose_basis_change(Basis_Change);
	Matrix<Integer> Inequ_on_Ker = Basis_Change.to_sublattice_dual(Inequalities);

	if(computation_type!="dual"){
		if (verbose) {
			cout <<endl<< "Computing extreme rays as support hyperplanes of the dual cone:";
		}
		Full_Cone<Integer> Dual_Cone(Inequ_on_Ker);
		Dual_Cone.support_hyperplanes();
		Matrix<Integer> Extreme_Rays=Dual_Cone.read_support_hyperplanes();
		run_mode_0(computation_type, Extreme_Rays, Out);
	}
	if(computation_type=="dual"){
		dim = Inequ_on_Ker.nr_of_columns();
		Integer norm;
		vector< Integer > hyperplane;
		multimap <Integer , vector <Integer> >  Help;
		typename multimap <Integer , vector <Integer> >::const_iterator ii;
		for (i = 1; i <= Inequ_on_Ker.nr_of_rows() ; i++) {
			hyperplane=Inequ_on_Ker.read(i);
			norm=0;
			for (j = 0; j <dim; j++) {
				norm+=Iabs(hyperplane[j]);
			}
			Help.insert(pair <Integer , vector <Integer> > (norm,hyperplane));
		}
		Matrix<Integer> Equations_Ordered(Inequ_on_Ker.nr_of_rows(),dim);
		i=1;
		for (ii=Help.begin(); ii != Help.end(); ii++) {
			Equations_Ordered.write(i,(*ii).second);
			i++;
		}
		Cone_Dual_Mode<Integer> ConeDM(Equations_Ordered);
		ConeDM.hilbert_basis_dual();
		//ConeDM zu einem Full_Cone<Integer> machen
		if ( ConeDM.Generators.rank() < ConeDM.dim ) {
			Sublattice_Representation<Integer> SR(ConeDM.Generators,true);
			ConeDM.to_sublattice(SR);
			Out.compose_basis_change(SR);
		}
		Full_Cone<Integer> Result(ConeDM);
		Result.dual_mode();
		Out.set_result(Result);
		Out.cone();
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void run_mode_10( string& computation_type,const Matrix<Integer>& Binomials, Output<Integer>& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for input type 10."<<endl;
		cerr<<"The program terminates."<<endl;
		global_error_handling();
	}

	int i,j, nr_of_monoid_generators=Binomials.nr_of_columns();
	Lineare_Transformation<Integer> Diagonalization=Transformation(Binomials);
	int rank=Diagonalization.get_rank();
	Matrix<Integer> Help=Diagonalization.get_right();
	Matrix<Integer> Generators(nr_of_monoid_generators,nr_of_monoid_generators-rank);
	for (i = 1; i <= nr_of_monoid_generators; i++) {
		for (j = rank+1; j <= nr_of_monoid_generators; j++) {
			Generators.write(i,j-rank,Help.read(i,j));
		}
	}
	Full_Cone<Integer> C(Generators);
	C.support_hyperplanes();
	Matrix<Integer> Supp_Hyp=C.read_support_hyperplanes();
	Matrix<Integer> Selected_Supp_Hyp_Trans=(Supp_Hyp.submatrix(Supp_Hyp.max_rank_submatrix_lex())).transpose();
	Matrix<Integer> Positive_Embedded_Generators=Generators.multiplication(Selected_Supp_Hyp_Trans);
	Out.set_original_generators(Positive_Embedded_Generators);
	run_mode_1( computation_type, Positive_Embedded_Generators, Out);
}

//---------------------------------------------------------------------------

template<typename Integer>
Full_Cone<Integer> make_computations(const string& computation_type, const Matrix<Integer>& Full_Cone_Generators){
	Full_Cone<Integer> C(Full_Cone_Generators);
	if (computation_type=="triangulation_hilbert_basis"){
		C.triangulation_hilbert_basis();
	}
	if (computation_type=="hilbert_basis"){
		C.hilbert_basis();
	}
	if (computation_type=="support_hyperplanes"){
		C.support_hyperplanes();
	}
	if (computation_type=="support_hyperplanes_pyramid"){
		C.support_hyperplanes_pyramid();
	}
	if (computation_type=="triangulation"){
		C.support_hyperplanes_triangulation();
	}
	if (computation_type=="triangulation_pyramid"){
		C.support_hyperplanes_triangulation_pyramid();
	}
	if (computation_type=="ht1_elements"){
		C.ht1_elements();
	}
	if (computation_type=="hilbert_polynomial"){
		C.hilbert_polynomial();
	}
	if (computation_type=="hilbert_basis_polynomial"){
		C.hilbert_basis_polynomial();
	}
	return C;
}

//---------------------------------------------------------------------------

}
