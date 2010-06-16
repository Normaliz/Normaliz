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
#include <list>
#include <iostream>
#include <algorithm>
#include <map>
using namespace std;

//---------------------------------------------------------------------------

#include "mode.h"
#include "integer.h"
#include "vector_operations.h"
#include "simplex.h"
#include "lineare_transformation.h"
#include "sublattice_representation.h"
#include "full_cone.h"
#include "cone_dual_mode.h"

//---------------------------------------------------------------------------

extern void  global_error_handling();

//---------------------------------------------------------------------------

void make_main_computation(const int& mode, string& computation_type,const Matrix& Input, Output& Out){
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
		case 4: run_mode_456(computation_type, Matrix(0,dim), Matrix(0,dim), Input, Out);break;
		case 5: run_mode_456(computation_type, Matrix(0,dim), Input, Matrix(0,dim), Out);break;
		case 6: run_mode_456(computation_type, Input, Matrix(0,dim), Matrix(0,dim), Out);break;
		case 10: run_mode_10(computation_type,Input, Out);break;
	}
}

//---------------------------------------------------------------------------

void run_mode_0( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in hilbert_basis mode."<<endl;
		computation_type="hilbert_basis";
	}
	Sublattice_Representation Basis_Change(Input,true);
	Matrix Full_Cone_Generators = Basis_Change.to_sublattice(Input);
	Full_Cone Result = make_computations(computation_type, Full_Cone_Generators);

	Out.set_result(Result);
	Out.compose_basis_change(Basis_Change);
	Out.cone();
}

//---------------------------------------------------------------------------

void run_mode_1( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in computation mode hilbert_basis."<<endl;
		computation_type="hilbert_basis";
	}
	Sublattice_Representation Basis_Change(Input,false);
	Matrix Full_Cone_Generators = Basis_Change.to_sublattice(Input);
	Full_Cone Result = make_computations(computation_type, Full_Cone_Generators);

	Out.set_result(Result);
	Out.compose_basis_change(Basis_Change);
	Out.cone();
}

//---------------------------------------------------------------------------

void run_mode_2( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in computation mode hilbert_basis."<<endl;
		computation_type="hilbert_basis";
	}
	int i,j,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
	Integer number;
	Matrix Generators(nr_rows,nr_columns+1,1);
	for(i=1; i<=nr_rows; i++){
		for(j=1; j<=nr_columns; j++) {
			number=Input.read(i,j);
			Generators.write(i,j,number);
		}
	}

	Sublattice_Representation Basis_Change(Generators,true);
	Matrix Full_Cone_Generators = Basis_Change.to_sublattice(Generators);
	Full_Cone Result = make_computations(computation_type, Full_Cone_Generators);

	Out.set_result(Result);
	Out.compose_basis_change(Basis_Change);
	Out.polytop();
}

//---------------------------------------------------------------------------

void run_mode_3( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for this input type."<<endl;
		cerr<<"The program will run in computation mode hilbert_basis."<<endl;
		computation_type="hilbert_basis";
	}
	int i,j,k,l,nr_rows=Input.nr_of_rows(), nr_columns=Input.nr_of_columns();
	bool primary=true;
	Integer number;
	Matrix Full_Cone_Generators(nr_rows+nr_columns,nr_columns+1,0);
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
	Matrix Prim_Test=Input;
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
	Full_Cone Result=make_computations(computation_type, Full_Cone_Generators);
	Out.set_result(Result);
	Out.rees(primary);
}

//---------------------------------------------------------------------------

void run_mode_4( string& computation_type,const Matrix& Input, const int& nr_equations, Output& Out){
	int i;
	int dim=Input.nr_of_columns();
	Matrix Equations(nr_equations, dim);
	int nr_inequalities=Input.nr_of_rows()-nr_equations;
	Matrix Inequalities(nr_inequalities, dim);
	for(i=1;i<=nr_inequalities;i++){
		Inequalities.write(i,Input.read(i));
	}
	for(i=1;i<=nr_equations;i++){
		Equations.write(i,Input.read(i+nr_inequalities));
	}
	run_mode_equ_inequ(computation_type,Equations,Inequalities,Out);
	
}

//---------------------------------------------------------------------------

void run_mode_5( string& computation_type,const Matrix& Input, Output& Out){
	int dim=Input.nr_of_columns();
	Matrix Inequalities(dim);
	run_mode_equ_inequ(computation_type,Input,Inequalities,Out);    
}

//---------------------------------------------------------------------------

void run_mode_456(string& computation_type, const Matrix& Congruences, Matrix Equations,  Matrix Inequalities, Output& Out) {
	

	int dim = 0;
	int nr_cong = Congruences.nr_of_rows();
	if (nr_cong > 0) {	
		dim = Congruences.nr_of_columns() -1;
		if (Equations.nr_of_rows() > 0 &&  Equations.nr_of_columns() != dim) {
			cerr << "Error: dimensions of input matricies do not match!";
			global_error_handling();
		}
	} else {
		dim = Equations.nr_of_columns();
	}
	// use positive orthant if no inequalities are given
	if (Inequalities.nr_of_rows() == 0) {
		Inequalities = Matrix(Equations.nr_of_columns()); 
	}
		
	// handle Congurences
	if (nr_cong > 0) {
		int i,j;
		
		//add slack variables
		Matrix Cong_Slack(nr_cong, dim+nr_cong);
		for (i = 1; i <= nr_cong; i++) {
			for (j = 1; j <= dim; j++) {
				Cong_Slack.write(i,j,Congruences.read(i,j));
			}
			Cong_Slack.write(i,dim+i,Congruences.read(i,dim+1));
		}

		//compute kernel
		Lineare_Transformation Diagonalization = Transformation(Cong_Slack);
		int rank = Diagonalization.get_rank();
		Matrix H = Diagonalization.get_right();
		Matrix Ker_Basis_Transpose(dim, dim+nr_cong-rank);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim+nr_cong; j++) {
				Ker_Basis_Transpose.write(i,j-rank,H.read(i,j));
			}
		}

		//TODO now a new linear transformation is computed, necessary??
		Sublattice_Representation Basis_Change(Ker_Basis_Transpose.transpose(),false);
		Out.compose_basis_change(Basis_Change);
		Equations = Equations.multiplication(Ker_Basis_Transpose);
		Inequalities = Inequalities.multiplication(Ker_Basis_Transpose);

	}

	run_mode_equ_inequ(computation_type, Equations, Inequalities, Out);
}

//---------------------------------------------------------------------------

void run_mode_equ_inequ( string& computation_type,const Matrix& Equations, const Matrix& Inequalities, Output& Out){
	int i,j,dim=Equations.nr_of_columns();
	Lineare_Transformation Diagonalization=Transformation(Equations);
	int rank=Diagonalization.get_rank();

	if(computation_type!="dual"){
		Matrix Help=Diagonalization.get_right();
		Matrix Ker_Basis_Transpose(dim,dim-rank);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim; j++) {
				Ker_Basis_Transpose.write(i,j-rank,Help.read(i,j));
			}
		}
		Matrix Inequ_on_Ker=Inequalities.multiplication(Ker_Basis_Transpose);

		Full_Cone Dual_Cone(Inequ_on_Ker);
		Dual_Cone.support_hyperplanes();
		Matrix Extreme_Rays=Dual_Cone.read_support_hyperplanes();
		Matrix Ker_Basis=Ker_Basis_Transpose.transpose();
		Matrix Generators=Extreme_Rays.multiplication(Ker_Basis);
		run_mode_0( computation_type ,Generators, Out);
	}
	if(computation_type=="dual"){
		Matrix H=Diagonalization.get_right();
		Matrix H_Inv=Diagonalization.get_right_inv();
		Matrix Ker_Basis_Transpose(dim,dim-rank);
		Matrix Ker_Basis_Transpose_Inv(dim-rank,dim);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim; j++) {
				Ker_Basis_Transpose.write(i,j-rank,H.read(i,j));
				Ker_Basis_Transpose_Inv.write(j-rank,i,H_Inv.read(j,i));
			}
		}
		Diagonalization.set_rank(dim-rank);
		Diagonalization.set_right(Ker_Basis_Transpose_Inv.transpose());
		Diagonalization.set_right_inv(Ker_Basis_Transpose.transpose());
		Diagonalization.set_center(Matrix(dim-rank));
		Out.compose_basis_change(Sublattice_Representation(Diagonalization,true));
		Matrix M=Inequalities.multiplication(Ker_Basis_Transpose);
		dim=M.nr_of_columns();
		Integer norm;
		vector< Integer > hyperplane;
		multimap <Integer , vector <Integer> >  Help;
		multimap <Integer , vector <Integer> >::const_iterator ii;
		for (i = 1; i <= M.nr_of_rows() ; i++) {
			hyperplane=M.read(i);
			norm=0;
			for (j = 0; j <dim; j++) {
				norm+=Iabs(hyperplane[j]);
			}
			Help.insert(pair <Integer , vector <Integer> > (norm,hyperplane));
		}
		Matrix Equations_Ordered(M.nr_of_rows(),dim);
		i=1;
		for (ii=Help.begin(); ii != Help.end(); ii++) {
			Equations_Ordered.write(i,(*ii).second);
			i++;
		}
		Cone_Dual_Mode Cone1(Equations_Ordered);
		Cone1.hilbert_basis_dual();
		//Cone1 zu einem Full_Cone machen
		Sublattice_Representation SR(Cone1.get_generators(),true);
		if ( SR.get_rank() < Cone1.dim ) {
			Cone1.to_sublattice(SR);
			Out.compose_basis_change(SR);
		}
		Full_Cone Result(Cone1);
		Result.dual_mode();
		Out.set_result(Result);
		Out.cone();
	}
}

//---------------------------------------------------------------------------

void run_mode_10( string& computation_type,const Matrix& Binomials, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Computation mode \"dual\" not implemented for input type 10."<<endl;
		cerr<<"The program terminates."<<endl;
		global_error_handling();
	}

	int i,j, nr_of_monoid_generators=Binomials.nr_of_columns();
	Lineare_Transformation Diagonalization=Transformation(Binomials);
	int rank=Diagonalization.get_rank();
	Matrix Help=Diagonalization.get_right();
	Matrix Generators(nr_of_monoid_generators,nr_of_monoid_generators-rank);
	for (i = 1; i <= nr_of_monoid_generators; i++) {
		for (j = rank+1; j <= nr_of_monoid_generators; j++) {
			Generators.write(i,j-rank,Help.read(i,j));
		}
	}
	Full_Cone C(Generators);
	C.support_hyperplanes();
	Matrix Supp_Hyp=C.read_support_hyperplanes();
	Matrix Selected_Supp_Hyp_Trans=(Supp_Hyp.submatrix(Supp_Hyp.max_rank_submatrix_lex())).transpose();
	Matrix Positive_Embedded_Generators=Generators.multiplication(Selected_Supp_Hyp_Trans);
	run_mode_1( computation_type, Positive_Embedded_Generators, Out);
}

//---------------------------------------------------------------------------

Full_Cone make_computations(const string& computation_type, const Matrix& Full_Cone_Generators){
	Full_Cone C(Full_Cone_Generators);
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

