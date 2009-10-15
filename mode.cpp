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

//---------------------------------------------------------------------------

extern void  global_error_handling();

//---------------------------------------------------------------------------

void make_main_computation(const int& mode,  string& computation_type,const Matrix& Input, Output& Out){
	if (mode<0 || mode>5) {
		cerr<<"warning: Unknown mode "<<mode<<". The program will run in mode 0."<<endl;
		run_mode_0(computation_type,Input, Out);
		return;
	}
	switch (mode){
		case 0: run_mode_0(computation_type,Input, Out);break;
		case 1: run_mode_1(computation_type,Input, Out);break;
		case 2: run_mode_2(computation_type,Input, Out);break;
		case 3: run_mode_3(computation_type,Input, Out);break;
		case 4: run_mode_4(computation_type,Input, Out);break;
		case 5: run_mode_5(computation_type,Input, Out);break;
	}
}

//---------------------------------------------------------------------------

void run_mode_0( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"computation type = dual not implemented in mode 0."<<endl;
		cerr<<"The program will run in normal mode."<<endl;
		computation_type="normal";
	}
	int i,j,rank;
	Lineare_Transformation Basis_Change=Transformation(Input);
	rank=Basis_Change.get_rank();
	if (rank==0) {
		cerr<<"error: Input matrix has rank 0. Plese check input data."<<endl;
		global_error_handling();
	}
	Matrix V=Basis_Change.get_right();
	Matrix V_Inv=Basis_Change.get_right_inv();
	Matrix Change_To_Full_Emb(Input.nr_of_columns(),rank);
	Matrix Change_To_Full_Emb_Inv(rank,Input.nr_of_columns());
	for (i = 1; i <=Input.nr_of_columns() ; i++) {
		for (j = 1; j <= rank; j++) {
			Change_To_Full_Emb.write(i,j,V.read(i,j));
			Change_To_Full_Emb_Inv.write(j,i,V_Inv.read(j,i));
		}
	}
	Matrix Diagonal(rank);
	Basis_Change.set_left(Input);
	Basis_Change.set_center(Diagonal);
	Basis_Change.set_right(Change_To_Full_Emb);
	Basis_Change.set_right_inv(Change_To_Full_Emb_Inv);
	Matrix Full_Cone_Generators=Input.multiplication(Change_To_Full_Emb);
	Full_Cone Result=make_computations(computation_type, Full_Cone_Generators);
	Out.set_result(Result);
	Out.set_basis_change(Basis_Change);
	Out.cone();
}

//---------------------------------------------------------------------------

void run_mode_1( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Computation type = dual not implemented in mode 1."<<endl;
		cerr<<"The program will run in Computation type normal."<<endl;
		computation_type="normal";
	}
	int i,j,rank;
	Lineare_Transformation Basis_Change=Transformation(Input);
	rank=Basis_Change.get_rank();
	if (rank==0) {
		cerr<<"error: Input matrix has rank 0. Please check input data."<<endl;
		global_error_handling();
	}
	Matrix V=Basis_Change.get_right();
	Matrix V_Inv=Basis_Change.get_right_inv();
	Matrix Change_To_Full_Emb(Input.nr_of_columns(),rank);
	Matrix Change_To_Full_Emb_Inv(rank,Input.nr_of_columns());
	for (i = 1; i <=Input.nr_of_columns() ; i++) {
		for (j = 1; j <= rank; j++) {
			Change_To_Full_Emb.write(i,j,V.read(i,j));
			Change_To_Full_Emb_Inv.write(j,i,V_Inv.read(j,i));
		}
	}
	Matrix D=Basis_Change.get_center();
	Matrix Diagonal(rank);
	for (i = 1; i <= rank; i++) {
		Diagonal.write(i,i,D.read(i,i));
	}
	Basis_Change.set_left(Input);
	Basis_Change.set_center(Diagonal);
	Basis_Change.set_right(Change_To_Full_Emb);
	Basis_Change.set_right_inv(Change_To_Full_Emb_Inv);
	Matrix Full_Cone_Generators=Input.multiplication(Change_To_Full_Emb);
	Full_Cone_Generators=Full_Cone_Generators.transpose();
	vector<Integer> v;
	for (i = 1; i <= rank; i++) {
		v=Full_Cone_Generators.read(i);
		v_scalar_division(v,Diagonal.read(i,i));
		Full_Cone_Generators.write(i,v);
	}
	Full_Cone_Generators=Full_Cone_Generators.transpose();
	Full_Cone Result=make_computations(computation_type, Full_Cone_Generators);
	Out.set_result(Result);
	Out.set_basis_change(Basis_Change);
	Out.cone();
}

//---------------------------------------------------------------------------

void run_mode_2( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Run mode type = dual not implemented in mode 2."<<endl;
		cerr<<"The program will run in normal mode."<<endl;
		computation_type="normal";
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
	Lineare_Transformation Basis_Change=Transformation(Generators);
	int rank=Basis_Change.get_rank();
	if (rank==0) {
		cerr<<"error: Input matrix has rank 0. Plese check input data."<<endl;
		global_error_handling();
	}
	Matrix V=Basis_Change.get_right();
	Matrix V_Inv=Basis_Change.get_right_inv();
	Matrix Change_To_Full_Emb(Generators.nr_of_columns(),rank);
	Matrix Change_To_Full_Emb_Inv(rank,Generators.nr_of_columns());
	for (i = 1; i <=Generators.nr_of_columns() ; i++) {
		for (j = 1; j <= rank; j++) {
			Change_To_Full_Emb.write(i,j,V.read(i,j));
			Change_To_Full_Emb_Inv.write(j,i,V_Inv.read(j,i));
		}
	}
	Matrix Diagonal(rank);
	Basis_Change.set_left(Input);
	Basis_Change.set_center(Diagonal);
	Basis_Change.set_right(Change_To_Full_Emb);
	Basis_Change.set_right_inv(Change_To_Full_Emb_Inv);
	Matrix Full_Cone_Generators=Generators.multiplication(Change_To_Full_Emb);
	Full_Cone Result=make_computations(computation_type, Full_Cone_Generators);
	Out.set_result(Result);
	Out.set_basis_change(Basis_Change);
	Out.polytop();
}

//---------------------------------------------------------------------------

void run_mode_3( string& computation_type,const Matrix& Input, Output& Out){
	if (computation_type=="dual") {
		cerr<<"Run mode type = dual not implemented in mode 3."<<endl;
		cerr<<"The program will run in normal mode."<<endl;
		computation_type="normal";
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

void run_mode_4( string& computation_type,const Matrix& Input, Output& Out){
	if (Input.rank()!=Input.nr_of_columns() ) {
		cerr<<"error: Rank is not maximal. In mode 4 the input matrix must be of maximal rank.";
		global_error_handling();
	}
	if(computation_type!="dual"){
		Full_Cone Help(Input);
		Help.support_hyperplanes();
		Matrix Generators=Help.read_support_hyperplanes();
		if (Generators.nr_of_rows()==0) {
			Matrix Trivial_Solution(1,Input.nr_of_columns(),0);
			Generators=Trivial_Solution;
			cerr<<"warning: The only solution of the system is 0.";
			global_error_handling();
		}
		run_mode_0( computation_type ,Generators, Out);
	}
	if(computation_type=="dual"){
		int i,j, dim=Input.nr_of_columns();
		Integer norm;
		vector< Integer > hyperplane;
		multimap <Integer , vector <Integer> >  Help;
		multimap <Integer , vector <Integer> >::const_iterator ii;
		for (i = 1; i <= Input.nr_of_rows() ; i++) {
			hyperplane=Input.read(i);
			norm=0;
			for (j = 0; j <dim; j++) {
				norm+=Iabs(hyperplane[j]);
			}
			Help.insert(pair <Integer , vector <Integer> > (norm,hyperplane));
		}
		Matrix Input_Ordered(Input.nr_of_rows(),dim);
		i=1;
		for (ii=Help.begin(); ii != Help.end(); ii++) {
			Input_Ordered.write(i,(*ii).second);
			i++;
		}
		Full_Cone Result(Input_Ordered); //in the mode dual the support hyperplanes are the generators
		//and as support hyperplanes we recover the generators
		Result.hilbert_basis_dual();
		Out.set_result(Result);
		Matrix I(dim); //identity matrix
		Lineare_Transformation Diagonalization=Transformation(I);
		Out.set_basis_change(Diagonalization);
		Out.dual();
	}
}

//---------------------------------------------------------------------------

void run_mode_5( string& computation_type,const Matrix& Input, Output& Out){
	if(computation_type!="dual"){
		int i,j,dim=Input.nr_of_columns();
		Lineare_Transformation Diagonalization=Transformation(Input);
		int rank=Diagonalization.get_rank();
		Matrix Help=Diagonalization.get_right();
		Matrix Ker_Basis_Transpose(dim,dim-rank);
		for (i = 1; i <= dim; i++) {
			for (j = rank+1; j <= dim; j++) {
				Ker_Basis_Transpose.write(i,j-rank,Help.read(i,j));
			}
		}
		Full_Cone Help_Cone(Ker_Basis_Transpose);
		Help_Cone.support_hyperplanes();
		Matrix Support_Hyperplanes=Help_Cone.read_support_hyperplanes();
		Matrix Ker_Basis=Ker_Basis_Transpose.transpose();
		Matrix Generators=Support_Hyperplanes.multiplication(Ker_Basis);
		if (Generators.nr_of_rows()==0) {
			Matrix Trivial_Solution(1,dim,0);
			Generators=Trivial_Solution;
			cerr<<"warning: The only solution of the system is 0.";
			global_error_handling();
		}
		run_mode_0( computation_type ,Generators, Out);
	}
	if(computation_type=="dual"){
		int i,j,dim=Input.nr_of_columns();
		Lineare_Transformation Diagonalization=Transformation(Input);
		int rank=Diagonalization.get_rank();
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
		Diagonalization.set_right(Ker_Basis_Transpose.transpose());
		Diagonalization.set_right_inv(Ker_Basis_Transpose_Inv.transpose());
		Out.set_basis_change(Diagonalization);
		Matrix M=Ker_Basis_Transpose;
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
		Matrix Input_Ordered(M.nr_of_rows(),dim);
		i=1;
		for (ii=Help.begin(); ii != Help.end(); ii++) {
			Input_Ordered.write(i,(*ii).second);
			i++;
		}
		Full_Cone Result(Input_Ordered); //in the mode dual the support hyperplanes are the generators
		//and as support hyperplanes we recover the generators
		Result.hilbert_basis_dual();
		Out.set_result(Result);
		Out.dual();
	}
}

//---------------------------------------------------------------------------

Full_Cone make_computations(const string& computation_type, const Matrix& Full_Cone_Generators){
	Full_Cone C(Full_Cone_Generators);
	if (computation_type=="normal"){
		C.hilbert_basis();
	}
	if (computation_type=="normal_compressed"){
		C.hilbert_basis(true);
	}
	if (computation_type=="support_hyperplanes"){
		C.support_hyperplanes();
	}
	if (computation_type=="triangulation"){
		C.support_hyperplanes_triangulation_multiplicity();
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



