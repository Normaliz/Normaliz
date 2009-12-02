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
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

//---------------------------------------------------------------------------

#include "output.h"
#include "integer.h"
#include "vector_operations.h"
#include "matrix.h"
#include "simplex.h"

//---------------------------------------------------------------------------

extern void global_error_handling();

//---------------------------------------------------------------------------

Output::Output(){
	out=true;
	inv=false;
	ext=false;
	esp=false;
	typ=false;
	egn=false;
	gen=false;
	sup=false;
	tri=false;
	ht1=false;
}

//---------------------------------------------------------------------------

Output::Output(const Output& Out){
	name=Out.name;
	out=Out.out;
	inv=Out.inv;
	ext=Out.ext;
	esp=Out.esp;
	typ=Out.typ;
	egn=Out.egn;
	gen=Out.gen;
	sup=Out.sup;
	tri=Out.tri;
	ht1=Out.ht1;
	Result=Out.Result;
	Basis_Change=Out.Basis_Change;
}

//---------------------------------------------------------------------------

Output::~Output(){
//automatic destructor
}

//---------------------------------------------------------------------------

void Output::read() const{
	cout<<"\nname="<<name<<"\n";
	cout<<"\nout="<<out<<"\n";
	cout<<"\ninv="<<inv<<"\n";
	cout<<"\next="<<ext<<"\n";
	cout<<"\nesp="<<esp<<"\n";
	cout<<"\ntyp="<<typ<<"\n";
	cout<<"\negn="<<egn<<"\n";
	cout<<"\ngen="<<gen<<"\n";
	cout<<"\nsup="<<sup<<"\n";
	cout<<"\ntri="<<tri<<"\n";
	cout<<"\nht1="<<ht1<<"\n";
	cout<<"\nResult is:\n";
	Result.read();
	cout<<"\nBasis Change is:\n";
	Basis_Change.read();
}

//---------------------------------------------------------------------------

void Output::set_name(const string& n){
	name=n;
}

//---------------------------------------------------------------------------

void Output::set_write_out(const bool& flag){
	out=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_inv(const bool& flag){
	inv=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_ext(const bool& flag){
	ext=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_esp(const bool& flag){
	esp=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_typ(const bool& flag){
	typ=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_egn(const bool& flag){
	egn=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_gen(const bool& flag){
	gen=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_sup(const bool& flag){
	sup=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_tri(const bool& flag) {
	tri=flag;
}

//---------------------------------------------------------------------------

void Output::set_write_ht1(const bool& flag) {
	ht1=flag;
}


//---------------------------------------------------------------------------

void Output::set_write_extra_files(){
	out=true;
	inv=true;
	ext=false;
	esp=false;
	typ=true;
	egn=false;
	gen=true;
	sup=true;
	tri=false;
	ht1=false;
}

//---------------------------------------------------------------------------

void Output::set_write_all_files(){
	out=true;
	inv=true;
	ext=true;
	esp=true;
	typ=true;
	egn=true;
	gen=true;
	sup=true;
	tri=true;
	ht1=true;
}
 
//---------------------------------------------------------------------------

void Output::set_result(const Full_Cone& C){
	Result=C;
}

//---------------------------------------------------------------------------

void Output::set_basis_change(const Lineare_Transformation& BC){
	Basis_Change=BC;
}

//---------------------------------------------------------------------------

void Output::write_matrix_ext(const Matrix& M) const{
	if (ext==true) {
	string file_name=name+".ext";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}



void Output::write_matrix_ext_1(const Matrix& M) const{
	if (ext==true) {
	string file_name=name+".ext";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<"1 "<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_esp(const Matrix& M) const{
	if (esp==true) {
	string file_name=name+".esp";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_typ(const Matrix& M) const{
	if (typ==true) {
	string file_name=name+".typ";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_egn(const Matrix& M) const {
	if (egn==true) {
	string file_name=name+".egn";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_gen(const Matrix& M) const {
	if (gen==true) {
	string file_name=name+".gen";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_sup(const Matrix& M) const{
	if (sup==true) {
	string file_name=name+".sup";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_tri(const Matrix& M) const{
	if (tri==true) {
	string file_name=name+".tri";
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	int	nr=M.nr_of_rows();
	int	nc=M.nr_of_columns();
	out<<nr<<endl<<nc<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			out<<M.read(i,j)<<" ";
			}
		out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::write_matrix_ht1(const Matrix& M) const{
	if (ht1==true) {
		string file_name=name+".ht1";
		const char* file=file_name.c_str();
		ofstream out(file);
		int i,j;
		int	nr=M.nr_of_rows();
		int	nc=M.nr_of_columns();
		out<<nr<<endl<<nc<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j <= nc; j++) {
				out<<M.read(i,j)<<" ";
			}
			out<<endl;
		}
		out.close();
	}
}

//---------------------------------------------------------------------------

void Output::cone()const{
	int i,j,k,nr,nc,rank=Basis_Change.get_rank(),max_decimal_length;    //read local data
	Integer buf;
	string status=Result.read_status();
	Matrix Generators=Basis_Change.get_left();
	Matrix Diagonal=Basis_Change.get_center();
	Matrix Change=Basis_Change.get_right();
	Matrix Change_Inv=Basis_Change.get_right_inv();
	Matrix Support_Hyperplanes_Full_Cone=Result.read_support_hyperplanes();
	
	write_matrix_esp(Support_Hyperplanes_Full_Cone);         //write the suport hyperplanes of the full dimensional cone
	if (status!="support hyperplanes"){           			 //write triangulation
		Matrix T=Result.read_triangulation_volume();
		write_matrix_tri(T);
	}
	
	
	if (out==true) {//printing .out file
		string name_open=name+".out"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream out(file);
		
		Matrix Hilbert_Basis;                                            //write Hilbert Basis
		if (status=="normal"||status=="hilbert basis polynomial") {
			Matrix Hilbert_Basis_Full_Cone=Result.read_hilbert_basis();
			write_matrix_egn(Hilbert_Basis_Full_Cone);
			Matrix V=Hilbert_Basis_Full_Cone.multiplication(Support_Hyperplanes_Full_Cone.transpose());
			write_matrix_typ(V);
			Hilbert_Basis_Full_Cone=Hilbert_Basis_Full_Cone.multiplication(Diagonal);
			Hilbert_Basis=Hilbert_Basis_Full_Cone.multiplication(Change_Inv);
			write_matrix_gen(Hilbert_Basis);
			nr=Hilbert_Basis.nr_of_rows();
			nc=Hilbert_Basis.nr_of_columns();
			max_decimal_length=Hilbert_Basis.maximal_decimal_length();
			out<<nr<<" generators of integral closure:"<<endl;
			for (i = 1; i <= nr; i++) {
				for (j = 1; j <= nc; j++) {
					buf = Hilbert_Basis.read(i,j);
					for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
						out<<" ";
					}
					out<<buf;
				}
				out<<endl;
			}
			out<<endl;
		}
		
		vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();          //write extreme rays
		int nr_ex_rays=0;
		for (i = 0; i <Ex_Rays_Marked.size(); i++) {
			 if (Ex_Rays_Marked[i]==true) {
				nr_ex_rays++;
			 }
		}
		vector<int> Ex_Rays_Position(nr_ex_rays);
		j=0;
		for (i = 0; i <Ex_Rays_Marked.size(); i++) {
			if (Ex_Rays_Marked[i]==true) {
				Ex_Rays_Position[j]=i+1;
				j++;
			}
		}
		Matrix Extreme_Rays=Generators.submatrix(Ex_Rays_Position);
		write_matrix_ext(Extreme_Rays);
		out<<nr_ex_rays<<" extreme rays:"<<endl;
		nr=Extreme_Rays.nr_of_rows();
		nc=Extreme_Rays.nr_of_columns();
		max_decimal_length=Extreme_Rays.maximal_decimal_length();
		for (i = 1; i <=nr; i++) {
			for (j = 1; j <=nc; j++) {
				buf = Extreme_Rays.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
			out<<endl;
		}
		out<<endl;
		
		if (rank==Change_Inv.nr_of_columns()){                   //write rank and index
			out<<"(original) semigroup has rank "<<rank<<" (maximal)"<<endl;
		}
		else {
			out<<"(original) semigroup has rank "<<rank<<endl;
		}
		Integer index=Basis_Change.get_index();
		Integer index_reduce=1;          //makes a diference between modes, in mode 1 index is 1 by default 
		for (i = 1; i <= rank; i++) {
			index_reduce*=Diagonal.read(i,i);
		}
		out<<"(original) semigroup is of index "<<index/index_reduce<<endl;
		out<<endl;
		
		Change=Change.transpose();  // preparing the matrix for transformationen on the dual space
		vector<Integer> v;
		Integer m;
		for (i = 1; i <= rank; i++) {
			m=index/Diagonal.read(i,i);
			v=Change.read(i);
			v_scalar_multiplication(v,m);
			Change.write(i,v);
		}
		
		if (rank==Change_Inv.nr_of_columns()) {                  //write the support hyperplanes
			Matrix Support_Hyperplanes=Support_Hyperplanes_Full_Cone.multiplication(Change);
			Support_Hyperplanes.make_prime();
			write_matrix_sup(Support_Hyperplanes);
			nr=Support_Hyperplanes.nr_of_rows();
			nc=Support_Hyperplanes.nr_of_columns();
			max_decimal_length=Support_Hyperplanes.maximal_decimal_length();
			out<<nr<<" support hyperplanes:"<<endl;
			for (i = 1; i <= nr; i++) {
				for (j = 1; j <= nc; j++) {
					buf = Support_Hyperplanes.read(i,j);
					for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
						out<<" ";
					}
					out<<buf;
				}
			out<<endl;
			}
		out<<endl;
		}
		
		if (Result.read_homogeneous()==false) {
			out<<"(original) semigroup is not homogeneous"<<endl;
		}
		else {
			if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
				Matrix Hom=Result.read_homogeneous_elements();
				Hom=Hom.multiplication(Diagonal);
				Hom=Hom.multiplication(Change_Inv);
				write_matrix_ht1(Hom);
				nr=Hom.nr_of_rows();
				nc=Hom.nr_of_columns();
				max_decimal_length=Hom.maximal_decimal_length();
				out<<nr<<" height 1 generators of integral closure:"<<endl;
				for (i = 1; i <= nr; i++) {
					for (j = 1; j <= nc; j++) {
						buf = Hom.read(i,j);
						for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
							out<<" ";
						}
						out<<buf;
					}
				out<<endl;
				}
				out<<endl;
		    }
			vector<Integer> Linear_Form=Result.read_linear_form();
			Linear_Form=Change.VxM(Linear_Form);
			Linear_Form=v_make_prime(Linear_Form);
			out<<"(original) semigroup is homogeneous via the linear form:"<<endl;
			for (i = 0; i < Linear_Form.size(); i++) {
				out<<Linear_Form[i]<<" ";
			}
			out<<endl<<endl;
			if (status!="support hyperplanes"){
				out<<"multiplicity = "<<Result.read_multiplicity()<<endl;
			}
			out<<endl;
			if (status=="hilbert polynomial"||status=="hilbert basis polynomial") {
				vector<Integer> h_vector=Result.read_h_vector();
				out<<"h-vector = ";
				for (i = 0; i < h_vector.size(); i++) {
					out<<h_vector[i]<<" ";
				}
				out<<endl<<endl;
				vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
				out<<"Hilbert polynomial : ";
				for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
					out<<hilbert_polynomial[i]<<"/"<<hilbert_polynomial[i+1]<<" ";
				}
			}
		}
		out.close();
	}
	
	
	
	if (inv==true) {//printing .inv file
	string name_open=name+".inv"; 							 //preparing output files
	const char* file=name_open.c_str();
	ofstream inv(file);
	
	Matrix Hilbert_Basis;                                            //write Hilbert Basis
	if (status=="normal"||status=="hilbert basis polynomial") {
		Matrix Hilbert_Basis_Full_Cone=Result.read_hilbert_basis();
		nr=Hilbert_Basis_Full_Cone.nr_of_rows();
		inv<<"integer hilbert_basis_elements = "<<nr<<endl;
	}
	
	vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();          //write extreme rays
	int nr_ex_rays=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			nr_ex_rays++;
		 }
	}
	inv<<"integer number_extreme_rays = "<<nr_ex_rays<<endl;
	inv<<"integer rank = "<<rank<<endl;
	Integer index=Basis_Change.get_index();
	Integer index_reduce=1;          //makes a diference between modes, in mode 1 index is 1 by default
	for (i = 1; i <= rank; i++) {
		index_reduce*=Diagonal.read(i,i);
	}
	inv<<"integer index = "<<index/index_reduce<<endl;
	inv<<"integer number_support_hyperplanes = "<<Support_Hyperplanes_Full_Cone.nr_of_rows()<<endl;
	Change=Basis_Change.get_right();
	Change=Change.transpose();  // preparing the matrix for transformationen on the dual space
	vector<Integer> v;
	Integer m;
	for (i = 1; i <= rank; i++) {
		m=index/Diagonal.read(i,i);
		v=Change.read(i);
		v_scalar_multiplication(v,m);
		Change.write(i,v);
	}
	
	if (Result.read_homogeneous()==false) {
		inv<<"boolean homogeneous = "<<"false"<<endl;
	}
	else {
		inv<<"boolean homogeneous = "<<"true"<<endl;
		if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
			Matrix Hom=Result.read_homogeneous_elements();
			nr=Hom.nr_of_rows();
			inv<<"integer height_1_elements = "<<nr<<endl;
		}
		vector<Integer> Linear_Form=Result.read_linear_form();
		Linear_Form=Change.VxM(Linear_Form);
		Linear_Form=v_make_prime(Linear_Form);
		inv<<"vector "<<Linear_Form.size()<<" homogeneous_weights = ";
		for (i = 0; i < Linear_Form.size(); i++) {
			inv<<Linear_Form[i]<<" ";
		}
		inv<<endl;
		if (status!="support hyperplanes"){
			inv<<"integer multiplicity = "<<Result.read_multiplicity()<<endl;
		}
		if (status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		vector<Integer> h_vector=Result.read_h_vector();
		inv<<"vector "<<h_vector.size()<<" h-vector = ";
		for (i = 0; i < h_vector.size(); i++) {
			inv<<h_vector[i]<<" ";
		}
		inv<<endl;
		Integer factorial=1;
		for (i = 2; i <rank; i++) {
			 factorial*=i;
		}
		vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
		inv<<"vector "<<h_vector.size()<<" hilbert_polynomial = ";
		for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
			inv<<hilbert_polynomial[i]*(factorial /hilbert_polynomial[i+1])<<" ";
		}
		inv<<endl;
		}
	}
	inv.close();
	}
}

//---------------------------------------------------------------------------

void Output::polytop()const{
	int i,j,k,nr,nc,rank=Basis_Change.get_rank(),max_decimal_length;    //read local data
	Integer buf;
	string status=Result.read_status();
	Matrix Generators=Basis_Change.get_left();
	Matrix Change=Basis_Change.get_right();
	Matrix Change_Inv=Basis_Change.get_right_inv();
	Matrix Support_Hyperplanes_Full_Cone=Result.read_support_hyperplanes();
	
	write_matrix_esp(Support_Hyperplanes_Full_Cone);         //write the suport hyperplanes of the full dimensional cone
	if (status!="support hyperplanes"){           			 //write triangulation
		Matrix T=Result.read_triangulation_volume();
		write_matrix_tri(T);
	}
	
	
	
	if (out==true ) {//printing .out file
	string name_open=name+".out"; 							 //preparing output files
	const char* file=name_open.c_str();
	ofstream out(file);
	
	Matrix Hilbert_Basis;                                            //write Hilbert Basis
	if (status=="normal"||status=="hilbert basis polynomial") {
		Matrix Hilbert_Basis_Full_Cone=Result.read_hilbert_basis();
		write_matrix_egn(Hilbert_Basis_Full_Cone);
		Matrix V=Hilbert_Basis_Full_Cone.multiplication(Support_Hyperplanes_Full_Cone.transpose());
		write_matrix_typ(V);
		Hilbert_Basis=Hilbert_Basis_Full_Cone.multiplication(Change_Inv);
		write_matrix_gen(Hilbert_Basis);
		nr=Hilbert_Basis.nr_of_rows();
		nc=Hilbert_Basis.nr_of_columns();
		max_decimal_length=Hilbert_Basis.maximal_decimal_length();
		out<<nr<<" generators of Ehrhart ring:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j <= nc; j++) {
				buf = Hilbert_Basis.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
		out<<endl;
		}
	out<<endl;
	}
	if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		Matrix Lattice_Points=Result.read_homogeneous_elements();
		Lattice_Points=Lattice_Points.multiplication(Change_Inv);
		nr=Lattice_Points.nr_of_rows();
		nc=Lattice_Points.nr_of_columns();
		max_decimal_length=Lattice_Points.maximal_decimal_length();
		out<<nr<<" lattice points in polytope:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j < nc; j++) {
				buf = Lattice_Points.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
		out<<endl;
		}
		out<<endl;
	}
	
	vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();          //write extreme rays
	int nr_ex_rays=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			nr_ex_rays++;
		 }
	}
	vector<int> Ex_Rays_Position(nr_ex_rays);
	j=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			Ex_Rays_Position[j]=i+1;
			j++;
		 }
	}
	Matrix Extreme_Rays=Generators.submatrix(Ex_Rays_Position);
	write_matrix_ext_1(Extreme_Rays);
	out<<nr_ex_rays<<" extreme points of polytope:"<<endl;
	nr=Extreme_Rays.nr_of_rows();
	nc=Extreme_Rays.nr_of_columns();
	max_decimal_length=Extreme_Rays.maximal_decimal_length();
	for (i = 1; i <=nr; i++) {
		for (j = 1; j <=nc; j++) {
				buf = Extreme_Rays.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
		}
		out<<endl;
	}
	out<<endl;
	
	Change=Change.transpose();  // preparing the matrix for transformationen on the dual space
	if (rank==Change_Inv.nr_of_columns()) {                  //write the support hyperplanes
		Matrix Support_Hyperplanes=Support_Hyperplanes_Full_Cone.multiplication(Change);
		Support_Hyperplanes.make_prime();
		write_matrix_sup(Support_Hyperplanes);
		nr=Support_Hyperplanes.nr_of_rows();
		nc=Support_Hyperplanes.nr_of_columns();
		for (i = 1; i <= nr; i++) {
	    	Support_Hyperplanes.write(i,nc,-Support_Hyperplanes.read(i,nc));
		}
		max_decimal_length=Support_Hyperplanes.maximal_decimal_length();
		out<<nr<<" support hyperplanes:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j < nc; j++) {
				buf = Support_Hyperplanes.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
			out<<" >=";
			buf = Support_Hyperplanes.read(i,j);
			for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
			}
			out<<buf;
			out<<endl;
		}
	out<<endl;
	}
	if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		Matrix Hom=Result.read_homogeneous_elements();
		Hom=Hom.multiplication(Change_Inv);
		write_matrix_ht1(Hom);
	}
	if (status!="support hyperplanes"){
		out<<"normalized volume = "<<Result.read_multiplicity()<<endl;
	}
	out<<endl;
	if (status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		vector<Integer> h_vector=Result.read_h_vector();
		out<<"h-vector = ";
		for (i = 0; i < h_vector.size(); i++) {
			out<<h_vector[i]<<" ";
		}
		out<<endl<<endl;
		vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
		out<<"Ehrhart polynomial : ";
		for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
			out<<hilbert_polynomial[i]<<"/"<<hilbert_polynomial[i+1]<<" ";
		}
	}
	out.close();
	}
	
	
	
	
	if (inv==true) {//printing .inv file
	string name_open=name+".inv"; 							 //preparing output files
	const char* file=name_open.c_str();
	ofstream inv(file);
	
	Matrix Hilbert_Basis;                                            //write Hilbert Basis
	if (status=="normal"||status=="hilbert basis polynomial") {
		Matrix Hilbert_Basis_Full_Cone=Result.read_hilbert_basis();
		nr=Hilbert_Basis_Full_Cone.nr_of_rows();
		inv<<"integer hilbert_basis_elements = "<<nr<<endl;
	}
	
	vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();          //write extreme rays
	int nr_ex_rays=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			nr_ex_rays++;
		 }
	}
	inv<<"integer number_extreme_rays = "<<nr_ex_rays<<endl;
	inv<<"integer rank = "<<rank<<endl;
	inv<<"integer index = "<<1<<endl;
	inv<<"integer number_support_hyperplanes = "<<Support_Hyperplanes_Full_Cone.nr_of_rows()<<endl;
	Change=Basis_Change.get_right();
	Change=Change.transpose();  // preparing the matrix for transformationen on the dual space
	
	if (Result.read_homogeneous()==false) {
		inv<<"boolean homogeneous = "<<"false"<<endl;
	}
	else {
		inv<<"boolean homogeneous = "<<"true"<<endl;
		if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		Matrix Hom=Result.read_homogeneous_elements();
		nr=Hom.nr_of_rows();
		inv<<"integer height_1_elements = "<<nr<<endl;
		}
		vector<Integer> Linear_Form=Result.read_linear_form();
		Linear_Form=Change.VxM(Linear_Form);
		Linear_Form=v_make_prime(Linear_Form);
		inv<<"vector "<<Linear_Form.size()<<" homogeneous_weights = ";
		for (i = 0; i < Linear_Form.size(); i++) {
			inv<<Linear_Form[i]<<" ";
		}
		inv<<endl;
		if (status!="support hyperplanes"){
			inv<<"integer multiplicity = "<<Result.read_multiplicity()<<endl;
		}
		if (status=="hilbert polynomial"||status=="hilbert basis polynomial") {
			vector<Integer> h_vector=Result.read_h_vector();
			inv<<"vector "<<h_vector.size()<<" h-vector = ";
			for (i = 0; i < h_vector.size(); i++) {
				inv<<h_vector[i]<<" ";
			}
			inv<<endl;
			Integer factorial=1;
			for (i = 2; i <rank; i++) {
				 factorial*=i;
			}
			vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
			inv<<"vector "<<h_vector.size()<<" hilbert_polynomial = ";
			for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
				inv<<hilbert_polynomial[i]*(factorial /hilbert_polynomial[i+1])<<" ";
			}
			inv<<endl;
		}
	}
	inv.close();
	}
}

//---------------------------------------------------------------------------

void Output::rees(const bool primary) const{
	int i,j,k,nr,nc,max_decimal_length;    //read local data
	Integer buf;
	string status=Result.read_status();
	Matrix Generators=Result.read_generators();
	Matrix Support_Hyperplanes=Result.read_support_hyperplanes();
	int nr_generators_ideal=0;
	
	write_matrix_esp(Support_Hyperplanes);         //write the suport hyperplanes of the full dimensional cone
	if (status!="support hyperplanes"){           			 //write triangulation
		Matrix T=Result.read_triangulation_volume();
		write_matrix_tri(T);
	}
	
	
	if (out==true) { //printing .out file
	string name_open=name+".out"; 							 //preparing output files
	const char* file=name_open.c_str();
	ofstream out(file);
	
	Matrix Hilbert_Basis;                                            //write Hilbert Basis
	if (status=="normal"||status=="hilbert basis polynomial") {
		Hilbert_Basis=Result.read_hilbert_basis();
		write_matrix_egn(Hilbert_Basis);
		Matrix V=Hilbert_Basis.multiplication(Support_Hyperplanes.transpose());
		write_matrix_typ(V);
		write_matrix_gen(Hilbert_Basis);
		nr=Hilbert_Basis.nr_of_rows();
		nc=Hilbert_Basis.nr_of_columns();
		max_decimal_length=Hilbert_Basis.maximal_decimal_length();
		out<<nr<<" generators of integral closure of the Rees algebra:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j <= nc; j++) {
				buf = Hilbert_Basis.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
				if (j==nc&&buf==1) {
					nr_generators_ideal++;
				}
			}
		out<<endl;
		}
	out<<endl;
	}
	
	vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();          //write extreme rays
	int nr_ex_rays=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			nr_ex_rays++;
		 }
	}
	vector<int> Ex_Rays_Position(nr_ex_rays);
	j=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			Ex_Rays_Position[j]=i+1;
			j++;
		 }
	}
	Matrix Extreme_Rays=Generators.submatrix(Ex_Rays_Position);
	write_matrix_ext(Extreme_Rays);
	out<<nr_ex_rays<<" extreme rays:"<<endl;
	nr=Extreme_Rays.nr_of_rows();
	nc=Extreme_Rays.nr_of_columns();
	max_decimal_length=Extreme_Rays.maximal_decimal_length();
	for (i = 1; i <=nr; i++) {
		for (j = 1; j <=nc; j++) {
				buf = Extreme_Rays.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
		}
		out<<endl;
	}
	out<<endl;
	
	if (status=="normal"||status=="hilbert basis polynomial") {
		nr=Hilbert_Basis.nr_of_rows();
		nc=Hilbert_Basis.nr_of_columns();
		max_decimal_length=Hilbert_Basis.maximal_decimal_length();
		out<<nr_generators_ideal<<" generators of integral closure of the ideal:"<<endl;
		for (i = 1; i <= nr; i++) {
			if (Hilbert_Basis.read(i,nc)==1) {
				for (j = 1; j < nc; j++) {
					buf = Hilbert_Basis.read(i,j);
					for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
						out<<" ";
					}
					out<<buf;
				}
			out<<endl;
			}
		}
	out<<endl;
	}
	
	write_matrix_sup(Support_Hyperplanes);    //write the support hyperplanes
	nr=Support_Hyperplanes.nr_of_rows();
	nc=Support_Hyperplanes.nr_of_columns();
	max_decimal_length=Support_Hyperplanes.maximal_decimal_length();
	out<<nr<<" support hyperplanes:"<<endl;
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			buf = Support_Hyperplanes.read(i,j);
			for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
				out<<" ";
			}
			out<<buf;
		}
	out<<endl;
	}
	out<<endl;
	
	if (Result.read_homogeneous()==false) {
		out<<"(original) semigroup is not homogeneous"<<endl;
	}
	else {
		if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		Matrix Hom=Result.read_homogeneous_elements();
		write_matrix_ht1(Hom);
		nr=Hom.nr_of_rows();
		nc=Hom.nr_of_columns();
		max_decimal_length=Hom.maximal_decimal_length();
		out<<nr<<" height 1 generators of integral closure of the Rees algebra:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j <= nc; j++) {
				buf = Hom.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
		out<<endl;
		}
		out<<endl;
		}
		vector<Integer> Linear_Form=Result.read_linear_form();
		out<<"(original) semigroup is homogeneous via the linear form:"<<endl;
		for (i = 0; i < Linear_Form.size(); i++) {
			out<<Linear_Form[i]<<" ";
		}
		out<<endl<<endl;
		if (status!="support hyperplanes"){
			out<<"multiplicity = "<<Result.read_multiplicity()<<endl;
		}
		out<<endl;
		if (status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		vector<Integer> h_vector=Result.read_h_vector();
		out<<"h-vector = ";
		for (i = 0; i < h_vector.size(); i++) {
			out<<h_vector[i]<<" ";
		}
		out<<endl<<endl;
		vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
		out<<"Hilbert polynomial : ";
		for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
			out<<hilbert_polynomial[i]<<"/"<<hilbert_polynomial[i+1]<<" ";
		}
		}
	out<<endl;
	}
	if (primary) {
		out<<"ideal is primary to the ideal generated by the indeterminates"<<endl;
		Integer primary_multiplicity=Result.primary_multiplicity();
		out<<"multiplicity of the ideal = "<<primary_multiplicity<<endl;
	}
	else {
		out<<"ideal is not primary to the ideal generated by the indeterminates"<<endl;
	}
	
	out.close();
	}
	
	if (inv==true) {//printing .inv file
	string name_open=name+".inv"; 							 //preparing output files
	const char* file=name_open.c_str();
	ofstream inv(file);
	
	Matrix Hilbert_Basis;                                            //write Hilbert Basis
	if (status=="normal"||status=="hilbert basis polynomial") {
		Matrix Hilbert_Basis_Full_Cone=Result.read_hilbert_basis();
		nr=Hilbert_Basis_Full_Cone.nr_of_rows();
		inv<<"integer hilbert_basis_elements = "<<nr<<endl;
	}
	
	vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();          //write extreme rays
	int nr_ex_rays=0;
	for (i = 0; i <Ex_Rays_Marked.size(); i++) {
		 if (Ex_Rays_Marked[i]==true) {
			nr_ex_rays++;
		 }
	}
	inv<<"integer number_extreme_rays = "<<nr_ex_rays<<endl;
	int rank=Result.read_dimension();
	inv<<"integer rank = "<<rank<<endl;
	inv<<"integer index = "<<1<<endl;
	Matrix Support_Hyperplanes_Full_Cone=Result.read_support_hyperplanes();
	inv<<"integer number_support_hyperplanes = "<<Support_Hyperplanes_Full_Cone.nr_of_rows()<<endl;
	
	if (Result.read_homogeneous()==false) {
		inv<<"boolean homogeneous = "<<"false"<<endl;
	}
	else {
		inv<<"boolean homogeneous = "<<"true"<<endl;
		if (status=="normal"||status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		Matrix Hom=Result.read_homogeneous_elements();
		nr=Hom.nr_of_rows();
		inv<<"integer height_1_elements = "<<nr<<endl;
		}
		vector<Integer> Linear_Form=Result.read_linear_form();
		inv<<"vector "<<Linear_Form.size()<<" homogeneous_weights = ";
		for (i = 0; i < Linear_Form.size(); i++) {
			inv<<Linear_Form[i]<<" ";
		}
		inv<<endl;
		if (status!="support hyperplanes"){
			inv<<"integer multiplicity = "<<Result.read_multiplicity()<<endl;
		}
		if (status=="hilbert polynomial"||status=="hilbert basis polynomial") {
		vector<Integer> h_vector=Result.read_h_vector();
		inv<<"vector "<<h_vector.size()<<" h-vector = ";
		for (i = 0; i < h_vector.size(); i++) {
			inv<<h_vector[i]<<" ";
		}
		inv<<endl;
		Integer factorial=1;
		for (i = 2; i <rank; i++) {
			 factorial*=i;
		}
		vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
		inv<<"vector "<<h_vector.size()<<" hilbert_polynomial = ";
		for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
			inv<<hilbert_polynomial[i]*(factorial /hilbert_polynomial[i+1])<<" ";
		}
		inv<<endl;
		}
	}
	if (primary) {
		inv<<"boolean primary = true"<<endl;
		Integer primary_multiplicity=Result.primary_multiplicity();
		inv<<"integer ideal_multiplicity = "<<primary_multiplicity<<endl;
	}
	else {
		inv<<"boolean primary = false"<<endl;
		}
	inv.close();
	}  
}

//--------------------------------------------------------------------------------

void Output::dual()const{
	int i,j,k,nr,nc,rank,max_decimal_length;    //read local data
	Integer buf;
	rank=Result.read_dimension();
	Matrix Change=Basis_Change.get_right();
	Matrix Extreme_Rays=Result.read_support_hyperplanes(); // in the dual mode the support hyperplanes are the generators
	Extreme_Rays=Extreme_Rays.multiplication(Change);
	write_matrix_ext(Extreme_Rays);
	Matrix Support_Hyperplanes=Result.read_generators();  // in the dual mode the extrem rays are the support hyperplanes
	write_matrix_esp(Support_Hyperplanes);         //write the suport hyperplanes of the full dimensional cone
	Matrix Hilbert_Basis=Result.read_hilbert_basis();
	write_matrix_egn(Hilbert_Basis);
	Hilbert_Basis=Hilbert_Basis.multiplication(Change);
	write_matrix_gen(Hilbert_Basis);
	
	
	
	
	if (out==true ) {//printing .out file
	string name_open=name+".out"; 							 //preparing output files
	const char* file=name_open.c_str();
	ofstream out(file);
	
	
	nr=Hilbert_Basis.nr_of_rows();
	nc=Hilbert_Basis.nr_of_columns();
	max_decimal_length=Hilbert_Basis.maximal_decimal_length();
	out<<nr<<" generators of integral closure:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j <= nc; j++) {
				buf = Hilbert_Basis.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
		out<<endl;
		}
	out<<endl;
	
	nr=Extreme_Rays.nr_of_rows();
	nc=Extreme_Rays.nr_of_columns();
	out<<nr<<" extreme rays:"<<endl;
	max_decimal_length=Extreme_Rays.maximal_decimal_length();
	for (i = 1; i <=nr; i++) {
		for (j = 1; j <=nc; j++) {
				buf = Extreme_Rays.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
		}
		out<<endl;
	}
	out<<endl;
	
	//write rank
	out<<"(original) semigroup has rank "<<rank<<endl<<endl;
	
	if (Result.read_homogeneous()==false) {
		out<<"(original) semigroup is not homogeneous"<<endl;
	}
	else {
		Matrix Hom=Result.read_homogeneous_elements();
		Hom=Hom.multiplication(Change);
		write_matrix_ht1(Hom);
		nr=Hom.nr_of_rows();
		nc=Hom.nr_of_columns();
		max_decimal_length=Hom.maximal_decimal_length();
		out<<nr<<" height 1 generators of integral closure:"<<endl;
		for (i = 1; i <= nr; i++) {
			for (j = 1; j <= nc; j++) {
				buf = Hom.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
		out<<endl;
		}
		out<<endl;
		vector<Integer> Linear_Form=Result.read_linear_form();
		Matrix Change_Inv=Basis_Change.get_right_inv();
		Change_Inv=Change_Inv.transpose();
		Linear_Form=Change_Inv.VxM(Linear_Form);
		out<<"(original) semigroup is homogeneous via the linear form:"<<endl;
		for (i = 0; i < Linear_Form.size(); i++) {
			out<<Linear_Form[i]<<" ";
		}
		out<<endl<<endl;
	    }
	
	out.close();
	}
	if (inv==true) {//printing .inv file
		string name_open=name+".inv"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream inv(file);
		
		//write Hilbert Basis
		inv<<"integer hilbert_basis_elements = "<<Hilbert_Basis.nr_of_rows()<<endl;
		//write extreme rays
		inv<<"integer number_extreme_rays = "<<Extreme_Rays.nr_of_rows()<<endl;
		inv<<"integer rank = "<<rank<<endl;
		
		if (Result.read_homogeneous()==false) {
			inv<<"boolean homogeneous = "<<"false"<<endl;
		}
		else {
			inv<<"boolean homogeneous = "<<"true"<<endl;
			Matrix Hom=Result.read_homogeneous_elements();
			nr=Hom.nr_of_rows();
			inv<<"integer height_1_elements = "<<nr<<endl;
			vector<Integer> Linear_Form=Result.read_linear_form();
			Matrix Change_Inv=Basis_Change.get_right_inv();
			Change_Inv=Change_Inv.transpose();
			Linear_Form=Change_Inv.VxM(Linear_Form);
			inv<<"vector "<<Linear_Form.size()<<" homogeneous_weights = ";
			for (i = 0; i < Linear_Form.size(); i++) {
				inv<<Linear_Form[i]<<" ";
			}
			inv<<endl;
		}
	inv.close();
	}
}


//---------------------------------------------------------------------------

void Output::error(string s) const{
	cerr <<"\nOutput "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

