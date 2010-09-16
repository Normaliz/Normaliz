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

template<typename Integer>
Output<Integer>::Output(){
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
	BC_set=false;
	Original_Generators = Matrix<Integer>(0);
}

//---------------------------------------------------------------------------

template<typename Integer>
Output<Integer>::Output(const Output<Integer>& Out){
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
	BC_set=Out.BC_set;
	Original_Generators = Out.Original_Generators;
}

//---------------------------------------------------------------------------

template<typename Integer>
Output<Integer>::~Output(){
	//automatic destructor
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::read() const{
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
	Result.print();
	cout<<"\nBasis Change is:\n";
//	Basis_Change.print();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_name(const string& n){
	name=n;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_out(const bool& flag){
	out=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_inv(const bool& flag){
	inv=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_ext(const bool& flag){
	ext=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_esp(const bool& flag){
	esp=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_typ(const bool& flag){
	typ=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_egn(const bool& flag){
	egn=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_gen(const bool& flag){
	gen=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_sup(const bool& flag){
	sup=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_tri(const bool& flag) {
	tri=flag;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_ht1(const bool& flag) {
	ht1=flag;
}


//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_write_extra_files(){
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

template<typename Integer>
void Output<Integer>::set_write_all_files(){
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

template<typename Integer>
void Output<Integer>::set_result(const Full_Cone<Integer>& C){
	Result=C;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::compose_basis_change(const Sublattice_Representation<Integer>& BC){
	if (BC_set) {
		Basis_Change.compose(BC);		
	} else {
		Basis_Change=BC;
		BC_set=true;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::set_original_generators(Matrix<Integer> OG) {
	Original_Generators = OG;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_ext(const Matrix<Integer>& M) const{
	if (ext==true) {
		M.print(name,"ext");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_esp(const Matrix<Integer>& M) const{
	if (esp==true) {
		M.print(name,"esp");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_typ(const Matrix<Integer>& M) const{
	if (typ==true) {
		M.print(name,"typ");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_egn(const Matrix<Integer>& M) const {
	if (egn==true) {
		M.print(name,"egn");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_gen(const Matrix<Integer>& M) const {
	if (gen==true) {
		M.print(name,"gen");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_sup(const Matrix<Integer>& M) const{
	if (sup==true) {
		M.print(name,"sup");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_tri(const Matrix<Integer>& M) const{
	if (tri==true) {
		M.print(name,"tri");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_matrix_ht1(const Matrix<Integer>& M) const{
	if (ht1==true) {
		M.print(name,"ht1");
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::write_inv_file() const{
	if (inv==true) {//printing .inv file
		int i,nr;
		int rank=Basis_Change.get_rank();
		string name_open=name+".inv"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream inv(file);

		Matrix<Integer> Hilbert_Basis;                                            //write Hilbert Basis
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			Matrix<Integer> Hilbert_Basis_Full_Cone=Result.read_hilbert_basis();
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
		inv<<"integer index = "<< Basis_Change.get_index() <<endl;
		inv<<"integer number_support_hyperplanes = "<<Result.read_support_hyperplanes().nr_of_rows()<<endl;

		if (Result.read_homogeneous()==false) {
			inv<<"boolean homogeneous = "<<"false"<<endl;
		}
		else {
			inv<<"boolean homogeneous = "<<"true"<<endl;
			if (Result.isComputed(ConeProperty::Ht1Elements)) {
				Matrix<Integer> Hom=Result.read_homogeneous_elements();
				nr=Hom.nr_of_rows();
				inv<<"integer height_1_elements = "<<nr<<endl;
			}
			vector<Integer> Linear_Form = Result.read_linear_form();
			Linear_Form = Basis_Change.from_sublattice_dual(Linear_Form);
			inv<<"vector "<<Linear_Form.size()<<" homogeneous_weights = ";
			for (i = 0; i < Linear_Form.size(); i++) {
				inv<<Linear_Form[i]<<" ";
			}
			inv<<endl;
			if (Result.isComputed(ConeProperty::Triangulation)){
				inv<<"integer multiplicity = "<<Result.read_multiplicity()<<endl;
			}
			if (Result.isComputed(ConeProperty::HVector)) {
				vector<Integer> h_vector=Result.read_h_vector();
				inv<<"vector "<<h_vector.size()<<" h-vector = ";
				for (i = 0; i < h_vector.size(); i++) {
					inv<<h_vector[i]<<" ";
				}
				inv<<endl;
			}
			if (Result.isComputed(ConeProperty::HilbertPolynomial)) {
				Integer factorial=1;
				for (i = 2; i <rank; i++) {
					factorial*=i;
				}
				vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
				inv<<"vector "<<hilbert_polynomial.size()/2<<" hilbert_polynomial = ";
				for (i = 0; i < hilbert_polynomial.size(); i+=2) {
					inv<<hilbert_polynomial[i]*(factorial /hilbert_polynomial[i+1])<<" ";
				}
				inv<<endl;
			}
		}
		inv.close();
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::cone() const{
cout<<"********************output::cone"<<endl;
	int i,j,nr,rank=Basis_Change.get_rank();    //read local data
	Matrix<Integer> Generators = Basis_Change.from_sublattice(Result.read_generators());
	Matrix<Integer> Support_Hyperplanes_Full_Cone = Result.read_support_hyperplanes();

	if (esp && Result.isComputed(ConeProperty::SupportHyperplanes)) {
		//write the suport hyperplanes of the full dimensional cone
		Support_Hyperplanes_Full_Cone.print(name,"esp");
	}
	if (tri && Result.isComputed(ConeProperty::Triangulation)) {     //write triangulation
		Result.read_triangulation_volume().print(name,"tri");
		Generators.print(name,"tgn");
	}

	if (out==true) {  //printing .out file
		string name_open=name+".out"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream out(file);

		// write "header" of the .out file
		int nr_orig_gens = Original_Generators.nr_of_rows();
		if (nr_orig_gens > 0) {
			out << nr_orig_gens <<" original generators"<<endl;
		}
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			out << Result.read_hilbert_basis().nr_of_rows() <<" Hilbert basis elements"<<endl;
		}
		if (Result.isComputed(ConeProperty::Ht1Elements)) {
			out << Result.read_homogeneous_elements().nr_of_rows() <<" height 1 Hilbert basis elements"<<endl;
		}
		if (Result.isComputed(ConeProperty::ExtremeRays)) {
			vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();
			int nr_ex_rays=0;
			for (i = 0; i <Ex_Rays_Marked.size(); i++) {
				if (Ex_Rays_Marked[i]==true) {
					nr_ex_rays++;
				}
			}
			out << nr_ex_rays <<" extreme rays"<<endl;
		}
		if (Result.isComputed(ConeProperty::SupportHyperplanes)) {
			out << Result.read_support_hyperplanes().nr_of_rows() <<" support hyperplanes"<<endl;
		}
		out<<endl;
		if (rank == Basis_Change.get_dim()){                   //write rank and index
			out<<"rank = "<<rank<<" (maximal)"<<endl;
		}
		else {
			out<<"rank = "<<rank<<endl;
		}
		out<<"index = "<< Basis_Change.get_index() <<endl;

		if (Result.isComputed(ConeProperty::IsIntegrallyClosed)) {
			if (Result.isIntegrallyClosed()) {
				out << "original monoid is integrally closed"<<endl;
			} else {
				out << "original monoid is not integrally closed"<<endl;
			}
		}
		out << endl;
		
		if (Result.read_homogeneous()==false) {
			out<<"extreme rays are not homogeneous"<<endl<<endl;
		} else {
			vector<Integer> Linear_Form=Result.read_linear_form();
			Linear_Form = Basis_Change.from_sublattice_dual(Linear_Form);
			out<<"extreme rays are homogeneous via the linear form:"<<endl;
			for (i = 0; i < Linear_Form.size(); i++) {
				out<<Linear_Form[i]<<" ";
			}
			out<<endl<<endl;
			if (Result.isComputed(ConeProperty::IsHt1HilbertBasis)) {
				if (Result.isHt1HilbertBasis()) {
					out << "Hilbert basis elements are homogeneous" << endl;
				} else {
					out << "Hilbert basis elements are not homogeneous" << endl;
				}
			}
			out<<endl;
			if (Result.isComputed(ConeProperty::Triangulation)){
				out<<"multiplicity = "<<Result.read_multiplicity()<<endl;
				out<<endl;
			}
			if (Result.isComputed(ConeProperty::HVector)) {
				vector<Integer> h_vector=Result.read_h_vector();
				out<<"h-vector:"<<endl;
				for (i = 0; i < h_vector.size(); i++) {
					out<<h_vector[i]<<" ";
				}
				out<<endl<<endl;
			}
			if (Result.isComputed(ConeProperty::HilbertPolynomial)) {
				vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
				out<<"Hilbert polynomial:"<<endl;
				for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
					out<<hilbert_polynomial[i]<<"/"<<hilbert_polynomial[i+1]<<" ";
				}
				out << endl<< endl;
			}
		}


		out << "***********************************************************************"
			 << endl << endl;


		if (nr_orig_gens > 0) {
			out << nr_orig_gens <<" original generators:"<<endl;
			Original_Generators.pretty_print(out);
		}
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			Matrix<Integer> Hilbert_Basis_Full_Cone = Result.read_hilbert_basis();
			write_matrix_egn(Hilbert_Basis_Full_Cone);
			if (typ) {
				Matrix<Integer> V=Hilbert_Basis_Full_Cone.multiplication(Support_Hyperplanes_Full_Cone.transpose());
				write_matrix_typ(V);
			}
			Matrix<Integer> Hilbert_Basis = Basis_Change.from_sublattice(Hilbert_Basis_Full_Cone);
			write_matrix_gen(Hilbert_Basis);
			nr=Hilbert_Basis.nr_of_rows();
			out<<nr<<" Hilbert basis elements:"<<endl;
			Hilbert_Basis.pretty_print(out);
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
			Matrix<Integer> Extreme_Rays=Generators.submatrix(Ex_Rays_Position);
			write_matrix_ext(Extreme_Rays);
			out<<nr_ex_rays<<" extreme rays:"<<endl;
			Extreme_Rays.pretty_print(out);

		{
			//write constrains (support hyperplanes, congruences, equations)
			Matrix<Integer> Support_Hyperplanes = Basis_Change.from_sublattice_dual(Support_Hyperplanes_Full_Cone);
			out << Support_Hyperplanes.nr_of_rows() <<" support hyperplanes:"<<endl;
			Support_Hyperplanes.pretty_print(out);
			
			//equations 
			int dim = Extreme_Rays.nr_of_columns();
			int nr_of_equ = dim-rank;
			Matrix<Integer> Equations(nr_of_equ,dim);
			if (nr_of_equ > 0) {
				Lineare_Transformation<Integer> NewLT = Transformation(Extreme_Rays);
				Matrix<Integer> Help = NewLT.get_right().transpose();
				for (i = 1+rank; i <= dim; i++) {
					Equations.write(i-rank,Help.read(i));
				}
	
				cout << nr_of_equ <<" equations:" <<endl;
				Equations.pretty_print(cout);
	
				out << nr_of_equ <<" equations:" <<endl;
				Equations.pretty_print(out);
			}
	
	
			//congruences
			Matrix<Integer> Congruences = Basis_Change.get_congruences();
			int nr_of_cong = Congruences.nr_of_rows();
			if (nr_of_cong > 0) {
				out << nr_of_cong <<" congruences:" <<endl;
				Congruences.pretty_print(out);
			}

			if(sup) {
				string cst_string = name+".cst";
				const char* cst_file = cst_string.c_str();
				ofstream cst_out(cst_file);
	
				Support_Hyperplanes.print(cst_out);
				cst_out<<"hyperplanes"<<endl;
				Equations.print(cst_out);
				cst_out<<"equations"<<endl;
				Congruences.print(cst_out);
				cst_out<<"congruences"<<endl;
	
				cst_out.close();
			}	
		}
		
		if (Result.read_homogeneous()) {
			if ( Result.isComputed(ConeProperty::Ht1Elements) ) {
				Matrix<Integer> Hom = Result.read_homogeneous_elements();
				Hom = Basis_Change.from_sublattice(Hom);
				write_matrix_ht1(Hom);
				nr=Hom.nr_of_rows();
				out<<nr<<" height 1 Hilbert basis elements:"<<endl;
				Hom.pretty_print(out);
			}
		}
		out.close();
	}

	write_inv_file();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::polytop() const{
	int i,j,k,nr,max_decimal_length;
	int dim = Basis_Change.get_dim(), rank=Basis_Change.get_rank();    //read local data
	Matrix<Integer> Generators = Basis_Change.from_sublattice(Result.read_generators());
	Matrix<Integer> Support_Hyperplanes_Full_Cone = Result.read_support_hyperplanes();

	if (esp && Result.isComputed(ConeProperty::SupportHyperplanes)) {
		//write the suport hyperplanes of the full dimensional cone
		Support_Hyperplanes_Full_Cone.print(name,"esp");
	}
	if (tri && Result.isComputed(ConeProperty::Triangulation)) {     //write triangulation
		Result.read_triangulation_volume().print(name,"tri");
		Generators.print(name,"tgn");
	}

	if (out==true) {  //printing .out file
		string name_open=name+".out"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream out(file);

		// write "header" of the .out file
		int nr_orig_gens = Original_Generators.nr_of_rows();
		if (nr_orig_gens > 0) {
			out << nr_orig_gens <<" original generators"<<endl;
		}
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			out << Result.read_hilbert_basis().nr_of_rows() <<" generators of Ehrhart ring"<<endl;
		}
		if (Result.isComputed(ConeProperty::Ht1Elements)) {
			out << Result.read_homogeneous_elements().nr_of_rows() <<" lattice points in polytope"<<endl;
		}
		if (Result.isComputed(ConeProperty::ExtremeRays)) {
			vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();
			int nr_ex_rays=0;
			for (i = 0; i <Ex_Rays_Marked.size(); i++) {
				if (Ex_Rays_Marked[i]==true) {
					nr_ex_rays++;
				}
			}
			out << nr_ex_rays <<" extreme points of polytope"<<endl;
		}
		if (Result.isComputed(ConeProperty::SupportHyperplanes)) {
			out << Result.read_support_hyperplanes().nr_of_rows() <<" support hyperplanes"<<endl;
		}
		out<<endl;

		if (Result.isComputed(ConeProperty::IsIntegrallyClosed)) {
			if (Result.isIntegrallyClosed()) {
				out << "polytope is integrally closed"<<endl;
			} else {
				out << "polytope is not integrally closed"<<endl;
			}
		}
		out << endl;
		out<<"dimension of the polytope = "<<rank-1<<endl;
		
		if (Result.read_homogeneous()) {
			if (Result.isComputed(ConeProperty::Triangulation)){
				out<<"normalized volume = " << Result.read_multiplicity()<<endl;
				out<<endl;
			}
			if (Result.isComputed(ConeProperty::HVector)) {
				vector<Integer> h_vector=Result.read_h_vector();
				out<<"h-vector:"<<endl;
				for (i = 0; i < h_vector.size(); i++) {
					out<<h_vector[i]<<" ";
				}
				out<<endl<<endl;
			}
			if (Result.isComputed(ConeProperty::HilbertPolynomial)) {
				vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
				out<<"Ehrhart polynomial:"<<endl;
				for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
					out<<hilbert_polynomial[i]<<"/"<<hilbert_polynomial[i+1]<<" ";
				}
				out << endl<< endl;
			}
		}


		out << "***********************************************************************"
			 << endl << endl;


		if (nr_orig_gens > 0) {
			out << nr_orig_gens <<" original generators:"<<endl;
			Original_Generators.pretty_print(out);
		}
		Matrix<Integer> Hilbert_Basis;                                            //write Hilbert Basis
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			Matrix<Integer> Hilbert_Basis_Full_Cone = Result.read_hilbert_basis();
			write_matrix_egn(Hilbert_Basis_Full_Cone);
			if (typ) {
				Matrix<Integer> V=Hilbert_Basis_Full_Cone.multiplication(Support_Hyperplanes_Full_Cone.transpose());
				write_matrix_typ(V);
			}
			Hilbert_Basis = Basis_Change.from_sublattice(Hilbert_Basis_Full_Cone);
			write_matrix_gen(Hilbert_Basis);
			int nr = Hilbert_Basis.nr_of_rows();
			out<<nr<<" generators of Ehrhart ring:"<<endl;
			Hilbert_Basis.pretty_print(out);
		}

		if ( Result.isComputed(ConeProperty::Ht1Elements) ) {
			Matrix<Integer> Hom = Result.read_homogeneous_elements();
			Hom = Basis_Change.from_sublattice(Hom);
			write_matrix_ht1(Hom);
			Hom.cut_columns(Hom.nr_of_columns()-1);
			nr=Hom.nr_of_rows();
			out<<nr<<" lattice points in polytope:"<<endl;
			Hom.pretty_print(out);
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
		Matrix<Integer> Extreme_Rays=Generators.submatrix(Ex_Rays_Position);
		write_matrix_ext(Extreme_Rays);
		out<<nr_ex_rays<<" extreme points of polytope:"<<endl;
		Matrix<Integer> Extreme_Rays_Cut(Extreme_Rays);
		Extreme_Rays_Cut.cut_columns(Extreme_Rays_Cut.nr_of_columns()-1);      //remove extra coordinate
		Extreme_Rays_Cut.pretty_print(out);


		//write constrains (support hyperplanes, congruences, equations)
		Matrix<Integer> Support_Hyperplanes = Basis_Change.from_sublattice_dual(Support_Hyperplanes_Full_Cone);
		Integer buf;
		int nr_sup = Support_Hyperplanes.nr_of_rows();
		max_decimal_length=Support_Hyperplanes.maximal_decimal_length();
		out<<nr_sup<<" support hyperplanes:"<<endl;
		for (i = 1; i <= nr_sup; i++) {
			for (j = 1; j < dim; j++) {
				buf = Support_Hyperplanes.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
			}
			out<<" >=";
			buf = - Support_Hyperplanes.read(i,j);
			for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
				out<<" ";
			}
			out<<buf;
			out<<endl;
		}
		out<<endl;
		

		//equations 
		Lineare_Transformation<Integer> NewLT = Transformation(Extreme_Rays);
		Matrix<Integer> Help = NewLT.get_right().transpose();
		Matrix<Integer> Equations(dim-rank,dim);
		for (i = 1+rank; i <= dim; i++) {
			Equations.write(i-rank,Help.read(i));
		}

		int nr_of_equ = Equations.nr_of_rows();
		if (nr_of_equ > 0) {
			max_decimal_length = Equations.maximal_decimal_length();
			out << nr_of_equ <<" equations:" <<endl;
			for (i = 1; i <= nr_of_equ; i++) {
				for (j = 1; j < dim; j++) {
					buf = Equations.read(i,j);
					for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
						out<<" ";
					}
					out<<buf;
				}
				out<<" = ";
				buf = - Equations.read(i,j);
				for (k= 0; k <= max_decimal_length-decimal_length(buf); k++) {
					out<<" ";
				}
				out<<buf;
				out<<endl;
			}
			out<<endl;
		}


		//congruences
		Matrix<Integer> Congruences = Basis_Change.get_congruences();
		int nr_of_cong = Congruences.nr_of_rows();
		if (nr_of_cong > 0) {
			out << nr_of_cong <<" congruences:" <<endl;
			Congruences.pretty_print(out);
		}


		if(sup) {
			string cst_string = name+".cst";
			const char* cst_file = cst_string.c_str();
			ofstream cst_out(cst_file);

			Support_Hyperplanes.print(cst_out);
			cst_out<<"hyperplanes"<<endl;
			Equations.print(cst_out);
			cst_out<<"equations"<<endl;
			Congruences.print(cst_out);
			cst_out<<"congruences"<<endl;

			cst_out.close();
		}	
		
		out.close();
	}

	write_inv_file();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::rees(const bool primary) const{
	int i,j,nr;
	int dim = Basis_Change.get_dim(); 
	int rank = Basis_Change.get_rank();
	int nr_generators_ideal=0;
	Matrix<Integer> Generators = Basis_Change.from_sublattice(Result.read_generators());
	Matrix<Integer> Support_Hyperplanes_Full_Cone = Result.read_support_hyperplanes();
	Matrix<Integer> Hilbert_Basis;
	vector<int> ideal_gen_key;

	if (esp && Result.isComputed(ConeProperty::SupportHyperplanes)) {
		//write the suport hyperplanes of the full dimensional cone
		Support_Hyperplanes_Full_Cone.print(name,"esp");
	}
	if (tri && Result.isComputed(ConeProperty::Triangulation)) {     //write triangulation
		Result.read_triangulation_volume().print(name,"tri");
		Generators.print(name,"tgn");
	}

	if (out==true) {  //printing .out file
		string name_open=name+".out"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream out(file);

		// write "header" of the .out file
		int nr_orig_gens = Original_Generators.nr_of_rows();
		if (nr_orig_gens > 0) {
			out << nr_orig_gens <<" original generators"<<endl;
		}
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			Hilbert_Basis = Basis_Change.from_sublattice(Result.read_hilbert_basis());
			nr = Hilbert_Basis.nr_of_rows();
			out << nr <<" generators of integral closure of the Rees algebra"<<endl;
			for (i = 1; i <= nr; i++) {
				if (Hilbert_Basis.read(i,dim)==1) {
					nr_generators_ideal++;
					ideal_gen_key.push_back(i);
				}
			}
			out << nr_generators_ideal <<" generators of integral closure of the ideal"<<endl;
		}
		if (Result.isComputed(ConeProperty::ExtremeRays)) {
			vector<bool> Ex_Rays_Marked=Result.read_extreme_rays();
			int nr_ex_rays=0;
			for (i = 0; i <Ex_Rays_Marked.size(); i++) {
				if (Ex_Rays_Marked[i]==true) {
					nr_ex_rays++;
				}
			}
			out << nr_ex_rays <<" extreme rays"<<endl;
		}
		if (Result.isComputed(ConeProperty::SupportHyperplanes)) {
			out << Result.read_support_hyperplanes().nr_of_rows() <<" support hyperplanes"<<endl;
		}
		out<<endl;
		if (rank == Basis_Change.get_dim()){                   //write rank and index
			out<<"rank = "<<rank<<" (maximal)"<<endl;
		}
		else {
			out<<"rank = "<<rank<<endl;
		}

		if (Result.isComputed(ConeProperty::IsIntegrallyClosed)) {
			if (Result.isIntegrallyClosed()) {
				out << "original monoid is integrally closed"<<endl;
			} else {
				out << "original monoid is not integrally closed"<<endl;
			}
		}
		out << endl;
		
		if (Result.read_homogeneous()==false) {
			out<<"extreme rays are not homogeneous"<<endl<<endl;
		} else {
			vector<Integer> Linear_Form=Result.read_linear_form();
			Linear_Form = Basis_Change.from_sublattice_dual(Linear_Form);
			out<<"extreme rays are homogeneous via the linear form:"<<endl;
			for (i = 0; i < Linear_Form.size(); i++) {
				out<<Linear_Form[i]<<" ";
			}
			out<<endl<<endl;
			if (Result.isComputed(ConeProperty::IsHt1HilbertBasis)) {
				if (Result.isHt1HilbertBasis()) {
					out << "generators of integral closure of the Rees algebra are homogeneous" << endl;
				} else {
					out << "generators of integral closure of the Rees algebra are not homogeneous" << endl;
				}
			}
			out<<endl;
			if (Result.isComputed(ConeProperty::Triangulation)){
				out<<"multiplicity = "<<Result.read_multiplicity()<<endl;
				out<<endl;
			}
			if (Result.isComputed(ConeProperty::HVector)) {
				vector<Integer> h_vector=Result.read_h_vector();
				out<<"h-vector:"<<endl;
				for (i = 0; i < h_vector.size(); i++) {
					out<<h_vector[i]<<" ";
				}
				out<<endl<<endl;
			}
			if (Result.isComputed(ConeProperty::HilbertPolynomial)) {
				vector<Integer> hilbert_polynomial=Result.read_hilbert_polynomial();
				out<<"Hilbert polynomial:"<<endl;
				for (i = 0; i < hilbert_polynomial.size(); i=i+2) {
					out<<hilbert_polynomial[i]<<"/"<<hilbert_polynomial[i+1]<<" ";
				}
				out << endl<< endl;
			}
		}

      if (primary) {
			out<<"ideal is primary to the ideal generated by the indeterminates"<<endl;
			Integer primary_multiplicity=Result.primary_multiplicity();
			out<<"multiplicity of the ideal = "<<primary_multiplicity<<endl;
		} else {
			out<<"ideal is not primary to the ideal generated by the indeterminates"<<endl;
		}
		out << endl;


		out << "***********************************************************************"
			 << endl << endl;


		if (nr_orig_gens > 0) {
			out << nr_orig_gens <<" original generators:"<<endl;
			Original_Generators.pretty_print(out);
		}
		if (Result.isComputed(ConeProperty::HilbertBasis)) {
			Matrix<Integer> Hilbert_Basis_Full_Cone = Result.read_hilbert_basis();
			write_matrix_egn(Hilbert_Basis_Full_Cone);
			if (typ) {
				Matrix<Integer> V=Hilbert_Basis_Full_Cone.multiplication(Support_Hyperplanes_Full_Cone.transpose());
				write_matrix_typ(V);
			}
			write_matrix_gen(Hilbert_Basis);
			nr=Hilbert_Basis.nr_of_rows();
			out<<nr<<" generators of integral closure of the Rees algebra:"<<endl;
			Hilbert_Basis.pretty_print(out);
			
			out << nr_generators_ideal <<" generators of integral closure of the ideal:"<<endl;
			Matrix<Integer> Ideal_Gens = Hilbert_Basis.submatrix(ideal_gen_key);
			Ideal_Gens.cut_columns(dim-1);
			Ideal_Gens.pretty_print(out);
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
		Matrix<Integer> Extreme_Rays=Generators.submatrix(Ex_Rays_Position);
		write_matrix_ext(Extreme_Rays);
		out<<nr_ex_rays<<" extreme rays:"<<endl;
		Extreme_Rays.pretty_print(out);


		//write constrains (support hyperplanes, congruences, equations)
		Matrix<Integer> Support_Hyperplanes = Basis_Change.from_sublattice_dual(Support_Hyperplanes_Full_Cone);
		out << Support_Hyperplanes.nr_of_rows() <<" support hyperplanes:"<<endl;
		Support_Hyperplanes.pretty_print(out);
		
		//equations 
		Lineare_Transformation<Integer> NewLT = Transformation(Extreme_Rays);
		Matrix<Integer> Help = NewLT.get_right().transpose();
		int dim = Extreme_Rays.nr_of_columns();
		Matrix<Integer> Equations(dim-rank,dim);
		for (i = 1+rank; i <= dim; i++) {
			Equations.write(i-rank,Help.read(i));
		}

		int nr_of_equ = Equations.nr_of_rows();
		if (nr_of_equ > 0) {
			out << nr_of_equ <<" equations:" <<endl;
			Equations.pretty_print(out);
		}


		//congruences
		Matrix<Integer> Congruences = Basis_Change.get_congruences();
		int nr_of_cong = Congruences.nr_of_rows();
		if (nr_of_cong > 0) {
			out << nr_of_cong <<" congruences:" <<endl;
			Congruences.pretty_print(out);
		}


		if(sup) {
			string cst_string = name+".cst";
			const char* cst_file = cst_string.c_str();
			ofstream cst_out(cst_file);

			Support_Hyperplanes.print(cst_out);
			cst_out<<"hyperplanes"<<endl;
			Equations.print(cst_out);
			cst_out<<"equations"<<endl;
			Congruences.print(cst_out);
			cst_out<<"congruences"<<endl;

			cst_out.close();
		}	
		
		
		if (Result.read_homogeneous()) {
			if ( Result.isComputed(ConeProperty::Ht1Elements) ) {
				Matrix<Integer> Hom = Result.read_homogeneous_elements();
				Hom = Basis_Change.from_sublattice(Hom);
				write_matrix_ht1(Hom);
			}
		}
		out.close();
	}



	if (inv==true) {//printing .inv file
		write_inv_file();

		
		string name_open=name+".inv"; 							 //preparing output files
		const char* file=name_open.c_str();
		ofstream inv(file,ios_base::app);

		if (primary) {
			inv<<"boolean primary = true"<<endl;
			Integer primary_multiplicity=Result.primary_multiplicity();
			inv<<"integer ideal_multiplicity = "<<primary_multiplicity<<endl;
		} else {
			inv<<"boolean primary = false"<<endl;
		}
		inv.close();
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Output<Integer>::error(string s) const{
	cerr <<"\nOutput<Integer> "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

