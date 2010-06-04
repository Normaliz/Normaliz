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
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

//---------------------------------------------------------------------------

#include "lineare_transformation.h"
#include "integer.h"

//---------------------------------------------------------------------------

extern bool test_arithmetic_overflow;
extern int overflow_test_modulus;
extern void global_error_handling();

//---------------------------------------------------------------------------

Lineare_Transformation::  Lineare_Transformation(){
	rk=0;
	status="non initialized";
	index=1;
}

//---------------------------------------------------------------------------

Lineare_Transformation::  Lineare_Transformation(const Matrix& M){
	rk=0;
	status="initialized, before transformation";
	index=1;
	Matrix Help1(M.nr_of_rows());
	Matrix Help2(M);
	Matrix Help3(M.nr_of_columns());
	Matrix Help4(M.nr_of_columns());
	Center=Help2;
	Right=Help3;
	Right_Inv=Help3;
}

//---------------------------------------------------------------------------

Lineare_Transformation::Lineare_Transformation(const Lineare_Transformation& LT){
	rk=LT.rk;
	status=LT.status;
	index=LT.index;
	Center=LT.Center;
	Right=LT.Right;
	Right_Inv=LT.Right_Inv;
}

//---------------------------------------------------------------------------

Lineare_Transformation::~Lineare_Transformation(){
	//automatic destructor
}

//---------------------------------------------------------------------------

void Lineare_Transformation::read() const{
	cout<<"\nRank="<<rk<<"\n";
	cout<<"\nStatus is "<<status<<".\n";
	cout<<"\nIndex="<<index<<"\n";
	cout<<"\nCenter matrix is:\n";
	Center.read();
	cout<<"\nRight matrix is:\n";
	Right.read();
	cout<<"\nRight_Inv matrix is:\n";
	Right_Inv.read();
}

//---------------------------------------------------------------------------

int Lineare_Transformation::get_rank() const{
	return rk;
}

//---------------------------------------------------------------------------

string Lineare_Transformation::get_status() const{
	return status;
}

//---------------------------------------------------------------------------

Integer Lineare_Transformation::get_index() const{
	return index;
}

//---------------------------------------------------------------------------

Matrix Lineare_Transformation::get_center()const{
	return Center;
}

//---------------------------------------------------------------------------

Matrix Lineare_Transformation::get_right() const{
	return Right;
}

//---------------------------------------------------------------------------

Matrix Lineare_Transformation::get_right_inv() const{
	return Right_Inv;
}

//---------------------------------------------------------------------------

void Lineare_Transformation::set_rank(const int rank) {
	rk = rank;
}

//---------------------------------------------------------------------------

void Lineare_Transformation::set_center(const Matrix& M){
	Center=M;
}

//---------------------------------------------------------------------------

void Lineare_Transformation::set_right(const Matrix& M){
	Right=M;
}

//---------------------------------------------------------------------------

void Lineare_Transformation::set_right_inv(const Matrix& M){
	Right_Inv=M;
}

//---------------------------------------------------------------------------

void Lineare_Transformation::exchange_rows(int row1, int row2){
	Center.exchange_rows(row1,row2);
}

//---------------------------------------------------------------------------

void Lineare_Transformation::exchange_columns(int col1, int col2){
	Center.exchange_columns(col1,col2);
	Right.exchange_columns(col1,col2);
	Right_Inv.exchange_rows(col1,col2);
}

//---------------------------------------------------------------------------

void Lineare_Transformation::reduce_row(int corner){
	Center.reduce_row(corner);
}

//---------------------------------------------------------------------------

void Lineare_Transformation::reduce_column(int corner){
	Center.reduce_column(corner, Right, Right_Inv);
}

//---------------------------------------------------------------------------

void Lineare_Transformation::transformation(){
	int r;
	int rk_max=min(Center.nr_of_rows(),Center.nr_of_columns());
	vector<int> piv(2,0);
	for (r = 1; r <= rk_max; r++) {
		piv=Center.pivot(r);
		if (piv[0]>0) {
			do {
				exchange_rows (r,piv[0]);
				exchange_columns (r,piv[1]);
				reduce_row (r);
				reduce_column (r);
				piv=Center.pivot(r);
			} while ((piv[0]>r)||(piv[1]>r));
		}
		else
			break;
	}
	rk=r-1;
	for (r = 1; r <= rk; r++) {
		index*=Center.read(r,r);
	}
	index=Iabs(index);
	status="initialized, after transformation";
}

//---------------------------------------------------------------------------

bool Lineare_Transformation::test_transformation(const Matrix& M,const int& m) const{
	int nc=Center.nr_of_columns();
	Matrix N=Right.multiplication(Right_Inv, m);
	Matrix I(nc);
	if (!I.equal(N,m)) {
		error("error: Lineare_Transformation::test_transformation failed.\nPossible arithmetic overflow in Lineare_transformation::transformation.");
		return false;
	}
	return true;
}

//---------------------------------------------------------------------------

void Lineare_Transformation::error(string s) const{
	cerr <<"\nLineare transformation "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

Lineare_Transformation Transformation(const Matrix& M) {
	Lineare_Transformation LT(M);
	LT.transformation();
	if (test_arithmetic_overflow==true) {
		bool testing=LT.test_transformation(M,overflow_test_modulus);
		if (testing==false) {
			cerr<<"\nThe linear transformation has failed.\n";
			global_error_handling();		// if test fails but global_error_handling does not exit the program
			return LT;						// return the erroneous Linear_Transformation LT
		}
	}
	return LT;
}

//---------------------------------------------------------------------------

