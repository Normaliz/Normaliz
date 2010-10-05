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
#include <iostream>
#include <string>
#include <algorithm>

#include "lineare_transformation.h"
#include "integer.h"

//---------------------------------------------------------------------------

namespace libnormaliz {

//---------------------------------------------------------------------------

template<typename Integer>
Lineare_Transformation<Integer>::Lineare_Transformation(){
	rk=0;
	status="non initialized";
	index=1;
}

//---------------------------------------------------------------------------

template<typename Integer>
Lineare_Transformation<Integer>::Lineare_Transformation(const Matrix<Integer>& M){
	rk=0;
	status="initialized, before transformation";
	index=1;
	Center    = Matrix<Integer>(M);
	Right     = Matrix<Integer>(M.nr_of_columns());
	Right_Inv = Matrix<Integer>(M.nr_of_columns());
}

//---------------------------------------------------------------------------

template<typename Integer>
Lineare_Transformation<Integer>::Lineare_Transformation(const Lineare_Transformation<Integer>& LT){
	rk=LT.rk;
	status=LT.status;
	index=LT.index;
	Center=LT.Center;
	Right=LT.Right;
	Right_Inv=LT.Right_Inv;
}

//---------------------------------------------------------------------------

template<typename Integer>
Lineare_Transformation<Integer>::~Lineare_Transformation(){
	//automatic destructor
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::read() const{
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

template<typename Integer>
int Lineare_Transformation<Integer>::get_rank() const{
	return rk;
}

//---------------------------------------------------------------------------

template<typename Integer>
string Lineare_Transformation<Integer>::get_status() const{
	return status;
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Lineare_Transformation<Integer>::get_index() const{
	return index;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Lineare_Transformation<Integer>::get_center()const{
	return Center;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Lineare_Transformation<Integer>::get_right() const{
	return Right;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Lineare_Transformation<Integer>::get_right_inv() const{
	return Right_Inv;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::set_rank(const int rank) {
	rk = rank;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::set_center(const Matrix<Integer>& M){
	Center=M;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::set_right(const Matrix<Integer>& M){
	Right=M;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::set_right_inv(const Matrix<Integer>& M){
	Right_Inv=M;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::exchange_rows(int row1, int row2){
	Center.exchange_rows(row1,row2);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::exchange_columns(int col1, int col2){
	Center.exchange_columns(col1,col2);
	Right.exchange_columns(col1,col2);
	Right_Inv.exchange_rows(col1,col2);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::reduce_row(int corner){
	Center.reduce_row(corner);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::reduce_column(int corner){
	Center.reduce_column(corner, Right, Right_Inv);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::transformation(){
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

template<typename Integer>
bool Lineare_Transformation<Integer>::test_transformation(const Matrix<Integer>& M,const int& m) const{
	int nc=Center.nr_of_columns();
	Matrix<Integer> N=Right.multiplication(Right_Inv, m);
	Matrix<Integer> I(nc);
	if (!I.equal(N,m)) {
		error_msg("error: Lineare_Transformation<Integer>::test_transformation failed.\nPossible arithmetic overflow in Lineare_transformation::transformation.");
		return false;
	}
	return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Lineare_Transformation<Integer>::error_msg(string s) const{
	cerr <<"\nLineare transformation "<< s<<"\n";
}

//---------------------------------------------------------------------------

template<typename Integer>
Lineare_Transformation<Integer> Transformation(const Matrix<Integer>& M) {
	Lineare_Transformation<Integer> LT(M);
	LT.transformation();
	if (test_arithmetic_overflow==true) {
		bool testing=LT.test_transformation(M,overflow_test_modulus);
		if (testing==false) {
			cerr<<"\nThe linear transformation has failed.\n";
			throw ArithmeticException();
		}
	}
	return LT;
}

//---------------------------------------------------------------------------

}
