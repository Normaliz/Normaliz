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
// This class implements the elementary divisor algoritm.
// An object Lineare_transformation consists of four matrices.
// Given a matrix M we want to obtain matrices U and V such that
// UMV=the diagonal form of M.
// In this class U is called Left, M Center and V Right.
// The inverse of V is also computed as Right_Inv.
// A Lineare_transormation should be initialized with a matrix M by the main
// constructor. After using the procedure transformation the result should be
// interpreted as follows:
// Left=U; Center=the diagonal form of M; Right=V;
// Right_Inv=the inverse of V; rk=rank of M.
// We recommend using the non-member function Transformation, which also test
// for possible arithmetic overflows
//---------------------------------------------------------------------------

#ifndef LINEARE_TRANSFORMATION_H
#define LINEARE_TRANSFORMATION_H

#include "matrix.h"
//---------------------------------------------------------------------------

class Lineare_Transformation {
  int rk;
  string status;
  Integer index;
  Matrix Left;
  Matrix Center;
  Matrix Right;
  Matrix Right_Inv;
//---------------------------------------------------------------------------
public:
//---------------------------------------------------------------------------
//						Construction and destruction
//---------------------------------------------------------------------------

  Lineare_Transformation();
  Lineare_Transformation(const Matrix& M);      //main constructor
  Lineare_Transformation(const Lineare_Transformation& LT);  //copy constructor
  ~Lineare_Transformation();            //destructor

//---------------------------------------------------------------------------
//						   Data acces
//---------------------------------------------------------------------------

  void read() const;                   // to be modified, just for tests
  int get_rank() const;              //returns rank if status is
							   // "initialized, after transformation"
  string get_status()const;       //returns status, may be:
							 // "non initialized"
							//  "initialized, before transformation"
						   //  "initialized, after transformation"
  Integer get_index() const;
  Matrix get_left() const;             //read left matrix
  Matrix get_center() const;          //read center matrix
  Matrix get_right() const;          //read right matrix
  Matrix get_right_inv() const;     //read the inverse of the right matrix
  void set_rank(const int rank);
  void set_left(const Matrix& M);             //write left matrix
  void set_center(const Matrix& M);          //write center matrix
  void set_right(const Matrix& M);          //write right matrix
  void set_right_inv(const Matrix& M);     //write the inverse of the right matrix

//---------------------------------------------------------------------------
//					  Rows and columns exchange
//---------------------------------------------------------------------------

  void exchange_rows(int row1, int row2);     //similar to Matrix::exchange_rows
  void exchange_columns(int col1, int col2); //similar to Matrix::exchange_columns

//---------------------------------------------------------------------------
//					Rows and columns reduction
//---------------------------------------------------------------------------

  void reduce_row(int corner);      //similar to Matrix::reduce_row
  void reduce_column(int corner);  //similar to Matrix::reduce_column

//---------------------------------------------------------------------------
//					 Algorithms
//---------------------------------------------------------------------------

  void transformation(); //makes the main computation
						//no tests for errors
//---------------------------------------------------------------------------
//Tests
//---------------------------------------------------------------------------

  bool test_transformation(const Matrix& M, const int& m) const;
									 // test the main computation
									// for arithmetic overflow
								   // uses multiplication mod m
//---------------------------------------------------------------------------
//					Error msg
//---------------------------------------------------------------------------

  void error(string s) const;
};
//class end *****************************************************************
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
Lineare_Transformation Transformation(const Matrix& M); //makes the main
										  //computation, test for errors
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
