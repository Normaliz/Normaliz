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
#ifndef MATRIX_HPP
#define MATRIX_HPP
//---------------------------------------------------------------------------
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

#include "integer.h"
//---------------------------------------------------------------------------

class Matrix {
  int nr;
  int nc;
  vector< vector<Integer> > elements;

//---------------------------------------------------------------------------
//				Private routines, used in the public routines
//---------------------------------------------------------------------------

  void max_rank_submatrix_lex(vector<int>& v, const int& rank) const;
  //v will be a vector with entries the indices of the first rows in lexicographic
  //order of this forming a submatrix of maximal rank.
  //v shoud be a vector of size 0 by call!!!
  
//---------------------------------------------------------------------------
public:
//---------------------------------------------------------------------------
//						Construction and destruction
//---------------------------------------------------------------------------

  Matrix();
  Matrix(int dim);                           //constructor of identity matrix
  Matrix(int row, int col);                 //main constructor, all entries 0
  Matrix(int row, int col, Integer value); //constructor, all entries are value
  Matrix(const vector< vector<Integer> >& elem);//constuctor, elements=elem
  Matrix(const Matrix& M);                //copy constructor
  ~Matrix();                             //destructor

//---------------------------------------------------------------------------
//							   Data acces
//---------------------------------------------------------------------------

  void write();                // to be modified, just for tests
  void write(int row, const vector<Integer>& data); //write  a row
  void write(int row, const vector<int>& data); //write  a row
  void write(int row, int col, Integer data);  // write data at (row,col)
  void print(const string& name, const string& suffix="in") const;         //  writes matrix into name.suffix
  void read() const;                 // to be modified, just for tests
  vector<Integer> read(int row) const;                   // read a row
  Integer read (int row, int col) const;         // read data at (row,col)
  int nr_of_rows() const;                       // returns nr
  int nr_of_columns() const;                   // returns nc
  void random();     // generates a pseudo random matrix for tests
  Matrix submatrix(const vector<int>& rows) const;  //returns a submatrix with rows
									  //corresponding to indices given by
									//the entries of rows, Numbering from 1 to n, not 0 to n-1 !
  vector<Integer> diagonale() const;     //returns the diagonale of this
								  //this should be a quadratic matrix
  int maximal_decimal_length() const;    //return the maximal number of decimals
									  //needed to write an entry

	void append(const Matrix& M); // appends the rows of M to this


	inline const Integer& get_elem(int row, int col) const {
		return elements[row-1][col-1];
	}
//---------------------------------------------------------------------------
//					Basic matrices operations
//---------------------------------------------------------------------------

  Matrix add(const Matrix& A) const;                       // returns this+A
  Matrix multiplication(const Matrix& A) const;          // returns this*A
  Matrix multiplication(const Matrix& A, int m) const;// returns this*A (mod m)
  bool equal(const Matrix& A) const;             // returns this==A
  bool equal(const Matrix& A, int m) const;     // returns this==A (mod m)
  Matrix transpose() const;                     // returns the transpose of this

//---------------------------------------------------------------------------
//							Scalar operations
//---------------------------------------------------------------------------

  void scalar_multiplication(const Integer& scalar);  //this=this*scalar
  void scalar_division(const Integer& scalar);
  //this=this div scalar, all the elements of this must be divisible with the scalar
  void reduction_modulo(const Integer& modulo);     //this=this mod scalar
  vector<Integer> make_prime();         //each row of this is reduced by its gcd
	//return a vector containing the gcd of the rows

//---------------------------------------------------------------------------
//							Vector operations
//---------------------------------------------------------------------------

   vector<Integer> MxV(const vector<Integer>& v) const;//returns this*V
   vector<Integer> VxM(const vector<Integer>& v) const;//returns V*this

//---------------------------------------------------------------------------
//						Rows and columns exchange
//---------------------------------------------------------------------------

  void exchange_rows(const int& row1, const int& row2);      //row1 is exchanged with row2
  void exchange_columns(const int& col1, const int& col2); // col1 is exchanged with col2

//---------------------------------------------------------------------------
//				Rows and columns reduction  in  respect to
//			the right-lower submatrix of this described by an int corner
//---------------------------------------------------------------------------

  void reduce_row(int corner);      //reduction by the corner-th row
  void reduce_row(int corner, Matrix& Left);//row reduction, Left used
  //for saving or copying the linear transformations
  void reduce_column(int corner);  //reduction by the corner-th column
  void reduce_column(int corner, Matrix& Right, Matrix& Right_Inv);
  //column reduction,  Right used for saving or copying the linear
  //transformations, Right_Inv used for saving the inverse linear transformations

//---------------------------------------------------------------------------
//						Pivots for rows/columns operations
//---------------------------------------------------------------------------

  vector<int> pivot(int corner); //Find the position of an element x with
  //0<abs(x)<=abs(y) for all y!=0 in the right-lower submatrix of this
  //described by an int corner
  int pivot_column(int col);  //Find the position of an element x with
  //0<abs(x)<=abs(y) for all y!=0 in the lower half of the column of this
  //described by an int col

//---------------------------------------------------------------------------
//							Matrices operations
//           --- this are more complicated algorithms ---
//---------------------------------------------------------------------------

  int diagonalize(); //computes rank and diagonalizes this, destructiv

  int rank() const; //returns rank, nondestructiv

  int rank_destructiv(); //returns rank, destructiv

  vector<int> max_rank_submatrix() const; //returns a vector with entries the
  //indices of the rows of this forming a submatrix of maximal rank

  vector<int>  max_rank_submatrix_lex() const; //returns a vector with entries
  //the indices of the first rows in lexicographic order of this forming
  //a submatrix of maximal rank.

  vector<int>  max_rank_submatrix_lex(const int& rank) const;
  //returns a vector with entries the indices of the first rows in lexicographic
  //order of this forming a submatrix of maximal rank, assuming that
  //the rank of this is known.


  Matrix solve(Matrix Right_side, Integer& det) const;// solves the system
  //this*Solution=|det(this)|*Right_side. this should be a quadratic
  //matrix with nonzero determinant. The determinat of this is saved in det.

  Matrix solve(Matrix Right_side, vector< Integer >& diagonal, Integer& det) const;// solves the system
  //this*Solution=|det(this)|*Right_side. this should be a quadratic
  //matrix with nonzero determinant. The determinat of this is saved in det,
  //and the diagonal of this after transformation into an upper triangular matrix
  //is saved in diagonal

  Matrix invert(vector< Integer >& diagonal, Integer& det) const;// solves the system
  //this*Solution=|det(this)|*I. this should be a quadratic
  //matrix with nonzero determinant. The determinat of this is saved in det,
  //and the diagonal of this after transformation into an upper triangular matrix
  //is saved in diagonal

  vector<Integer> homogeneous (bool& homogeneous) const;// solves the system
  //this*Solution=(1,1,...). this should be a   m x n , m>=n,
  //matrix of maxinal rank. The existence of solution  is marked in homogeneous
  
  vector<Integer> homogeneous_low_dim (bool& homogeneous) const;
  //same as homogeneous but also works with not maximal rank
  //uses a linear transformation to get a full rank matrix
//---------------------------------------------------------------------------
//								Tests
//---------------------------------------------------------------------------

  bool test_solve(const Matrix& Solution, const Matrix& Right_side, const Integer& det,const int& m) const;
  // test the main computation for arithmetic overflow
  // uses multiplication mod m

  bool test_invert(const Matrix& Solution, const Integer& det,const int& m) const;
  // test the main computation for arithmetic overflow
  // uses multiplication mod m

//---------------------------------------------------------------------------
//							Error msg
//---------------------------------------------------------------------------

  void error(string s) const;
};
//class end *****************************************************************
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

Matrix Solve(const Matrix& Left_side, const Matrix& Right_side,Integer& det);
// solves the system Left_side*Solution=det(Left_side)*Right_side and tests for
//errors. Left_side should be a quadratic matrix with nonzero determinant.
//The determinat of Left_side is saved in det.

Matrix Invert(const Matrix& Left_side, vector< Integer >& diagonal ,Integer& det);
// solves the system Left_side*Solution=det(Left_side)*Right_side and tests for
//errors. Left_side should be a quadratic matrix with nonzero determinant.
//The determinat of Left_side is saved in det.

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
