/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

//---------------------------------------------------------------------------
#ifndef MATRIX_HPP
#define MATRIX_HPP
//---------------------------------------------------------------------------


#include <vector>
#include <list>
#include <iostream>
#include <string>

#include "libnormaliz.h"
#include "integer.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using std::list;
using std::vector;
using std::string;

template<typename Integer> class Matrix {

    template<typename> friend class Matrix;
    template<typename> friend class Lineare_Transformation;
    
    // public:

    size_t nr;
    size_t nc;
    vector< vector<Integer> > elem;

//---------------------------------------------------------------------------
//              Private routines, used in the public routines
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//                      Rows and columns exchange
//---------------------------------------------------------------------------

    void exchange_rows(const size_t& row1, const size_t& row2);      //row1 is exchanged with row2
    void exchange_columns(const size_t& col1, const size_t& col2); // col1 is exchanged with col2

//---------------------------------------------------------------------------
//              Row and column reduction
//---------------------------------------------------------------------------
    // return value false undicates failure because of overflow
    // for all the routines below
    
    // reduction via integer division and elemntary transformations
    bool reduce_row(size_t corner);      //reduction by the corner-th row
    bool reduce_row (size_t row, size_t col); // corner at position (row,col)
    bool reduce_row(size_t corner, Matrix& LeftRHS);//row reduction, Left used
    //for saving or linear systems with right hand side RHS
    // bool reduce_column(size_t corner);  //reduction by the corner-th column
    
    // bool reduce_column(size_t corner, Matrix& Right, Matrix& Right_Inv); --- not in use presently
    //column reduction,  Right used for saving or copying the linear
    //transformations, Right_Inv used for saving the inverse linear transformations
    
    // replaces two rows by linear combinations of them
    // bool linear_comb_rows(const size_t& row,const size_t& i,const size_t& col,  
    //        const Integer& u,const Integer& v,const Integer& w,const Integer&  z);
            
    // replaces two rows by linear combinations of them
    bool linear_comb_columns(const size_t& col,const size_t& j,
            const Integer& u,const Integer& w,const Integer& v,const Integer& z);
                       
    // use the extended Euclidean algorithm for row reduction instead of elemntary transformations
    // bool gcd_reduce_row (size_t row, size_t col);
    // bool gcd_reduce_row (size_t corner);
    
    // the same for column
    bool gcd_reduce_column (size_t corner, Matrix<Integer>& Right);
    
//---------------------------------------------------------------------------
//                      Work horses
//---------------------------------------------------------------------------

    // takes product of the diagonal elem
    void do_compute_vol(bool& success);  
        
    // Solve system with coefficients and right hand side in one matrix, using elementary transformations
    // solution replaces right hand side
    bool solve_destructive_elem_inner(Integer& denom);
    void solve_destructive_elem(Integer& denom);
                    
    size_t row_echelon_inner_elem(bool& success); // does the work and checks for overflows
    size_t row_echelon_inner_bareiss(bool& success);
    // size_t row_echelon_inner_gcd(bool& success); 
    
    size_t row_echelon(bool& success); // transforms this into row echelon form and returns rank
    // reduces the rows a matrix in row echelon form upwards, from left to right
    bool reduce_rows_upwards();
    size_t row_echelon_reduce(bool& success); // combines row_echelon and reduce_rows_upwards
    
    // The Bareiss routine does NOT use Z-invertible transformations
    size_t row_echelon_bareiss(bool& success);
    

//---------------------------------------------------------------------------
//                      Pivots for rows/columns operations
//---------------------------------------------------------------------------

    vector<long> pivot(size_t corner); //Find the position of an element x with
    //0<abs(x)<=abs(y) for all y!=0 in the right-lower submatrix of this
    //described by an int corner
    // vector<long> max_pivot(size_t corner);
    long pivot_column(size_t col);  //Find the position of an element x with
    //0<abs(x)<=abs(y) for all y!=0 in the lower half of the column of this
    //described by an int col
    
    long pivot_column(size_t row,size_t col); //in column col starting from row
    
                    
public:    
  
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//                      Construction and destruction
//---------------------------------------------------------------------------

    Matrix();
    Matrix(size_t dim);                           //constructor of unit matrix
    Matrix(size_t row, size_t col);                 //main constructor, all entries 0
    Matrix(size_t row, size_t col, Integer value); //constructor, all entries set to value
    Matrix(const vector< vector<Integer> >& elem); //constuctor, elem=elem
    Matrix(const list< vector<Integer> >& elems);

//---------------------------------------------------------------------------
//                             Data access
//---------------------------------------------------------------------------

    void write(std::istream& in = std::cin);                // to be modified, just for tests
    void write(size_t row, const vector<Integer>& data); //write a row
    void write(size_t row, const vector<int>& data); //write a row
    void write_column(size_t col, const vector<Integer>& data); //write a column
    void write(size_t row, size_t col, Integer data);  // write data at (row,col)
    void print(const string& name, const string& suffix) const;         //  writes matrix into name.suffix
    void print(std::ostream& out) const;          // writes matrix to the stream
    void pretty_print(std::ostream& out, bool with_row_nr=false) const;  // writes matrix in a nice format to the stream
    void read() const;                 // to be modified, just for tests
    vector<Integer> read(size_t row) const;                   // read a row
    Integer read (size_t row, size_t col) const;         // read data at (row,col)
    size_t nr_of_rows() const;                       // returns nr
    size_t nr_of_columns() const;                   // returns nc
    /* generates a pseudo random matrix for tests, entries form 0 to mod-1 */
    void random(int mod=3);

    /* returns a submatrix with rows corresponding to indices given by
     * the entries of rows, Numbering from 0 to n-1 ! */
    Matrix submatrix(const vector<key_t>& rows) const;
    Matrix submatrix(const vector<int>& rows) const;
    Matrix submatrix(const vector<bool>& rows) const;

    Matrix& remove_zero_rows(); // remove zero rows, modifies this

    vector<Integer> diagonal() const;     //returns the diagonale of this
                                  //this should be a quadratic matrix
    size_t maximal_decimal_length() const;    //return the maximal number of decimals
                                      //needed to write an entry

    void append(const Matrix& M); // appends the rows of M to this
    void append(const vector<vector<Integer> >& M); // the same, but for another type of matrix
    void append(const vector<Integer>& v); // append the row v to this
    void cut_columns(size_t c); // remove columns, only the first c columns will survive

    inline const Integer& get_elem(size_t row, size_t col) const {
        return elem[row][col];
    }
    inline const vector< vector<Integer> >& get_elements() const {
        return elem;
    }
    inline vector<Integer> const& operator[] (size_t row) const {
        return elem[row];
    }
    inline vector<Integer>& operator[] (size_t row) { 
        return elem[row];
    }
    void set_nc(size_t col){
        nc=col;
    }

//---------------------------------------------------------------------------
//                  Basic matrices operations
//---------------------------------------------------------------------------

    Matrix add(const Matrix& A) const;                       // returns this+A
    Matrix multiplication(const Matrix& A) const;          // returns this*A
    Matrix multiplication(const Matrix& A, long m) const;// returns this*A (mod m)
    Matrix<Integer> multiplication_cut(const Matrix<Integer>& A, const size_t& c) const; // returns 
    // this*(first c columns of A)
    bool equal(const Matrix& A) const;             // returns this==A
    bool equal(const Matrix& A, long m) const;     // returns this==A (mod m)
    Matrix transpose() const;                     // returns the transpose of this
    
    bool is_diagonal() const;

//---------------------------------------------------------------------------
//                          Scalar operations
//---------------------------------------------------------------------------

    void scalar_multiplication(const Integer& scalar);  //this=this*scalar
    void scalar_division(const Integer& scalar);
    //this=this div scalar, all the elem of this must be divisible with the scalar
    void reduction_modulo(const Integer& modulo);     //this=this mod scalar
    Integer matrix_gcd() const; //returns the gcd of all elem
    vector<Integer> make_prime();         //each row of this is reduced by its gcd
    //return a vector containing the gcd of the rows

    Matrix multiply_rows(const vector<Integer>& m) const;  //returns matrix were row i is multiplied by m[i]

//---------------------------------------------------------------------------
//                          Vector operations
//---------------------------------------------------------------------------

   vector<Integer> MxV(const vector<Integer>& v) const;//returns this*V
   vector<Integer> VxM(const vector<Integer>& v) const;//returns V*this



//---------------------------------------------------------------------------
//                          Matrices operations
//           --- this are more complicated algorithms ---
//---------------------------------------------------------------------------

// Normal forms

    // transforms matrix in lower triangular form via column transformations
    // assumes that the rk is the rank and that the matrix is zero after the first rk rows
    // Right = Right*(column transformation of this call)
    bool column_trigonalize(size_t rk, Matrix<Integer>& Right);
    // combines row_echelon_reduce and column_trigonalize
    // returns column transformation matrix
    Matrix<Integer> row_column_trigonalize(size_t& rk, bool& success);
    
// rank and determinant

    size_t rank() const; //returns rank, nondestructive
    size_t rank_destructive(); //returns rank, destructive
    
    Integer vol_destructive();
    Integer vol() const;
    
// find linearly indepenpendent submatrix of maximal rank

    vector<key_t>  max_rank_submatrix_lex() const; //returns a vector with entries
    //the indices of the first rows in lexicographic order of this forming
    //a submatrix of maximal rank.
    
// Solution of linear systems with square matrix
  
    // In the following routines denom is the absolute value of the determinant of the
    // left side matrix ( =this).
    
    // Solve system with coefficients and right hand side in one matrix
    // solution replaces right hand side
    // using elementary transformations:
    void solve_destructive(vector< Integer>& diagonal, Integer& denom);
    // allowing arbitrary transformations:
    void solve_destructive(Integer& denom);

    Matrix solve(const Matrix& Right_side, vector< Integer >& diagonal, Integer& denom) const;// solves the system
    //this*Solution=denom*Right_side. this should be a quadratic /matrix with nonzero determinant.
    //The diagonal of this after transformation into an upper triangular matrix
    //is saved in diagonal
    // Variant:
    Matrix solve(const Matrix& Right_side, Integer& denom) const;

    // "This" and Right_side get destroyed!
    Matrix solve_destructive(Matrix& Right_side, vector< Integer >& diagonal, Integer& denom);

    // Returns the solution of the system in Solution (for efficiency)
    void solve_destructive_Sol(Matrix<Integer>& Right_side, vector< Integer >& diagonal, 
                    Integer& denom, Matrix<Integer>& Solution);
                    
// For non-square matrices
                    
    // The next two solve routines do not require the matrix to be square.
    // However, we want rank = number of columns, ensuring unique solvability
    
    vector<Integer> solve_rectangular(const vector<Integer>& v, Integer& denom) const;
    // computes solution vector for right side v, solution over the rationals
    // matrix needs not be quadratic, but must have rank = number of columns
    // with denominator denom. 
    // gcd of denom and solution is extracted !!!!!
    
    vector<Integer> solve_ZZ(const vector<Integer>& v) const;
    // computes solution vector for right side v
    // insists on integrality of the solution
                    
// homogenous linear systems

    Matrix<Integer> kernel () const;
    // computes a ZZ-basis of the solutions of (*this)x=0
    // the basis is formed by the ROWS of the returned matrix
                    
// inverse matrix
                    
    Matrix invert(vector< Integer >& diagonal, Integer& denom) const;// solves the system
    //this*Solution=denom*I. "this" should be a quadratic matrix with nonzero determinant. 
    //The diagonal of this after transformation into an upper triangular matrix
    //is saved in diagonal
                    
// find linear form that is constant on the rows 

    vector<Integer> find_linear_form () const;
    // Tries to find a linear form which gives the same value an all rows of this
    // this should be a m x n matrix (m>=n) of maxinal rank
    // returns an empty vector if there does not exist such a linear form
  
    vector<Integer> find_linear_form_low_dim () const;
    //same as find_linear_form but also works with not maximal rank
    //uses a linear transformation to get a full rank matrix   

};
//class end *****************************************************************

//---------------------------------------------------------------------------
//                  Conversion between integer types
//---------------------------------------------------------------------------

template<typename Integer>
void mat_to_mpz(const Matrix<Integer>& mat, Matrix<mpz_class>& mpz_mat);

template<typename Integer>
void mat_to_Int(const Matrix<mpz_class>& mpz_mat, Matrix<Integer>& mat);

} // namespace

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
