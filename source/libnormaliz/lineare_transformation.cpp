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

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>

#include "lineare_transformation.h"
#include "integer.h"
#include "matrix.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

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
    
    bool success=transformation();
    // cout << "Success " << success << endl;
    if(success && !using_GMP<Integer>()){
        success=test_transformation(overflow_test_modulus);
        //  if(!success)
        //  cout << endl << "******** Test daneben!!" << endl;
    }
    if(!success){
        Matrix<mpz_class> mpz_M(M.nr_of_rows(),M.nr_of_columns());
        mat_to_mpz(M,mpz_M);
        Lineare_Transformation<mpz_class> mpz_LT(mpz_M);
        // mpz_LT.transformation();
        // mpz_LT.read();
        mat_to_Int(mpz_LT.Center,Center);
        mat_to_Int(mpz_LT.Right,Right);
        mat_to_Int(mpz_LT.Right_Inv,Right_Inv);
        rk=mpz_LT.rk;
        index=to_Int<Integer>(mpz_LT.index);
        status=mpz_LT.status;       
    }
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
long Lineare_Transformation<Integer>::get_rank() const{
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
void Lineare_Transformation<Integer>::set_rank(const size_t rank) {
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
/*
template<typename Integer>
void Lineare_Transformation<Integer>::exchange_rows(size_t row1, size_t row2){
    Center.exchange_rows(row1,row2);
}
*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
void Lineare_Transformation<Integer>::exchange_columns(size_t col1, size_t col2){
    Center.exchange_columns(col1,col2);
    Right.exchange_columns(col1,col2);
    Right_Inv.exchange_rows(col1,col2);
}
*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Lineare_Transformation<Integer>::reduce_row(size_t corner){
    return Center.reduce_row(corner);
}
*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Lineare_Transformation<Integer>::gcd_reduce_row(size_t corner){
    return Center.gcd_reduce_row(corner);
}

*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Lineare_Transformation<Integer>::reduce_column(size_t corner){
    return Center.reduce_column(corner, Right,Right_Inv);
}
*/

//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Lineare_Transformation<Integer>::gcd_reduce_column(size_t corner){
    return Center.gcd_reduce_column(corner, Right);
}
*/

//---------------------------------------------------------------------------


// new
template<typename Integer>
bool Lineare_Transformation<Integer>::transformation(){
    long r;
    // long rk_max=min(Center.nr_of_rows(),Center.nr_of_columns());
    bool success=true;
    
    while(true){
        rk=Center.row_echelon_inner(success);
        if(!success)
            return false;
        if(rk==0)
            break;
            
        /* cout << "----------------------" << endl;
        cout << "Nach Rows " << endl << endl;            
        Center.pretty_print(cout); */
            
        bool is_diagonal=true;
        for(size_t i=0;i<rk;++i){
            for(size_t j=0;j<Center.nc;++j){
                if(i!=j && Center[i][j]!=0){
                    is_diagonal=false;
                    break;
                }
            }
            if(!is_diagonal)
                break;
        }
        
        if(is_diagonal)
            break;
            
        success=Center.reduce_rows_upwards();
        if(!success)
            return false;
        
        /* cout << "----------------------" << endl;
        cout << "Nach Upwards " << endl << endl;             
        Center.pretty_print(cout); */
        
        success=Center.column_triangulate(rk,Right);
        if(!success)
            return false;
        
        /*cout << "----------------------" << endl;
        cout << "Mach Columns " << endl << endl;        
        Center.pretty_print(cout);*/
                                
    }
    for (r = 0; r < rk; r++) {
        index*=Center[r][r];
    }
    index=Iabs(index);
    
    vector<Integer> Diag(Right.nr_of_columns());
    Integer denom;    
    Right_Inv=Right.invert(Diag,denom);
    
    status="initialized, after transformation";
    return true;
}


//---------------------------------------------------------------------------

template<typename Integer>
bool Lineare_Transformation<Integer>::test_transformation(const size_t& m) const{
    size_t nc=Center.nr_of_columns();
    Matrix<Integer> N=Right.multiplication(Right_Inv, m);
    Matrix<Integer> I(nc);
    return I.equal(N,m);
}

//---------------------------------------------------------------------------

}



// basement filled with old routines
/*
                if (test==false) {
                errorOutput()<<"Arithmetic failure in linear transformation. Most likely overflow.\n";
                throw ArithmeticException();
*/


//---------------------------------------------------------------------------
/*

// middle

template<typename Integer>
bool Lineare_Transformation<Integer>::transformation(){
    long r;
    long rk_max=min(Center.nr_of_rows(),Center.nr_of_columns());
    vector<long> piv(2,0);
    for (r = 0; r < rk_max; r++) {
        piv=Center.pivot(r);
        
        cout << "piv " << piv;
        if(piv[0]<0)
            break;
                                
        exchange_rows (r,piv[0]);
        exchange_columns (r,piv[1]);
        
        bool done;
        
        do{
        piv=Center.pivot(r);
                   exchange_rows (r,piv[0]);
        exchange_columns (r,piv[1]);
            cout << "r " << r << endl; 
            Center.pretty_print(cout);
            cout << "----------" << endl; 
                    piv=Center.pivot(r);
                   exchange_rows (r,piv[0]);
        exchange_columns (r,piv[1]);     
            gcd_reduce_column(r); 
            Center.pretty_print(cout);
            cout << "=========" << endl;           
            gcd_reduce_row(r);
            Center.pretty_print(cout);
            cout << "+++++++++" << endl;
            done=true;                
            for(size_t j=r+1;j<Center.nr_of_columns();++j)
                if(Center[r][j]!=0)
                    done=false;                    
        } while(!done);
    }
    rk=r;
    for (r = 0; r < rk; r++) {
        index*=Center[r][r];
    }
    index=Iabs(index);
    status="initialized, after transformation";
    return true;
}
*/

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

/*

// old 
template<typename Integer>
bool Lineare_Transformation<Integer>::transformation(){
    long r;
    long rk_max=min(Center.nr_of_rows(),Center.nr_of_columns());
    vector<long> piv(2,0);
    for (r = 0; r < rk_max; r++) {
        piv=Center.pivot(r);
        if (piv[0]>=0) {
            do {
                exchange_rows (r,piv[0]);
                exchange_columns (r,piv[1]);
                if(!reduce_row (r))
                    return false;
                if(!reduce_column (r))
                    return false;
                Center.pretty_print(cout);
            cout << "+++++++++" << endl;
                piv=Center.pivot(r);
            } while ((piv[0]>r)||(piv[1]>r));
        }
        else
            break;
    }
    rk=r;
    for (r = 0; r < rk; r++) {
        index*=Center.read(r,r);
    }
    index=Iabs(index);
    status="initialized, after transformation";
    return true;
}
*/
