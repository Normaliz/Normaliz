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

#include <fstream>
#include <algorithm>
#include <math.h>

#include "matrix.h"
#include "vector_operations.h"
#include "lineare_transformation.h"
#include "normaliz_exception.h"
#include "sublattice_representation.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;


//---------------------------------------------------------------------------
//Public
//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(){
    nr=0;
    nc=0;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(size_t dim){
    nr=dim;
    nc=dim;
    elem = vector< vector<Integer> >(dim, vector<Integer>(dim));
    for (size_t i = 0; i < dim; i++) {
        elem[i][i]=1;
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(size_t row, size_t col){
    nr=row;
    nc=col;
    elem = vector< vector<Integer> >(row, vector<Integer>(col));
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(size_t row, size_t col, Integer value){
    nr=row;
    nc=col;
    elem = vector< vector<Integer> > (row, vector<Integer>(col,value));
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(const vector< vector<Integer> >& new_elem){
    nr=new_elem.size();
    if (nr>0) {
        nc=new_elem[0].size();
        elem=new_elem;
        //check if all rows have the same length
        for (size_t i=1; i<nr; i++) {
            if (elem[i].size() != nc) {
                errorOutput() << "Inconsistent lengths of rows in matrix!" << endl;
                throw BadInputException();
            }
        }
    } else {
        nc=0;
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(const list< vector<Integer> >& new_elem){
    nr = new_elem.size();
    elem = vector< vector<Integer> > (nr);
    nc = 0;
    size_t i=0;
    typename list< vector<Integer> >::const_iterator it=new_elem.begin();
    for(; it!=new_elem.end(); ++it, ++i) {
        if(i == 0) {
            nc = (*it).size();
        } else {
            if ((*it).size() != nc) {
                errorOutput() << "Inconsistent lengths of rows in matrix!" << endl;
                throw BadInputException();
            }
        }
        elem[i]=(*it);
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(istream& in){
    size_t i,j;
    for(i=0; i<nr; i++){
        for(j=0; j<nc; j++) {
            in >> elem[i][j];
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(size_t row, const vector<Integer>& data){
    assert(row >= 0);
    assert(row < nr); 
    assert(nc == data.size());
    
    elem[row]=data;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(size_t row, const vector<int>& data){
    assert(row >= 0);
    assert(row < nr); 
    assert(nc == data.size());

    for (size_t i = 0; i < nc; i++) {
        elem[row][i]=data[i];
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write_column(size_t col, const vector<Integer>& data){
    assert(col >= 0);
    assert(col < nc); 
    assert(nr == data.size());

    for (size_t i = 0; i < nr; i++) {
        elem[i][col]=data[i];
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(size_t row, size_t col, Integer data){
    assert(row >= 0);
    assert(row < nr); 
    assert(col >= 0);
    assert(col < nc); 

    elem[row][col]=data;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::print(const string& name,const string& suffix) const{
    string file_name = name+"."+suffix;
    const char* file = file_name.c_str();
    ofstream out(file);
    print(out);
    out.close();
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::print(ostream& out) const{
    size_t i,j;
    out<<nr<<endl<<nc<<endl;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            out<<elem[i][j]<<" ";
        }
        out<<endl;
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::pretty_print(ostream& out, bool with_row_nr) const{
    size_t i,j,k;
    size_t max_length = maximal_decimal_length();
    size_t max_index_length = decimal_length(nr);
    for (i = 0; i < nr; i++) {
        if (with_row_nr) {
            for (k= 0; k <= max_index_length - decimal_length(i); k++) {
                out<<" ";
            }
            out << i << ": ";
        }
        for (j = 0; j < nc; j++) {
            for (k= 0; k <= max_length - decimal_length(elem[i][j]); k++) {
                out<<" ";
            }
            out<<elem[i][j];
        }
        out<<endl;
    }
}
//---------------------------------------------------------------------------


template<typename Integer>
void Matrix<Integer>::read() const{      //to overload for files
    size_t i,j;
    for(i=0; i<nr; i++){
        cout << "\n" ;
        for(j=0; j<nc; j++) {
            cout << elem[i][j] << " ";
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::read(size_t row) const{
    assert(row >= 0);
    assert(row < nr);

    return elem[row];
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::read (size_t row, size_t col) const{
    assert(row >= 0);
    assert(row < nr); 
    assert(col >= 0);
    assert(col < nc); 

    return elem[row][col];
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::nr_of_rows () const{
    return nr;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::nr_of_columns () const{
    return nc;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::random (int mod) {
    size_t i,j;
    int k;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            k = rand();
            elem[i][j]=k % mod;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::submatrix(const vector<key_t>& rows) const{
    size_t size=rows.size(), j;
    Matrix<Integer> M(size, nc);
    for (size_t i=0; i < size; i++) {
        j=rows[i];
        assert(j >= 0);
        assert(j < nr);
        M.elem[i]=elem[j];
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::submatrix(const vector<int>& rows) const{
    size_t size=rows.size(), j;
    Matrix<Integer> M(size, nc);
    for (size_t i=0; i < size; i++) {
        j=rows[i];
        assert(j >= 0);
        assert(j < nr);
        M.elem[i]=elem[j];
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::submatrix(const vector<bool>& rows) const{
    assert(rows.size() == nr);
    size_t size=0;
    for (size_t i = 0; i <rows.size(); i++) {
        if (rows[i]) {
            size++;
        }
    }
    Matrix<Integer> M(size, nc);
    size_t j = 0;
    for (size_t i = 0; i < nr; i++) {
        if (rows[i]) {
            M.elem[j++] = elem[i];
        }
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>& Matrix<Integer>::remove_zero_rows() {
    size_t from = 0, to = 0; // maintain to <= from
    while (from < nr && v_is_zero(elem[from])) from++; //skip zero rows
    while (from < nr) {  // go over matrix
        // now from is a non-zero row
        if (to != from) swap(elem[to],elem[from]);
        ++to; ++from;
        while (from < nr && v_is_zero(elem[from])) from++; //skip zero rows
    }
    nr = to;
    elem.resize(nr);
    return *this;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::diagonal() const{
    assert(nr == nc); 
    vector<Integer> diag(nr);
    for(size_t i=0; i<nr;i++){
        diag[i]=elem[i][i];
    }
    return diag;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::maximal_decimal_length() const{
    size_t i,j,maxim=0;
    for (i = 0; i <nr; i++) {
        for (j = 0; j <nc; j++) {
            maxim=max(maxim,decimal_length(elem[i][j]));
        }
    }
    return maxim;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::append(const Matrix<Integer>& M) {
    assert (nc == M.nc);
    elem.reserve(nr+M.nr);
    for (size_t i=0; i<M.nr; i++) {
        elem.push_back(M.elem[i]);
    }
    nr += M.nr;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::append(const vector<vector<Integer> >& M) {
    if(M.size()==0)
        return;
    assert (nc == M[0].size());
    elem.reserve(nr+M.size());
    for (size_t i=0; i<M.size(); i++) {
        elem.push_back(M[i]);
    }
    nr += M.size();
}

//---------------------------------------------------------------------------


template<typename Integer>
void Matrix<Integer>::append(const vector<Integer>& V) {
    assert (nc == V.size());
    elem.push_back(V);
    nr++;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::cut_columns(size_t c) {
    assert (c >= 0);
    assert (c <= nc);
    for (size_t i=0; i<nr; i++) {
        elem[i].resize(c);
    }
    nc = c;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::add(const Matrix<Integer>& A) const{
    assert (nr == A.nr);
    assert (nc == A.nc);
    
    Matrix<Integer> B(nr,nc);
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            B.elem[i][j]=elem[i][j]+A.elem[i][j];
        }
    }
    return B;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::multiplication(const Matrix<Integer>& A) const{
    assert (nc == A.nr);

    Matrix<Integer> B(nr,A.nc,0);  //initialized with 0
    size_t i,j,k;
    for(i=0; i<B.nr;i++){
        for(j=0; j<B.nc; j++){
            for(k=0; k<nc; k++){
                B.elem[i][j]=B.elem[i][j]+elem[i][k]*A.elem[k][j];
            }
        }
    }
    return B;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::multiplication_cut(const Matrix<Integer>& A, const size_t& c) const{
    assert (nc == A.nr);
    assert(c<= A.nc);

    Matrix<Integer> B(nr,c,0);  //initialized with 0
    size_t i,j,k;
    for(i=0; i<B.nr;i++){
        for(j=0; j<c; j++){
            for(k=0; k<nc; k++){
                B.elem[i][j]=B.elem[i][j]+elem[i][k]*A.elem[k][j];
            }
        }
    }
    return B;
}


//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::multiplication(const Matrix<Integer>& A, long m) const{
    assert (nc == A.nr);

    Matrix<Integer> B(nr,A.nc,0);  //initialized with 0
    size_t i,j,k;
    for(i=0; i<B.nr;i++){
        for(j=0; j<B.nc; j++){
                for(k=0; k<nc; k++){
                B.elem[i][j]=(B.elem[i][j]+elem[i][k]*A.elem[k][j])%m;
                if (B.elem[i][j]<0) {
                    B.elem[i][j]=B.elem[i][j]+m;
                }
            }
        }
    }
    return B;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::equal(const Matrix<Integer>& A) const{
    if ((nr!=A.nr)||(nc!=A.nc)){  return false; }
    size_t i,j;
    for (i=0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            if (elem[i][j]!=A.elem[i][j]) {
                return false;
            }
        }
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::equal(const Matrix<Integer>& A, long m) const{
    if ((nr!=A.nr)||(nc!=A.nc)){  return false; }
    size_t i,j;
    for (i=0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            if (((elem[i][j]-A.elem[i][j]) % m)!=0) {
                return false;
            }
        }
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::transpose()const{
    Matrix<Integer> B(nc,nr);
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            B.elem[j][i]=elem[i][j];
        }
    }
    return B;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::scalar_multiplication(const Integer& scalar){
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            elem[i][j] *= scalar;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::scalar_division(const Integer& scalar){
    size_t i,j;
    assert(scalar != 0);
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            assert (elem[i][j]%scalar == 0);
            elem[i][j] /= scalar;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::reduction_modulo(const Integer& modulo){
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            elem[i][j] %= modulo;
            if (elem[i][j] < 0) {
                elem[i][j] += modulo;
            }
        }
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::matrix_gcd() const{
    Integer g=0,h;
    for (size_t i = 0; i <nr; i++) {
        h = v_gcd(elem[i]);
        g = gcd<Integer>(g, h);
        if (g==1) return g;
    }
    return g;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::make_prime() {
    vector<Integer> g(nr);
    for (size_t i = 0; i <nr; i++) {
        g[i] = v_make_prime(elem[i]);
    }
    return g;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::multiply_rows(const vector<Integer>& m) const{  //row i is multiplied by m[i]
  Matrix M = Matrix(nr,nc);
  size_t i,j;
  for (i = 0; i<nr; i++) {
     for (j = 0; j<nc; j++) {
        M.elem[i][j] = elem[i][j]*m[i];
     }
  }
  return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::MxV(const vector<Integer>& v) const{
    assert (nc == v.size());
    vector<Integer> w(nr);
    for(size_t i=0; i<nr;i++){
        w[i]=v_scalar_product(elem[i],v);
    }
    return w;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::VxM(const vector<Integer>& v) const{
    assert (nr == v.size());
    vector<Integer> w(nc,0);
    size_t i,j;
    for (i=0; i<nc; i++){
        for (j=0; j<nr; j++){
            w[i] += v[j]*elem[j][i];
        }
    }
    return w;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::is_diagonal() const{

    for(size_t i=0;i<nr;++i)
        for(size_t j=0;j<nc;++j)
            if(i!=j && elem[i][j]!=0)
                return false;
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<long> Matrix<Integer>::pivot(size_t corner){
    assert(corner < nc);
    assert(corner < nr);
    size_t i,j;
    Integer help=0;
    vector<long> v(2,-1);

    for (i = corner; i < nr; i++) {
        for (j = corner; j < nc; j++) {
            if (elem[i][j]!=0) {
                if ((help==0)||(Iabs(elem[i][j])<help)) {
                    help=Iabs(elem[i][j]);
                    v[0]=i;
                    v[1]=j;
                    if (help == 1) return v;
                }
            }
        }
    }
    
    return v;
}

//---------------------------------------------------------------------------

template<typename Integer>
long Matrix<Integer>::pivot_column(size_t row,size_t col){
    assert(col < nc);
    assert(row < nr);
    size_t i;
    long j=-1;
    Integer help=0;

    for (i = row; i < nr; i++) {
        if (elem[i][col]!=0) {
            if ((help==0)||(Iabs(elem[i][col])<help)) {
                help=Iabs(elem[i][col]);
                j=i;
                if (help == 1) return j;
            }
        }
    }

    return j;
}

//---------------------------------------------------------------------------

template<typename Integer>
long Matrix<Integer>::pivot_column(size_t col){
    return pivot_column(col,col);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::exchange_rows(const size_t& row1, const size_t& row2){
    if (row1 == row2) return;
    assert(row1 < nr);
    assert(row2 < nr);
    elem[row1].swap(elem[row2]);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::exchange_columns(const size_t& col1, const size_t& col2){
    if (col1 == col2) return;
    assert(col1 < nc);
    assert(col2 < nc);
    for(size_t i=0; i<nr;i++){
        std::swap(elem[i][col1], elem[i][col2]);
    }
}

//---------------------------------------------------------------------------
 
template<typename Integer>
bool Matrix<Integer>::reduce_row (size_t row, size_t col) {
    assert(col < nc);
    assert(row < nr);
    size_t i,j;
    Integer help;
    for (i =row+1; i < nr; i++) {
        if (elem[i][col]!=0) {
            help=elem[i][col] / elem[row][col];
            for (j = col; j < nc; j++) {
                elem[i][j] -= help*elem[row][j];
                if (!check_range(elem[i][j]) ) {
                    return false;
                }
            }
            // v_el_trans<Integer>(elem[row],elem[i],-help,col);
        }
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool  Matrix<Integer>::reduce_row (size_t corner) {
    return reduce_row(corner,corner);
}

//---------------------------------------------------------------------------
 
template<typename Integer>
bool Matrix<Integer>::reduce_rows_upwards () {
// assumes that "this" is in row echelon form
// and reduces eevery column in which the rank jumps 
// by its lowest element

    for(size_t row=1;row<nr;++row){
        size_t col;
        for(col=0;col<nc;++col)
            if(elem[row][col]!=0)
                break;
        if(col==nc)
            continue;
        if(elem[row][col]<0)
            v_scalar_multiplication<Integer>(elem[row],-1);
        
        for(long i=row-1;i>=0;--i){
            Integer quot, rem;
            
            minimal_remainder(elem[i][col],elem[row][col],quot,rem);
            // cout << "Zahlen " << elem[row][col] << " " << elem[i][col] << " " << quot << " " << rem << endl;
            elem[i][col]=rem;
            for(size_t j=col+1;j<nc;++j){
                elem[i][j]-=quot* elem[row][j];
                if ( !check_range(elem[i][j]) ) {
                    return false;
                }
            }                            
            // cout << "Nach Abzug " << elem[i];                  
        }
    }
    return true;
}

//---------------------------------------------------------------------------
 
template<typename Integer>
bool Matrix<Integer>::linear_comb_columns(const size_t& col,const size_t& j,
            const Integer& u,const Integer& w,const Integer& v,const Integer& z){
                       
    for(size_t i=0;i<nr;++i){
        Integer rescue=elem[i][col];
        elem[i][col]=u*elem[i][col]+v*elem[i][j];
        elem[i][j]=w*rescue+z*elem[i][j];
        if ( (!check_range(elem[i][col])  || !check_range(elem[i][j]) )) {
            return false;
        }        
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::gcd_reduce_column (size_t corner, Matrix<Integer>& Right){
    assert(corner < nc);
    assert(corner < nr);
    Integer d,u,w,z,v;
    for(size_t j=corner+1;j<nc;++j){
       d=ext_gcd(elem[corner][corner],elem[corner][j],u,v);
       w=-elem[corner][j]/d;
       z=elem[corner][corner]/d;
       // Now we multiply the submatrix formed by columns "corner" and "j" 
       // and rows corner,...,nr from the right by the 2x2 matrix
       // | u w |
       // | v z |              
       if(!linear_comb_columns(corner,j,u,w,v,z))
           return false; 
       if(!Right.linear_comb_columns(corner,j,u,w,v,z))
           return false;  
    }   
    return true;
}


//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::column_trigonalize(size_t rk, Matrix<Integer>& Right) { 
    assert(Right.nr == nc);
    assert(Right.nc == nc);
    vector<long> piv(2,0);       
    for(size_t j=0;j<rk;++j){
            piv=pivot(j);
            assert(piv[0]>=0); // pritect against wrong rank
            exchange_rows (j,piv[0]);
            exchange_columns (j,piv[1]);
            Right.exchange_columns(j,piv[1]);
            if(!gcd_reduce_column(j, Right))
                return false;
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::compute_vol(bool& success){
        
    assert(nr<=nc);
    
    Integer det, test_det = 1;
    if(do_arithmetic_check<Integer>())
        for (size_t i=0; i<nr; i++){
            test_det=(test_det*elem[i][i]%overflow_test_modulus)%overflow_test_modulus;
      test_det=Iabs(test_det);  
    }
    
    det=1;
    for(size_t i=0;i<nr;++i)
        det*=elem[i][i];           
    det=Iabs(det);
    
    if(do_arithmetic_check<Integer>() && test_det!=det%overflow_test_modulus){
        success=false;
    }
    return det;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::row_echelon_inner_elem(bool& success){

    size_t pc=0;
    long piv=0, rk=0;
    success=true;

    if(nr==0)
        return 0;
    
    for (rk = 0; rk < (long) nr; rk++){
        for(;pc<nc;pc++){
            piv=pivot_column(rk,pc);
            if(piv>=0)
                break;
        }
        if(pc==nc)
            break;
        do{
            exchange_rows (rk,piv);
            if(!reduce_row(rk,pc)){
                success=false;
                return rk;
            }
            piv=pivot_column(rk,pc);
        }while (piv>rk);
    }
                
    return rk;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::row_echelon_inner_bareiss(bool& success, Integer& det){
// no overflow checks since this is supposed to be only used with GMP

    success=true;
    if(nr==0)
        return 0;
    assert(using_GMP<Integer>());

    size_t pc=0;
    long piv=0, rk=0;
    vector<bool> last_time_mult(nr,false),this_time_mult(nr,false);
    Integer last_div=1,this_div=1;
    size_t this_time_exp=0,last_time_exp=0;
    Integer det_factor=1;
    
    for (rk = 0; rk < (long) nr; rk++){
    
        for(;pc<nc;pc++){
            piv=pivot_column(rk,pc);
            if(piv>=0)
                break;
        }
        if(pc==nc)
            break;

        exchange_rows (rk,piv);
        swap(last_time_mult[rk],last_time_mult[piv]);
        
        if(!last_time_mult[rk])
            for(size_t i=0;i<nr;++i)
                last_time_mult[i]=false;
         
        Integer a=elem[rk][pc];
        this_div=Iabs(a);
        this_time_exp=0;
        
        for(size_t i=rk+1;i<nr;++i){
            if(elem[i][pc]==0){
                this_time_mult[i]=false;
                continue;
            }
            this_time_exp++;
            this_time_mult[i]=true;
            bool divide=last_time_mult[i] && (last_div!=1);
            if(divide)
                last_time_exp--;
            Integer b=elem[i][pc];
            elem[i][pc]=0;
            if(a==1){
                for(size_t j=pc+1;j<nc;++j){
                    elem[i][j]=elem[i][j]-b*elem[rk][j]; 
                    if(divide){
                        elem[i][j]/=last_div;
                    }
                }            
            }
            else{
                if(a==-1){
                    for(size_t j=pc+1;j<nc;++j){
                        elem[i][j]=-elem[i][j]-b*elem[rk][j]; 
                        if(divide){
                            elem[i][j]/=last_div;
                        }            
                    }
                }
                else{            
                    for(size_t j=pc+1;j<nc;++j){
                        elem[i][j]=a*elem[i][j]-b*elem[rk][j]; 
                       if(divide){
                            elem[i][j]/=last_div;
                        }                   
                    }
                }
            }
        }
        // pretty_print(cout);
        // cout << "-------- " << last_time_mult[rk] << " " << last_div << " " << this_time_exp << " " << last_time_exp << endl;
        for(size_t i=0;i<last_time_exp;++i)
            det_factor*=last_div;
        last_time_mult=this_time_mult;
        last_div=this_div;
        last_time_exp=this_time_exp;

    }
    
    det=0;
    if(nr<=nc && rk==nr){ // must allow nonsquare matrices
        det=1;
        for(size_t i=0;i<nr;++i)
            det*=elem[i][i];            
        det=Iabs<Integer>(det/det_factor);        
    }
    return rk;
}


//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::row_echelon_reduce(bool& success){

    size_t rk=row_echelon_inner_elem(success);
    if(success)
        reduce_rows_upwards();
    return rk;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::row_column_trigonalize(size_t& rk, bool& success) {

    Matrix<Integer> Right(nc);
    rk=row_echelon_reduce(success);
    if(success)
        success=column_trigonalize(rk,Right); 
    return Right; 
} 

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::row_echelon(bool& success, bool do_compute_vol, Integer& det){
    
    if(using_GMP<Integer>())
        return row_echelon_inner_bareiss(success,det);
    else{ 
        size_t rk=row_echelon_inner_elem(success);
        if(do_compute_vol)
            det=compute_vol(success);
        return rk;
    }
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::row_echelon(bool& success){
    
    Integer dummy;
    return row_echelon(success,false,dummy);
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::row_echelon(bool& success, Integer& det){
    
    return row_echelon(success,true,det);
}


//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::rank_destructive(){

    bool success;
    if(!do_arithmetic_check<Integer>()){
        size_t rk=row_echelon(success);
        if(!success){
            errorOutput()<<"Arithmetic failure in matrix operation. Most likely overflow.\n";
            throw ArithmeticException();
        }
        return rk;       
    }
    
    Matrix<Integer> Copy=*this;
    size_t rk=row_echelon(success);
    if(!success){
        Matrix<mpz_class> mpz_this(nr,nc);
        mat_to_mpz(Copy,mpz_this);
        rk=mpz_this.row_echelon(success);
    }
    return rk;                               
}
//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::rank() const{
    Matrix<Integer> N(*this);
    return N.rank_destructive();
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::vol_destructive(){

    assert(nr==nc);
    if(nr==0)
        return 1;

    bool success; 
    Integer det;  
        
    if(!do_arithmetic_check<Integer>()){
        row_echelon(success,det);
        if(!success){
            errorOutput()<<"Arithmetic failure in matrix operation. Most likely overflow.\n";
            throw ArithmeticException();
        }
        return det;     
    }
    
    Matrix<Integer> Copy=*this;
    row_echelon(success,det);
    if(success)
        return det;
        
    Matrix<mpz_class> mpz_this(nr,nc);
    mat_to_mpz(Copy,mpz_this);
    mpz_class vol=mpz_this.vol_destructive();
    return to_Int<Integer>(vol);
}
//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::vol() const{
    Matrix<Integer> N(*this);
    return N.vol_destructive();
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<key_t>  Matrix<Integer>::max_rank_submatrix_lex_inner(bool& success) const{

    success=true;
    size_t max_rank=min(nr,nc);
    Matrix<Integer> Test(max_rank,nc);
    Test.nr=0;
    vector<key_t> col;
    col.reserve(max_rank);
    vector<key_t> key;
    key.reserve(max_rank);
    
    vector<vector<bool> > col_done(max_rank,vector<bool>(nc,false));
    
    vector<Integer> Test_vec(nc);
     
    for(size_t i=0;i<nr;++i){    
        Test_vec=elem[i];            
        for(size_t k=0;k<Test.nr;++k){
            if(Test_vec[col[k]]==0)
                continue;
            Integer a=Test[k][col[k]];
            Integer b=Test_vec[col[k]];
            for(size_t j=0;j<nc;++j)
                if(!col_done[k][j]){
                Test_vec[j]=a*Test_vec[j]-b*Test[k][j];
                if (!check_range(Test_vec[j]) ) {
                    success=false;
                    return key;
                }
            }
        }
        
        size_t j=0;
        for(;j<nc;++j)
            if(Test_vec[j]!=0)
                break;
        if(j==nc)
            continue;
            
        col.push_back(j); 
        
        for(size_t k=0;k<Test.nr;++k)
            col_done[Test.nr][col[k]]=true;    
        
        v_make_prime(Test_vec);
        key.push_back(i);

        /* size_t rk=submatrix(key).rank();
        if(rk!=key.size()){
         cout << "ALARM" << endl;
         exit(0);
        }*/
        Test.nr++;
        Test[Test.nr-1]=Test_vec;
            
        if(key.size()==max_rank)
            break;   
    }    
    return key;                
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<key_t>  Matrix<Integer>::max_rank_submatrix_lex() const{

    bool success;
    vector<key_t> key=max_rank_submatrix_lex_inner(success);
    if(!success){
        Matrix<mpz_class> mpz_this(nr,nc);
        mat_to_mpz(*this,mpz_this);
        key=mpz_this.max_rank_submatrix_lex_inner(success);    
    }
    return key;
}


//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::solve_destructive_inner(bool ZZinvertible,Integer& denom) {

    assert(nc>=nr);
    size_t dim=nr;
    bool success;  
    
    size_t rk;
    
    if(ZZinvertible){
        rk=row_echelon_inner_elem(success); 
        if(!success)
            return false;        
        assert(rk==nr);
        denom=compute_vol(success);
    }
    else{
        rk=row_echelon(success,denom);
        if(!success)
            return false;    
    }

    // cout << "denom " << denom << endl<< "------------" << endl;
    if (denom==0) { 
        if(!do_arithmetic_check<Integer>()){
            errorOutput() << "Cannot solve system (denom=0)!" << endl;
            throw ArithmeticException();
        }
        else
            return false;            
    }

    Integer S;
    size_t i;
    long j;
    size_t k;
    for (i = nr; i < nc; i++) {
        for (j = dim-1; j >= 0; j--) {
            S=denom*elem[j][i];
            for (k = j+1; k < dim; k++) {
                S-=elem[j][k]*elem[k][i];
            }
            if(!check_range(S))
                return false;
            elem[j][i]=S/elem[j][j];
        }
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::solve_destructive_outer(bool ZZinvertible, Integer& denom) {

    if(!do_arithmetic_check<Integer>()){
        if(!solve_destructive_inner(ZZinvertible,denom)){
            errorOutput()<<"Arithmetic failure in matrix operation. Most likely overflow.\n";
            throw ArithmeticException();        
        }
        return;        
    }
    
    Matrix<Integer> Copy=*this;
    if(!solve_destructive_inner(ZZinvertible,denom)){  // <---------------------- eventuell Probe einbauen
        Matrix<mpz_class> mpz_this(nr,nc);
        mat_to_mpz(Copy,mpz_this);
        mpz_class mpz_denom;
        mpz_this.solve_destructive_inner(ZZinvertible, mpz_denom);
        mat_to_Int(mpz_this,*this);
        denom=to_Int<Integer>(mpz_denom);
    }    
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::solve_destructive(vector< Integer >& diagonal,Integer& denom) {

    solve_destructive_outer(true,denom);
    assert(diagonal.size()==nr);
    for(size_t i=0;i<nr;++i)
        diagonal[i]=elem[i][i];
}


//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::solve_destructive(Integer& denom){

    solve_destructive_outer(false,denom);
}


//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve_destructive(Matrix<Integer>& Right_side, vector< Integer >& diagonal, Integer& denom) {
    size_t dim=Right_side.nr;
    // cout << endl << "Sol.nc " << Solution.nc << " Sol.nr " << Solution.nr << " " << nr_sys << endl;
    assert(nr == nc);
    assert(nc == dim);
    assert(dim == diagonal.size());
    Matrix<Integer> M(nr,nc+Right_side.nc);
    for(size_t i=0;i<nr;++i){
        for(size_t j=0;j<nc;++j)
            M[i][j]=elem[i][j];
        for(size_t j=nc;j<M.nc;++j)
            M[i][j]=Right_side[i][j-nc];
    }
    M.solve_destructive(diagonal,denom);
    Matrix<Integer> Solution(Right_side.nr,Right_side.nc); 
    for(size_t i=0;i<nr;++i){
        for(size_t j=0;j<Right_side.nc;++j)
            Solution[i][j]=M[i][j+nc];    
    }
    return Solution;   
}
    

//--------------------------------------------------------------------------
/*
template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve_destructive(Matrix<Integer>& Right_side, vector< Integer >& diagonal, Integer& denom) {

    Matrix<Integer> Solution(Right_side.nr,Right_side.nc);  
    solve_destructive_Sol(Right_side,diagonal,denom,Solution);
    return Solution;
}
*/
//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve(const Matrix<Integer>& Right_side, vector< Integer >& diagonal, Integer& denom) const {
    Matrix<Integer> Left_side(*this);
    Matrix<Integer> Copy_Right_Side=Right_side;
    return Left_side.solve_destructive(Copy_Right_Side, diagonal, denom);
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve(const Matrix<Integer>& Right_side, Integer& denom) const {
    // Matrix<Integer> Left_side(*this);
    vector<Integer> dummy_diag(nr);
    return solve(Right_side, dummy_diag, denom);
}   

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::invert(vector< Integer >& diagonal, Integer& denom) const{
    assert(nr == nc);
    assert(nr == diagonal.size());
    Matrix<Integer> Left_side(*this);
    Matrix<Integer> Right_side(nr);

    return Left_side.solve_destructive(Right_side,diagonal,denom);
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::solve_rectangular(const vector<Integer>& v, Integer& denom) const {
    if (nc == 0 || nr == 0) { //return zero-vector as solution
        return vector<Integer>(nc,0);
    }
    size_t i;
    vector<key_t>  rows=max_rank_submatrix_lex();
    Matrix<Integer> Left_Side=submatrix(rows);
    assert(nc == Left_Side.nr); //otherwise input hadn't full rank //TODO 
    Matrix<Integer> Right_Side(v.size(),1);
    Right_Side.write_column(0,v);
    Right_Side = Right_Side.submatrix(rows);
    Matrix<Integer> Solution=Left_Side.solve(Right_Side, denom);
    vector<Integer> Linear_Form(nc);
    for (i = 0; i <nc; i++) {
        Linear_Form[i] = Solution[i][0];  // the solution vector is called Linear_Form
    }
    // cout << denom << " v= " << v;
    // denom/=v_make_prime(Linear_Form);
    vector<Integer> test = MxV(Linear_Form); // we have solved the system by taking a square submatrix
    // cout << denom << " v= " << v;         // now we must test whether the solution satisfies the full system
    // cout << denom << " t= " << test; 
    for (i = 0; i <nr; i++) {
        if (test[i] != denom * v[i]){
            return vector<Integer>();
        }
    }
    Integer total_gcd =gcd(denom,v_gcd(Linear_Form)); // extract the gcd of denom and solution
    denom/=total_gcd;
    v_scalar_division(Linear_Form,total_gcd);
    return Linear_Form;
}
//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::solve_ZZ(const vector<Integer>& v) const {

    Integer denom;
    vector<Integer> result=solve_rectangular(v,denom);
    if(denom!=1)
        result.clear();
    return result;
}
//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::find_linear_form() const {

    Integer denom;
    vector<Integer> result=solve_rectangular(vector<Integer>(nr,1),denom);
    v_make_prime(result);
    return result;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::find_linear_form_low_dim () const{
    size_t rank=(*this).rank();
    if (rank == 0) { //return zero-vector as linear form
        return vector<Integer>(nc,0);
    }
    if (rank == nc) { // basis change not necessary
        return (*this).find_linear_form();
    }

    Sublattice_Representation<Integer> Basis_Change(*this,true);
    vector<Integer> Linear_Form=Basis_Change.to_sublattice(*this).find_linear_form();
    if(Linear_Form.size()!=0)
        Linear_Form=Basis_Change.from_sublattice_dual(Linear_Form);

    return Linear_Form;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::kernel () const{
// computes a ZZ-basis of the solutions of (*this)x=0
// the basis is formed by the ROWS of the returned matrix

    size_t dim=nc;
    if(nr==0)
        return(Matrix<Integer>(dim));

    Matrix<Integer> Copy(*this);
    size_t rank;
    bool success;
    Matrix<Integer> Transf=Copy.row_column_trigonalize(rank,success);
    if(!success){
        Matrix<mpz_class> mpz_Copy(nr,nc);
        mat_to_mpz(*this,mpz_Copy);
        Matrix<mpz_class> mpz_Transf=mpz_Copy.row_column_trigonalize(rank,success);
        mat_to_Int(mpz_Transf,Transf);    
    }
    
    Matrix<Integer> ker_basis(dim-rank,dim);
    Matrix<Integer> Help =Transf.transpose();
    for (size_t i = rank; i < dim; i++) 
            ker_basis[i-rank]=Help[i];
    return(ker_basis);
}


//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Classless conversion routines
//---------------------------------------------------------------------------

template<typename Integer>
void mat_to_mpz(const Matrix<Integer>& mat, Matrix<mpz_class>& mpz_mat){
    size_t nrows=min(mat.nr_of_rows(),mpz_mat.nr_of_rows()); // we allow the matrices to have different sizes
    size_t ncols=min(mat.nr_of_columns(),mpz_mat.nr_of_columns());
    for(size_t i=0; i<nrows;++i)
        for(size_t j=0; j<ncols;++j)
            mpz_mat[i][j]=to_mpz(mat[i][j]);
}

template<typename Integer>
void mat_to_Int(const Matrix<mpz_class>& mpz_mat, Matrix<Integer>& mat){
    size_t nrows=min(mat.nr_of_rows(),mpz_mat.nr_of_rows()); // we allow the matrices to have different sizes
    size_t ncols=min(mat.nr_of_columns(),mpz_mat.nr_of_columns());
    for(size_t i=0; i<nrows;++i)
        for(size_t j=0; j<ncols;++j)
            mat[i][j]=to_Int<Integer>(mpz_mat[i][j]);
}
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Basement (routines currently not in use)
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
/*
// version with minimal remainder
template<typename Integer>
void Matrix<Integer>::reduce_row (size_t corner, Matrix<Integer>& Left) {
    assert(corner < nc);
    assert(corner < nr);
    assert(Left.nr == nr);
    size_t i,j;
    Integer help2=elem[corner][corner],quot,rem;
    for ( i = corner+1; i < nr; i++) {
        if (elem[i][corner]!=0) {
            minimal_remainder(elem[i][corner],help2,quot,rem);
            elem[i][corner]=rem;
            for (j = corner+1; j < nc; j++) {
                elem[i][j] -= quot*elem[corner][j];
                if ( !check_range(elem[i][j]) ) {
                    errorOutput()<<"Arithmetic failure in reduce_row. Most likely overflow.\n";
                    throw ArithmeticException();
                }
            }
            for (j = 0; j < Left.nc; j++) {
                Left.elem[i][j] -= quot*Left.elem[corner][j];
                if ( !check_range(Left.elem[i][j]) ) {
                    errorOutput()<<"Arithmetic failure in reduce_row. Most likely overflow.\n";
                    throw ArithmeticException();
                }
            }
        }
    }
}
*/

//---------------------------------------------------------------------------

/*
// version with minimal remainder
template<typename Integer>
void Matrix<Integer>::reduce_row (size_t row, size_t col) {
    assert(col < nc);
    assert(row < nr);
    size_t i,j;
    Integer quot, rem;
    for (i =row+1; i < nr; i++) {
        if (elem[i][col]!=0) {        
            minimal_remainder(elem[i][col], elem[row][col], quot, rem);
            elem[i][col]=rem;            
            for (j = col+1; j < nc; j++) {
                elem[i][j] -= quot*elem[row][j];
            }
        }
    }
}
*/



//---------------------------------------------------------------------------
/*
template<typename Integer>
Matrix<Integer> solve(const Matrix<Integer>& Left_side, const Matrix<Integer>& Right_side,Integer& denom){
    return Left_side.solve(Right_side,denom);
}
*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
Matrix<Integer> invert(const Matrix<Integer>& Left_side, vector< Integer >& diagonal, Integer& denom){
    return Left_side.invert(diagonal,denom);
}
*/

//---------------------------------------------------------------------------

/*
template<typename Integer>
bool Matrix<Integer>::reduce_column (size_t corner) {
    assert(corner < nc);
    assert(corner < nr);
    size_t i,j;
    Integer help1, help2=elem[corner][corner];
    for ( j = corner+1; j < nc; j++) {
        help1=elem[corner][j] / help2;
        if (help1!=0) {
            for (i = corner; i < nr; i++) {
                elem[i][j] -= help1*elem[i][corner];
                if ( !check_range(elem[i][j]) ) {
                    return false; 
                }
            }
        }
    }
    return true;
}
*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Matrix<Integer>::reduce_column (size_t corner, Matrix<Integer>& Right, Matrix<Integer>& Right_Inv) {
    assert(corner < nc);
    assert(corner < nr);
    assert(Right.nr == nc);
    assert(Right.nc == nc);
    assert(Right_Inv.nr == nc);
    assert(Right_Inv.nc ==nc);
    size_t i,j;
    Integer help1, help2=elem[corner][corner];
    for ( j = corner+1; j < nc; j++) {
        help1=elem[corner][j] / help2;
        if (help1!=0) {
            for (i = corner; i < nr; i++) {
                elem[i][j] -= help1*elem[i][corner];
                    if ( !check_range(elem[i][j]) ) {
                    return false;
                }
            }
            for (i = 0; i < nc; i++) {
                Right.elem[i][j] -= help1*Right.elem[i][corner];
                Right_Inv.elem[corner][i] += help1*Right_Inv.elem[j][i];
                if ( (!check_range(Right.elem[i][j])  ||
                        !check_range(Right_Inv.elem[corner][i]) ) ) {
                    return false;
               } 
            }
        }
    }
    return true;
}
*/

//---------------------------------------------------------------------------
/*
template<typename Integer>
vector<long> Matrix<Integer>::max_pivot(size_t corner){
    assert(corner < nc);
    assert(corner < nr);
    size_t i,j;
    Integer help=0;
    vector<long> v(2,-1);

    for (i = corner; i < nr; i++) {
        for (j = corner; j < nc; j++) {
            if (elem[i][j]!=0) {
                if ((help==0)||(Iabs(elem[i][j])>help)) {
                    help=Iabs(elem[i][j]);
                    v[0]=i;
                    v[1]=j;
                    if (help == 1) return v;
                }
            }
        }
    }
    
    return v;
}
*/

//---------------------------------------------------------------------------
/*
template<typename Integer>
size_t Matrix<Integer>::row_echelon_inner_gcd(bool compute_vol, bool& success){

    size_t pc=0;
    long piv=0, rk=0;
    success=true;
    
    for (rk = 0; rk < (long) nr; rk++){
        for(;pc<nc;pc++){
            piv=pivot_column(rk,pc);
            if(piv>=0)
                break;
        }
        if(pc==nc)
            break;
        exchange_rows (rk,piv);
        gcd_reduce_row(rk,pc);        
    }
    
    if(compute_vol)    
        do_compute_vol(success);        
    
    return rk;
}

*/

//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Matrix<Integer>::linear_comb_rows(const size_t& row,const size_t& i,const size_t& col,
            const Integer& u,const Integer& v,const Integer& w,const Integer&  z){
              
    for(size_t j=col;j<nc;++j){
        Integer rescue=elem[row][j];
        elem[row][j]=u*elem[row][j]+v*elem[i][j];
        elem[i][j]=w*rescue+z*elem[i][j];
        if ( (!check_range(elem[row][j])  || !check_range(elem[i][j]) )) {
            return false;
        }    
    }
    return true;
}
*/
//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Matrix<Integer>::gcd_reduce_row (size_t row, size_t col) {
    assert(row<nr);
    assert(col<nc);

    Integer u,v,w,z,d;      
    for(size_t i=row+1;i<nr;++i){
        if(elem[i][col]==0)
            continue;
        if(Iabs(elem[row][col])==1){
            if(elem[row][col]==-1)            
                v_scalar_multiplication<Integer>(elem[row],-1);
            Integer b=elem[i][col];
            for(size_t j=col;j<nc;++j)
                elem[i][j]-=b*elem[row][j];
            continue;
        }
        d=ext_gcd(elem[row][col],elem[i][col],u,v);
        w=-elem[i][col]/d;
        z=elem[row][col]/d; 
        // Now we multiply the submatrix formed by rows "row" and "i" and columns
        // col,...,nc from the left  by the 2x2 matrix
        // | u v |
        // | w z |      
        if(!linear_comb_rows(row,i,col,u,v,w,z))
            return false;
    }
    return true;
}
*/

/*
template<typename Integer>
vector<key_t>  Matrix<Integer>::max_rank_submatrix_lex() const{
    
    Matrix<Integer> T(nc,nr);
    for(size_t i=0;i<nr;++i)
        for(size_t j=0;j<nc;++j)
            T[j][i]=elem[i][j];
    size_t rk=T.row_echelon(false);
    vector<key_t> key(rk);
    for(size_t i=0;i<rk;++i){
        size_t j=0;
        for(;j<T.nc;++j)
            if(T[i][j]!=0)
                break;
        key[i]=j;    
    }
    return key;
    
}
*/

//---------------------------------------------------------------------------
/*
template<typename Integer>
bool Matrix<Integer>::reduce_row (size_t corner, Matrix<Integer>& RHS) {
    assert(corner < nc);
    assert(corner < nr);
    assert(RHS.nr == nr);
    size_t i,j;
    Integer help1, help2=elem[corner][corner];
    for ( i = corner+1; i < nr; i++) {
        if (elem[i][corner]!=0) {
            help1=elem[i][corner] / help2;
            for (j = corner; j < nc; j++) {
                elem[i][j] -= help1*elem[corner][j];
                if (!check_range(elem[i][j]) ){ 
                    return false;
                } 
            }
            for (j = 0; j < RHS.nc; j++) {
                RHS.elem[i][j] -= help1*RHS.elem[corner][j];
                if ( !check_range(RHS.elem[i][j]) ) {
                }
            }
        }
    }
    return true;
}
*/
}  // namespace
