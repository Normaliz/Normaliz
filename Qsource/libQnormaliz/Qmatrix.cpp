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
#include <set>
#include <algorithm>
#include <math.h>

#include "libQnormaliz/Qmatrix.h"
#include "libQnormaliz/Qvector_operations.h"
#include "libQnormaliz/Qnormaliz_exception.h"
#include "libQnormaliz/Qsublattice_representation.h"

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//Public
//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>::Matrix(){
    nr=0;
    nc=0;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>::Matrix(size_t dim){
    nr=dim;
    nc=dim;
    elem = vector< vector<Number> >(dim, vector<Number>(dim));
    for (size_t i = 0; i < dim; i++) {
        elem[i][i]=1;
    }
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>::Matrix(size_t row, size_t col){
    nr=row;
    nc=col;
    elem = vector< vector<Number> >(row, vector<Number>(col));
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>::Matrix(size_t row, size_t col, Number value){
    nr=row;
    nc=col;
    elem = vector< vector<Number> > (row, vector<Number>(col,value));
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>::Matrix(const vector< vector<Number> >& new_elem){
    nr=new_elem.size();
    if (nr>0) {
        nc=new_elem[0].size();
        elem=new_elem;
        //check if all rows have the same length
        for (size_t i=1; i<nr; i++) {
            if (elem[i].size() != nc) {
                throw BadInputException("Inconsistent lengths of rows in matrix!");
            }
        }
    } else {
        nc=0;
    }
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>::Matrix(const list< vector<Number> >& new_elem){
    nr = new_elem.size();
    elem = vector< vector<Number> > (nr);
    nc = 0;
    size_t i=0;
    typename list< vector<Number> >::const_iterator it=new_elem.begin();
    for(; it!=new_elem.end(); ++it, ++i) {
        if(i == 0) {
            nc = (*it).size();
        } else {
            if ((*it).size() != nc) {
                throw BadInputException("Inconsistent lengths of rows in matrix!");
            }
        }
        elem[i]=(*it);
    }
}

//---------------------------------------------------------------------------
/*
template<typename Number>
void Matrix<Number>::write(istream& in){
    size_t i,j;
    for(i=0; i<nr; i++){
        for(j=0; j<nc; j++) {
            in >> elem[i][j];
        }
    }
}
*/
//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::write_column(size_t col, const vector<Number>& data){
    assert(col >= 0);
    assert(col < nc); 
    assert(nr == data.size());

    for (size_t i = 0; i < nr; i++) {
        elem[i][col]=data[i];
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::print(const string& name,const string& suffix) const{
    string file_name = name+"."+suffix;
    const char* file = file_name.c_str();
    ofstream out(file);
    print(out);
    out.close();
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::print_append(const string& name,const string& suffix) const{
    string file_name = name+"."+suffix;
    const char* file = file_name.c_str();
    ofstream out(file,ios_base::app);
    print(out);
    out.close();
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::print(ostream& out) const{
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

template<typename Number>
void Matrix<Number>::pretty_print(ostream& out, bool with_row_nr) const{
    if(nr>1000000 && !with_row_nr){
        print(out);
        return;
    }
    size_t i,j,k;
    vector<size_t> max_length = maximal_decimal_length_columnwise();
    size_t max_index_length = decimal_length(nr);
    for (i = 0; i < nr; i++) {
        if (with_row_nr) {
            for (k= 0; k <= max_index_length - decimal_length(i); k++) {
                out<<" ";
            }
            out << i << ": ";
        }
        for (j = 0; j < nc; j++) {
            ostringstream to_print;
            to_print << elem[i][j];
            for (k= 0; k <= max_length[j] - to_print.str().size(); k++) {
                out<<" ";
            }
            out<< to_print.str();
        }
        out<<endl;
    }
}

/*
 * string to_print;
            ostringstream(to_print) << elem[i][j];
            cout << elem[i][j] << " S " << to_print << " L " << decimal_length(elem[i][j]) << endl;
            for (k= 0; k <= max_length[j] - to_print.size(); k++) {
                out<<" ";
            }
            out << to_print;
*/
//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::nr_of_rows () const{
    return nr;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::nr_of_columns () const{
    return nc;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::set_nr_of_columns(size_t c){
    nc=c;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::random (int mod) {
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
  
template<typename Number>
void Matrix<Number>::set_zero() {
    size_t i,j;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            elem[i][j] = 0;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::select_submatrix(const Matrix<Number>& mother, const vector<key_t>& rows){

    assert(nr>=rows.size());
    assert(nc>=mother.nc);
    
    size_t size=rows.size(), j;
    for (size_t i=0; i < size; i++) {
        j=rows[i];
        for(size_t k=0;k<mother.nc;++k)
            elem[i][k]=mother[j][k];
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::select_submatrix_trans(const Matrix<Number>& mother, const vector<key_t>& rows){

    assert(nc>=rows.size());
    assert(nr>=mother.nc);
    
    size_t size=rows.size(), j;
    for (size_t i=0; i < size; i++) {
        j=rows[i];
        for(size_t k=0;k<mother.nc;++k)
            elem[k][i]=mother[j][k];
    }
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::submatrix(const vector<key_t>& rows) const{
    size_t size=rows.size(), j;
    Matrix<Number> M(size, nc);
    for (size_t i=0; i < size; i++) {
        j=rows[i];
        assert(j >= 0);
        assert(j < nr);
        M.elem[i]=elem[j];
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::submatrix(const vector<int>& rows) const{
    size_t size=rows.size(), j;
    Matrix<Number> M(size, nc);
    for (size_t i=0; i < size; i++) {
        j=rows[i];
        assert(j >= 0);
        assert(j < nr);
        M.elem[i]=elem[j];
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::submatrix(const vector<bool>& rows) const{
    assert(rows.size() == nr);
    size_t size=0;
    for (size_t i = 0; i <rows.size(); i++) {
        if (rows[i]) {
            size++;
        }
    }
    Matrix<Number> M(size, nc);
    size_t j = 0;
    for (size_t i = 0; i < nr; i++) {
        if (rows[i]) {
            M.elem[j++] = elem[i];
        }
    }
    return M;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number>& Matrix<Number>::remove_zero_rows() {
    size_t from = 0, to = 0; // maintain to <= from
    while (from < nr && v_is_zero(elem[from])) from++; //skip zero rows
    while (from < nr) {  // go over matrix
        // now from is a non-zero row
        if (to != from) elem[to].swap(elem[from]);
        ++to; ++from;
        while (from < nr && v_is_zero(elem[from])) from++; //skip zero rows
    }
    nr = to;
    elem.resize(nr);
    return *this;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::swap(Matrix<Number>& x) {
    size_t tmp = nr; nr = x.nr; x.nr = tmp;
    tmp = nc; nc = x.nc; x.nc = tmp;
    elem.swap(x.elem);
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::resize(size_t nr_rows, size_t nr_cols) {
    nc = nr_cols; //for adding new rows with the right length
    resize(nr_rows);
    resize_columns(nr_cols);
}

template<typename Number>
void Matrix<Number>::resize(size_t nr_rows) {
    if (nr_rows > elem.size()) {
        elem.resize(nr_rows, vector<Number>(nc));
    }
    nr = nr_rows;
}

template<typename Number>
void Matrix<Number>::resize_columns(size_t nr_cols) {
    for (size_t i=0; i<nr; i++) {
        elem[i].resize(nr_cols);
    }
    nc = nr_cols;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::diagonal() const{
    assert(nr == nc); 
    vector<Number> diag(nr);
    for(size_t i=0; i<nr;i++){
        diag[i]=elem[i][i];
    }
    return diag;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::maximal_decimal_length() const{
    size_t i,maxim=0;
    vector<size_t> maxim_col;
    maxim_col=maximal_decimal_length_columnwise();
    for (i = 0; i <nr; i++)
        maxim=max(maxim,maxim_col[i]);
    return maxim;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<size_t> Matrix<Number>::maximal_decimal_length_columnwise() const{
    size_t i,j=0;
    vector<size_t> maxim(nc,0);
    for (i = 0; i <nr; i++) {
        for (j = 0; j <nc; j++) {
            maxim[j]=max(maxim[j],decimal_length(elem[i][j]));
/*            if(elem[i][j]<0){
                if(elem[i][j]<neg_max[j])
                    neg_max[j]=elem[i][j];
                continue;
            }
            if(elem[i][j]>pos_max[j])
                pos_max[j]=elem[i][j];
*/
        }
    }
    /* for(size_t j=0;j<nc;++j)
        maxim[j]=max(decimal_length(neg_max[j]),decimal_length(pos_max[j])); */
    return maxim;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::append(const Matrix<Number>& M) {
    assert (nc == M.nc);
    elem.reserve(nr+M.nr);
    for (size_t i=0; i<M.nr; i++) {
        elem.push_back(M.elem[i]);
    }
    nr += M.nr;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::append(const vector<vector<Number> >& M) {
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

template<typename Number>
void Matrix<Number>::append(const vector<Number>& V) {
    assert (nc == V.size());
    elem.push_back(V);
    nr++;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::append_column(const vector<Number>& v) {
    assert (nr == v.size());
    for (size_t i=0; i<nr; i++) {
        elem[i].resize(nc+1);
        elem[i][nc] = v[i];
    }
    nc++;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::remove_row(const vector<Number>& row) {
    size_t tmp_nr = nr;
    for (size_t i = 1; i <= tmp_nr; ++i) {
        if (elem[tmp_nr-i] == row) {
            elem.erase(elem.begin()+(tmp_nr-i));
            nr--;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::remove_duplicate_and_zero_rows() {
    bool remove_some = false;
    vector<bool> key(nr, true);

    set<vector<Number> > SortedRows;
    SortedRows.insert( vector<Number>(nc,0) );
    typename set<vector<Number> >::iterator found;
    for (size_t i = 0; i<nr; i++) {
        found = SortedRows.find(elem[i]);
        if (found != SortedRows.end()) {
            key[i] = false;
            remove_some = true;
        }
        else
            SortedRows.insert(found,elem[i]);
    }

    if (remove_some) {
        *this = submatrix(key);
    }
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::add(const Matrix<Number>& A) const{
    assert (nr == A.nr);
    assert (nc == A.nc);
    
    Matrix<Number> B(nr,nc);
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            B.elem[i][j]=elem[i][j]+A.elem[i][j];
        }
    }
    return B;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::multiplication(const Matrix<Number>& A) const{
    assert (nc == A.nr);

    Matrix<Number> B(nr,A.nc,0);  //initialized with 0
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

template<typename Number>
Matrix<Number> Matrix<Number>::multiplication_cut(const Matrix<Number>& A, const size_t& c) const{
    assert (nc == A.nr);
    assert(c<= A.nc);

    Matrix<Number> B(nr,c,0);  //initialized with 0
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
/*
template<typename Number>
Matrix<Number> Matrix<Number>::multiplication(const Matrix<Number>& A, long m) const{
    assert (nc == A.nr);

    Matrix<Number> B(nr,A.nc,0);  //initialized with 0
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
*/

//---------------------------------------------------------------------------

template<typename Number>
bool Matrix<Number>::equal(const Matrix<Number>& A) const{
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
/*
template<typename Number>
bool Matrix<Number>::equal(const Matrix<Number>& A, long m) const{
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
} */

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::transpose()const{
    Matrix<Number> B(nc,nr);
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            B.elem[j][i]=elem[i][j];
        }
    }
    return B;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::scalar_multiplication(const Number& scalar){
    size_t i,j;
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            elem[i][j] *= scalar;
        }
    }
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::scalar_division(const Number& scalar){
    size_t i,j;
    assert(scalar != 0);
    for(i=0; i<nr;i++){
        for(j=0; j<nc; j++){
            // assert (elem[i][j]%scalar == 0);
            elem[i][j] /= scalar;
        }
    }
}

//---------------------------------------------------------------------------

/*
template<typename Number>
void Matrix<Number>::reduction_modulo(const Number& modulo){
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
*/


//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::simplify_rows() {
    // vector<Number> g(nr);
    vector<Number> dummy;
    for (size_t i = 0; i <nr; i++) {
        v_simplify(elem[i],dummy);
    }
    // return g;
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::multiply_rows(const vector<Number>& m) const{  //row i is multiplied by m[i]
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

template<typename Number>
void Matrix<Number>::MxV(vector<Number>& result, const vector<Number>& v) const{
    assert (nc == v.size());
    result.resize(nr);
    for(size_t i=0; i<nr;i++){
        result[i]=v_scalar_product(elem[i],v);
    }
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::MxV(const vector<Number>& v) const{
    vector<Number> w(nr);
    MxV(w, v);
    return w;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::VxM(const vector<Number>& v) const{
    assert (nr == v.size());
    vector<Number> w(nc,0);
    size_t i,j;
    for (i=0; i<nc; i++){
        for (j=0; j<nr; j++){
            w[i] += v[j]*elem[j][i];
        }
        if(!check_range(w[i]))
            break;
    }

        return w;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::VxM_div(const vector<Number>& v, const Number& divisor, bool& success) const{
    assert (nr == v.size());
    vector<Number> w(nc,0);
    success=true;
    size_t i,j;
    for (i=0; i<nc; i++){
        for (j=0; j<nr; j++){
            w[i] += v[j]*elem[j][i];
        }
        if(!check_range(w[i])){
            success=false;
            break;
        }
    }

    if(success)      
        v_scalar_division(w,divisor);  
        
    return w;
}

//---------------------------------------------------------------------------

template<typename Number>
bool Matrix<Number>::is_diagonal() const{

    for(size_t i=0;i<nr;++i)
        for(size_t j=0;j<nc;++j)
            if(i!=j && elem[i][j]!=0)
                return false;
    return true;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<long> Matrix<Number>::pivot(size_t corner){
    assert(corner < nc);
    assert(corner < nr);
    size_t i,j;
    Number help=0;
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

template<typename Number>
long Matrix<Number>::pivot_column(size_t row,size_t col){
    assert(col < nc);
    assert(row < nr);
    size_t i;
    long j=-1;
    Number help=0;

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

template<typename Number>
long Matrix<Number>::pivot_column(size_t col){
    return pivot_column(col,col);
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::exchange_rows(const size_t& row1, const size_t& row2){
    if (row1 == row2) return;
    assert(row1 < nr);
    assert(row2 < nr);
    elem[row1].swap(elem[row2]);
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::exchange_columns(const size_t& col1, const size_t& col2){
    if (col1 == col2) return;
    assert(col1 < nc);
    assert(col2 < nc);
    for(size_t i=0; i<nr;i++){
        std::swap(elem[i][col1], elem[i][col2]);
    }
}

//---------------------------------------------------------------------------
 
template<typename Number>
bool Matrix<Number>::reduce_row (size_t row, size_t col) {
    assert(col < nc);
    assert(row < nr);
    size_t i,j;
    Number help;
    for (i =row+1; i < nr; i++) {
        if (elem[i][col]!=0) {
            help=elem[i][col] / elem[row][col];
            for (j = col; j < nc; j++) {
                elem[i][j] -= help*elem[row][j];
                if (!check_range(elem[i][j]) ) {
                    return false;
                }
            }
            // v_el_trans<Number>(elem[row],elem[i],-help,col);
        }
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Number>
bool  Matrix<Number>::reduce_row (size_t corner) {
    return reduce_row(corner,corner);
}

//---------------------------------------------------------------------------
 
template<typename Number>
bool Matrix<Number>::reduce_rows_upwards () {
// assumes that "this" is in row echelon form
// and reduces eevery column in which the rank jumps 
// by its lowest element
//
// Aplies v_simplify to make rows "nice"
    if(nr==0)
        return true;

    for(size_t row=0;row<nr;++row){
        size_t col;
        for(col=0;col<nc;++col)
            if(elem[row][col]!=0)
                break;
        if(col==nc) // zero row
            continue;
        if(elem[row][col]<0)
            v_scalar_multiplication<Number>(elem[row],-1); // make corner posizive
        
        for(long i=row-1;i>=0;--i){
            Number quot;            
            //minimal_remainder(elem[i][col],elem[row][col],quot,rem);
            quot=elem[i][col]/elem[row][col];
            elem[i][col]=0; // rem
            for(size_t j=col+1;j<nc;++j){
                elem[i][j]-=quot* elem[row][j];
            }                                           
        }
    }
    
    simplify_rows();
           
    return true;
}

//---------------------------------------------------------------------------
 
template<typename Number>
bool Matrix<Number>::linear_comb_columns(const size_t& col,const size_t& j,
            const Number& u,const Number& w,const Number& v,const Number& z){
                       
    for(size_t i=0;i<nr;++i){
        Number rescue=elem[i][col];
        elem[i][col]=u*elem[i][col]+v*elem[i][j];
        elem[i][j]=w*rescue+z*elem[i][j];
        if ( (!check_range(elem[i][col])  || !check_range(elem[i][j]) )) {
            return false;
        }        
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Number>
bool Matrix<Number>::gcd_reduce_column (size_t corner, Matrix<Number>& Right){
    assert(corner < nc);
    assert(corner < nr);
    Number d,u,w,z,v;
    for(size_t j=corner+1;j<nc;++j){
       d =elem[corner][corner],elem[corner]; // ext_gcd(elem[corner][corner],elem[corner][j],u,v);
       u=1;
       v=0;
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

template<typename Number>
bool Matrix<Number>::column_trigonalize(size_t rk, Matrix<Number>& Right) { 
    assert(Right.nr == nc);
    assert(Right.nc == nc);
    vector<long> piv(2,0);       
    for(size_t j=0;j<rk;++j){
            piv=pivot(j);
            assert(piv[0]>=0); // protect against wrong rank
            exchange_rows (j,piv[0]);
            exchange_columns (j,piv[1]);
            Right.exchange_columns(j,piv[1]);
            if(!gcd_reduce_column(j, Right))
                return false;
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Number>
Number Matrix<Number>::compute_vol(bool& success){
        
    assert(nr<=nc);
    
    Number det=1;
    for(size_t i=0;i<nr;++i){
        det*=elem[i][i]; 
        if(!check_range(det)){
            success=false;
            return 0;
        }
    }
            
    det=Iabs(det);
    success=true;
    return det;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon_inner_elem(bool& success){

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

/*
template<typename Number>
size_t Matrix<Number>::row_echelon_inner_bareiss(bool& success, Number& det){
// no overflow checks since this is supposed to be only used with GMP

    success=true;
    if(nr==0)
        return 0;
    assert(using_GMP<Number>());

    size_t pc=0;
    long piv=0, rk=0;
    vector<bool> last_time_mult(nr,false),this_time_mult(nr,false);
    Number last_div=1,this_div=1;
    size_t this_time_exp=0,last_time_exp=0;
    Number det_factor=1;
    
    for (rk = 0; rk < (long) nr; rk++){

        for(;pc<nc;pc++){
            piv=pivot_column(rk,pc);
            if(piv>=0)
                break;
        }
        if(pc==nc)
            break;
                        
        if(!last_time_mult[piv]){
            for(size_t k=rk;k<nr;++k)
                if(elem[k][pc]!=0 && last_time_mult[k]){
                    piv=k;
                    break;                
                }        
        }        
        
        exchange_rows (rk,piv);
        v_bool_entry_swap(last_time_mult,rk,piv);
        
        if(!last_time_mult[rk])
            for(size_t i=0;i<nr;++i)
                last_time_mult[i]=false;
                     
        Number a=elem[rk][pc];
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
            Number b=elem[i][pc];
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
        
        for(size_t i=0;i<last_time_exp;++i)
            det_factor*=last_div;
        last_time_mult=this_time_mult;
        last_div=this_div;
        last_time_exp=this_time_exp;
    }
    
    det=0;
    if (nr <= nc && rk == (long) nr) { // must allow nonsquare matrices
        det=1;
        for(size_t i=0;i<nr;++i)
            det*=elem[i][i];            
        det=Iabs<Number>(det/det_factor);        
    }
    
    return rk;
}
*/

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon_reduce(bool& success){

    size_t rk=row_echelon_inner_elem(success);
    if(success)
        success=reduce_rows_upwards();
    return rk;
}

//---------------------------------------------------------------------------

template<typename Number>
Number Matrix<Number>::full_rank_index(bool& success){

    size_t rk=row_echelon_inner_elem(success);
    Number index=1;
    if(success){
        for(size_t i=0;i<rk;++i){
            index*=elem[i][i];
            if(!check_range(index)){
                success=false;
                index=0;
                return index;
            }
        }
    }
    assert(rk==nc); // must have full rank
    index=Iabs(index);
    return index;
}
//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::row_column_trigonalize(size_t& rk, bool& success) {

    Matrix<Number> Right(nc);
    rk=row_echelon_reduce(success);
    if(success)
        success=column_trigonalize(rk,Right); 
    return Right; 
} 

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon(bool& success, bool do_compute_vol, Number& det){
    
/*    if(using_GMP<Number>()){
        return row_echelon_inner_bareiss(success,det);;
    }
    else{ */
        size_t rk=row_echelon_inner_elem(success);
        if(do_compute_vol)
            det=compute_vol(success);
        return rk;
//    }
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon(bool& success){
    
    Number dummy;
    return row_echelon(success,false,dummy);
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon(bool& success, Number& det){
    
    return row_echelon(success,true,det);
}



//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::rank_submatrix(const Matrix<Number>& mother, const vector<key_t>& key){

    assert(nc>=mother.nc);
    if(nr<key.size()){
        elem.resize(key.size(),vector<Number>(nc,0));
        nr=key.size();    
    }
    size_t save_nr=nr;
    size_t save_nc=nc;
    nr=key.size();
    nc=mother.nc;

    select_submatrix(mother,key);

    bool success;
    size_t rk=row_echelon(success);
    
    nr=save_nr;
    nc=save_nc;
    return rk;                               
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::rank_submatrix(const vector<key_t>& key) const{

    Matrix<Number> work(key.size(),nc);
    return work.rank_submatrix(*this,key);              
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::rank() const{
    vector<key_t> key(nr);
    for(size_t i=0;i<nr;++i)
        key[i]=i;
    return rank_submatrix(key);
}

//---------------------------------------------------------------------------

template<typename Number>
Number Matrix<Number>::vol_submatrix(const Matrix<Number>& mother, const vector<key_t>& key){

    assert(nc>=mother.nc);
    if(nr<key.size()){
        elem.resize(key.size(),vector<Number>(nc,0));
        nr=key.size();    
    }
    size_t save_nr=nr;
    size_t save_nc=nc;
    nr=key.size();
    nc=mother.nc;
    
    select_submatrix(mother,key);

    bool success;
    Number det;
    row_echelon(success,det);
    
    nr=save_nr;
    nc=save_nc;
    return det;                               
}
//---------------------------------------------------------------------------

template<typename Number>
Number Matrix<Number>::vol_submatrix(const vector<key_t>& key) const{

    Matrix<Number> work(key.size(),nc);
    return work.vol_submatrix(*this,key);              
}

//---------------------------------------------------------------------------

template<typename Number>
Number Matrix<Number>::vol() const{
    vector<key_t> key(nr);
    for(size_t i=0;i<nr;++i)
        key[i]=i;
    return vol_submatrix(key);
}

//---------------------------------------------------------------------------

template<typename Number>
vector<key_t>  Matrix<Number>::max_rank_submatrix_lex_inner(bool& success) const{
    
    vector<Number> dummy;

    success=true;
    size_t max_rank=min(nr,nc);
    Matrix<Number> Test(max_rank,nc);
    Test.nr=0;
    vector<key_t> col;
    col.reserve(max_rank);
    vector<key_t> key;
    key.reserve(max_rank);
    size_t rk=0;
    
    vector<vector<bool> > col_done(max_rank,vector<bool>(nc,false));
    
    vector<Number> Test_vec(nc);
     
    for(size_t i=0;i<nr;++i){    
        Test_vec=elem[i];            
        for(size_t k=0;k<rk;++k){
            if(Test_vec[col[k]]==0)
                continue;
            Number a=Test[k][col[k]];
            Number b=Test_vec[col[k]];
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
        if(j==nc)     // Test_vec=0
            continue;
            
        col.push_back(j);
        key.push_back(i);
        
        if(rk>0){
            col_done[rk]=col_done[rk-1];
            col_done[rk][col[rk-1]]=true;
        }

        Test.nr++;
        rk++;
        v_simplify(Test_vec,dummy);
        Test[rk-1]=Test_vec;
            
        if(rk==max_rank)
            break;   
    }    
    return key;                
}

//---------------------------------------------------------------------------

template<typename Number>
vector<key_t>  Matrix<Number>::max_rank_submatrix_lex() const{
    bool success;
    vector<key_t> key=max_rank_submatrix_lex_inner(success);
    return key;
}

//---------------------------------------------------------------------------

template<typename Number>
bool Matrix<Number>::solve_destructive_inner(bool ZZinvertible,Number& denom) {

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

    if (denom==0) { 
        if(using_GMP<Number>()){
            errorOutput() << "Cannot solve system (denom=0)!" << endl;
            throw ArithmeticException();
        }
        else
            return false;            
    }

    Number S;
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

template<typename Number>
void Matrix<Number>::customize_solution(size_t dim, Number& denom, size_t red_col, 
                     size_t sign_col, bool make_sol_prime) {
    
    return;
                         
 /*   assert(!(make_sol_prime && (sign_col>0 || red_col>0)));

    for(size_t j=0;j<red_col;++j){  // reduce first red_col columns of solution mod denom
        for(size_t k=0;k<dim;++k){
          elem[k][dim+j]%=denom;
          if(elem[k][dim+j]<0)
              elem[k][dim+j]+=Iabs(denom);
       }
    }
        
    for(size_t j=0;j<sign_col;++j)   // replace entries in the next sign_col columns by their signs
      for(size_t k=0;k<dim;++k){
        if(elem[k][dim+red_col+j]>0){
            elem[k][dim+red_col+j]=1;
            continue;
        } 
        if(elem[k][dim+red_col+j]<0){
            elem[k][dim+red_col+j]=-1;
            continue;
        }       
      }
      
    if(make_sol_prime) // make columns of solution coprime if wanted
        make_cols_prime(dim,nc-1); */
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::solve_system_submatrix_outer(const Matrix<Number>& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
        Number& denom, bool ZZ_invertible, bool transpose, size_t red_col, size_t sign_col, 
        bool compute_denom, bool make_sol_prime) {
     
    size_t dim=mother.nc;
    assert(key.size()==dim);
    assert(nr==dim);
    assert(dim+RS.size()<=nc);
    size_t save_nc=nc;
    nc=dim+RS.size();
    
    if(transpose)
       select_submatrix_trans(mother,key);           
    else
       select_submatrix(mother,key);
               
    for(size_t i=0;i<dim;++i)
       for(size_t k=0;k<RS.size();++k)
           elem[i][k+dim]= (*RS[k])[i];
    
    if(solve_destructive_inner(ZZ_invertible,denom)){
        customize_solution(dim, denom,red_col,sign_col,make_sol_prime);        
    } /*
    else{          
       #pragma omp atomic
       GMP_mat++;
    
       Matrix<mpz_class> mpz_this(nr,nc);
       mpz_class mpz_denom;
       if(transpose)
           mpz_submatrix_trans(mpz_this,mother,key);
       else            
           mpz_submatrix(mpz_this,mother,key);
           
       for(size_t i=0;i<dim;++i)
           for(size_t k=0;k<RS.size();++k)
               convert(mpz_this[i][k+dim], (*RS[k])[i]);
       mpz_this.solve_destructive_inner(ZZ_invertible,mpz_denom);
       mpz_this.customize_solution(dim, mpz_denom,red_col,sign_col,make_sol_prime);           
          
       for(size_t i=0;i<dim;++i)  // replace left side by 0, except diagonal if ZZ_invetible
          for(size_t j=0;j<dim;++j){
            if(i!=j || !ZZ_invertible)
                mpz_this[i][j]=0;              
          }
              
       mat_to_Int(mpz_this,*this);
       if(compute_denom)
           convert(denom, mpz_denom);                
    }*/
    nc=save_nc;      
}


//---------------------------------------------------------------------------


template<typename Number>
void Matrix<Number>::solve_system_submatrix(const Matrix<Number>& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         vector< Number >& diagonal, Number& denom, size_t red_col, size_t sign_col) {

    solve_system_submatrix_outer(mother,key,RS,denom,true,false,red_col,sign_col);
    assert(diagonal.size()==nr);
    for(size_t i=0;i<nr;++i)
        diagonal[i]=elem[i][i];
                 
}



//---------------------------------------------------------------------------
// the same without diagonal
template<typename Number>
void Matrix<Number>::solve_system_submatrix(const Matrix<Number>& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         Number& denom, size_t red_col, size_t sign_col, bool compute_denom, bool make_sol_prime) {

    solve_system_submatrix_outer(mother,key,RS,denom,false,false,red_col,sign_col, 
                compute_denom, make_sol_prime);
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::solve_system_submatrix_trans(const Matrix<Number>& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         Number& denom, size_t red_col, size_t sign_col) {
         
    solve_system_submatrix_outer(mother,key,RS,denom,false,true,red_col,sign_col);
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::extract_solution() const {
    assert(nc>=nr);
    Matrix<Number> Solution(nr,nc-nr); 
    for(size_t i=0;i<nr;++i){
        for(size_t j=0;j<Solution.nc;++j)
            Solution[i][j]=elem[i][j+nr];    
    }
    return Solution;  
}

//---------------------------------------------------------------------------

template<typename Number>
vector<vector<Number>* > Matrix<Number>::row_pointers(){

    vector<vector<Number>* > pointers(nr);
    for(size_t i=0;i<nr;++i)
        pointers[i]=&(elem[i]);
    return pointers;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<vector<Number>* > Matrix<Number>::submatrix_pointers(const vector<key_t>& key){

    vector<vector<Number>* > pointers(key.size());
    for(size_t i=0;i<key.size();++i)
        pointers[i]=&(elem[key[i]]);
    return pointers;
}
//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::solve(const Matrix<Number>& Right_side,vector<Number>& diagonal,Number& denom) const {

    Matrix<Number> M(nr,nc+Right_side.nc);
    vector<key_t> key=identity_key(nr);
    Matrix<Number> RS_trans=Right_side.transpose();
    vector<vector<Number>* > RS=RS_trans.row_pointers();
    M.solve_system_submatrix(*this,key,RS,diagonal,denom,0,0);
    return M.extract_solution(); 
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::solve(const Matrix<Number>& Right_side, Number& denom) const {

    Matrix<Number> M(nr,nc+Right_side.nc);
    vector<key_t> key=identity_key(nr);
    Matrix<Number> RS_trans=Right_side.transpose();
    vector<vector<Number>* > RS=RS_trans.row_pointers();
    M.solve_system_submatrix(*this,key,RS,denom,0,0);
    return M.extract_solution(); 
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::invert(Number& denom) const{
    assert(nr == nc);
    Matrix<Number> Right_side(nr);

    return solve(Right_side,denom);
}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::bundle_matrices(const Matrix<Number>& Right_side) const {

    assert(nr == nc);
    assert(nc == Right_side.nr);
    Matrix<Number> M(nr,nc+Right_side.nc);
    for(size_t i=0;i<nr;++i){
        for(size_t j=0;j<nc;++j)
            M[i][j]=elem[i][j];
        for(size_t j=nc;j<M.nc;++j)
            M[i][j]=Right_side[i][j-nc];
    }
    return M;
}
//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::invert_unprotected(Number& denom, bool& success) const{
    assert(nr == nc);
    Matrix<Number> Right_side(nr);
    Matrix<Number> M=bundle_matrices(Right_side);
    success=M.solve_destructive_inner(false,denom);
    return M.extract_solution();;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::invert_submatrix(const vector<key_t>& key, Number& denom, Matrix<Number>& Inv, bool compute_denom, bool make_sol_prime) const{
    assert(key.size() == nc);
    Matrix<Number> unit_mat(key.size());
    Matrix<Number> M(key.size(),2*key.size());        
    vector<vector<Number>* > RS_pointers=unit_mat.row_pointers();
    M.solve_system_submatrix(*this,key,RS_pointers,denom,0,0, compute_denom, make_sol_prime);
    Inv=M.extract_solution();;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::simplex_data(const vector<key_t>& key, Matrix<Number>& Supp, Number& vol, bool compute_vol) const{
    assert(key.size() == nc);
    invert_submatrix(key,vol,Supp,compute_vol,true);
    Supp=Supp.transpose();
    // Supp.make_prime(); now done internally -- but not in Q !! Therefore
    Supp.simplify_rows();
}
//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::solve_rectangular(const vector<Number>& v, Number& denom) const {
    if (nc == 0 || nr == 0) { //return zero-vector as solution
        return vector<Number>(nc,0);
    }
    size_t i;
    vector<key_t>  rows=max_rank_submatrix_lex();
    Matrix<Number> Left_Side=submatrix(rows);
    assert(nc == Left_Side.nr); //otherwise input hadn't full rank //TODO 
    Matrix<Number> Right_Side(v.size(),1);
    Right_Side.write_column(0,v);
    Right_Side = Right_Side.submatrix(rows);
    Matrix<Number> Solution=Left_Side.solve(Right_Side, denom);
    vector<Number> Linear_Form(nc);
    for (i = 0; i <nc; i++) {
        Linear_Form[i] = Solution[i][0];  // the solution vector is called Linear_Form
    }
    vector<Number> test = MxV(Linear_Form); // we have solved the system by taking a square submatrix
                        // now we must test whether the solution satisfies the full system
    for (i = 0; i <nr; i++) {
        if (test[i] != denom * v[i]){
            return vector<Number>();
        }
    }
    Number total_gcd = 1; // libnormaliz::gcd(denom,v_gcd(Linear_Form)); // extract the gcd of denom and solution
    denom/=total_gcd;
    v_scalar_division(Linear_Form,total_gcd);
    return Linear_Form;
}
//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::solve_ZZ(const vector<Number>& v) const {

    Number denom;
    vector<Number> result=solve_rectangular(v,denom);
    if(denom!=1)
        result.clear();
    return result;
}
//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::find_linear_form() const {

    Number denom;
    vector<Number> result=solve_rectangular(vector<Number>(nr,1),denom);
    return result;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> Matrix<Number>::find_linear_form_low_dim () const{
    size_t rank=(*this).rank();
    if (rank == 0) { //return zero-vector as linear form
        return vector<Number>(nc,0);
    }
    if (rank == nc) { // basis change not necessary
        return (*this).find_linear_form();
    }

    Sublattice_Representation<Number> Basis_Change(*this,true);
    vector<Number> Linear_Form=Basis_Change.to_sublattice(*this).find_linear_form();
    if(Linear_Form.size()!=0)
        Linear_Form=Basis_Change.from_sublattice_dual(Linear_Form);

    return Linear_Form;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon_reduce(){

    size_t rk;
    Matrix<Number> Copy(*this);
    bool success;
    rk=row_echelon_reduce(success);

        Shrink_nr_rows(rk);
        return rk;

}
//---------------------------------------------------------------------------

template<typename Number>
Number Matrix<Number>::full_rank_index() const{
    
    Matrix<Number> Copy(*this);
    Number index;
    bool success;
    index=Copy.full_rank_index(success);

        return index;
}

//---------------------------------------------------------------------------

template<typename Number>
size_t Matrix<Number>::row_echelon(){
    
    Matrix<Number> Copy(*this);
    bool success;
    size_t rk;
    rk=row_echelon(success);

        Shrink_nr_rows(rk);
        return rk;

}

//---------------------------------------------------------------------------

template<typename Number>
Matrix<Number> Matrix<Number>::kernel () const{
// computes a ZZ-basis of the solutions of (*this)x=0
// the basis is formed by the rOWS of the returned matrix

    size_t dim=nc;
    if(nr==0)
        return(Matrix<Number>(dim));

    Matrix<Number> Copy(*this);
    size_t rank;
    bool success;
    Matrix<Number> Transf=Copy.row_column_trigonalize(rank,success);
    
    Matrix<Number> ker_basis(dim-rank,dim);
    Matrix<Number> Help =Transf.transpose();
    for (size_t i = rank; i < dim; i++) 
            ker_basis[i-rank]=Help[i];
    ker_basis.row_echelon_reduce();
    return(ker_basis);
}

//---------------------------------------------------------------------------
// Converts "this" into (column almost) Hermite normal form, returns column transformation matrix
/*template<typename Number>
Matrix<Number> Matrix<Number>::AlmostHermite(size_t& rk){

    Matrix<Number> Copy=*this;
    Matrix<Number> Transf;
    bool success;
    Transf=row_column_trigonalize(rk,success);

        return Transf;
} */

//---------------------------------------------------------------------------
// Classless conversion routines
//---------------------------------------------------------------------------

template<typename ToType, typename FromType>
void convert(Matrix<ToType>& to_mat, const Matrix<FromType>& from_mat){
    size_t nrows = from_mat.nr_of_rows();
    size_t ncols = from_mat.nr_of_columns();
    to_mat.resize(nrows, ncols);
    for(size_t i=0; i<nrows; ++i)
        for(size_t j=0; j<ncols; ++j)
            convert(to_mat[i][j], from_mat[i][j]);
}

//---------------------------------------------------------------------------


template<typename Number>
bool weight_lex(const order_helper<Number>& a, const order_helper<Number>& b){
    
        if(a.weight < b.weight)
            return true;
        if(a.weight==b.weight)
            if(*(a.v)< *(b.v))
                return true;
        return false;
}

//---------------------------------------------------------------------------

template<typename Number>
void Matrix<Number>::order_rows_by_perm(const vector<key_t>& perm){
    order_by_perm(elem,perm);    
}

template<typename Number>
Matrix<Number>& Matrix<Number>::sort_by_weights(const Matrix<Number>& Weights, vector<bool> absolute){
    if(nr<=1)
        return *this;
    vector<key_t> perm=perm_by_weights(Weights,absolute);
    order_by_perm(elem,perm);
    return *this;   
}

template<typename Number>
Matrix<Number>& Matrix<Number>::sort_lex(){
    if(nr<=1)
        return *this;
    vector<key_t> perm=perm_by_weights(Matrix<Number>(0,nc),vector<bool>(0));
    order_by_perm(elem,perm);
    return *this;    
}

template<typename Number>
vector<key_t> Matrix<Number>::perm_by_weights(const Matrix<Number>& Weights, vector<bool> absolute){
// the smallest entry is the row with index perm[0], then perm[1] etc.
    
    assert(Weights.nc==nc);
    assert(absolute.size()==Weights.nr);

    list<order_helper<Number> > order;
    order_helper<Number> entry;
    entry.weight.resize(Weights.nr);
    
    for(key_t i=0;i<nr; ++i){
        for(size_t j=0;j<Weights.nr;++j){
            if(absolute[j])
                entry.weight[j]=v_scalar_product(Weights[j],v_abs_value(elem[i]));
            else
                entry.weight[j]=v_scalar_product(Weights[j],elem[i]);                
        }
        entry.index=i;
        entry.v=&(elem[i]);
        order.push_back(entry);        
    }
    order.sort(weight_lex<Number>);
    vector<key_t> perm(nr);
    typename list<order_helper<Number> >::const_iterator ord=order.begin();
    for(key_t i=0;i<nr;++i, ++ord)
        perm[i]=ord->index; 
    
    return perm;
}

//---------------------------------------------------

/* template<typename Number>
Matrix<Number> Matrix<Number>::solve_congruences(bool& zero_modulus) const{
 
    
    zero_modulus=false;
    size_t i,j;
    size_t nr_cong=nr, dim=nc-1;
    if(nr_cong==0)
        return Matrix<Number>(dim); // give back unit matrix
    
    //add slack variables to convert congruences into equaitions
    Matrix<Number> Cong_Slack(nr_cong, dim+nr_cong);
    for (i = 0; i < nr_cong; i++) {
        for (j = 0; j < dim; j++) {
            Cong_Slack[i][j]=elem[i][j];
        }
        Cong_Slack[i][dim+i]=elem[i][dim];
        if(elem[i][dim]==0){
            zero_modulus=true;
            return Matrix<Number>(0,dim);
        }
    }
    
    //compute kernel
    
    Matrix<Number> Help=Cong_Slack.kernel(); // gives the solutions to the the system with slack variables
    Matrix<Number> Ker_Basis(dim,dim);   // must now project to first dim coordinates to get rid of them
    for(size_t i=0;i<dim;++i)
        for(size_t j=0;j<dim;++j)
            Ker_Basis[i][j]=Help[i][j];
    return Ker_Basis;
        
} */

//---------------------------------------------------

template<typename Number>
void Matrix<Number>::saturate(){
    
    // *this=kernel().kernel();
    return;    // no saturation necessary over a field
}

//---------------------------------------------------

template<typename Number>
vector<key_t> Matrix<Number>::max_and_min(const vector<Number>& L, const vector<Number>& norm) const{

    vector<key_t> result(2,0);
    if(nr==0)
        return result;
    key_t maxind=0,minind=0;
    Number maxval=v_scalar_product(L,elem[0]);
    Number maxnorm=1,minnorm=1;
    if(norm.size()>0){
        maxnorm=v_scalar_product(norm,elem[0]);
        minnorm=maxnorm;              
    }
    Number minval=maxval;
    for(key_t i=0;i<nr;++i){
        Number val=v_scalar_product(L,elem[i]);
        if(norm.size()==0){
            if(val>maxval){
                maxind=i;
                maxval=val;            
            }
            if(val<minval){
                minind=i;
                minval=val;            
            }
        }
        else{
            Number nm=v_scalar_product(norm,elem[i]);
            if(maxnorm*val>nm*maxval){
                maxind=i;
                maxval=val;            
            }
            if(minnorm*val<nm*minval){
                minind=i;
                minval=val;            
            }
        }            
    }
    result[0]=maxind;
    result[1]=minind;
    return result;
}

/*
template<typename Number>
size_t Matrix<Number>::extreme_points_first(const vector<Number> norm){
    
    if(nr==0)
        return 1;
    
    vector<long long> norm_copy;
    
    size_t nr_extr=0;
    Matrix<long long> HelpMat(nr,nc);
    try{
        convert(HelpMat,*this);
        convert(norm_copy,norm);
    }
    catch(ArithmeticException){
        return nr_extr;        
    }

    HelpMat.sort_lex();
    
    vector<bool> marked(nr,false);
    size_t no_success=0;
    // size_t nr_attempt=0;
    while(true){
        // nr_attempt++; cout << nr_attempt << endl;
        vector<long long> L=v_random<long long>(nc,10);
        vector<key_t> max_min_ind;
        max_min_ind=HelpMat.max_and_min(L,norm_copy);
            
        if(marked[max_min_ind[0]] && marked[max_min_ind[1]])
            no_success++;
        else
            no_success=0;
        if(no_success > 1000)
            break;
        marked[max_min_ind[0]]=true;
        marked[max_min_ind[1]]=true;
    }
    Matrix<long long> Extr(nr_extr,nc);  // the recognized extreme rays
    Matrix<long long> NonExtr(nr_extr,nc); // the other generators
    size_t j=0;
    vector<key_t> perm(nr);
    for(size_t i=0;i<nr;++i) {
        if(marked[i]){
            perm[j]=i;;
            j++;
        }
    }
    nr_extr=j;
    for(size_t i=0;i<nr;++i) {
        if(!marked[i]){
            perm[j]=i;;
            j++;
        }
    }
    order_rows_by_perm(perm);    
    // cout << nr_extr << "extreme points found"  << endl;
    return nr_extr;
    // exit(0);
}

template<typename Number>
vector<Number> Matrix<Number>::find_inner_point(){
    vector<key_t> simplex=max_rank_submatrix_lex();
    vector<Number> point(nc);
    for(size_t i=0;i<simplex.size();++i)
        point=v_add(point,elem[simplex[i]]);
   return point;    
}
*/

template<typename Number>
void Matrix<Number>::Shrink_nr_rows(size_t new_nr_rows){

    if(new_nr_rows>=nr)
        return;
    nr=new_nr_rows;
    elem.resize(nr);
}

template<typename Number>
Matrix<Number>  readMatrix(const string project){
// reads one matrix from file with name project
// format: nr of rows, nr of colimns, entries
// all separated by white space
    
    string name_in=project;
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false){
        cerr << "Cannot find input file" << endl;
        exit(1);
    }
    
    int nrows,ncols;
    in >> nrows;
    in >> ncols;
    
    if(nrows==0 || ncols==0){
        cerr << "Matrix empty" << endl;
        exit(1);
    }
    
    
    int i,j,entry;
    Matrix<Number> result(nrows,ncols);
    
    for(i=0;i<nrows;++i)
        for(j=0;j<ncols;++j){
            in >> entry;
            result[i][j]=entry;
        }
    return result;
}

/*
#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template Matrix<long>  readMatrix(const string project);
#endif // NMZ_MIC_OFFLOAD
template Matrix<long long>  readMatrix(const string project);
template Matrix<mpz_class>  readMatrix(const string project);
*/

}  // namespace
