/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#include <fstream>
#include <algorithm>

#include "matrix.h"
#include "vector_operations.h"
#include "lineare_transformation.h"
#include "normaliz_exception.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//Private
//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::max_rank_submatrix_lex(vector<size_t>& v, const size_t& rank) const{
	size_t level=v.size();
	if (level==rank) {
		return;
	}
	if (level==0) {
		v.push_back(1);
	}
	else{
		v.push_back(v[level-1]+1);
	}
	for (; v[level] <= nr; v[level]++) {
		Matrix<Integer> S=submatrix(v);
		if (S.rank()==S.nr_of_rows()) {
			max_rank_submatrix_lex(v,rank);
			return;
		}
	}
}

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
	assert(dim>=0);
	nr=dim;
	nc=dim;
	elements = vector< vector<Integer> >(dim, vector<Integer>(dim));
	for (size_t i = 0; i < dim; i++) {
		elements[i][i]=1;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(size_t row, size_t col){
	assert(row>=0);
	assert(col>=0);
	nr=row;
	nc=col;
	elements = vector< vector<Integer> >(row, vector<Integer>(col));
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(size_t row, size_t col, Integer value){
	assert(row>=0);
	assert(col>=0);
	nr=row;
	nc=col;
	elements = vector< vector<Integer> > (row, vector<Integer>(col,value));
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(const vector< vector<Integer> >& elem){
	nr=elem.size();
	nc=elem[0].size();
	elements=elem;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(const list< vector<Integer> >& elem){
	nr = elem.size();
	elements = vector< vector<Integer> > (nr);
	nc = 0;
	size_t i=0;
	for(typename list< vector<Integer> >::const_iterator it=elem.begin(); it!=elem.end(); ++it, ++i) {
		if(i == 0) {
			nc = (*it).size();
		} else {
			assert((*it).size() == nc);
		}
		elements[i]=(*it);
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::Matrix(const Matrix<Integer>& M){
	nr=M.nr;
	nc=M.nc;
	elements=M.elements;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer>::~Matrix(){
	//automatic destructor
}


//---------------------------------------------------------------------------

template<typename Integer>
list< vector<Integer> > Matrix<Integer>::to_list(){
	//list< vector<Integer> > elemlist(nr, vector<Integer>());
	list< vector<Integer> > elemlist = list< vector<Integer> >();
	//typename list< vector<Integer> >::iterator it=elemlist.begin();
	for (size_t i=0; i<nr; ++i) {
		//(*(it++)).swap(elements[i]);
		elemlist.push_back(elements[i]);
	}
	return elemlist;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(){
	size_t i,j;
	for(i=0; i<nr; i++){
		for(j=0; j<nc; j++) {
			cin >> elements[i][j];
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(size_t row, const vector<Integer>& data){
	assert(row >= 1);
	assert(row <= nr); 
	assert(nc == data.size());
	
	elements[row-1]=data;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(size_t row, const vector<int>& data){
	assert(row >= 1);
	assert(row <= nr); 
	assert(nc == data.size());

	for (size_t i = 0; i <nc ; i++) {
		elements[row-1][i]=data[i];
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::write(size_t row, size_t col, Integer data){
	assert(row >= 1);
	assert(row <= nr); 
	assert(col >= 1);
	assert(col <= nc); 

	elements[row-1][col-1]=data;
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
			out<<elements[i][j]<<" ";
		}
		out<<endl;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::pretty_print(ostream& out) const{
	size_t i,j,k;
	size_t max_length = maximal_decimal_length();
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nc; j++) {
			for (k= 0; k <= max_length - decimal_length(elements[i][j]); k++) {
				out<<" ";
			}
			out<<elements[i][j];
		}
		out<<endl;
	}
	out<<endl;
}
//---------------------------------------------------------------------------


template<typename Integer>
void Matrix<Integer>::read() const{      //to overload for files
	size_t i,j;
	for(i=0; i<nr; i++){
		cout << "\n" ;
		for(j=0; j<nc; j++) {
			cout << elements[i][j] << " ";
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::read(size_t row) const{
	assert(row >= 1);
	assert(row <= nr); 

	return elements[row-1];
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::read (size_t row, size_t col) const{
	assert(row >= 1);
	assert(row <= nr); 
	assert(col >= 1);
	assert(col <= nc); 

	return elements[row-1][col-1];
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
	size_t i,j,k;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nc; j++) {
			k = rand();
			elements[i][j]=k % mod;
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::submatrix(const vector<size_t>& rows) const{
	size_t size=rows.size(), j;
	Matrix<Integer> M(size, nc);
	for (size_t i=0; i < size; i++) {
		j=rows[i]-1;
		assert(j >= 0);
		assert(j < nr);
		M.elements[i]=elements[j];
	}
	return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::submatrix(const vector<int>& rows) const{
	size_t size=rows.size(), j;
	Matrix<Integer> M(size, nc);
	for (size_t i=0; i < size; i++) {
		j=rows[i]-1;
		assert(j >= 0);
		assert(j < nr);
		M.elements[i]=elements[j];
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
			M.elements[j++] = elements[i];
		}
	}
	return M;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::diagonale() const{
	assert(nr == nc); 
	vector<Integer> diag(nr);
	for(size_t i=0; i<nr;i++){
		diag[i]=elements[i][i];
	}
	return diag;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::maximal_decimal_length() const{
	size_t i,j,maxim=0;
	for (i = 0; i <nr; i++) {
		for (j = 0; j <nc; j++) {
			maxim=max(maxim,decimal_length(elements[i][j]));
		}
	}
	return maxim;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::append(const Matrix<Integer>& M) {
	assert (nc == M.nc);
	elements.reserve(nr+M.nr);
	for (size_t i=0; i<M.nr; i++) {
		elements.push_back(M.elements[i]);
	}
	nr += M.nr;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::append(const vector<Integer>& V) {
	assert (nc == V.size());
	elements.push_back(V);
	nr++;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::cut_columns(size_t c) {
	assert (c >= 0);
	assert (c <= nc);
	for (size_t i=0; i<nr; i++) {
		elements[i].resize(c);
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
			B.elements[i][j]=elements[i][j]+A.elements[i][j];
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
				B.elements[i][j]=B.elements[i][j]+elements[i][k]*A.elements[k][j];
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
				B.elements[i][j]=(B.elements[i][j]+elements[i][k]*A.elements[k][j])%m;
				if (B.elements[i][j]<0) {
					B.elements[i][j]=B.elements[i][j]+m;
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
			if (elements[i][j]!=A.elements[i][j]) {
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
			if (((elements[i][j]-A.elements[i][j]) % m)!=0) {
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
			B.elements[j][i]=elements[i][j];
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
			elements[i][j] *= scalar;
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
			assert (elements[i][j]%scalar == 0);
			elements[i][j] /= scalar;
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::reduction_modulo(const Integer& modulo){
	size_t i,j;
	for(i=0; i<nr;i++){
		for(j=0; j<nc; j++){
			elements[i][j] %= modulo;
			if (elements[i][j] < 0) {
				elements[i][j] += modulo;
			}
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Matrix<Integer>::matrix_gcd() const{
	Integer g=0,h;
	for (size_t i = 0; i <nr; i++) {
		h = v_gcd(elements[i]);
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
		elements[i]=v_make_prime(elements[i],g[i]);
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
        M.elements[i][j] = elements[i][j]*m[i];
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
		w[i]=v_scalar_product(elements[i],v);
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
			w[i] += v[j]*elements[j][i];
 		}
	}
	return w;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::exchange_rows(const size_t& row1, const size_t& row2){
	if (row1 == row2) return;
	assert(row1 > 0);
	assert(row1 <= nr);
	assert(row2 > 0);
	assert(row2 <= nr);
	elements[row1-1].swap(elements[row2-1]);
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::exchange_columns(const size_t& col1, const size_t& col2){
	if (col1 == col2) return;
	assert(col1 > 0);
	assert(col1 <= nc);
	assert(col2 > 0);
	assert(col2 <= nc);
	register const size_t c1=col1-1;
	register const size_t c2=col2-1;
	register size_t i;
	Integer help;
	for(i=0; i<nr;i++){
		help=elements[i][c1];
		elements[i][c1]= elements[i][c2];
		elements[i][c2]=help;
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::reduce_row (size_t corner) {
	assert(corner > 0);
	assert(corner <= nc);
	assert(corner <= nr);
	register size_t i,j;
	Integer help;
	for ( i = corner; i < nr; i++) {
		if (elements[i][corner-1]!=0) {
			help=elements[i][corner-1] / elements[corner-1][corner-1];
			for (j = corner-1; j < nc; j++) {
				elements[i][j] -= help*elements[corner-1][j];
			}
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::reduce_row (size_t corner, Matrix<Integer>& Left) {
	assert(corner > 0);
	assert(corner <= nc);
	assert(corner <= nr);
	assert(Left.nr == nr);
	size_t i,j;
	Integer help1, help2=elements[corner-1][corner-1];
	for ( i = corner; i < nr; i++) {
		help1=elements[i][corner-1] / help2;
		if (help1!=0) {
			for (j = corner-1; j < nc; j++) {
				elements[i][j] -= help1*elements[corner-1][j];
			}
			for (j = 0; j < Left.nc; j++) {
				Left.elements[i][j] -= help1*Left.elements[corner-1][j];
			}
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::reduce_column (size_t corner) {
	assert(corner > 0);
	assert(corner <= nc);
	assert(corner <= nr);
	size_t i,j;
	Integer help1, help2=elements[corner-1][corner-1];
	for ( j = corner; j < nc; j++) {
		help1=elements[corner-1][j] / help2;
		if (help1!=0) {
			for (i = corner-1; i < nr; i++) {
				elements[i][j] -= help1*elements[i][corner-1];
			}
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
void Matrix<Integer>::reduce_column (size_t corner, Matrix<Integer>& Right, Matrix<Integer>& Right_Inv) {
	assert(corner > 0);
	assert(corner <= nc);
	assert(corner <= nr);
	assert(Right.nr == nc);
	assert(Right.nc == nc);
	assert(Right_Inv.nr == nc);
	assert(Right_Inv.nc ==nc);
	size_t i,j;
	Integer help1, help2=elements[corner-1][corner-1];
	for ( j = corner; j < nc; j++) {
		help1=elements[corner-1][j] / help2;
		if (help1!=0) {
			for (i = corner-1; i < nr; i++) {
				elements[i][j] -= help1*elements[i][corner-1];
			}
			for (i = 0; i < nc; i++) {
				Right.elements[i][j] -= help1*Right.elements[i][corner-1];
				Right_Inv.elements[corner-1][i] += help1*Right_Inv.elements[j][i];
			}
		}
	}
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t> Matrix<Integer>::pivot(size_t corner){
	assert(corner > 0);
	assert(corner <= nc);
	assert(corner <= nr);
	size_t i,j;
	Integer help=0;
	vector<size_t> v(2,0);

	for (i = corner-1; i < nr; i++) {
		for (j = corner-1; j < nc; j++) {
			if (elements[i][j]!=0) {
				if ((help==0)||(Iabs(elements[i][j])<help)) {
					help=Iabs(elements[i][j]);
					v[0]=i+1;
					v[1]=j+1;
					if (help == 1) return v;
				}
			}
		}
	}

	return v;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::pivot_column(size_t col){
	assert(col > 0);
	assert(col <= nc);
	assert(col <= nr);
	size_t i,j=0;
	Integer help=0;

	for (i = col-1; i < nr; i++) {
		if (elements[i][col-1]!=0) {
			if ((help==0)||(Iabs(elements[i][col-1])<help)) {
				help=Iabs(elements[i][col-1]);
				j=i+1;
			}
			if (help == 1) break;
		}
	}

	return j;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::diagonalize(){
	size_t rk;
	size_t rk_max=min(nr,nc);
	vector<size_t> piv(2,0);
	for (rk = 1; rk <= rk_max; rk++) {
		piv=pivot(rk);
		if (piv[0]>0) {
			do {
				exchange_rows (rk,piv[0]);
				exchange_columns (rk,piv[1]);
				reduce_row (rk);
				reduce_column (rk);
				piv=pivot(rk);
			} while ((piv[0]>rk)||(piv[1]>rk));
		}
		else
			break;
	}
	return rk-1;
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::rank() const{
    Matrix<Integer> N=*this;
    return N.rank_destructiv();
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Matrix<Integer>::rank_destructiv(){
	register size_t rk,i,j,Min_Row, rk_max=min(nr,nc);
	register bool empty;
	Integer Test, Min;
	for (rk = 1; rk <= rk_max; rk++) {
		for (i = rk; i <= nr; i++) {
			for (j = rk; j <= nc; j++)
				if (elements[i-1][j-1]!=0)
					break;
			if (j<=nc)
				break;
		}
		if (i>nr)
			break;
		if (rk!=i)
			exchange_rows (rk,i);
		if (rk!=j)
			exchange_columns (rk,j);
		do {
			Min=Iabs(elements[rk-1][rk-1]);
			Min_Row=rk;
			empty=true;
			for (i = rk+1; i <= nr; i++) {
				Test=Iabs(elements[i-1][rk-1]);
				empty=empty && (Test==0);
				if (Test!=0&& (Test<Min)) {
					Min=Test;
					Min_Row=i;
				}
			}
			if (Min_Row!=rk) {
				exchange_rows (rk,Min_Row);    
			}
			reduce_row (rk);
		} while (!empty);
	}
	return rk-1;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t> Matrix<Integer>::max_rank_submatrix() const{
	//may be optimized in two ways
	//first only a triangular matrix is realy needed, no full diagonalization is necesary
	//second the matrix Rows_Exchanges may be computed by Lineare_transformation::transformation
	size_t rk,i,j,k;
	size_t rk_max=min(nr,nc);
	vector<size_t> piv(2,0);
	Matrix<Integer> M(*this);
	Matrix<Integer> Rows_Exchanges(nr);
	for (rk = 1; rk <= rk_max; rk++) {
		piv=M.pivot(rk);
		if (piv[0]>0) {
			do {
				M.exchange_rows (rk,piv[0]);
				Rows_Exchanges.exchange_columns(rk,piv[0]);
				M.exchange_columns (rk,piv[1]);
				M.reduce_row (rk);
				M.reduce_column (rk);  //optimization posible here
				piv=M.pivot(rk);
			} while ((piv[0]>rk)||(piv[1]>rk));
		}
		else
			break;
	}
	rk=rk-1;
	M=Rows_Exchanges.multiplication(M);
	vector<size_t> simplex(rk);
	k=0;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nc; j++) {
			if (M.elements[i][j]!=0) {
				simplex[k]=i+1;
				k++;
			}
		}
	}
	return simplex;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t>  Matrix<Integer>::max_rank_submatrix_lex() const{
	size_t rk=rank();
	vector<size_t> v(0);
	max_rank_submatrix_lex(v,rk);
	return v;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t>  Matrix<Integer>::max_rank_submatrix_lex(const size_t& rank) const {
	vector<size_t> v(0);
	max_rank_submatrix_lex(v,rank);
	return v;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve(Matrix<Integer> Right_side, Integer& denom) const {
	Matrix<Integer> Left_side(*this);
	vector<Integer> dummy_diag(nr);
	return Left_side.solve_destructiv(Right_side, dummy_diag, denom);
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve(Matrix<Integer> Right_side, vector< Integer >& diagonal, Integer& denom) const {
	Matrix<Integer> Left_side(*this);
	return Left_side.solve_destructiv(Right_side, diagonal, denom);
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::solve_destructiv(Matrix<Integer>& Right_side, vector< Integer >& diagonal, Integer& denom) {
	size_t dim=Right_side.nr;
	size_t nr_sys=Right_side.nc;
	assert(nr == nc);
	assert(nc == dim);
	assert(dim == diagonal.size());

	Matrix<Integer> Solution(dim,nr_sys);
	Integer S;
	size_t piv,rk,i;

	for (rk = 1; rk <= dim; rk++) {
		piv=(*this).pivot_column(rk);
		if (piv>0) {
			do {
				(*this).exchange_rows (rk,piv);
				Right_side.exchange_rows (rk,piv);
				(*this).reduce_row(rk, Right_side);
				piv=(*this).pivot_column(rk);
			} while (piv>rk);
		}
	}
	denom=(*this).elements[0][0];
	diagonal[0]= (*this).elements[0][0];
	for (i = 1; i < dim; i++) {
		denom*=(*this).elements[i][i];
		diagonal[i]= (*this).elements[i][i];
	}

	if (denom==0) { 
		throw NormalizException(); //TODO welche Exception?
	}

	denom=Iabs(denom);
	int j;
	size_t k;
	for (i = 0; i < nr_sys; i++) {
		for (j = dim-1; j >= 0; j--) {
			S=denom*Right_side.elements[j][i];
			for (k = j+1; k < dim; k++) {
				S-=(*this).elements[j][k]*Solution.elements[k][i];
			}
			Solution.elements[j][i]=S/(*this).elements[j][j];
		}
	}
	return Solution;

}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Matrix<Integer>::invert(vector< Integer >& diagonal, Integer& denom) const{
	assert(nr == nc);
	assert(nr == diagonal.size());
	Matrix<Integer> Left_side(*this);
	Matrix<Integer> Right_side(nr);

	return Left_side.solve_destructiv(Right_side,diagonal,denom);
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::homogeneous (bool& homogeneous) const{
	if (nc == 0 || nr == 0) { //return zero-vector as linear form
		homogeneous=true;
		return vector<Integer>(nc,0);
	}
	size_t i;
	Integer denom,buffer;
	vector<size_t>  rows=max_rank_submatrix_lex();
	Matrix<Integer> Left_Side=submatrix(rows);
	assert(nc == Left_Side.nr); //otherwise input hadn't full rank //TODO 
	Matrix<Integer> Right_Side(nc,1,1);
	Matrix<Integer> Solution=Solve(Left_Side, Right_Side, denom);
	vector<Integer> Linear_Form(nc);
	for (i = 0; i <nc; i++) {
		buffer=Solution.read(i+1,1);
		Linear_Form[i]=buffer/denom;
	}
	vector<Integer> test_homogeneous=MxV(Linear_Form);
	for (i = 0; i <nr; i++) {
		if (test_homogeneous[i]!=1) {
			homogeneous=false;
			vector<Integer> F;
			return F;
		}
	}
	homogeneous=true;
	return Linear_Form;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Matrix<Integer>::homogeneous_low_dim (bool& homogeneous) const{
	size_t rank=(*this).rank();
	if (rank == 0) { //return zero-vector as linear form
		homogeneous=true;
		return vector<Integer>(nc,0);
	}
	if (rank == nc) { // basis change not necessary
		return (*this).homogeneous(homogeneous);
	}

	// prepare basis change
	vector <size_t> key = max_rank_submatrix_lex(rank);
	Matrix<Integer> Full_Rank_Matrix = submatrix(key);  // has maximal number of linear independent lines
	Lineare_Transformation<Integer> Basis_Change = Transformation(Full_Rank_Matrix);
	rank=Basis_Change.get_rank();
	Matrix<Integer> V=Basis_Change.get_right();
	Matrix<Integer> Change_To_Full_Emb(nc,rank);
	size_t i,j;
	for (i = 1; i <=nc; i++) {
		for (j = 1; j <= rank; j++) {
			Change_To_Full_Emb.write(i,j,V.read(i,j));
		}
	}
	
	//apply basis change
	Matrix<Integer> Full_Cone_Generators = Full_Rank_Matrix.multiplication(Change_To_Full_Emb);
	//compute linear form
	vector<Integer> Linear_Form = Full_Cone_Generators.homogeneous(homogeneous);
	if (homogeneous) {
		//lift linear form back
		Change_To_Full_Emb = Change_To_Full_Emb.transpose();  // preparing the matrix for transformation on the dual space
		vector<Integer> v;
		Integer index=Basis_Change.get_index();
		if (index != 1) { 
			homogeneous = false;
			return Linear_Form;
		}
		Linear_Form=Change_To_Full_Emb.VxM(Linear_Form);
		Linear_Form=v_make_prime(Linear_Form);

		//check if all rows are in height 1
		for (i=0; i<nr; i++) {
			if (v_scalar_product(read(i+1),Linear_Form) != 1) {
				homogeneous=false;
				return Linear_Form;
			}
		}
	}
	return Linear_Form;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::test_solve(const Matrix<Integer>& Solution, const Matrix<Integer>& Right_side,
		const Integer& denom,const long& m) const{
	Matrix<Integer> LS=multiplication(Solution,m);
	Matrix<Integer> RS=Right_side;
	RS.scalar_multiplication(denom);
	if (LS.equal(RS,m)!=true) {
		throw ArithmeticException();
		return false;
	}
	return true;

}

//---------------------------------------------------------------------------

template<typename Integer>
bool Matrix<Integer>::test_invert(const Matrix<Integer>& Solution, const Integer& denom,const long& m) const{
	Matrix<Integer> LS=multiplication(Solution,m);
	Matrix<Integer> RS(nr);
	RS.scalar_multiplication(denom);
	if (LS.equal(RS,m)!=true) {
		throw ArithmeticException();
		return false;
	}
	return true;

}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Solve(const Matrix<Integer>& Left_side, const Matrix<Integer>& Right_side,Integer& denom){
	Matrix<Integer> S=Left_side.solve(Right_side,denom);
	if (test_arithmetic_overflow) {
		bool testing=Left_side.test_solve(S,Right_side,denom,overflow_test_modulus);
		if (testing==false) {
			throw ArithmeticException();
		}
	}
	return S;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Invert(const Matrix<Integer>& Left_side, vector< Integer >& diagonal, Integer& denom){
	Matrix<Integer> S=Left_side.invert(diagonal,denom);
	if (test_arithmetic_overflow) {
		bool testing=Left_side.test_invert(S,denom,overflow_test_modulus);
		if (testing==false) {
			throw ArithmeticException();
		}
	}
	return S;
}

//---------------------------------------------------------------------------

}
