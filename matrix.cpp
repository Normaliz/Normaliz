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

#include <fstream>
using namespace std;

//---------------------------------------------------------------------------

#include "matrix.h"
#include "vector_operations.h"
#include "lineare_transformation.h"

//---------------------------------------------------------------------------

extern bool test_arithmetic_overflow;
extern int overflow_test_modulus;
extern void global_error_handling();

//---------------------------------------------------------------------------
//Private
//---------------------------------------------------------------------------

void Matrix::max_rank_submatrix_lex(vector<int>& v, const int& rank) const{
	int level=v.size();
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
		Matrix S=submatrix(v);
		if (S.rank()==S.nr_of_rows()) {
			max_rank_submatrix_lex(v,rank);
			return;
		}
	}
}

//---------------------------------------------------------------------------
//Public
//---------------------------------------------------------------------------

Matrix::Matrix(){
	nr=0;
	nc=0;
}

//---------------------------------------------------------------------------

Matrix::Matrix(int dim){
	nr=dim;
	nc=dim;
	elements = vector< vector<Integer> >(dim, vector<Integer>(dim));
	for (int i = 0; i < dim; i++) {
		elements[i][i]=1;
	}
}

//---------------------------------------------------------------------------

Matrix::Matrix(int row, int col){
	nr=row;
	nc=col;
	elements = vector< vector<Integer> >(row, vector<Integer>(col));
}

//---------------------------------------------------------------------------

Matrix::Matrix(int row, int col, Integer value){
	nr=row;
	nc=col;
	elements = vector< vector<Integer> > (row, vector<Integer>(col,value));
}

//---------------------------------------------------------------------------

Matrix::Matrix(const vector< vector<Integer> >& elem){
	nr=elem.size();
	nc=elem[0].size();
	elements=elem;
}

//---------------------------------------------------------------------------

Matrix::Matrix(const Matrix& M){
	nr=M.nr;
	nc=M.nc;
	elements=M.elements;
}

//---------------------------------------------------------------------------

Matrix::~Matrix(){
	//automatic destructor
}

//---------------------------------------------------------------------------

void Matrix::write(){      //to overload for files
	int i,j;
	for(i=0; i<nr; i++){
		for(j=0; j<nc; j++) {
			cin >> elements[i][j];
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::write(int row, const vector<Integer>& data){
	if ((row<1)||(row>nr)|| (nc!=data.size())) {
		error("error:Bad argument passed to Matrix::write.");
	}
	else
		elements[row-1]=data;
}

//---------------------------------------------------------------------------

void Matrix::write(int row, const vector<int>& data){
	if ((row<1)||(row>nr)|| (nc!=data.size())) {
		error("error:Bad argument passed to Matrix::write.");
	}
	else
		for (int i = 0; i <nc ; i++) {
			elements[row-1][i]=data[i];
		}
}

//---------------------------------------------------------------------------

void Matrix::write(int row, int col, Integer data){
	if ((row<1)||(row>nr) || (col<1)||(col>nc)) {
		error("error:Bad argument passed to Matrix::write.");
	}
	else
		elements[row-1][col-1]=data;
}

//---------------------------------------------------------------------------

void Matrix::print(const string& name,const string& suffix) const{
	string file_name=name+"."+suffix;
	const char* file=file_name.c_str();
	ofstream out(file);
	int i,j;
	out<<nr<<endl<<nc<<endl;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nc; j++) {
			out<<elements[i][j]<<" ";
		}
		out<<endl;
	}
	out.close();
}

//---------------------------------------------------------------------------


void Matrix::read() const{      //to overload for files
	int i,j;
	for(i=0; i<nr; i++){
		cout << "\n" ;
		for(j=0; j<nc; j++) {
			cout << elements[i][j] << " ";
		}
	}
}

//---------------------------------------------------------------------------

vector<Integer> Matrix::read(int row) const{
	if (row<1||row>nr)  {
		error("error: Bad argument passed to Matrix::read.");
		vector<Integer> v(0);
		return v;
	}
	else
		return elements[row-1];
}

//---------------------------------------------------------------------------

Integer Matrix::read (int row, int col) const{
	if ((row>nr) || (col>nc)) {
		error("error: Bad argument passed to Matrix::read.");
		return 0;
	}
	else
		return elements[row-1][col-1];
}

//---------------------------------------------------------------------------

int Matrix::nr_of_rows () const{
	return nr;
}

//---------------------------------------------------------------------------

int Matrix::nr_of_columns () const{
	return nc;
}

//---------------------------------------------------------------------------

void Matrix::random () {
	int i,j,k;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nc; j++) {
			k = rand();
			elements[i][j]=k % 3;
		}
	}
}

//---------------------------------------------------------------------------

Matrix Matrix::submatrix(const vector<int>& rows) const{
	int size=rows.size(), j;
	Matrix M(size, nc);
	for (int i=0; i < size; i++) {
		j=rows[i]-1;
		if (nr-1<j) {
			error("error: Bad argument passed to Matrix::submatrix.");
			return M;
		}
		else {
			M.elements[i]=elements[j];
		}
	}
	return M;
}

//---------------------------------------------------------------------------

vector<Integer> Matrix::diagonale() const{
	if (nr!= nc) {
		error("error: Bad argument passed to Matrix::diagonale.");
		vector<Integer> diag(0);
		return diag;
	}
	else  {
		vector<Integer> diag(nr);
		for(int i=0; i<nr;i++){
			diag[i]=elements[i][i];
		}
		return diag;
	}
}

//---------------------------------------------------------------------------

int Matrix::maximal_decimal_length() const{
	int i,j,maxim=0;
	for (i = 0; i <nr; i++) {
		for (j = 0; j <nc; j++) {
			maxim=max(maxim,decimal_length(elements[i][j]));
		}
	}
	return maxim;
}

//---------------------------------------------------------------------------

Matrix Matrix::add(const Matrix& A) const{
	if ((nr!= A.nr)||(nc!=A.nc)) {
		error("error: Bad argument passed to Matrix::add.");
		Matrix M;
		return M;
	}
	else  {
		Matrix B(nr,nc);
		int i,j;
		for(i=0; i<nr;i++){
			for(j=0; j<nc; j++){
				B.elements[i][j]=elements[i][j]+A.elements[i][j];
			}
		}
		return B;
	}
}

//---------------------------------------------------------------------------

Matrix Matrix::multiplication(const Matrix& A) const{
	if ((nc!=A.nr)) {
		error("error: Bad argument passed to Matrix::multiplication.");
		Matrix M;
		return M;
	}
	else  {
		Matrix B(nr,A.nc,0);  //initialized with 0
		int i,j,k;
		for(i=0; i<B.nr;i++){
			for(j=0; j<B.nc; j++){
				for(k=0; k<nc; k++){
					B.elements[i][j]=B.elements[i][j]+elements[i][k]*A.elements[k][j];
				}
			}
		}
		return B;
	}
}

//---------------------------------------------------------------------------

Matrix Matrix::multiplication(const Matrix& A, int m) const{
	if ((nc!=A.nr)) {
		error("error: Bad argument passed to Matrix::multiplication.");
		Matrix M;
		return M;
	}
	else  {
		Matrix B(nr,A.nc,0);  //initialized with 0
		int i,j,k;
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
}

//---------------------------------------------------------------------------

bool Matrix::equal(const Matrix& A) const{
	if ((nr!=A.nr)||(nc!=A.nc)){  return false; }
	int i,j;
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

bool Matrix::equal(const Matrix& A, int m) const{
	if ((nr!=A.nr)||(nc!=A.nc)){  return false; }
	int i,j;
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

Matrix Matrix::transpose()const{
	Matrix B(nc,nr);
	int i,j;
	for(i=0; i<nr;i++){
		for(j=0; j<nc; j++){
			B.elements[j][i]=elements[i][j];
		}
	}
	return B;
}

//---------------------------------------------------------------------------

void Matrix::scalar_multiplication(const Integer& scalar){
	int i,j;
	for(i=0; i<nr;i++){
		for(j=0; j<nc; j++){
			elements[i][j]=elements[i][j]*scalar;
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::scalar_division(const Integer& scalar){
	int i,j;
	for(i=0; i<nr;i++){
		for(j=0; j<nc; j++){
			if (elements[i][j]%scalar!=0) {
				error("error: Bad argument passed to Matrix::scalar_division.");
				return;
			}
			else {
				elements[i][j]=elements[i][j] / scalar;
			}
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::reduction_modulo(const Integer& modulo){
	int i,j;
	for(i=0; i<nr;i++){
		for(j=0; j<nc; j++){
			elements[i][j]=elements[i][j] % modulo;
			if (elements[i][j]<0) {
				elements[i][j]=elements[i][j]+modulo;
			}
		}
	}
}

//---------------------------------------------------------------------------

vector<Integer> Matrix::make_prime(){
	vector<Integer> g(nr);
	for (int i = 0; i <nr; i++) {
		elements[i]=v_make_prime(elements[i],g[i]);
	}
	return g;
}

//---------------------------------------------------------------------------

vector<Integer> Matrix::MxV(const vector<Integer>& v) const{
	if ((nc!=v.size())) {
		error("error: Bad argument passed to Matrix::MxV.");
		vector<Integer> w;
		return w;
	}
	else  {
		vector<Integer> w(nr);
		for(int i=0; i<nr;i++){
			w[i]=v_scalar_product(elements[i],v);
		}
		return w;
	}
}

//---------------------------------------------------------------------------

vector<Integer> Matrix::VxM(const vector<Integer>& v) const{
	if ((nr!=v.size())) {
		error("error: Bad argument passed to Matrix::VxM.");
		vector<Integer> w;
		return w;
	}
	else  {
		vector<Integer> w(nc,0);
		int i,j;
		for (i=0; i<nc; i++){
			for (j=0; j<nr; j++){
				w[i]=w[i]+v[j]*elements[j][i];
			}
		}
		return w;
	}
}

//---------------------------------------------------------------------------

void Matrix::exchange_rows(const int& row1, const int& row2){
	if ((row1>nr) || (row2>nr)) {
		error("error: Bad argument passed to Matrix::exchange_rows.");
	}
	else {
		elements[row1-1].swap(elements[row2-1]);
	}
}

//---------------------------------------------------------------------------

void Matrix::exchange_columns(const int& col1, const int& col2){
	if ((col1>nc) || (col2>nc)) {
		error("error: Bad argument passed to Matrix::exchange_columns.");
	}
	else {
		register const int c1=col1-1;
		register const int c2=col2-1;
		register int i;
		Integer help;
		for(i=0; i<nr;i++){
			help=elements[i][c1];
			elements[i][c1]= elements[i][c2];
			elements[i][c2]=help;
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::reduce_row (int corner) {
	if ((corner>nr)||(corner>nc)) {
		error("error: Bad argument passed to Matrix::reduce_row.");
	}
	else{
		register int i,j;
		Integer help;
		for ( i = corner; i < nr; i++) {
			if (elements[i][corner-1]!=0) {
				help=elements[i][corner-1] / elements[corner-1][corner-1];
				for (j = corner-1; j < nc; j++) {
					elements[i][j]=elements[i][j]-help*elements[corner-1][j];
				}
			}
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::reduce_row (int corner, Matrix& Left) {
	if ((corner>nr)||(corner>nc)||(Left.nr!=nr)) {
		error("error: Bad argument passed to Matrix::reduce_row.");
	}
	else {
		int i,j;
		Integer help1, help2=elements[corner-1][corner-1];
		for ( i = corner; i < nr; i++) {
			help1=elements[i][corner-1] / help2;
			if (help1!=0) {
				for (j = corner-1; j < nc; j++) {
					elements[i][j]=elements[i][j]-help1*elements[corner-1][j];
				}
				for (j = 0; j < Left.nc; j++) {
					Left.elements[i][j]=Left.elements[i][j]-help1*Left.elements[corner-1][j];
				}
			}
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::reduce_column (int corner) {
	if ((corner>nr)||(corner>nc)) {
		error("error: Bad argument passed to Matrix::reduce_column.");
	}
	else{
		int i,j;
		Integer help1, help2=elements[corner-1][corner-1];
		for ( j = corner; j < nc; j++) {
			help1=elements[corner-1][j] / help2;
			if (help1!=0) {
				for (i = corner-1; i < nr; i++) {
					elements[i][j]=elements[i][j]-help1*elements[i][corner-1];
				}
			}
		}
	}
}

//---------------------------------------------------------------------------

void Matrix::reduce_column (int corner, Matrix& Right, Matrix& Right_Inv) {
	if ((corner>nr)||(corner>nc)||(Right.nr!=nc)||(Right.nc!=nc)||(Right_Inv.nr!=nc)||(Right_Inv.nc!=nc)) {
		error("error: Bad argument passed to Matrix::reduce_columen.");
	}
	else {
		int i,j;
		Integer help1, help2=elements[corner-1][corner-1];
		for ( j = corner; j < nc; j++) {
			help1=elements[corner-1][j] / help2;
			if (help1!=0) {
				for (i = corner-1; i < nr; i++) {
					elements[i][j]=elements[i][j]-help1*elements[i][corner-1];
				}
				for (i = 0; i < nc; i++) {
					Right.elements[i][j]=Right.elements[i][j]-help1*Right.elements[i][corner-1];
					Right_Inv.elements[corner-1][i]=Right_Inv.elements[corner-1][i]+help1*Right_Inv.elements[j][i];
				}
			}
		}
	}
}

//---------------------------------------------------------------------------

vector<int> Matrix::pivot(int corner){
	int i,j;
	Integer help=0;
	vector<int> v(2,0);
	if ((corner>nr)||(corner>nc)) {
		error("error: Bad argument passed to Matrix::pivot.");
		return v;
	}
	else{
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
	}
	return v;
}

//---------------------------------------------------------------------------

int Matrix::pivot_column(int col){
	int i,j=0;
	Integer help=0;
	if ((col>nr)||(col>nc)) {
		error("error: Bad argument passed to Matrix::pivot_column.");
		return 0;
	}
	else{
		for (i = col-1; i < nr; i++) {
			if (elements[i][col-1]!=0) {
				if ((help==0)||(Iabs(elements[i][col-1])<help)) {
					help=Iabs(elements[i][col-1]);
					j=i+1;
				}
			}
		}
	}
	return j;
}

//---------------------------------------------------------------------------

int Matrix::diagonalize(){
	int rk;
	int rk_max=min(nr,nc);
	vector<int> piv(2,0);
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

int Matrix::rank() const{
	int rk;
	int rk_max=min(nr,nc);
	vector<int> piv(2,0);
	Matrix M(*this);
	for (rk = 1; rk <= rk_max; rk++) {
		piv=M.pivot(rk);
		if (piv[0]>0) {
			do {
				M.exchange_rows (rk,piv[0]);
				M.exchange_columns (rk,piv[1]);
				M.reduce_row (rk);
				M.reduce_column (rk);
				piv=M.pivot(rk);
			} while ((piv[0]>rk)||(piv[1]>rk));
		}
		else
			break;
	}
	return rk-1;
}

//---------------------------------------------------------------------------

int Matrix::rank_destructiv(){
	register int rk,i,j,Min_Row, rk_max=min(nr,nc);
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

vector<int> Matrix::max_rank_submatrix() const{
	//may be optimized in two ways
	//first only a triangular matrix is realy neaded, no full diagonalization is necesary
	//second the matrix Rows_Exchanges may be computed by Lineare_transformation::transformation
	int rk,i,j,k;
	int rk_max=min(nr,nc);
	vector<int> piv(2,0);
	Matrix M(*this);
	Matrix Rows_Exchanges(nr);
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
	vector<int> simplex(rk);
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

vector<int>  Matrix::max_rank_submatrix_lex() const{
	int rk=rank();
	vector<int> v(0);
	max_rank_submatrix_lex(v,rk);
	return v;
}

//---------------------------------------------------------------------------

vector<int>  Matrix::max_rank_submatrix_lex(const int& rank) const {
	vector<int> v(0);
	max_rank_submatrix_lex(v,rank);
	return v;
}

//---------------------------------------------------------------------------

Matrix Matrix::solve(Matrix Right_side, Integer& det) const {
	int dim=Right_side.nr;
	int nr_sys=Right_side.nc;
	if ((nr!=nc)||(nc!=dim)) {
		error("error: Bad argument passed to Matrix::solve(Matrix, Integer).");
		Matrix Solution;
		return Solution;
	}
	else {
		Matrix Left_side(*this);
		Matrix Solution(dim,nr_sys);
		Integer S;
		int piv,rk,i,j,k;
		for (rk = 1; rk <= dim; rk++) {
			piv=Left_side.pivot_column(rk);
			if (piv>0) {
				do {
					Left_side.exchange_rows (rk,piv);
					Right_side.exchange_rows (rk,piv);
					Left_side.reduce_row(rk, Right_side);
					piv=Left_side.pivot_column(rk);
				} while (piv>rk);
			}
		}
		det=Left_side.elements[0][0];
		for (i = 1; i < dim; i++) {
			det*=Left_side.elements[i][i];
		}
		if (det==0) {
			error("warning: Determinant=0 in Matrix::solve.");
			return Solution;
		}
		else {
			Integer d=Iabs(det);
			for (i = 0; i < nr_sys; i++) {
				for (j = dim-1; j >= 0; j--) {
					S=Iabs(d)*Right_side.elements[j][i];
					for (k = j+1; k < dim; k++) {
						S-=Left_side.elements[j][k]*Solution.elements[k][i];
					}
					Solution.elements[j][i]=S/Left_side.elements[j][j];
				}
			}
			return Solution;
		}
	}
}

//---------------------------------------------------------------------------

Matrix Matrix::solve(Matrix Right_side, vector< Integer >& diagonal, Integer& det) const {
	int dim=Right_side.nr;
	int nr_sys=Right_side.nc;
	if ((nr!=nc)||(nc!=dim)||(dim!=diagonal.size())) {
		error("error: Bad argument passed to Matrix::solve(Matrix, vector, Integer).");
		Matrix Solution;
		return Solution;
	}
	else {
		Matrix Left_side(*this);
		Matrix Solution(dim,nr_sys);
		Integer S;
		int piv,rk,i,j,k;
		for (rk = 1; rk <= dim; rk++) {
			piv=Left_side.pivot_column(rk);
			if (piv>0) {
				do {
					Left_side.exchange_rows (rk,piv);
					Right_side.exchange_rows (rk,piv);
					Left_side.reduce_row(rk, Right_side);
					piv=Left_side.pivot_column(rk);
				} while (piv>rk);
			}
		}
		det=Left_side.elements[0][0];
		diagonal[0]= Left_side.elements[0][0];
		for (i = 1; i < dim; i++) {
			det*=Left_side.elements[i][i];
			diagonal[i]= Left_side.elements[i][i];
		}
		if (det==0) {
			error("warning: Determinant=0 in Matrix::solve.");
			return Solution;
		}
		else {
			Integer d=Iabs(det);
			for (i = 0; i < nr_sys; i++) {
				for (j = dim-1; j >= 0; j--) {
					S=Iabs(d)*Right_side.elements[j][i];
					for (k = j+1; k < dim; k++) {
						S-=Left_side.elements[j][k]*Solution.elements[k][i];
					}
					Solution.elements[j][i]=S/Left_side.elements[j][j];
				}
			}
			return Solution;
		}
	}
}

//---------------------------------------------------------------------------

Matrix Matrix::invert(vector< Integer >& diagonal, Integer& det) const{
	if ((nr!=nc)||(nr!=diagonal.size())) {
		error("error: Bad argument passed to Matrix::invert.");
		Matrix Solution;
		return Solution;
	}
	else{
		Matrix Left_side(*this);
		Matrix Right_side(nr);
		Matrix Solution(nr,nr);
		Integer S;
		int piv,rk,i,j,k;
		for (rk = 1; rk <= nr; rk++) {
			piv=Left_side.pivot_column(rk);
			if (piv>0) {
				do {
					Left_side.exchange_rows (rk,piv);
					Right_side.exchange_rows (rk,piv);
					Left_side.reduce_row(rk, Right_side);
					piv=Left_side.pivot_column(rk);
				} while (piv>rk);
			}
		}
		det=Left_side.elements[0][0];
		diagonal[0]= Left_side.elements[0][0];
		for (i = 1; i < nr; i++) {
			det*=Left_side.elements[i][i];
			diagonal[i]= Left_side.elements[i][i];
		}
		if (det==0) {
			error("error: Determinant=0 in Matrix::invert. Non-invertible Matrix.");
			return Solution;
		}
		else {
			Integer d=Iabs(det);
			for (i = 0; i < nr; i++) {
				for (j = nr-1; j >= 0; j--) {
					S=Iabs(d)*Right_side.elements[j][i];
					for (k = j+1; k < nr; k++) {
						S-=Left_side.elements[j][k]*Solution.elements[k][i];
					}
					Solution.elements[j][i]=S/Left_side.elements[j][j];
				}
			}
			return Solution;
		}
	}
}

//---------------------------------------------------------------------------

vector<Integer> Matrix::homogeneous (bool& homogeneous) const{
	if (nc == 0 || nr == 0) { //return zero-vector as linear form
		homogeneous=true;
		return vector<Integer>(nc,0);
	}
	int i;
	Integer det,buffer;
	vector<int>  rows=max_rank_submatrix_lex();
	Matrix Left_Side=submatrix(rows);
	Matrix Right_Side(nc,1,1);
	Matrix Solution=Solve(Left_Side, Right_Side, det);
	det=Iabs(det);
	vector<Integer> Linear_Form(nc);
	for (i = 0; i <nc; i++) {
		buffer=Solution.read(i+1,1);
		Linear_Form[i]=buffer/det;
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

vector<Integer> Matrix::homogeneous_low_dim (bool& homogeneous) const{
	int rank=(*this).rank();
	if (rank == 0) { //return zero-vector as linear form
		homogeneous=true;
		return vector<Integer>(nc,0);
	}
	if (rank == nc) { // basis change not necessary
		return (*this).homogeneous(homogeneous);
	}

	// prepare basis change
	vector <int> key = max_rank_submatrix_lex(rank);
	Matrix full_rank_matrix = submatrix(key);  // has maximal number of linear independent lines
	Lineare_Transformation Basis_Change = Transformation(full_rank_matrix);
	rank=Basis_Change.get_rank();
	Matrix V=Basis_Change.get_right();
	Matrix Change_To_Full_Emb(nc,rank);
	int i,j;
	for (i = 1; i <=nc; i++) {
		for (j = 1; j <= rank; j++) {
			Change_To_Full_Emb.write(i,j,V.read(i,j));
		}
	}
	
	//apply basis change
	Matrix Full_Cone_Generators = full_rank_matrix.multiplication(Change_To_Full_Emb);
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

bool Matrix::test_solve(const Matrix& Solution, const Matrix& Right_side,
		const Integer& det,const int& m) const{
	Matrix LS=multiplication(Solution,m);
	Matrix RS=Right_side;
	RS.scalar_multiplication(Iabs(det));
	if (LS.equal(RS,m)!=true) {
		error("error: Matrix::test_solve failed.\nPossible arithmetic overflow in Matrix::solve.\n");
		return false;
	}
	return true;

}

//---------------------------------------------------------------------------

bool Matrix::test_invert(const Matrix& Solution,
		const Integer& det,const int& m) const{
	Matrix LS=multiplication(Solution,m);
	Matrix RS(nr);
	RS.scalar_multiplication(Iabs(det));
	if (LS.equal(RS,m)!=true) {
		error("error: Matrix::test_invert failed.\nPossible arithmetic overflow in Matrix::invert.\n");
		return false;
	}
	return true;

}

//---------------------------------------------------------------------------

void Matrix::error(string s) const{
	cerr <<"\nMatrix "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

Matrix Solve(const Matrix& Left_side, const Matrix& Right_side,Integer& det){
	Matrix S=Left_side.solve(Right_side,det);
	if (test_arithmetic_overflow==true) {
		bool testing=Left_side.test_solve(S,Right_side,det,overflow_test_modulus);
		if (testing==false) {
			cout<<"\nSolving the linear system of equations has failed.\n";
			global_error_handling();
		}
	}
	return S;
}

//---------------------------------------------------------------------------

Matrix Invert(const Matrix& Left_side,  vector< Integer >& diagonal ,Integer& det){
	Matrix S=Left_side.invert(diagonal,det);
	if (test_arithmetic_overflow==true) {
		bool testing=Left_side.test_invert(S,det,overflow_test_modulus);
		if (testing==false) {
			cout<<"\nInverting the matrix has failed.\n";
			global_error_handling();
		}
	}
	return S;
}

//---------------------------------------------------------------------------

