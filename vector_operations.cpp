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
#include <set>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

//---------------------------------------------------------------------------

#include "integer.h"
#include "vector_operations.h"

//---------------------------------------------------------------------------

extern void global_error_handling();

//---------------------------------------------------------------------------

void v_write(vector<Integer>& v){
	int i,s=v.size();
	for (i=0; i <s; i++) {
		cin>>v[i];
	}
}

//---------------------------------------------------------------------------

int v_read(const vector<Integer>& v){
	int i,s=v.size();
	cout<<"\n";
	for (i=0; i <s; i++) {
		cout<<v[i]<<" ";
	}
	return s;
}

//---------------------------------------------------------------------------

int v_read(const vector<int>& v){
	int i,s=v.size();
	cout<<"\n";
	for (i=0; i <s; i++) {
		cout<<v[i]<<" ";
	}
	return s;
}

//---------------------------------------------------------------------------

int v_read(const vector<bool>& v){
	int i,s=v.size();
	cout<<"\n";
	for (i=0; i <s; i++) {
		cout<<v[i]<<" ";
	}
	return s;
}

//---------------------------------------------------------------------------

Integer v_scalar_product(const vector<Integer>& av,const vector<Integer>& bv){

	//loop stretching ; brings some small speed improvement
	Integer ans = 0;
	register int i,n=av.size();
	vector<Integer>::const_iterator    a=av.begin(),b=bv.begin();

	if( n >= 64 )
	{
		for( i = 0; i < ( n >> 6 ); ++i, a += 64, b += 64 ){
			ans += a[0] * b[0];
			ans += a[1] * b[1];
			ans += a[2] * b[2];
			ans += a[3] * b[3];
			ans += a[4] * b[4];
			ans += a[5] * b[5];
			ans += a[6] * b[6];
			ans += a[7] * b[7];
			ans += a[8] * b[8];
			ans += a[9] * b[9];
			ans += a[10] * b[10];
			ans += a[11] * b[11];
			ans += a[12] * b[12];
			ans += a[13] * b[13];
			ans += a[14] * b[14];
			ans += a[15] * b[15];
			ans += a[16] * b[16];
			ans += a[17] * b[17];
			ans += a[18] * b[18];
			ans += a[19] * b[19];
			ans += a[20] * b[20];
			ans += a[21] * b[21];
			ans += a[22] * b[22];
			ans += a[23] * b[23];
			ans += a[24] * b[24];
			ans += a[25] * b[25];
			ans += a[26] * b[26];
			ans += a[27] * b[27];
			ans += a[28] * b[28];
			ans += a[29] * b[29];
			ans += a[30] * b[30];
			ans += a[31] * b[31];
			ans += a[32] * b[32];
			ans += a[33] * b[33];
			ans += a[34] * b[34];
			ans += a[35] * b[35];
			ans += a[36] * b[36];
			ans += a[37] * b[37];
			ans += a[38] * b[38];
			ans += a[39] * b[39];
			ans += a[40] * b[40];
			ans += a[41] * b[41];
			ans += a[42] * b[42];
			ans += a[43] * b[43];
			ans += a[44] * b[44];
			ans += a[45] * b[45];
			ans += a[46] * b[46];
			ans += a[47] * b[47];
			ans += a[48] * b[48];
			ans += a[49] * b[49];
			ans += a[50] * b[50];
			ans += a[51] * b[51];
			ans += a[52] * b[52];
			ans += a[53] * b[53];
			ans += a[54] * b[54];
			ans += a[55] * b[55];
			ans += a[56] * b[56];
			ans += a[57] * b[57];
			ans += a[58] * b[58];
			ans += a[59] * b[59];
			ans += a[60] * b[60];
			ans += a[61] * b[61];
			ans += a[62] * b[62];
			ans += a[63] * b[63];
		}
		n -= i << 6;
	}
	for( i = 0; i < n; ++i )
		ans += a[i] * b[i];
	return( ans );
}

//---------------------------------------------------------------------------

Integer v_scalar_product_unequal_vectors_end(const vector<Integer>& a,const vector<Integer>& b){
	Integer ans = 0;
	register int i,n=a.size(),m=b.size();
	for (i = 1; i <= n; i++) {
		ans+=a[n-i]*b[m-i];
	}
	return ans;
}


//---------------------------------------------------------------------------

vector<Integer> v_merge(const vector<Integer>& a,const vector<Integer>& b){
	int s1=a.size(), s2=b.size(), i;
	vector<Integer> c(s1+s2);
	for (i = 0; i < s1; i++) {
		c[i]=a[i];
	}
	for (i = 0; i < s2; i++) {
		c[s1+i]=b[i];
	}
	return c;
}

//---------------------------------------------------------------------------

vector<Integer> v_add(const vector<Integer>& a,const vector<Integer>& b){
	register int i,s=a.size();
	vector<Integer> d(s);
	for (i = 0; i <s; i++) {
		d[i]=a[i]+b[i];
	}
	return d;
}

//---------------------------------------------------------------------------

vector<int> v_complement(const int& a, const vector<int>& v){
	int i,s=v.size();
	vector<int> u(s-1);
	for (i = 0; i <a; i++) {
		u[i]=v[i];
	}
	for (; i <s-1; i++) {
		u[i]=v[i+1];
	}
	return u;
}

//---------------------------------------------------------------------------

vector<Integer> v_abs(const vector<Integer>& v){
	int i,j=v.size();
	vector<Integer> w(j,0);
	if (j==0) {
		v_error("error: Bad argument passed to v_abs.\n");
		return w;
	}
	else{
		for (i = 0; i < j; i++) {
			w[i]=Iabs(v[i]);
		}
	}
	return w;
}

//---------------------------------------------------------------------------

Integer v_gcd(const vector<Integer>& v){
	int i,j=v.size();
	if (j==0) {
		v_error("error: Bad argument passed to v_gcd.");
		return 0;
	}
	else{
		if (v[0]==1) {
			return 1;
		}
		Integer g=v[0];
		for (i = 1; i < j; i++) {
			g=gcd(g,v[i]);
			if (g==1) {
				return 1;
			}
		}
		return g;
	}
}

//---------------------------------------------------------------------------

Integer v_lcm(const vector<Integer>& v){
	int i,j=v.size();
	if (j==0) {
		v_error("error: Bad argument passed to v_lcm.");
		return 0;
	}
	else{
		Integer g=v[0];
		for (i = 1; i < j; i++) {
			g=lcm(g,v[i]);
			if (g==0) {
				return 0;
			}
		}
		return g;
	}
}

//---------------------------------------------------------------------------

vector<Integer> v_make_prime(const vector<Integer>& v){
	int i,j=v.size();
	vector<Integer> w(j,0);
	if (j==0) {
		v_error("error: Bad argument passed to v_make_prime.");
		return w;
	}
	else{
		Integer g=v_gcd(v);
		if (g==0) {
			return w;
		}
		else {
			for (i = 0; i < j; i++) {
				w[i]=v[i]/g;
			}
		}
		return w;
	}
}

//---------------------------------------------------------------------------

vector<Integer> v_make_prime(const vector<Integer>& v,Integer& g){
	int i,j=v.size();
	vector<Integer> w(j,0);
	if (j==0) {
		v_error("error: Bad argument passed to v_make_prime.");
		return w;
	}
	else{
		g=v_gcd(v);
		if (g==0) {
			return w;
		}
		else {
			for (i = 0; i < j; i++) {
				w[i]=v[i]/g;
			}
		}
		return w;
	}
}

//---------------------------------------------------------------------------

void v_scalar_multiplication(vector<Integer>& v, const Integer& scalar){
	int i,j=v.size();
	for (i = 0; i <j; i++) {
		v[i]=v[i]*scalar;
	}
}

vector<Integer> v_scalar_multiplication_two(const vector<Integer>& v, const Integer& scalar){
	int i,j=v.size();
	vector<Integer> w(j);
	for (i = 0; i <j; i++) {
		w[i]=v[i]*scalar;
	}
	return w;
}

//---------------------------------------------------------------------------

void v_scalar_division(vector<Integer>& v, const Integer& scalar){
	int i,j=v.size();
	for (i = 0; i <j; i++) {
		if (v[i]%scalar!=0) {
			v_error("error: Bad argument passed to vector_operations::v_scalar_division.");
			return;
		}
		else{
			v[i]=v[i]/scalar;
		}
	}
}

//---------------------------------------------------------------------------

void v_reduction_modulo(vector<Integer>& v, const Integer& modulo){
	int i,j=v.size();
	for (i = 0; i <j; i++) {
		v[i]=v[i]%modulo;
		if (v[i]<0) {
			v[i]=v[i]+modulo;
		}
	}
}

//---------------------------------------------------------------------------

int v_diference_ordered_fast(const vector<int>& u,const vector<int>& v){
	int i,j,k, s=u.size();
	i=0;
	// erste Stelle suchen, an der Unterschied
	while (u[i]==v[i]){
		i++ ;
	}
	j=s-1;
	// letzte Stelle suchen, an der Unterschied
	while (u[j]==v[j]){
		j--;
	}
	if(i==j)
		return u[i];

	if(u[i]<v[i])    // bei "true" kommt u[i] nicht in v vor, sonst v[i] nicht in u
	{
		if(u[j]>v[j]) // bei "true" kommt auch u[j] nicht in v vor
			return 0;  // jetzt u[i] und u[j] nicht in v
		else
		{
			for(k=i+1;k<=j;k++)    // u[i] nicht in v und wir pruefen
				if(u[k]!=v[k-1])  // ob u[i+1],...,u[j] in v
					return 0;  // nein
			return u[i];           // es fehlt wirklich nur u[i]
		}
	}
	else  // v[i] nicht in u
	{
		if(u[j]<v[j]) // bei true kommt auch v[j] nicht in u vor
			return 0;   // jetzt v[i] und v[j] nicht in u
		else
		{
			for(k=i;k<=j-1;k++)        // v[i] nicht in u
				if(u[k]!=v[k+1]) // und wir pruefen, ob v[i+1],...v[j] in u
					return 0;  // nein
			return u[j];         // es fehlt nur u[j] in v
		}
	}
}

//---------------------------------------------------------------------------

bool v_test_scalar_product(const vector<Integer>& av,const vector<Integer>& bv, const Integer& result, const int& m){
	Integer ans = 0;
	register int i,n=av.size();
	vector<Integer>::const_iterator    a=av.begin(),b=bv.begin();

	if( n >= 64 )
	{
		for( i = 0; i < ( n >> 6 ); ++i, a += 64, b += 64 ){
			ans += a[0] * b[0];
			ans += a[1] * b[1];
			ans += a[2] * b[2];
			ans += a[3] * b[3];
			ans %=m;
			ans += a[4] * b[4];
			ans += a[5] * b[5];
			ans += a[6] * b[6];
			ans += a[7] * b[7];
			ans %=m;
			ans += a[8] * b[8];
			ans += a[9] * b[9];
			ans += a[10] * b[10];
			ans += a[11] * b[11];
			ans %=m;
			ans += a[12] * b[12];
			ans += a[13] * b[13];
			ans += a[14] * b[14];
			ans += a[15] * b[15];
			ans %=m;
			ans += a[16] * b[16];
			ans += a[17] * b[17];
			ans += a[18] * b[18];
			ans += a[19] * b[19];
			ans %=m;
			ans += a[20] * b[20];
			ans += a[21] * b[21];
			ans += a[22] * b[22];
			ans += a[23] * b[23];
			ans %=m;
			ans += a[24] * b[24];
			ans += a[25] * b[25];
			ans += a[26] * b[26];
			ans += a[27] * b[27];
			ans %=m;
			ans += a[28] * b[28];
			ans += a[29] * b[29];
			ans += a[30] * b[30];
			ans += a[31] * b[31];
			ans %=m;
			ans += a[32] * b[32];
			ans += a[33] * b[33];
			ans += a[34] * b[34];
			ans += a[35] * b[35];
			ans %=m;
			ans += a[36] * b[36];
			ans += a[37] * b[37];
			ans += a[38] * b[38];
			ans += a[39] * b[39];
			ans %=m;
			ans += a[40] * b[40];
			ans += a[41] * b[41];
			ans += a[42] * b[42];
			ans += a[43] * b[43];
			ans %=m;
			ans += a[44] * b[44];
			ans += a[45] * b[45];
			ans += a[46] * b[46];
			ans += a[47] * b[47];
			ans %=m;
			ans += a[48] * b[48];
			ans += a[49] * b[49];
			ans += a[50] * b[50];
			ans += a[51] * b[51];
			ans %=m;
			ans += a[52] * b[52];
			ans += a[53] * b[53];
			ans += a[54] * b[54];
			ans += a[55] * b[55];
			ans %=m;
			ans += a[56] * b[56];
			ans += a[57] * b[57];
			ans += a[58] * b[58];
			ans += a[59] * b[59];
			ans %=m;
			ans += a[60] * b[60];
			ans += a[61] * b[61];
			ans += a[62] * b[62];
			ans += a[63] * b[63];
			ans %=m;
		}
		n -= i << 6;
	}
	for( i = 0; i < n; ++i ){
		ans += a[i] * b[i];
		ans %=m;
	}
	if (((result-ans) % m)!=0) {
		return false;
	}
	return true;
}

//---------------------------------------------------------------------------

void v_error(string s){
	cerr <<"\nVector "<< s<<"\n";
	global_error_handling();
}

//---------------------------------------------------------------------------

