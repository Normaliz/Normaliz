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

#include <iostream>
#include <string>
#include <algorithm>
#include <list>

#include "libQnormaliz/Qinteger.h"
#include "libQnormaliz/Qvector_operations.h"
#include "libQnormaliz/Qmatrix.h"

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

//---------------------------------------------------------------------------

template<typename Number>
Number v_scalar_product(const vector<Number>& av,const vector<Number>& bv){
    //loop stretching ; brings some small speed improvement

    Number ans = 0;
    size_t i,n=av.size();

    typename vector<Number>::const_iterator a=av.begin(), b=bv.begin();

    if( n >= 16 )
    {
        for( i = 0; i < ( n >> 4 ); ++i, a += 16, b +=16 ){
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
        }

        n -= i<<4;
    }

    if( n >= 8)
    {
        ans += a[0] * b[0];
        ans += a[1] * b[1];
        ans += a[2] * b[2];
        ans += a[3] * b[3];
        ans += a[4] * b[4];
        ans += a[5] * b[5];
        ans += a[6] * b[6];
        ans += a[7] * b[7];

        n -= 8;
        a += 8;
        b += 8;
    }

    if( n >= 4)
    {
        ans += a[0] * b[0];
        ans += a[1] * b[1];
        ans += a[2] * b[2];
        ans += a[3] * b[3];

        n -= 4;
        a += 4;
        b += 4;
    }

    if( n >= 2)
    {
        ans += a[0] * b[0];
        ans += a[1] * b[1];

        n -= 2;
        a += 2;
        b += 2;
    }

    if(n>0)
        ans += a[0]*b[0];
        
    return ans;
}

//---------------------------------------------------------------------------

template<typename Number>
Number v_scalar_product_unequal_vectors_end(const vector<Number>& a,const vector<Number>& b){
    Number ans = 0;
    size_t i,n=a.size(),m=b.size();
    for (i = 1; i <= n; i++) {
        ans+=a[n-i]*b[m-i];
    }
    return ans;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> v_add(const vector<Number>& a,const vector<Number>& b){
   assert(a.size() == b.size());
    size_t i,s=a.size();
    vector<Number> d(s);
    for (i = 0; i <s; i++) {
        d[i]=a[i]+b[i];
    }
    return d;
}

//---------------------------------------------------------------------------

template<typename Number>
void v_add_result(vector<Number>& result, const size_t s, const vector<Number>& a,const vector<Number>& b){
   assert(a.size() == b.size() && a.size() == result.size());
    size_t i;
    // vector<Number> d(s);
    for (i = 0; i <s; i++) {
        result[i]=a[i]+b[i];
    }
    // return d;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number>& v_abs(vector<Number>& v){
    size_t i, size=v.size();
    for (i = 0; i < size; i++) {
        if (v[i]<0) v[i] = Iabs(v[i]);
    }
    return v;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<Number> v_abs_value(vector<Number>& v){
    size_t i, size=v.size();
    vector<Number> w=v;
    for (i = 0; i < size; i++) {
        if (v[i]<0) w[i] = Iabs(v[i]);
    }
    return w;
}

//---------------------------------------------------------------------------

// the following function removes the denominators and then extracts the Gcd of the numerators
mpq_class v_simplify(vector<mpq_class>& v){
    size_t size=v.size();
    mpz_class d=1;
    for (size_t i = 0; i < size; i++)
        d=lcm(d,v[i].get_den());  // we use the GMP function
    for (size_t i = 0; i < size; i++)
        v[i]*=d;
    mpz_class g=0;
    for (size_t i = 0; i < size; i++)
        g=gcd(g,v[i].get_num());  //  we use the GMP function
    if (g==0)
        return 0;
    for (size_t i = 0; i < size; i++)
        v[i]/=g;
    return 1;
}
//---------------------------------------------------------------------------

template<typename Number>
void v_scalar_division(vector<Number>& v, const Number& scalar){
    size_t i,size=v.size();
    for (i = 0; i <size; i++) {
        v[i] /= scalar;
    }
}

//---------------------------------------------------------------------------

template<typename T>
vector<T> v_merge(const vector<T>& a, const T& b) {
    size_t s=a.size();
    vector<T> c(s+1);
    for (size_t i = 0; i < s; i++) {
        c[i]=a[i];
    }
    c[s] = b;
    return c;
}

//---------------------------------------------------------------------------

template<typename T>
vector<T> v_merge(const vector<T>& a,const vector<T>& b){
    size_t s1=a.size(), s2=b.size(), i;
    vector<T> c(s1+s2);
    for (i = 0; i < s1; i++) {
        c[i]=a[i];
    }
    for (i = 0; i < s2; i++) {
        c[s1+i]=b[i];
    }
    return c;
}
//---------------------------------------------------------------------------

template<typename T>
vector<T> v_cut_front(const vector<T>& v, size_t size){
    size_t s,k;
    vector<T> tmp(size);
    s=v.size()-size;
    for (k = 0; k < size; k++) {
        tmp[k]=v[s+k];
    }
    return tmp;
}

//---------------------------------------------------------------------------

template<typename Number>
vector<key_t> v_non_zero_pos(const vector<Number>& v){
    vector<key_t> key;
    size_t size=v.size();
    key.reserve(size);
    for (key_t i = 0; i <size; i++) {
        if (v[i]!=0) {
            key.push_back(i);
        }
    }
    return key;
}

//---------------------------------------------------------------------------

template<typename Number>
bool v_is_zero(const vector<Number>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] != 0) return false;
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Number>
bool v_is_symmetric(const vector<Number>& v) {
    for (size_t i = 0; i < v.size()/2; ++i) {
        if (v[i] != v[v.size()-1-i]) return false;
    }
    return true;
}

//---------------------------------------------------------------------------

template<typename Number>
bool v_is_nonnegative(const vector<Number>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] <0) return false;
    }
    return true;
}


//---------------------------------------------------------------------------

template<typename Number>
void v_el_trans(const vector<Number>& av,vector<Number>& bv, const Number& F, const size_t& start){

    size_t i,n=av.size();

    typename vector<Number>::const_iterator a=av.begin();
    typename vector<Number>::iterator b=bv.begin();

    a += start;
    b += start;
    n -= start;


    if( n >= 8 )
    {
        for( i = 0; i < ( n >> 3 ); ++i, a += 8, b += 8 ){
            b[0] += F*a[0];
            b[1] += F*a[1];
            b[2] += F*a[2];
            b[3] += F*a[3];
            b[4] += F*a[4];
            b[5] += F*a[5];
            b[6] += F*a[6];
            b[7] += F*a[7];
        }
        n -= i << 3;
    }

    if( n >= 4)
    {
        b[0] += F*a[0];
        b[1] += F*a[1];
        b[2] += F*a[2];
        b[3] += F*a[3];

        n -=4;
        a +=4;
        b +=4;
    }

    if( n >= 2)
    {
        b[0] += F*a[0];
        b[1] += F*a[1];

        n -=2;
        a +=2;
        b +=2;
    }

    if(n>0)
        b[0] += F*a[0];
}

//---------------------------------------------------------------

vector<bool> v_bool_andnot(const vector<bool>& a, const vector<bool>& b) {
    assert(a.size() == b.size());
    vector<bool> result(a);
    for (size_t i=0; i<b.size(); ++i) {
        if (b[i])
            result[i]=false;
    }
    return result;
}

// swaps entry i and j of the vector<bool> v
void v_bool_entry_swap(vector<bool>& v, size_t i, size_t j) {
    if (v[i] != v[j]) {
        v[i].flip();
        v[j].flip();
    }
}


vector<key_t> identity_key(size_t n){
    vector<key_t> key(n);
    for(size_t k=0;k<n;++k)
        key[k]=k;
    return key;
}

//---------------------------------------------------------------
// Sorting

template <typename T>
void order_by_perm(vector<T>& v, const vector<key_t>& permfix){
    
    vector<key_t> perm=permfix; // we may want to use permfix a second time
    vector<key_t> inv(perm.size());
    for(key_t i=0;i<perm.size();++i)
        inv[perm[i]]=i;
    for(key_t i=0;i<perm.size();++i){
        key_t j=perm[i];
        swap(v[i],v[perm[i]]);        
        swap(perm[i],perm[inv[i]]);        
        swap(inv[i],inv[j]);                
    }
}

// vector<bool> is special
template <>
void order_by_perm(vector<bool>& v, const vector<key_t>& permfix){
    
    vector<key_t> perm=permfix; // we may want to use permfix a second time
    vector<key_t> inv(perm.size());
    for(key_t i=0;i<perm.size();++i)
        inv[perm[i]]=i;
    for(key_t i=0;i<perm.size();++i){
        key_t j=perm[i];
        // v.swap(v[i],v[perm[i]]);
        v_bool_entry_swap(v,i,perm[i]);
        swap(perm[i],perm[inv[i]]);        
        swap(inv[i],inv[j]);                
    }
}

// make random vector of length n with entries between -m and m
template <typename Number>
vector<Number> v_random(size_t n, long m){
    vector<Number> result(n);
    for(size_t i=0;i<n;++i)
        result[i]=rand()%(2*m+1)-m;
    return result;    
}

template bool v_is_nonnegative<long>(const vector<long>&);
template bool v_is_nonnegative<long long>(const vector<long long>&);
template bool v_is_nonnegative<mpz_class>(const vector<mpz_class>&);

template bool v_is_symmetric<long>(const vector<long>&);
template bool v_is_symmetric<long long>(const vector<long long>&);
template bool v_is_symmetric<mpz_class>(const vector<mpz_class>&);
// 

template void v_add_result<long     >(vector<long     >&, size_t, const vector<long     >&, const vector<long     >&);
template void v_add_result<long long>(vector<long long>&, size_t, const vector<long long>&, const vector<long long>&);
template void v_add_result<mpz_class>(vector<mpz_class>&, size_t, const vector<mpz_class>&, const vector<mpz_class>&);

template long v_scalar_product(const vector<long>& a,const vector<long>& b);
template long long v_scalar_product(const vector<long long>& a,const vector<long long>& b);
template mpz_class v_scalar_product(const vector<mpz_class>& a,const vector<mpz_class>& b);

} // end namespace libQnormaliz
