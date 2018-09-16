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

#ifndef QINTEGER_H_
#define QINTEGER_H_

#include <libQnormaliz/Qgeneral.h>

#ifdef ENFNORMALIZ
#include <e-antic/renfxx.h>
#endif

#include <list>
#include <vector>
#include <string>
#include <limits.h>

// #include "libnormaliz/integer.h"

//---------------------------------------------------------------------------

namespace libQnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//                     Basic functions
//---------------------------------------------------------------------------

// returns the absolute value of a
template<typename Number> inline Number Iabs(const Number& a) {
    return (a>=0) ? (a) : Number(-a);
}

//---------------------------------------------------------------------------
//                     Conversions and checks
//---------------------------------------------------------------------------

// convert val to ret
// does the conversion and returns false if it fails
bool try_convert(long& ret, const long long& val);
inline bool try_convert(long long& ret, const long& val) {assert(false); return true;}
bool try_convert(long& ret, const mpz_class& val);
bool try_convert(long long& ret, const mpz_class& val);
inline bool try_convert(mpz_class& ret, const long& val) {assert(false); return true;}
bool try_convert(mpz_class& ret, const long long& val);

// template for same typ "conversion"
template<typename Type>
inline bool try_convert(Type& ret, const Type& val) {ret = val; return true;}


bool fits_long_range(long long a);

template<typename Number>
inline bool using_GMP() {
  return false;
}

template<>
inline bool using_GMP<mpq_class>() {
  return true;
}

template<typename Number>
inline bool using_renf() {
  return false;
}

#ifdef ENFNORMALIZ
template<>
inline bool using_renf<renf_elem_class>() {
  return true;
}
#endif

//---------------------------------------------------------------------------

// Should be completely remoced:
template<typename Number>
inline bool check_range(const Number& m) {
    return true;
}

//---------------------------------------------------------------------------
//                     Special functions
//---------------------------------------------------------------------------

//return the number of decimals, needed to write the Number a
template<typename Number> size_t decimal_length(Number a);

//returns b!/a!
template<typename Number> Number permutations(const size_t& a, const size_t& b);
template<typename Number> Number permutations_modulo(const size_t& a, const size_t& b, long m);

//---------------------------------------------------------------------------
//                     String conversion functions
//---------------------------------------------------------------------------

// forward declaration to silence clang error:
// 'operator<<' should be declared prior to the call site or in the global namespace
template <typename T> std::ostream& operator<< (std::ostream& out, const vector<T>& vec);

template<typename Number> string toString(Number a) {
    ostringstream ostream;
    ostream << a;
    return ostream.str();
}
/* template<> inline string toString(mpz_class a) {
    return a.get_str();
}*/
template<> inline string toString(mpq_class a) {
    return a.get_str();
}

//-------------------------------------------------------

double MPQ_to_nmz_float(const mpq_class& val){
    mpz_class bound=1;
    for(size_t i=0;i<60; ++i)
        bound*=10;
    mpz_class gmp_num=val.get_num(),gmp_den=val.get_den();
    while(Iabs(gmp_num) > bound && Iabs(gmp_den) > bound){
            gmp_num/=10;
            gmp_den/=10;
    }
    double num,den;
    // convert(num,gmp_num);
    num=gmp_num.get_d();
    // convert(den,gmp_den);
    den=gmp_den.get_d();
    return num/den;
}

mpq_class DEC_fraction_to_mpq(string s){

    size_t skip=0;               // skip leading spaces
    for(;skip<s.size();++skip){
            if(!isspace(s[skip]))
                break;
    }
    s=s.substr(skip);

    mpz_class sign=1;
    if(s[0]=='+')
        s=s.substr(1);
    else
        if(s[0]=='-'){
            s=s.substr(1);
            sign=-1;
    }
    
    if(s[0]=='+' || s[0]=='-')
            throw BadInputException("Error in decimal fraction "+s);
    
    string int_string,frac_string,exp_string;
    size_t frac_part_length=0;
    size_t pos_point=s.find(".");
    size_t pos_E=s.find("e");
    if(pos_point!=string::npos){
        int_string=s.substr(0,pos_point);
        if(pos_E!=string::npos){
            frac_part_length=pos_E-(pos_point+1);
        }
        else
            frac_part_length=s.size()-(pos_point+1);
        frac_string=s.substr(pos_point+1,frac_part_length);
        if(frac_string[0]=='+' || frac_string[0]=='-')
            throw BadInputException("Error in decimal fraction "+s);
    }
    else
        int_string=s.substr(0,pos_E);
    if(pos_E!=string::npos)
        exp_string=s.substr(pos_E+1,s.size()-(pos_E+1));
    
    /* cout << "int  " << int_string << endl;
    cout << "frac " << frac_string << endl;
    cout << "exp  " << exp_string << endl; */
    
    // remove leading 0
    while(int_string[0]=='0')
        int_string=int_string.substr(1);
    while(frac_string[0]=='0')
        frac_string=frac_string.substr(1);
    
    mpq_class int_part, frac_part, exp_part;
    if(!int_string.empty())
        int_part=mpz_class(int_string);

    if(pos_E==0)
        int_part=1;
    
    // cout << "int_part " << int_part << endl;
    
    mpz_class den=1;
    if(!frac_string.empty()){
        frac_part=mpz_class(frac_string);
        for(size_t i=0;i<frac_part_length;++i)
            den*=10;        
    }
    // cout << "frac_part " << frac_part << endl;
    mpq_class result=int_part;
    if(frac_part!=0)
        result+=frac_part/den;
    if(!exp_string.empty()){
        mpz_class expo(exp_string); // we take mpz_class because it has better error checking
        // long expo=stol(exp_string);        
        mpz_class abs_expo=Iabs(expo);
        // cout << "expo " << expo << endl;
        mpz_class factor=1;
        for(long i=0;i< abs_expo;++i)
            factor*=10;
        if(expo>=0)
            result*=factor;
        else
            result/=factor;
    }
    /* cout <<" result " << sign*result << endl;
    cout << "==========" << endl; */
    return sign*result;
}

#ifdef ENFNORMALIZ
inline mpq_class approx_to_mpq(const renf_elem_class& x){

    stringstream str_str;
    str_str << x;
    string str=str_str.str();

    string nf_str, approx_str;
    bool rational=true;
    bool nf_finished=false;
    for(size_t i=0;i<str.size();++i){
        if(str[i]=='a')
            rational=false;
        if(str[i]=='(' || str[i]==')')
            continue;
        if(str[i]=='~' || str[i]=='='){
            nf_finished=true;
            continue;
        }
        if(nf_finished)
            approx_str+=str[i];
        else
            nf_str+=str[i];
        
    }
    if(rational)
        return mpq_class(nf_str);
    else{
        return DEC_fraction_to_mpq(approx_str);        
    }
}
#endif

inline mpq_class approx_to_mpq(const mpq_class& x){
        return x;
}

template<typename Number>
double approx_to_double(const Number& x){
        return MPQ_to_nmz_float(approx_to_mpq(x));
}

template<typename Number>
vector<mpq_class> approx_to_mpq(const vector<Number>& ori){
    
    vector<mpq_class> res(ori.size());
    for(size_t i=0;i<ori.size();++i)
        res[i]=approx_to_mpq(ori[i]);
    return res;
}



} // end libnormaliz

//---------------------------------------------------------------------------
#endif /* INTEGER_H_ */
//---------------------------------------------------------------------------
