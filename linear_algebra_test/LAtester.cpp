/*
 * Test for the linear algebra of Normaliz 2.5
 * Copyright (C) 2010-2014 Christof Soeger
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

#include <stdlib.h>
#include <time.h>
#include <vector>

#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

#include "../source/libnormaliz/matrix.h"
#include "../source/libnormaliz/libnormaliz.h"
#include "../source/libnormaliz/integer.h"
#include "../source/libnormaliz/vector_operations.h"
#include "gmpxx.h"
using namespace libnormaliz;


//---------------------------------------------------------------------------

mpz_class binary_length(mpz_class a) {
    mpz_class log2 = 1; //for the sign
    while (a > 0) {
        a /= 2;
        ++log2;
    }
    return log2;
}

//---------------------------------------------------------------------------

int main(int argc, char* argv[]) {

#ifdef USE_MPZ
    typedef mpz_class I;
#else
    typedef long long I;
#endif
    int N = 100000; // iterations
    if (argc >= 3) istringstream(argv[2]) >> N;

    libnormaliz::test_arithmetic_overflow = false;  //default: false
    
    ifstream fin(argv[1],ifstream::in);
    if (!fin.is_open()) {
        cerr << "Cannot open file \"" << argv[1] << "\"."<< endl;
        exit(-1);
    }
    int nr, nc;
    fin >> nr;
    fin >> nc;
    Matrix<I> M(nr,nc);

    for(int i=0; i<nr; i++) {
        for(int j=0; j<nc; j++) {
            fin >> M[i][j];
        }
    }
    fin.close();
    if (M.rank()!=nr) {
        vector<unsigned int> key = M.max_rank_submatrix_lex();
        M = M.submatrix(key);
        nr = M.nr_of_rows();
    }

    if (nr<nc) {
        cerr << "Need a matrix with rank = nr columns!" << endl;
        exit(-2);
    }
    //M.print(cout);

    vector<I> diag(nc);
    I denom;
    size_t rank_res;

    Matrix<I> Result;
    Matrix<I> Copy(M);
    Matrix<I> En(nc);
    Matrix<I> Empty(nc,0);
    Matrix<I> RS_copy(En);

    // compute norms
    mpz_class m=0, a;
    long d = nc; 
    for (size_t i=0; i<nc; ++i) {
        for (size_t j=0; j<nc; ++j) {
            a = abs(to_mpz(M[i][j]));
            if (a>m) m=a;
        }
    }
    mpz_class lb =1;
    for (long i=1; i<=d; ++i) lb *= to_mpz(i)*to_mpz(d);
    cout << "max="<< m ;//<< "  d! * m^d = " << lb; 
    cout << "   " << binary_length(lb) << " bits" << endl;

    mpz_class hwb = 1;
    for (long i=0; i<d; ++i) hwb *= to_mpz(m);
    for (long i=0; i<d/2; ++i) hwb *= to_mpz(d);
    if (d%2 == 1) hwb *= sqrt((double)d);
    cout << "weak hadamard bound " << binary_length(hwb) << " bits"<< endl;

    cout << "squared row norms: ";
    mpz_class sq_v_norm, sq_norm=1, v_M=0;
    for (size_t i=0; i<nc; ++i) {
        sq_v_norm=0;
        for (size_t j=0; j<nc; ++j) 
            sq_v_norm += to_mpz(M[i][j])*to_mpz(M[i][j]);
        cout << sq_v_norm << " ";
        sq_norm *= sq_v_norm;
        if (sq_v_norm > v_M) v_M = sq_v_norm;
    }
    cout << endl;
//  cout << "squared norm: " << sq_norm << endl;
    mpz_class s;
    mpz_sqrt (s.get_mpz_t(), sq_norm.get_mpz_t());
    cout << "norm: " << s << "   " << binary_length(s) << " bits"<< endl;

    cout << "M = "<<v_M<<endl;
    cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    // prep run
    Result = M.invert(diag, denom);
    Result = M.solve(RS_copy, diag, denom);

    for (size_t i=0; i<nc; ++i) if (diag[i]<0) diag[i] = - diag[i];
    cout << "diag " << diag;
    int r=0;
    for (size_t i=0; i<nc; ++i) if (diag[i]!=1) r++;
    Matrix<I> RS2(nc,r);
    size_t j = 0;
    for (size_t i=0; i<nc; ++i) if (diag[i]>1) RS2[i,j++];
    Matrix<I> RS2_copy(RS2);

    RS_copy = En;
    Result = Copy.solve_destructive(RS_copy, diag, denom);
    Copy = M;
    denom = Copy.vol_destructive();
    rank_res = M.rank();

    clock_t t;


/*    cout << "Inverting Matrix              " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Result = M.invert(diag, denom);
    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;


    cout << "Inverting by solve            " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        RS_copy = En;
        Result = M.solve(RS_copy, diag, denom);

    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;


    cout << "Inverting by solve_dest       " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        RS_copy = En;
        Result = Copy.solve_destructive(RS_copy, diag, denom);

    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;

*/
    cout << "Inverting by solve_dest_Sol   " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        RS_copy = En;
        Copy.solve_destructive_Sol(RS_copy, diag, denom, Result);

    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;


    cout << "Solving for "<<r<<" right sides     " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        RS2_copy = RS2;
        Copy.solve_destructive_Sol(RS2_copy, diag, denom, Result);

    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;

    cout << "Trigonalize                   " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        Copy.solve_destructive_Sol(Empty, diag, denom, Result);

    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;


    cout << "vol_destructive               " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        denom = Copy.vol_destructive();
    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;

/*
    cout << "vol_bareiss_destructive       " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        denom = Copy.vol_bareiss_destructive();
    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl ;
*/


    cout << "rank                          " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        rank_res = M.rank();
    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl;


    cout << "rank_destructive              " << N << " times..." << flush;
    t = clock();
    for(int i=0; i< N; i++) {
        Copy = M;
        rank_res = Copy.rank_destructive();
    }
    t = clock() - t;
    cout << " took "<< ((float)t)/CLOCKS_PER_SEC <<" sec." << endl;

    return 0;
}
