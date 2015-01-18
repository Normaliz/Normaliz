#include <stdlib.h>
#include <vector>
#include <map>

#include <list>

#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <omp.h>
using namespace std;

#include "libnormaliz.h" 
#include "cone.h"
#include "vector_operations.h"
#include "lineare_transformation.h"
#include "cone_property.h"
#include "integer.h"
#include "matrix.h"
//#include "libnormaliz/libnormaliz.cpp"
using namespace libnormaliz;

template<typename Integer>
Matrix<Integer>  ReadMat(string project){
// reads one matrix from .in file

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
    Matrix<Integer> result(nrows,ncols);

    for(i=0;i<nrows;++i)
        for(j=0;j<ncols;++j){
            in >> entry;
            result[i][j]=entry;
        }
    
    
    return result;
}




int main(int argc, char* argv[])
{

    // typedef long long Integer;
    typedef mpz_class Integer;
    bool success;
    
    string project=string(argv[1]);
    
    Matrix<Integer> M=ReadMat<Integer>(project);
    M.pretty_print(cout);
    Matrix<Integer> N=M;
    
    size_t rank;
    
    Matrix<Integer> Transf=N.SmithNormalForm(rank);
    
    Integer denom;
    
    Transf.pretty_print(cout); 
    
    cout <<"--------------------" << endl;
    
    N.pretty_print(cout);
    
    cout <<"--------------------" << endl;
    
    Transf.invert(denom).pretty_print(cout);
    
     cout <<"--------------------" << endl;
     
     Transf.invert(denom).multiplication(Transf).pretty_print(cout);
    
    

    
        /*
    
    for(size_t i=0;i<100000;++i){
    
        N=M;
        N.max_rank_submatrix_lex();
    }    

    
    // cout << "Vol " << N.row_echelon_inner_gcd(success) << endl;
    // N.pretty_print(cout);
    
        //cout << "Rang " << N.row_echelon_bareiss(success) << endl;
        
    // N.rank_destructive();
    
    for(size_t i=0;i<1000000;++i){
    
        N=M;
        // N.row_echelon_inner_gcd(success);
        // N.rank_destructive();
        N.row_echelon_inner_elem(success);
        // N.row_echelon_inner_bareiss(success);
        
        
    }
        

    N.pretty_print(cout);
    
    */
    
    
    exit(0);
    
}
