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
    cout << "---------" << endl;
    Matrix<Integer> N=M;
    
    size_t dim=M.nr_of_columns();
    
    vector<Integer> Grading(dim,0);
    Grading [Grading.size()-1]=1;
    
    Matrix<Integer> WeightsGradL1(0,dim);  // weight matrix for ordering, first row Grading, then L1
    // if(isComputed(ConeProperty::Grading))
        WeightsGradL1.append(Grading);
    WeightsGradL1.append(vector<Integer>(dim,1));
    vector<bool> AbsWGL(WeightsGradL1.nr_of_rows(),false);
    AbsWGL[AbsWGL.size()-1]=true;
    
    Matrix<Integer> WeightsL1(0,dim);  // only L1
    WeightsL1.append(vector<Integer>(dim,1));
    vector<bool> AbsWL(WeightsL1.nr_of_rows(),false);
    AbsWL[AbsWL.size()-1]=true;
 
    WeightsGradL1.pretty_print(cout);
    cout << "---------" << endl;
    WeightsL1.pretty_print(cout);
    cout << "---------" << endl;
    
    M.sort_by_weights(WeightsL1,AbsWL);
    M.pretty_print(cout);
    
    
    exit(0);
    
}
