#include <stdlib.h>
#include <vector>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

#include "libnormaliz/libnormaliz.h"

using namespace libnormaliz;

typedef long long Integer;


int main(int argc, char* argv[]){

    Matrix<Integer> First(24); // =readMatrix<Integer>("first.mat");
    Cone<Integer> MyCone(Type::inequalities, First);
    MyCone.setVerbose(true);
    MyCone.setFaceCodimBound(2);
    MyCone.compute(ConeProperty::Dynamic, ConeProperty::FaceLattice);    
    Matrix<Integer> Second=readMatrix<Integer>("second.mat");
    Second.pretty_print(cout);
    MyCone.addInput(Type::inequalities,Second);
    MyCone.compute(ConeProperty::FaceLattice);    
    Matrix<Integer> Third=readMatrix<Integer>("third.mat");
    MyCone.addInput(Type::equations,Third);
    Third.pretty_print(cout);
    MyCone.setFaceCodimBound(1);
    MyCone.compute(ConeProperty::FaceLattice);
    
    MyCone.write_cone_output("MyCone");
}  //end main
