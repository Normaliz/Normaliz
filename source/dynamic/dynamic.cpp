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

    Matrix<Integer> First=readMat("first.mat");
    Cone<Integer> MyCone(Type::inequalities, First);
    MyCone.setFaceCodimBound(2);
    MyCone.compute(ConeProperty::Dynamic, ConeProperty::FaceLattice);    
    Matrix<Integer> Second("second.mat");
    MyCone.addInput(Type::inequalities,Second);
    MyCone.compute(ConeProperty::FaceLattice);    
    Matrix<Integer> Third("third.mat");
    MyCone.addInput(Type::equation,Third);
    MyCone.setFaceCodimBound(1);
    MyCone.compute(ConeProperty::FaceLattice); 
}  //end main
