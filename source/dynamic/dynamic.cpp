#include <cstdlib>
#include <vector>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

#include "libnormaliz/libnormaliz.h"

using namespace libnormaliz;

typedef long long Integer;

int main(int argc, char* argv[]) {
    Matrix<Integer> First(24);  // =readMatrix<Integer>("first.mat");
    Cone<Integer> MyCone(Type::inequalities, First);
    MyCone.setVerbose(true);
    MyCone.setFaceCodimBound(2);
    MyCone.compute(ConeProperty::Dynamic, ConeProperty::FaceLattice);
    Matrix<Integer> Second = readMatrix<Integer>("second.mat");
    Second.pretty_print(cout);
    MyCone.modifyCone(Type::inequalities, Second);
    MyCone.compute(ConeProperty::FaceLattice);
    Matrix<Integer> Facets = MyCone.getSupportHyperplanesMatrix();
    MyCone.write_cone_output("MyConeAfterSecond");
    map<dynamic_bitset, int> FL = MyCone.getFaceLattice();
    auto FaceIt = FL.end();
    FaceIt--;
    dynamic_bitset Indicator = FaceIt->first;
    cout << "Codim of last face " << FaceIt->second << endl;
    cout << "Indicator of last face " << Indicator << endl;
    size_t dim = MyCone.getEmbeddingDim();
    Matrix<Integer> FaceEq(0, dim);
    for (size_t i = 0; i < Indicator.size(); ++i) {
        if (Indicator[i])
            FaceEq.append(Facets[i]);
    }
    cout << "Equations of last face" << endl;
    FaceEq.pretty_print(cout);
    Cone<Integer> FaceCone(Type::inequalities, Facets, Type::equations, FaceEq);
    FaceCone.compute(ConeProperty::ExtremeRays);
    cout << " Extreme rays of last face " << endl;
    FaceCone.getExtremeRaysMatrix().pretty_print(cout);
    Matrix<Integer> Third = readMatrix<Integer>("third.mat");
    MyCone.modifyCone(Type::equations, Third);
    MyCone.setFaceCodimBound(1);
    MyCone.compute(ConeProperty::FaceLattice);

    MyCone.write_cone_output("MyConeAfterThird");
}  // end main
