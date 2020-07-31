#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

#include "libnormaliz/libnormaliz.h"

using namespace libnormaliz;

typedef long long Integer;

int main(int argc, char* argv[]) {
    Matrix<Integer> Gens = readMatrix<Integer>(string("small_gens.mat"));

    /*Cone<Integer> C(Type::cone,Gens);
    C. setVerbose(true);

    C.compute(ConeProperty::DefaultMode);

    C.getVerticesFloatMatrix().pretty_print(cout);

    C.setExpansionDegree(25);
    cout << C.getHilbertSeries().getExpansion();
    cout << "========================" << endl;

    C.setExpansionDegree(50);
    cout << C.getHilbertSeries().getExpansion();
    cout << "========================" << endl;

    C.compute(ConeProperty::HSOP);

    cout << "HSOP " << C.getHilbertSeries().getHSOPNum();

    C.setNrCoeffQuasiPol(2);
    cout << "========================" << endl;
    Matrix<mpz_class> Q(C.getHilbertSeries().getHilbertQuasiPolynomial());
    Q.pretty_print(cout);
    cout << "========================" << endl;
    C.setNrCoeffQuasiPol(5);
    Q=C.getHilbertSeries().getHilbertQuasiPolynomial();
    Q.pretty_print(cout);*/

    vector<Cone<Integer> > ParCones(16);
#pragma omp parallel for
    for (size_t i = 0; i < ParCones.size(); ++i) {
        ParCones[i] = Cone<Integer>(Type::cone, Gens);
        // ParCones[i].setVerbose(true);
        switch (i % 8) {
            case 0:
                ParCones[i].compute(ConeProperty::DefaultMode);
                break;
            case 1:
                ParCones[i].compute(ConeProperty::DualMode, ConeProperty::Deg1Elements);
                break;
            case 2:
                ParCones[i].compute(ConeProperty::Projection, ConeProperty::Deg1Elements);
                break;
            case 3:
                ParCones[i].compute(ConeProperty::ProjectionFloat, ConeProperty::Deg1Elements);
                break;
            case 4:
                ParCones[i].compute(ConeProperty::Approximate, ConeProperty::Deg1Elements, ConeProperty::IsGorenstein);
                break;
            case 5:
                ParCones[i].compute(ConeProperty::SupportHyperplanes);
                break;
            case 6:
                ParCones[i].compute(ConeProperty::IntegerHull);
                break;
            case 7:
                ParCones[i].compute(ConeProperty::IsIntegrallyClosed);
                break;
            default:
                break;
        }
    }

}  // end main
