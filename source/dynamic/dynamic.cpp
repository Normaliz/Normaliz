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

int main(int argc, char* argv[]) {

    vector<vector<Integer> > inhom_ineq={{2,0,1}};
    Cone<Integer> C(Type::inhom_inequalities, inhom_ineq);
    C.compute(ConeProperty::IntegerHull);
    Cone<Integer> D=C.getIntegerHullCone();
    vector<vector<Integer> > supp_hyps=D.getSupportHyperplanes();
    vector<vector<Integer> > latt_points=D.getLatticePoints();
    cout << "Lattice points" << endl;
    cout << latt_points;
    cout << "Inequalities" << endl;
    cout << supp_hyps;
    vector<vector<Integer> > equations=D.getSublattice().getEquations();
    cout << "Number of equations is " << equations.size() << endl;
}
