#include <cstdlib>
#include <vector>
#include <fstream>
#include <ctime>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

#include "libnormaliz/libnormaliz.h"

using namespace libnormaliz;

typedef long long Integer;

Cone<Integer> rand_simplex(size_t dim, long bound) {
    vector<vector<Integer> > vertices(dim + 1, vector<Integer>(dim));
    while (true) {  // an eternal loop ...
        for (size_t i = 0; i <= dim; ++i) {
            for (size_t j = 0; j < dim; ++j)
                vertices[i][j] = rand() % (bound + 1);
        }
        Cone<Integer> Simplex(Type::polytope, vertices);
        // we must check the rank and normality
        if (Simplex.getRank() == dim + 1 && Simplex.isDeg1HilbertBasis())
            return Simplex;
    }
    vector<vector<Integer> > dummy_gen(1, vector<Integer>(1, 1));  // to make the compiler happy
    return Cone<Integer>(Type::cone, dummy_gen);
}

bool exists_jump_over(Cone<Integer>& Polytope, const vector<vector<Integer> >& jump_cands) {
    vector<vector<Integer> > test_polytope = Polytope.getExtremeRays();
    test_polytope.resize(test_polytope.size() + 1);
    for (const auto& jump_cand : jump_cands) {
        test_polytope[test_polytope.size() - 1] = jump_cand;
        Cone<Integer> TestCone(Type::cone, test_polytope);
        if (TestCone.getNrDeg1Elements() != Polytope.getNrDeg1Elements() + 1)
            continue;
        if (TestCone.isDeg1HilbertBasis())
            return true;
    }
    return false;
}

vector<Integer> lattice_widths(Cone<Integer>& Polytope) {
    if (!Polytope.isDeg1ExtremeRays()) {
        cerr << "Cone in lattice_widths is not defined by lattice polytope" << endl;
        exit(1);
    }
    vector<Integer> widths(Polytope.getNrExtremeRays());
    for (size_t i = 0; i < Polytope.getNrSupportHyperplanes(); ++i) {
        widths[i] = 0;
        for (size_t j = 0; j < Polytope.getNrExtremeRays(); ++j) {
            // v_scalar_product is a useful function from vector_operations.h
            Integer test = v_scalar_product(Polytope.getSupportHyperplanes()[i], Polytope.getExtremeRays()[j]);
            if (test > widths[i])
                widths[i] = test;
        }
    }
    return widths;
}

int main(int argc, char* argv[]) {
    time_t ticks;
    srand(time(&ticks));
    cout << "Seed " << ticks << endl;  // we may want to reproduce the run

    size_t polytope_dim = 4;
    size_t cone_dim = polytope_dim + 1;
    long bound = 6;
    vector<Integer> grading(cone_dim, 0);  // at some points we need the explicit grading
    grading[polytope_dim] = 1;

    size_t nr_simplex = 0;  // for the progress report

    while (true) {
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        Cone<Integer> Candidate = rand_simplex(polytope_dim, bound);
        nr_simplex++;
        if (nr_simplex % 1000 == 0)
            cout << "simplex " << nr_simplex << endl;
        vector<vector<Integer> > supp_hyps_moved = Candidate.getSupportHyperplanes();
        for (auto& i : supp_hyps_moved)
            i[polytope_dim] += 1;
        Cone<Integer> Candidate1(Type::inequalities, supp_hyps_moved, Type::grading, to_matrix(grading));
        if (Candidate1.getNrDeg1Elements() > Candidate.getNrDeg1Elements())
            continue;  // there exists a point of height 1
        cout << "No ht 1 jump"
             << " #latt " << Candidate.getNrDeg1Elements() << endl;
        // move the hyperplanes further outward
        for (auto& i : supp_hyps_moved)
            i[polytope_dim] += polytope_dim;
        Cone<Integer> Candidate2(Type::inequalities, supp_hyps_moved, Type::grading, to_matrix(grading));
        cout << "Testing " << Candidate2.getNrDeg1Elements() << " jump candidates" << endl;
        // including the lattice points in P
        if (exists_jump_over(Candidate, Candidate2.getDeg1Elements()))
            continue;
        cout << "No ht <= 1+dim jump" << endl;
        vector<Integer> widths = lattice_widths(Candidate);
        for (size_t i = 0; i < supp_hyps_moved.size(); ++i)
            supp_hyps_moved[i][polytope_dim] += -polytope_dim + (widths[i]) * (polytope_dim - 2);
        vector<vector<mpz_class> > mpz_supp_hyps;
        convert(mpz_supp_hyps, supp_hyps_moved);
        vector<mpz_class> mpz_grading = convertTo<vector<mpz_class> >(grading);
#ifdef _OPENMP
        omp_set_num_threads(4);
#endif
        Cone<mpz_class> Candidate3(Type::inequalities, mpz_supp_hyps, Type::grading, to_matrix(mpz_grading));
        Candidate3.compute(ConeProperty::Deg1Elements, ConeProperty::Approximate);
        vector<vector<Integer> > jumps_cand;  // for conversion from mpz_class
        convert(jumps_cand, Candidate3.getDeg1Elements());
        cout << "Testing " << jumps_cand.size() << " jump candidates" << endl;
        if (exists_jump_over(Candidate, jumps_cand))
            continue;
        cout << "Maximal simplex found" << endl;
        cout << "Vertices" << endl;
        Candidate.getExtremeRaysMatrix().pretty_print(cout);
        cout << "Number of lattice points = " << Candidate.getNrDeg1Elements();
        cout << " Multiplicity = " << Candidate.getMultiplicity() << endl;
    }  // end while
}  // end main
