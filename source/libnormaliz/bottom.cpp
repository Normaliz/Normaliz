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

#ifdef NMZ_SCIP
#include <stdlib.h>
#include <math.h>

#include <iostream>
//#include <sstream>
#include <algorithm>

#include "bottom.h"
#include "libnormaliz.h"
#include "vector_operations.h"
#include "integer.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"  //TODO needed?
#include "scip/cons_linear.h"
//#include "objscip/objscip.h"


namespace libnormaliz {
using namespace std;

const long ScipBound = 1000000;

template<typename Integer>
vector<Integer> opt_sol(SCIP* scip, const Matrix<Integer>& gens, const Matrix<Integer>& SuppHyp, const vector<Integer>& grading);

template<typename Integer>
void bottom_points_inner(SCIP* scip, list< vector<Integer> >& new_points, Matrix<Integer> gens);


// don't do it with mpz_class TODO
//template<> void bottom_points(list< vector<mpz_class> >& new_points, Matrix<mpz_class> gens) {}

long long stellar_det_sum;

template<typename Integer>
void bottom_points(list< vector<Integer> >& new_points, Matrix<Integer> gens) { //TODO lieber Referenz fuer Matrix?

    // setup scip enviorenment
    SCIP* scip;
    SCIPcreate(& scip);
    SCIPincludeDefaultPlugins(scip);
//    SCIPsetMessagehdlr(scip,NULL);  // deactivate scip output
    
	SCIPsetIntParam(scip, "display/verblevel", 0); 

	SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", -1); 
	SCIPsetIntParam(scip, "branching/pscost/priority", 1000000); 
//	SCIPsetIntParam(scip, "nodeselection/uct/stdpriority", 1000000); 

	SCIPsetRealParam(scip, "numerics/epsilon", 1e-10); 
	SCIPsetRealParam(scip, "numerics/feastol", 1e-9); 

    stellar_det_sum = 0;
    bottom_points_inner(scip, new_points, gens);

    cout << "stellar_det_sum = " << stellar_det_sum << endl;

    SCIPfree(& scip);

    

}


template<typename Integer>
void bottom_points_inner(SCIP* scip, list< vector<Integer> >& new_points, Matrix<Integer> gens) { //TODO lieber Referenz fuer Matrix?

    vector<Integer>  grading = gens.find_linear_form();
    Integer volume;
    int dim = gens[0].size();
    Matrix<Integer> Support_Hyperplanes = gens.invert(volume);

    if (volume < ScipBound) {
        stellar_det_sum += explicit_cast_to_long(volume);
        return;
    }

    Support_Hyperplanes = Support_Hyperplanes.transpose();
    Support_Hyperplanes.make_prime();
    /*
       cout << "gens " << endl;
       gens.pretty_print(cout);
       cout << "supp hyps " << endl;
       Support_Hyperplanes.pretty_print(cout);
       cout << "grading " << grading;
       */
    // call scip
    vector<Integer> new_point = opt_sol(scip, gens, Support_Hyperplanes, grading);
    if ( !new_point.empty() ){

        if (find(new_points.begin(), new_points.end(),new_point) == new_points.end())
            new_points.push_back(new_point);
        Matrix<Integer> stellar_gens(gens);

        int nr_hyps = 0;
        for (int i=0; i<dim; ++i) {
            if (v_scalar_product(Support_Hyperplanes[i], new_point) != 0) {
                stellar_gens[i] = new_point;
                bottom_points_inner(scip, new_points, stellar_gens);

                stellar_gens[i] = gens[i];
            } else nr_hyps++;

        }
//       cout << new_point << " liegt in " << nr_hyps << endl;

    }
    else {
//      cout << "Not using " << new_point;
        stellar_det_sum += explicit_cast_to_long(volume);
    }
    return;
}


template<typename Integer>
Integer max_in_col(const Matrix<Integer>& M, size_t j) {
    Integer max = 0;
    for (size_t i=0; i<M.nr_of_rows(); ++i) {
        if (M[i][j] > max) max = M[i][j];
    }
    return max;
}

template<typename Integer>
Integer min_in_col(const Matrix<Integer>& M, size_t j) {
    Integer min = 0;
    for (size_t i=0; i<M.nr_of_rows(); ++i) {
        if (M[i][j] < min) min = M[i][j];
    }
    return min;
}

double convert_to_double(mpz_class a) {
    return a.get_d();
}

double convert_to_double(long a) {
    return a;
}

double convert_to_double(long long a) {
    return a;
}

template<typename Integer>
vector<Integer> opt_sol(SCIP* scip,
                        const Matrix<Integer>& gens, const Matrix<Integer>& SuppHyp,
                        const vector<Integer>& grading) {
    //    double SCIP_scalar =1.0;
    double upper_bound = convert_to_double(v_scalar_product(grading,gens[0]) -1);
    // TODO make the test more strict
    long dim = grading.size();
    // create variables
    SCIP_VAR** x = new SCIP_VAR*[dim];
    char name[SCIP_MAXSTRLEN];
    SCIPcreateProbBasic(scip, "extra_points");
    for (long i=0; i<dim; i++) {
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
        SCIPcreateVarBasic(scip, &x[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), convert_to_double(grading[i]), SCIP_VARTYPE_INTEGER);
        //SCIPcreateVarBasic(scip, &x[i], name, min_in_col(gens,i), max_in_col(gens, i), grading[i], SCIP_VARTYPE_INTEGER);
        SCIPaddVar(scip, x[i]);
    }

    // create constraints
    // vector< vector<Integer> > SuppHyp(MyCone.getSupportHyperplanes());
    double* ineq = new double[dim];
    long nrSuppHyp = SuppHyp.nr_of_rows();
    for( long i = 0; i < nrSuppHyp; ++i )
    {
        SCIP_CONS* cons;
        for (long j=0; j<dim; j++)
            ineq[j] = convert_to_double(SuppHyp[i][j]);
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ineq_%d", i);
        SCIPcreateConsBasicLinear(scip, &cons, name, dim, x, ineq, 0.0, SCIPinfinity(scip));
        SCIPaddCons(scip, cons);
        SCIPreleaseCons(scip, &cons);
    }

    SCIP_CONS* cons;
    for (long j=0; j<dim; j++)
        ineq[j] = convert_to_double(grading[j]);
//    SCIPcreateConsBasicLinear(scip, &cons, "non_zero", dim, x, ineq, 1.0, upper_bound);
    SCIPcreateConsBasicLinear(scip, &cons, "non_zero", dim, x, ineq, 1.0, SCIPinfinity(scip));
    SCIPsetObjlimit(scip,upper_bound);
    SCIPaddCons(scip, cons);
    SCIPreleaseCons(scip, &cons);

    // give original generators as hints to scip
    SCIP_SOL* input_sol;
    SCIP_Bool stored;
    SCIPcreateOrigSol(scip, &input_sol, NULL);
    for (long i=0; i<dim; i++) {
        for (long j=0; j<dim; j++) {
            SCIPsetSolVal(scip, input_sol, x[j], convert_to_double(gens[i][j]));
        }
        //SCIPprintSol(scip, input_sol, NULL, TRUE);
        SCIPaddSol(scip, input_sol, &stored);
    }
    SCIPfreeSol(scip, &input_sol);
    
//    SCIPinfoMessage(scip, NULL, "Original problem:\n");
//    SCIPprintOrigProblem(scip, NULL, NULL, FALSE);
//    SCIPinfoMessage(scip, NULL, "\nSolving...\n");
#ifndef NDEBUG 
        FILE* file = fopen("mostrecent.lp","w");
        assert (file != NULL);
        SCIPprintOrigProblem(scip, file, "lp", FALSE);
        SCIPwriteParams(scip, "mostrecent.set", TRUE, TRUE);
        fclose(file);
#endif

    SCIPsolve(scip);

    //    SCIPprintStatistics(scip, NULL);
    vector<Integer> sol_vec(dim);
    if( SCIPgetNLimSolsFound(scip) > 0 ) // solutions respecting objective limit (ie not our input solutions)
    {
        SCIP_SOL* sol = SCIPgetBestSol(scip);
//        SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) ;

        for (int i=0;i<dim;i++) {
            sol_vec[i] = explicit_cast_to_long(SCIPconvertRealToLongint(scip,SCIPgetSolVal(scip,sol,x[i])));
        }
        Integer sc = v_scalar_product(sol_vec,grading);
cout << "solution " << sol_vec << " | " << sc << endl;
#ifndef NDEBUG 
        FILE* file = fopen("mostrecent.lp","w");
        assert (file != NULL);
        SCIPprintOrigProblem(scip, file, "lp", FALSE);
        SCIPwriteParams(scip, "mostrecent.set", TRUE, TRUE);
        fclose(file);
#endif

    } else {
        //		cout << "ich komme nieeeeeeeeemals vor!!!" << endl;
        //		exit(5); //TODO kann jetzt doch vorkommen
        return vector<Integer>();
    }
    
    for (int j=0;j<dim;j++) SCIPreleaseVar(scip, &x[j]);
    SCIPfreeProb(scip);
    return sol_vec; 
}

template void bottom_points(list< vector<long> >& new_points, Matrix<long> gens);
template void bottom_points(list< vector<long long> >& new_points, Matrix<long long> gens);
template void bottom_points(list< vector<mpz_class> >& new_points, Matrix<mpz_class> gens);

} // namespace

#endif // NMZ_SCIP
