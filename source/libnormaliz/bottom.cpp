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
#include <omp.h>

#include <iostream>
//#include <sstream>
#include <algorithm>
#include <queue>

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

long ScipBound = 1000000;

template<typename Integer>
vector<Integer> opt_sol(SCIP* scip, const Matrix<Integer>& gens, const Matrix<Integer>& SuppHyp, const vector<Integer>& grading);

template<typename Integer>
void bottom_points_inner(SCIP* scip, Matrix<Integer>& gens, list< vector<Integer> >& new_points, vector< Matrix<Integer> >& q_gens);

double convert_to_double(mpz_class a) {
    return a.get_d();
}

double convert_to_double(long a) {
    return a;
}

double convert_to_double(long long a) {
    return a;
}

// TODO do not use global variables
long long stellar_det_sum;

template<typename Integer>
void bottom_points(list< vector<Integer> >& new_points, Matrix<Integer> gens) { //TODO lieber Referenz fuer Matrix?
	
	cout << "Using SCIP to compute points from bottom. BOUND: " <<ScipBound<< endl;
    stellar_det_sum = 0;
    vector< Matrix<Integer> > q_gens;
    q_gens.push_back(gens);
    int level = 0;

    #pragma omp parallel reduction(+:stellar_det_sum)
    {
    // setup scip enviorenment
    SCIP* scip;
    SCIPcreate(& scip);
    SCIPincludeDefaultPlugins(scip);
//    SCIPsetMessagehdlr(scip,NULL);  // deactivate scip output
    
	SCIPsetIntParam(scip, "display/verblevel", 0); 
	
	// modify timing for better parallelization
//	SCIPsetBoolParam(scip, "timing/enabled", FALSE);
	SCIPsetBoolParam(scip, "timing/statistictiming", FALSE);
	SCIPsetBoolParam(scip, "timing/rareclockcheck", TRUE);


	SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", -1); 
	SCIPsetIntParam(scip, "branching/pscost/priority", 1000000); 
//	SCIPsetIntParam(scip, "nodeselection/uct/stdpriority", 1000000); 


    vector< Matrix<Integer> > local_q_gens;
    list< vector<Integer> > local_new_points;

    while (!q_gens.empty()) {
        #pragma omp single
        cout << q_gens.size() << " simplices on level " << level++ << endl;
        #pragma omp for schedule(static)
        for (size_t i = 0; i < q_gens.size(); ++i) {
            bottom_points_inner(scip, q_gens[i], local_new_points, local_q_gens);
        }
        #pragma omp single
        {
			q_gens.clear();
		}
        #pragma omp critical
        {
            q_gens.insert(q_gens.end(),local_q_gens.begin(),local_q_gens.end());
        }
        local_q_gens.clear();
        #pragma omp barrier
    }
    #pragma omp critical
    {
        new_points.splice(new_points.end(), local_new_points, local_new_points.begin(), local_new_points.end());
    }

    SCIPfree(& scip);
    } // end parallel

    cout  << new_points.size() << " new points accumulated" << endl;
    new_points.sort();
    new_points.unique();
    cout << "of which " << new_points.size() << " are unique" << endl;

    cout << "stellar_det_sum = " << stellar_det_sum << endl;

}


template<typename Integer>
void bottom_points_inner(SCIP* scip, Matrix<Integer>& gens, list< vector<Integer> >& local_new_points,
                         vector< Matrix<Integer> >& local_q_gens) {

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
    // set time limit according to volume   
    double time_limit = pow(log10(convert_to_double(volume)),2);
    SCIPsetRealParam(scip, "limits/time", time_limit);
    // call scip
    vector<Integer> new_point = opt_sol(scip, gens, Support_Hyperplanes, grading);
    if ( !new_point.empty() ){

        //if (find(local_new_points.begin(), local_new_points.end(),new_point) == local_new_points.end())
        local_new_points.push_back(new_point);
        Matrix<Integer> stellar_gens(gens);

        int nr_hyps = 0;
        for (int i=0; i<dim; ++i) {
            if (v_scalar_product(Support_Hyperplanes[i], new_point) != 0) {
                stellar_gens[i] = new_point;
                local_q_gens.push_back(stellar_gens);

                stellar_gens[i] = gens[i];
            } else nr_hyps++;

        }
        //#pragma omp critical(VERBOSE)
       //cout << new_point << " liegt in " << nr_hyps <<" hyperebenen" << endl;

    }
    else {
//      cout << "Not using " << new_point;
        stellar_det_sum += explicit_cast_to_long(volume);
    }
    return;
}

template<typename Integer>
double max_in_col(const Matrix<Integer>& M, size_t j) {
	// warum setzen wir das auf 0???
    Integer max = 0;
    //Integer max = M[1][j];
    for (size_t i=0; i<M.nr_of_rows(); ++i) {
        if (M[i][j] > max) max = M[i][j];
    }
    return convert_to_double(max);
}

template<typename Integer>
double max_in_col_no_zero(const Matrix<Integer>& M, size_t j) {
	// warum setzen wir das auf 0???
    //Integer max = 0;
    Integer max = M[1][j];
    for (size_t i=0; i<M.nr_of_rows(); ++i) {
        if (M[i][j] > max) max = M[i][j];
    }
    return convert_to_double(max);
}

template<typename Integer>
double min_in_col(const Matrix<Integer>& M, size_t j) {
	// warum setzen wir das auf 0???
    Integer min = 0;
    //Integer min = M[1][j];
    for (size_t i=0; i<M.nr_of_rows(); ++i) {
        if (M[i][j] < min) min = M[i][j];
    }
    return convert_to_double(min);
}

template<typename Integer>
double min_in_col_no_zero(const Matrix<Integer>& M, size_t j) {
	// warum setzen wir das auf 0???
    //Integer min = 0;
    Integer min = M[1][j];
    for (size_t i=0; i<M.nr_of_rows(); ++i) {
        if (M[i][j] < min) min = M[i][j];
    }
    return convert_to_double(min);
}

template<typename Integer>
vector<Integer> opt_sol(SCIP* scip,
                        const Matrix<Integer>& gens, const Matrix<Integer>& SuppHyp,
                        const vector<Integer>& grading) {
    double SCIPFactor =1.0;
    double upper_bound = convert_to_double(v_scalar_product(grading,gens[0]))-1;
    upper_bound*=SCIPFactor;
    // TODO make the test more strict
    long dim = grading.size();
    // create variables
    SCIP_VAR** x = new SCIP_VAR*[dim];
    char name[SCIP_MAXSTRLEN];
    SCIPcreateProbBasic(scip, "extra_points");
    for (long i=0; i<dim; i++) {
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
//        SCIPcreateVarBasic(scip, &x[i], name, -SCIPinfinity(scip), SCIPinfinity(scip),
//                           convert_to_double(grading[i]), SCIP_VARTYPE_INTEGER);
        SCIPcreateVarBasic(scip, &x[i], name, min_in_col(gens,i), max_in_col(gens, i),
                             convert_to_double(grading[i]), SCIP_VARTYPE_INTEGER);
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
    // setup non-zero constraints
    // if all extreme rays have the same sign in one dimension, add the x_i>=1 or x_i<=-1 constraint
    
    for (long i=0; i<dim; i++){
		double min = min_in_col_no_zero(gens,i);
		double max = max_in_col_no_zero(gens,i);
		if (min*max>0){
			cout << "same sign in variable " << i << endl;
			for (int j=0;j<dim;j++) ineq[j]=0;
			ineq[i]=1;
			if (min>0){
				SCIPcreateConsBasicLinear(scip, &cons, "non_zero", dim, x, ineq, 1.0, SCIPinfinity(scip));
			}
			else {

				SCIPcreateConsBasicLinear(scip, &cons, "non_zero", dim, x, ineq, -SCIPinfinity(scip), -1.0);
			}
			break;
		}
		if (i==dim-1){
			cout << "no same sign. using bound disjunction" << endl;
			// set bound disjunction
			
			SCIP_VAR** double_x = new SCIP_VAR*[2*dim];
			SCIP_BOUNDTYPE* boundtypes = new SCIP_BOUNDTYPE[2*dim];
			SCIP_Real* bounds = new SCIP_Real[2*dim];
			for (long i=0; i<dim;i++) {
				double_x[2*i] = x[i];
				double_x[2*i+1] = x[i];
				boundtypes[2*i]= SCIP_BOUNDTYPE_LOWER;
				boundtypes[2*i+1] = SCIP_BOUNDTYPE_UPPER;
				bounds[2*i] = 1.0;
				bounds[2*i+1] = -1.0;
			}
			SCIPcreateConsBasicBounddisjunction	(scip, &cons,"non_zero",2*dim,double_x,boundtypes,bounds);		
		
		/*
		// this type of constraints procedures numerical problems:
		for (long j=0; j<dim; j++)
        ineq[j] = convert_to_double(grading[j]);
		SCIPcreateConsBasicLinear(scip, &cons, "non_zero", dim, x, ineq, 1.0, SCIPinfinity(scip));
		*/
		}
	}
	
	

	// set objective limit. feasible solution has to have at least this objective value
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

//#ifndef NDEBUG_BLA 
        //FILE* file = fopen("mostrecent.lp","w");
        //assert (file != NULL);
        //SCIPprintOrigProblem(scip, file, "lp", FALSE);
        //SCIPwriteParams(scip, "mostrecent.set", TRUE, TRUE);
        //fclose(file);
//#endif

	// set numerics
	Integer maxabs = v_max_abs(grading);
	double epsilon = max(1e-20,min(1/(convert_to_double(maxabs)*10),1e-10));
	//cout << "epsilon is in region " << log10(epsilon) << endl;
	double feastol = max(1e-17,epsilon*10);
	SCIPsetRealParam(scip, "numerics/epsilon", epsilon); 
	SCIPsetRealParam(scip, "numerics/feastol", feastol); 

    SCIPsolve(scip);
    //SCIPprintStatistics(scip, NULL);
    vector<Integer> sol_vec(dim);
	if(SCIPgetStatus(scip) == SCIP_STATUS_TIMELIMIT) cout << "time limit reached!" << endl;
    if( SCIPgetNLimSolsFound(scip) > 0 ) // solutions respecting objective limit (ie not our input solutions)
    {
        SCIP_SOL* sol = SCIPgetBestSol(scip);
        //SCIPprintOrigProblem(scip, NULL, NULL, FALSE);
        //SCIPprintSol(scip, sol, NULL, FALSE) ;

        for (int i=0;i<dim;i++) {
            sol_vec[i] = explicit_cast_to_long(SCIPconvertRealToLongint(scip,SCIPgetSolVal(scip,sol,x[i])));
        }

        // HOTFIX to avoid pseudo solution
        // if(v_scalar_product(grading,sol_vec)>upper_bound) return vector<Integer>();
        //for (int i=0;i<nrSuppHyp;i++){
        //    if((v_scalar_product(SuppHyp[i],sol_vec))<0) return vector<Integer>();
        //    }
        // if((v_scalar_product(grading,sol_vec))<1) return vector<Integer>();

        if(v_scalar_product(grading,sol_vec)>upper_bound){
			Integer sc = v_scalar_product(sol_vec,grading);
				#pragma omp critical(VERBOSE)
				{
					cout << "solution does not respect upper bound! here's the data:" << endl;
					cout << "upper bound: " << upper_bound << endl;
					cout << "grading: " << grading;
					cout << "hyperplanes:" << endl;
					SuppHyp.pretty_print(cout);
					cout << "generators:" << endl;
					gens.pretty_print(cout);
					cout << sc << " | solution " << sol_vec;
					cout << "epsilon: " << epsilon << endl;
					SCIPprintOrigProblem(scip, NULL, NULL, FALSE);
					SCIPprintSol(scip, sol, NULL, FALSE) ;
					cout << "write files... " << endl;
					FILE* file = fopen("mostrecent.lp","w");
					assert (file != NULL);
					SCIPprintOrigProblem(scip, file, "lp", FALSE);
					SCIPwriteParams(scip, "mostrecent.set", TRUE, TRUE);
					fclose(file);
					assert(v_scalar_product(grading,sol_vec)<=upper_bound);
			
				}
			return vector<Integer>();
			}
			
		
        for (int i=0;i<nrSuppHyp;i++) {
			if((v_scalar_product(SuppHyp[i],sol_vec))<0) {
				Integer sc = v_scalar_product(sol_vec,grading);
				#pragma omp critical(VERBOSE)
				{
					cout << "solution does not respect hyperplanes! here's the data:" << endl;
					cout << "the hyperplane: " << SuppHyp[i];
					cout << "grading: " << grading;
					cout << "hyperplanes:" << endl;
					SuppHyp.pretty_print(cout);
					cout << "generators:" << endl;
					gens.pretty_print(cout);
					cout << sc << " | solution " << sol_vec;
					cout << "epsilon: " << epsilon << endl;
					SCIPprintOrigProblem(scip, NULL, NULL, FALSE);
					SCIPprintSol(scip, sol, NULL, FALSE) ;
					cout << "write files... " << endl;
					FILE* file = fopen("mostrecent.lp","w");
					assert (file != NULL);
					SCIPprintOrigProblem(scip, file, "lp", FALSE);
					SCIPwriteParams(scip, "mostrecent.set", TRUE, TRUE);
					fclose(file);
					assert((v_scalar_product(SuppHyp[i],sol_vec))>=0);
			
				}
			return vector<Integer>();
			}
		}
        if((v_scalar_product(grading,sol_vec))<1) {
		Integer sc = v_scalar_product(sol_vec,grading);
		#pragma omp critical(VERBOSE)
		{
			cout << "the solution does not respect the nonzero condition! here's the data:" << endl;
			cout << "grading: " << grading;
			cout << "hyperplanes:" << endl;
			SuppHyp.pretty_print(cout);
			cout << "generators:" << endl;
			gens.pretty_print(cout);
			cout << sc << " | solution " << sol_vec;
			cout << "epsilon: " << epsilon << endl;
			SCIPprintOrigProblem(scip, NULL, NULL, FALSE);
			SCIPprintSol(scip, sol, NULL, FALSE) ;
			cout << "write files... " << endl;
			FILE* file = fopen("mostrecent.lp","w");
			assert (file != NULL);
			SCIPprintOrigProblem(scip, file, "lp", FALSE);
			SCIPwriteParams(scip, "mostrecent.set", TRUE, TRUE);
			fclose(file);
			assert((v_scalar_product(grading,sol_vec))>=1);
			
		}
			return vector<Integer>();
        }
        assert(v_scalar_product(grading,sol_vec)<=upper_bound);
        for (int i=0;i<nrSuppHyp;i++) assert((v_scalar_product(SuppHyp[i],sol_vec))>=0);
        assert((v_scalar_product(grading,sol_vec))>=1);
        //Integer sc = v_scalar_product(sol_vec,grading);
		//#pragma omp critical(VERBOSE)
		//cout << sc << " | solution " << sol_vec;

    } else {
        return vector<Integer>();
    }
    
    for (int j=0;j<dim;j++) SCIPreleaseVar(scip, &x[j]);
    SCIPfreeProb(scip);
    return sol_vec; 
}

template void bottom_points(list< vector<long> >& new_points, Matrix<long> gens);
template void bottom_points(list< vector<long long> >& new_points, Matrix<long long> gense);
template void bottom_points(list< vector<mpz_class> >& new_points, Matrix<mpz_class> gens);

} // namespace

#endif // NMZ_SCIP
