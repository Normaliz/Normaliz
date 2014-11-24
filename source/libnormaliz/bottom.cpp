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
#include <vector>
#include <list>

#include <sstream>
#include <algorithm>

#include "Normaliz.h"
#include "libnormaliz.h"
#include "cone.h"
//#include "libnormaliz.cpp"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"  //TODO needed?
#include "scip/cons_linear.h"
//#include "objscip/objscip.h"


namespace libnormaliz {
using namespace std;

void bottom_points(list< vector<Integer> >& new_points, Matrix<Integer> gens) {

    vector<Integer>  grading = gens.find_linear_form();
    Integer volume;
    vector<Integer> diagonal = vector< Integer >(dim);
    Matrix<Integer> Support_Hyperplanes=invert(gens, diagonal, volume);

    if (volume < 100000000) // 100 mio
        return;

    Support_Hyperplanes = Support_Hyperplanes.transpose();
    Support_Hyperplanes.make_prime();

    // call scip
    vector<Integer> new_point = opt_sol(Support_Hyperplanes, grading);

    if ( !new_point.empty()
       && v_scalar_product(grading,new_point) < v_scalar_product(grading,gens[0])) {

        new_points.push_back(new_point);
        Matrix<Integer> stellar_gens(gens);

        for (i=0; i<dim; ++i) {
            if (v_scalar_product(Support_Hyperplanes[i], new_point) != 0) {
                stellar_gens[i] = new_point;
                bottom_points(new_points, stellar_gens);
                
                stellar_gens[i] = gens[i];
            }
            
        }

    }
  return;
}

template<typename Integer> vector<Integer> opt_sol(Matrix<Integer> SuppHyp, vector<Integer> grading){
	long dim = grading.size();
	cout << "grading " << grading;
	// setup scip enviorenment
    SCIP* scip;
    SCIPcreate(& scip);
    SCIPincludeDefaultPlugins(scip);
    SCIPsetMessagehdlr(scip,NULL);
    // create variables
    SCIP_VAR** x = new SCIP_VAR*[dim];
    char name[SCIP_MAXSTRLEN];
    SCIPcreateProbBasic(scip, "extra_points");
    for (long i=0; i<dim; i++) {
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
        SCIPcreateVarBasic(scip, &x[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), grading[i], SCIP_VARTYPE_INTEGER) ;
        SCIPaddVar(scip, x[i]);
    }
	// create constraints
    // vector< vector<Integer> > SuppHyp(MyCone.getSupportHyperplanes());
    double* ineq = new double[dim];
    long nrSuppHyp = SuppHyp.size();
    for( long i = 0; i < nrSuppHyp; ++i )
    {
        SCIP_CONS* cons;
        for (long j=0; j<dim; j++)
            ineq[j] = SuppHyp[i][j];
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ineq_%d", i);
        SCIPcreateConsBasicLinear(scip, &cons, name, dim, x, ineq, 0.0, SCIPinfinity(scip));
        SCIPaddCons(scip, cons);
        SCIPreleaseCons(scip, &cons);
    }

    SCIP_CONS* cons;
    for (long j=0; j<dim; j++)
        ineq[j] = grading[j];
    SCIPcreateConsBasicLinear(scip, &cons, "non_zero", dim, x, ineq, 1.0, SCIPinfinity(scip));
    SCIPaddCons(scip, cons);
    SCIPreleaseCons(scip, &cons);
        
//   SCIPinfoMessage(scip, NULL, "Original problem:\n");
//   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
//   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
    SCIPsolve(scip);
    SCIPfreeTransform(scip);

    vector<Integer> sol_vec(dim);
	if( SCIPgetNSols(scip) > 0 )
    {
      SCIP_SOL* sol = SCIPgetBestSol(scip);
//   SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
	  
	  for (int i=0;i<dim;i++) sol_vec[i] = SCIPgetSolVal(scip,sol,x[i]);
	  
	  cout << "solution " << sol_vec;
	  
    } else {
	    return vector<Integer>();
    }
    
	for (int j=0;j<dim;j++) SCIPreleaseVar(scip, &x[j]);
	SCIPfree(&scip);
	return sol_vec; 
}

} // namespace

#endif // NMZ_SCIP
