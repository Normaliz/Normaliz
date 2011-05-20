/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 */

#ifndef LIBNORMALIZ_H_
#define LIBNORMALIZ_H_

#include <iostream>

namespace libnormaliz {

enum InputType {
	integral_closure,
	normalization,
	polytope,
	rees_algebra,
//	constraints, //TODO what to do with a matrix of type constraints? 
	lattice_ideal
};
//it would be nice if ConstraintType is a subset of InputType

//TODO add homogeneous for first 3?
enum ConstraintType {
	hyperplanes,
	equations,
	congruences,
	inhomogeneous_hyperplanes,
	inhomogeneous_equations,
	inhomogeneous_congruences
};

//TODO entweder in namespace oder Namen eindeutig machen
//namespace Mode {
enum ComputationMode {
	supportHyperplanes,
	volumeTriangulation,
	volumeLarge,
	height1Elements,
	hilbertBasisTriangulation,
	hilbertBasisLarge,
	hilbertPolynomial,
	hilbertPolynomialLarge,
	hilbertBasisPolynomial,
	hilbertBasisPolynomialLarge,
	dual
};
//} //end namespace mode

extern bool verbose;

/* if test_arithmetic_overflow is true, many operations are also done
 * modulo overflow_test_modulus to ensure the correctness of the calculations
 */
extern bool test_arithmetic_overflow;
extern long overflow_test_modulus;

/* methods to set and use the output streams */
void setVerboseOutput(std::ostream&);
void setErrorOutput(std::ostream&);

std::ostream& verboseOutput();
std::ostream& errorOutput();


} /* end namespace libnormaliz */

#endif /* LIBNORMALIZ_H_ */
