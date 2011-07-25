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

#include "libnormaliz.h"
#include "general.h"

namespace libnormaliz {

bool verbose = false;

bool test_arithmetic_overflow = false;
long overflow_test_modulus = 15401;

namespace {
	std::ostream* verbose_ostream_ptr = &std::cout;
	std::ostream* error_ostream_ptr = &std::cerr;
} // end anonymous namespace, only accessible in this file (and when it is included)

void setVerboseOutput(std::ostream& v_out) {
	verbose_ostream_ptr = &v_out;
}

void setErrorOutput(std::ostream& e_out) {
	error_ostream_ptr = &e_out;
}

std::ostream& verboseOutput() {
	return *verbose_ostream_ptr;
}

std::ostream& errorOutput() {
	return *error_ostream_ptr;
}

} /* end namespace libnormaliz */

#include "integer.cpp"
#include "vector_operations.cpp"
#include "matrix.cpp"
#include "simplex.cpp"
#include "list_operations.cpp"
#include "lineare_transformation.cpp"
#include "sublattice_representation.cpp"
#include "full_cone.cpp"
#include "cone_dual_mode.cpp"
#include "cone.cpp"


namespace libnormaliz {

template class Cone<long>;
template class Cone<long long int>;
template class Cone<mpz_class>;

template class Matrix<long>;
template class Matrix<long long int>;
template class Matrix<mpz_class>;

template class Sublattice_Representation<long>;
template class Sublattice_Representation<long long int>;
template class Sublattice_Representation<mpz_class>;

template class Lineare_Transformation<long>;
template class Lineare_Transformation<long long int>;
template class Lineare_Transformation<mpz_class>;

template Lineare_Transformation<long> Transformation(const Matrix<long>& M);
template Lineare_Transformation<long long int> Transformation(const Matrix<long long int>& M);
template Lineare_Transformation<mpz_class> Transformation(const Matrix<mpz_class>& M);

template size_t decimal_length<long>(long);
template size_t decimal_length<long long int>(long long int);
template size_t decimal_length<mpz_class>(mpz_class);

}
