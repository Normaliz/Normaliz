/*
 * Normaliz 2.5
 * Copyright (C) 2007-2010  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

namespace libnormaliz {

bool verbose = false;

bool test_arithmetic_overflow = false;
int overflow_test_modulus = 15401;

std::ostream& verbose_ostream = std::cout;
std::ostream& error_ostream = std::cerr;

// this function determinates if and how the program will be terminated in case of errors
void global_error_handling(){
	cout<<"Some error detected. The program will throw exception."<<endl;
	throw 1;
}

}

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

template class Cone<long long>;
template class Cone<mpz_class>;

template class Matrix<long long>;
template class Matrix<mpz_class>;

template class Sublattice_Representation<long long>;
template class Sublattice_Representation<mpz_class>;

template class Lineare_Transformation<long long>;
template class Lineare_Transformation<mpz_class>;

template int decimal_length<long long>(long long);
template int decimal_length<mpz_class>(mpz_class);

}
