/*
 * libnormaliz.cpp
 *
 *  Created on: Aug 24, 2010
 *      Author: csoeger
 */

#include "libnormaliz.h"
#include <stdlib.h> //for exit()

namespace libnormaliz {

bool verbose = false;

bool test_arithmetic_overflow = false;
int overflow_test_modulus = 15401;

std::ostream& verbose_ostream = std::cout;
std::ostream& error_ostream = std::cerr;

// this function determinates if and how the program will be terminated in case of errors
void global_error_handling(){
	cout<<"Some error detected. The program will be terminated."<<endl;
	throw 1;
	exit(1);
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
#include "output.cpp"
#include "mode.cpp"

namespace libnormaliz {

template class Matrix<long long>;
template class Matrix<mpz_class>;
template class Output<long long>;
template class Output<mpz_class>;

template void make_main_computation<long long>(const int&, string&, const Matrix<long long>&, Output<long long>&);
template void run_mode_456<long long>(string&, const Matrix<long long>&, Matrix<long long>, Matrix<long long>, Output<long long>&);
template void make_main_computation<mpz_class>(const int&, string&, const Matrix<mpz_class>&, Output<mpz_class>&);
template void run_mode_456<mpz_class>(string&, const Matrix<mpz_class>&, Matrix<mpz_class>, Matrix<mpz_class>, Output<mpz_class>&);

}
