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

namespace libnormaliz {

bool verbose = false;

bool test_arithmetic_overflow = false;
long overflow_test_modulus = 15401;

size_t RecBoundFactor = 5000000;

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
#include "cone_property.cpp"

