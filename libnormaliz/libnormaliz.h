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

#ifndef LIBNORMALIZ_H_
#define LIBNORMALIZ_H_

#include <iostream>
#include <assert.h>

#ifdef _WIN32 //for 32 and 64 bit windows
	#include <mpirxx.h>
#else
	#include <gmpxx.h>
#endif


namespace libnormaliz {

extern bool verbose;

/* if test_arithmetic_overflow is true, many operations are also done
 * modulo overflow_test_modulus to ensure the correctness of the calculations
 */
extern bool test_arithmetic_overflow;
extern int overflow_test_modulus;

extern std::ostream& verbose_ostream;
extern std::ostream& error_ostream;

//extern const unsigned int major_version, minor_version;  //TODO version


/* An enumeration of things, that can be computed for a cone.
 * The namespace prevents interfering with other names.
 */
namespace ConeProperty {
	enum Enum {
		Generators,
		ExtremeRays,
		SupportHyperplanes,
		Triangulation,
		Multiplicity,
		HilbertBasis,
		Ht1Elements,
		HVector,
		HilbertPolynomial,
		LinearForm,
		IsPointed,
		IsHt1Generated,
		IsHt1ExtremeRays,
		IsHt1HilbertBasis,
		IsIntegrallyClosed,
		EnumSize //this has to be the last entry, to get the number of entries in the enum
	};
}

}

#endif /* LIBNORMALIZ_H_ */
