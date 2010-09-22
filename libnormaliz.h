/*
 * libnormaliz.h
 *
 *  Created on: Aug 23, 2010
 *      Author: csoeger
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

using namespace std;

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
 * The namespace prevents interferring with other names.
 */
namespace ConeProperty {
	enum Enum {
		ExtremeRays,
		SupportHyperplanes,
		Triangulation,
		HilbertBasis,
		Ht1Elements,
		HVector,
		HilbertPolynomial,
		IsPointed,
		IsHt1Generated,
		IsHt1ExtremeRays,
		IsHt1HilbertBasis,
		IsIntegrallyClosed,
		EnumSize // this has to be the last entry, to get the number of entries in the enum
	};
}

#endif /* LIBNORMALIZ_H_ */
