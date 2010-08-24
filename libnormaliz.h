/*
 * libnormaliz.h
 *
 *  Created on: Aug 23, 2010
 *      Author: csoeger
 */

#ifndef LIBNORMALIZ_H_
#define LIBNORMALIZ_H_

#include <iostream>


extern bool verbose;
extern bool test_arithmetic_overflow;
extern ostream verbose_ostream;
extern ostream error_ostream;
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

//#include "integer.h"
//#include "cone.h"

#endif /* LIBNORMALIZ_H_ */
