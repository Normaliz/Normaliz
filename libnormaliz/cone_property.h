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

#ifndef CONE_PROPERTY_H_
#define CONE_PROPERTY_H_

#include <bitset>

namespace libnormaliz {

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

typedef std::bitset<ConeProperty::EnumSize> ConeProperties;
/*class ConeProperties : public virtual std::bitset<ConeProperty::EnumSize> {
	std::bitset<N>& set (ConeProperty::Enum Property) {

	}
};*/

}

#endif /* CONE_PROPERTY_H_ */
