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

#ifndef CONE_PROPERTY_H_
#define CONE_PROPERTY_H_

#include <bitset>
#include <iostream>

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
		GeneratorsOfToricRing,
		ReesPrimary,
		EnumSize //this has to be the last entry, to get the number of entries in the enum
	};
}

template<typename Integer> class Cone;
template<typename Integer> class Full_Cone;

class ConeProperties {
public:
	/* Constructors */
	ConeProperties();
	ConeProperties(ConeProperty::Enum);
	ConeProperties(ConeProperty::Enum, ConeProperty::Enum);
	ConeProperties(const std::bitset<ConeProperty::EnumSize>&);

	/* set properties */
	ConeProperties& set(ConeProperty::Enum, bool value=true);
	ConeProperties& set(ConeProperty::Enum, ConeProperty::Enum);
	ConeProperties& set(const ConeProperties&);

	/* reset (=unset) properties */
	ConeProperties& reset(ConeProperty::Enum Property);
	ConeProperties& reset(const ConeProperties&);

	/* test which/how many properties are set */
	bool test(ConeProperty::Enum Property) const;
	bool any() const;
	bool none() const;
	size_t count () const;


	/* print it in a nice way */
	void print(std::ostream& out);


private:
	std::bitset<ConeProperty::EnumSize> CPs;

	template<typename Integer> friend class Cone;
	template<typename Integer> friend class Full_Cone;
};

}

#endif /* CONE_PROPERTY_H_ */
