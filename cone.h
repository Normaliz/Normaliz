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

#ifndef CONE_H_
#define CONE_H_

#include <string>
#include "libnormaliz.h"
#include "full_cone.h"
#include "sublattice_representation.h"
#include "matrix.h"

namespace libnormaliz {

template<typename Integer>
class Cone {
	int dim;
//	Full_Cone<Integer> FullDimCone;
	Sublattice_Representation<Integer> ChangeToFullDim;  //always use compose_basis_change() !
	bool BC_set, OrigGens_set;
	bitset<ConeProperty::EnumSize> is_Computed;
	list< vector<Integer> > OriginalGenerators;
	list< vector<Integer> > Generators;
	list< vector<Integer> > SupportHyperplanes;

	void compose_basis_change(const Sublattice_Representation<Integer>& SR); // composes SR

//---------------------------------------------------------------------------
//                  Progress input, depending on input_type
//---------------------------------------------------------------------------

	void prepare_input_type_0(const list< vector<Integer> >& Input);
	void prepare_input_type_1(const list< vector<Integer> >& Input);
	void prepare_input_type_2(const list< vector<Integer> >& Input);
	void prepare_input_type_3(const list< vector<Integer> >& Input);
	void prepare_input_type_10(const list< vector<Integer> >& Binomials);
//TODO reihenfolge anpassen
	void prepare_input_type_456(const list< vector<Integer> >& Congruences, const list< vector<Integer> >& Equations, const list< vector<Integer> >& Inequalities);
	void prepare_input_type_45(const Matrix<Integer>& Equations, const Matrix<Integer>& Inequalities);

	/* only used by the constructors */
	void initialize();
//---------------------------------------------------------------------------
//                          public methods
//---------------------------------------------------------------------------
public:
	/* Constructors, they preprocess the input */
	Cone(const list< vector<Integer> >& input, int type);
	Cone(const list< vector<Integer> >& Inequalities, const list< vector<Integer> >& Equations, const list< vector<Integer> >& Congruences);

	/* do computation */
	void compute(const string& mode);

	/* check what is computed */
	bool isComputed(ConeProperty::Enum prop) const;

	/* getter */
	list< vector<Integer> > const& getExtremeRays() const;
	list< vector<Integer> > const& getSupportHyperplanes() const;
	list< vector<Integer> > const& getTriangulation() const;
	list< vector<Integer> > const& getHilbertBasis() const;
	list< vector<Integer> > const& getHt1Elements() const;
	list< vector<Integer> > const& getHVector() const;
	list< vector<Integer> > const& getHilbertPolynomial() const;
	      vector<Integer>   const& getLinearFunction() const;
	bool isPointed() const;
	bool isHt1Generated() const;
	bool isHt1ExtremeRays() const;
	bool isHt1HilbertBasis() const;
	bool isIntegrallyClosed() const;
};

}  //end namespace libnormaliz

#endif /* CONE_H_ */
