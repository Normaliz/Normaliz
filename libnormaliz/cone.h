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

#ifndef CONE_H_
#define CONE_H_

#include <string>
#include "libnormaliz.h"
#include "cone_property.h"
#include "full_cone.h"
#include "sublattice_representation.h"
#include "matrix.h"

namespace libnormaliz {
using std::list;
using std::vector;


template<typename Integer>
class Cone {
	size_t dim;

//	Full_Cone<Integer> FullDimCone;
	Sublattice_Representation<Integer> BasisChange;  //always use compose_basis_change() !
	bool BC_set;
	ConeProperties is_Computed;
	list< vector<Integer> > GeneratorsOfToricRing;
	list< vector<Integer> > Generators;
	vector<bool> ExtremeRays;
	list< vector<Integer> > SupportHyperplanes;
	list< pair<vector<size_t>, Integer> > Triangulation;
	Integer multiplicity;
	list< vector<Integer> > HilbertBasis;
	list< vector<Integer> > Ht1Elements;
	vector<Integer> HVector;
	vector<Integer> HilbertPolynomial;
	vector<Integer> LinearForm;
	bool pointed;
	bool ht1_extreme_rays;
	bool ht1_hilbert_basis;
	bool integrally_closed;
	bool rees_primary;
	Integer ReesPrimaryMultiplicity;

	void compose_basis_change(const Sublattice_Representation<Integer>& SR); // composes SR


	/* Progress input, depending on input_type */
	void prepare_input_type_0(const list< vector<Integer> >& Input);
	void prepare_input_type_1(const list< vector<Integer> >& Input);
	void prepare_input_type_2(const list< vector<Integer> >& Input);
	void prepare_input_type_3(const list< vector<Integer> >& Input);
	void prepare_input_type_10(const list< vector<Integer> >& Binomials);
//TODO reihenfolge anpassen
	void prepare_input_type_456(const list< vector<Integer> >& Congruences, const list< vector<Integer> >& Equations, const list< vector<Integer> >& Inequalities);
	void prepare_input_type_45(const list< vector<Integer> >& Equations, const list< vector<Integer> >& Inequalities);

	/* only used by the constructors */
	void initialize();

	/* compute method for the dual_mode, used in compute(string) */
	void compute_dual();

	/* extract the data from Full_Cone, this may remove data from Full_Cone!*/
	void extract_data(Full_Cone<Integer>& FC);

//---------------------------------------------------------------------------
//                          public methods
//---------------------------------------------------------------------------
public:
	/* Constructors, they preprocess the input */
	Cone(const list< vector<Integer> >& input, int type);

	Cone(const list< vector<Integer> >& Inequalities,
	     const list< vector<Integer> >& Equations,
	     const list< vector<Integer> >& Congruences);

	/* do computation */
	void compute(ConeProperties to_compute);
	void compute(const string& mode);

	/* check what is computed */
	bool isComputed(ConeProperty::Enum prop) const;

	/* getter */
	Sublattice_Representation<Integer> const& getBasisChange() const;
	list< vector<Integer> > const& getGeneratorsOfToricRing() const;
	list< vector<Integer> > const& getGenerators() const;
	vector<bool> const& getExtremeRays() const;
	list< vector<Integer> > const& getSupportHyperplanes() const;
	list< pair<vector<size_t>, Integer> > const& getTriangulation() const;
	list< vector<Integer> > const& getHilbertBasis() const;
	list< vector<Integer> > const& getHt1Elements() const;
	vector<Integer> const& getHVector() const;
	vector<Integer> const& getHilbertPolynomial() const;
	vector<Integer> const& getLinearForm() const;
	Integer const& getMultiplicity() const;
	bool isPointed() const;
	bool isHt1ExtremeRays() const;
	bool isHt1HilbertBasis() const;
	bool isIntegrallyClosed() const;
	bool isReesPrimary() const;
	Integer const& getReesPrimaryMultiplicity() const;
};

}  //end namespace libnormaliz

#endif /* CONE_H_ */
