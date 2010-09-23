/*
 * cone.h
 *
 *  Created on: Aug 24, 2010
 *      Author: csoeger
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
	Full_Cone<Integer> FullDimCone;
	Sublattice_Representation<Integer> ChangeToFullDim;
	list< vector<Integer> > OriginalGenerators;

public:
	Cone(list< vector<Integer> > input, int type);
	Cone(list< vector<Integer> > Inequalities, list< vector<Integer> > Equations, list< vector<Integer> > Congruences=list< vector<Interger> >() );

	void compute(string mode);
	bool isComputed(ConeProperty::Enum prop) const;

	//getter
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

}

#endif /* CONE_H_ */
