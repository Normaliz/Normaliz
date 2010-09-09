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


template<typename Integer>
class Cone {
	Full_Cone<Integer> FullDimCone;
	Sublattice_Representation<Integer> ChangeToFullDim;
	list< vector<Integer> > OriginalGenerators;

public:
	Cone(list< vector<Integer> > input1, int type1);
	Cone(list< vector<Integer> > input1, int type1, list< vector<Integer> > input2, int type2, list< vector<Integer> > input3, int type3);

	void compute(string mode);
	bool isComputed(ConeProperty::Enum prop) const;

	//getter
	list< vector<Integer> > getExtremeRays() const;
	list< vector<Integer> > getSupportHyperplanes() const;
	list< vector<Integer> > getTriangulation() const;
	list< vector<Integer> > getHilbertBasis() const;
	list< vector<Integer> > getHt1Elements() const;
	list< vector<Integer> > getHVector() const;
	list< vector<Integer> > getHilbertPolynomial() const;
	bool isPointed() const;
	bool isHt1Generated() const;
	bool isHt1ExtremeRays() const;
	bool isHt1HilbertBasis() const;
	bool isIntegrallyClosed() const;
};


#endif /* CONE_H_ */
