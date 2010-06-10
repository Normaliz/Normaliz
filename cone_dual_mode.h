/*
* Normaliz 2.2
* Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim
* With contributions by Christof Soeger
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

#ifndef CONE_DUAL_MODE_H
#define CONE_DUAL_MODE_H

#include <set>
#include <list>
#include "integer.h"
#include "matrix.h"
#include "sublattice_representation.h"

class Cone_Dual_Mode {
public:
	int dim;
	int nr_sh;
	int hyp_size;
	
	Matrix SupportHyperplanes;
	list<vector<Integer> > Generators;
	list<vector<Integer> > Hilbert_Basis;

/* ---------------------------------------------------------------------------
 *				Private routines, used in the public routines
 * ---------------------------------------------------------------------------
 */

	/* Returns true if new_element is reducible versus the elements in Ired
	 * used for dual algorithm */
	bool reduce(list<vector<Integer> > & Ired, const vector<Integer> & new_element, const int & size);
	bool reduce(list<vector<Integer> *> & Ired, const vector<Integer> & new_element, const int & size);

	/* reduce Red versus Ired */
	void reduce(list<vector<Integer> > & Ired, list<vector<Integer> > & Red, const int & size);

	/* adds a new element irreducible to the Hilbert basis
	 * the new elements must come from a  structure sorted by total degree
	 * used for dual algorithm */
	void reduce_and_insert(const vector<Integer> & new_element, const int & size);
	/* select extreme rays  by reduction
	 * used for the dual algorithm */
	void reduce_and_insert_extreme(const vector<Integer> & new_element);

	/* computes the Hilbert basis after adding a support hyperplane with the dual algorithm */
	void cut_with_halfspace_hilbert_basis(const int & hyp_counter, const bool & lifting, vector<Integer> & halfspace);
	/* computes the Hilbert basis after adding a support hyperplane with the dual algorithm , general case */
	Matrix cut_with_halfspace(const int & hyp_counter, const Matrix & Basis_Max_Subspace);

	/* computes the extreme rays using reduction, used for the dual algorithm */
	void extreme_rays_reduction();
	/* computes the extreme rays using rank test, used for the dual algorithm */
	void extreme_rays_rank();

	void relevant_support_hyperplanes();

	Cone_Dual_Mode();
	Cone_Dual_Mode(Matrix M);            //main constructor
	Cone_Dual_Mode(const Cone_Dual_Mode & C); //copy constructor
	~Cone_Dual_Mode();                   //destructor

/*---------------------------------------------------------------------------
 *						Data access
 *---------------------------------------------------------------------------
 */
	void print() const;                //to be modified, just for tests
	Matrix get_support_hyperplanes() const;
	Matrix get_generators() const;
	Matrix read_hilbert_basis() const;



/*---------------------------------------------------------------------------
 *				Computation Methods
 *---------------------------------------------------------------------------
 */
	void hilbert_basis_dual();

	/* transforms all data to the sublattice */
	void to_sublattice(Sublattice_Representation SR);

	void error(string s) const;
};
//class end *****************************************************************

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
