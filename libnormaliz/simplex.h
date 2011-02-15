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

//---------------------------------------------------------------------------

// The matrix Map is assumed to be such that the rank of Map equals
// the number of columns of Map and no test is performed in the constructor

//---------------------------------------------------------------------------
#ifndef SIMPLEX_H
#define SIMPLEX_H

//---------------------------------------------------------------------------

#include <vector>
#include <list>

#include "libnormaliz.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using std::list;
using std::vector;

template<typename Integer>
class Simplex {
  size_t dim;
  string status;
  Integer volume;
  vector<size_t> key;
  Matrix<Integer> Generators;
  vector< Integer > diagonal;
  vector< Integer > multiplicators;
  vector<size_t> New_Face; // to use with a shelling
  Matrix<Integer> Support_Hyperplanes;
  list< vector<Integer> > Hilbert_Basis;
  list< vector<Integer> > Ht1_Elements;
  vector<Integer> H_Vector;

//---------------------------------------------------------------------------
//                 Private routines, used in the public routines
//---------------------------------------------------------------------------

  void reduce_and_insert_interior(const vector< Integer >& new_element);
  //adds a new element to the Hilbert basis

//---------------------------------------------------------------------------
public:
//---------------------------------------------------------------------------
//                      Construction and destruction
//---------------------------------------------------------------------------

  Simplex();
  Simplex(const vector<size_t>& k);  //constructor, only key is initialized
  Simplex(const Matrix<Integer>& Map);  //contructor of the first in lexicographic
  //order simplex inside Map, the rank of Map is assumed to equal the number of
  //columns of Map
  Simplex(const vector<size_t>& k, const Matrix<Integer>& Map); //main constuctor
  //the rank of M is assumed to equal the number of columns of M
  Simplex(const Simplex<Integer>& S);   //copy constructor
  ~Simplex();                 //destructor

//---------------------------------------------------------------------------
//                          Data acces
//---------------------------------------------------------------------------

  void write_new_face(const vector<size_t>& Face);//writes the new face in case of a shelling
  void read() const;                        // to be modified, just for tests
  void read_k() const;                        // to be modified, just for tests
  size_t read_dimension() const;              // returns dim
  string read_status() const;              // returns status
  void write_volume(const Integer& vol);  // writes volume
  Integer read_volume() const;            // returns volume
  vector<size_t> read_key() const;          // returns key
  Matrix<Integer> read_generators() const;        // returns generators
  vector<Integer> read_diagonal() const;    // returns diagonal
  vector<Integer> read_multiplicators() const;    // returns multiplicators
  vector<size_t> read_new_face() const;    // returns new face
  size_t read_new_face_size() const;    // returns new face size
  Matrix<Integer> read_support_hyperplanes() const;  // returns the support hyperplanes
  Matrix<Integer> read_hilbert_basis()const; //read the Hilbert basis
  list< vector<Integer> > read_ht1_elements()const; //read the ht1 elements
  const list< vector<Integer> >& acces_hilbert_basis()const; //read the Hilbert basis
  vector<Integer> read_h_vector() const; //returns the h-vector
  size_t read_hilbert_basis_size() const; //returns the size of the Hilbert basis

//---------------------------------------------------------------------------
//                          Algoritms
//---------------------------------------------------------------------------

  int compare(const Simplex<Integer>& S) const; //compare the key of this with the key of S
  void initialize(const Matrix<Integer>& Map); //change status from "key initialized" to "initialized"
  void hilbert_basis_interior(); // computes the Hilbert basis,
  //the generators are not considered !!!, status must be "initialized"
  void hilbert_basis_interior(const Matrix<Integer>& Map); // computes the Hilbert basis,
  //the generators are not considered !!!, status may be "key initialized" or "initialized"
  void ht1_elements(const vector<Integer>& Form); //computes the volume and ht1 elements (unique without generators, stored in Hilbert_Basis)
  void h_vector(const vector<Integer>& Form);
  // computes the new elements of h vector in case of a shelling,
  // status must be "initialized"
  void h_vector(const Matrix<Integer>& Map, const vector<Integer>& Form); // computes
  // computes the new elements of h vector in case of a shelling,
  // status must be "key initialized"
  void hilbert_basis_interior_h_vector(const vector<Integer>& Form); // computes the Hilbert basis
  //and the new elements of h vector in case of a shelling,
  //the generators are not considered !!!, status must be "initialized"
  void hilbert_basis_interior_h_vector(const Matrix<Integer>& Map, const vector<Integer>& Form); // computes the Hilbert basis
  //and the new elements of h vector in case of a shelling,
  //the generators are not considered !!!, status must be "key initialized"
  
  Integer evaluate(Full_Cone<Integer>& C, const Integer& height);
  
};
//class end *****************************************************************

}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------

