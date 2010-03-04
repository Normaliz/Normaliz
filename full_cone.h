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

#ifndef FULL_CONE_H
#define FULL_CONE_H

#include <set>
#include <list>
#include "integer.h"
#include "matrix.h"
#include "simplex.h"


class Full_Cone{
  int dim;
  int nr_gen;
  int hyp_size;
  string status;
  bool homogeneous;
  vector<Integer> Linear_Form;
  Integer multiplicity;
  Matrix Generators;
  vector<bool> Extreme_Rays;
  list< vector<Integer> > Support_Hyperplanes;
  list< Simplex > Triangulation;
  list< vector<Integer> > Hilbert_Basis;
  list< vector<Integer> > Homogeneous_Elements;
  vector<Integer> H_Vector;
  vector<Integer> Hilbert_Polynomial;

  
//---------------------------------------------------------------------------
//			Private routines, used in the public routines
//---------------------------------------------------------------------------

  void add_hyperplane(const int& size,const vector<Integer>& positive_gen,const vector<Integer>& negative_gen);
  void transform_values(const int& size, const vector <int> & test_key);
  void add_simplex(const int& new_generator,const int& size,const vector<int>& col,const vector<int>& col_inv);
  void reduce_and_insert(const vector< Integer >& new_element);
  //adds a new element to the Hilbert basis
  void reduce_and_insert_speed(const vector< Integer >& new_element);
  //adds a new element to the Hilbert basis
  //faster as above, provided enough memory is available
  
   //retuns true if new_element is reducible versus the elements in Ired
  bool is_reducible(list< vector<Integer> >& Ired, const vector<Integer>& new_element);
  bool is_reducible(list< vector<Integer>* >& Ired, const vector<Integer>& new_element);
  
  bool reduce ( list < vector < Integer > > &  Ired , const vector< Integer >& new_element , const int& size );
  bool reduce ( list < vector < Integer >* > &  Ired , const vector< Integer >& new_element , const int& size );
  //retuns true if new_element is reducible versus the elements in Ired
  // used for dual algorithm
  void reduce (list < vector < Integer > > & Ired, list < vector< Integer > >& Red, const int& size );
 //reduce Red versus Ired
  void reduce_and_insert(const vector< Integer > & new_element, const int& size);
  //adds a new element ireducible to the  Hilbert basis
  // used for dual algorithm
  //the new elements must come from a  structure sorted by total degree
  void reduce_and_insert(const Matrix& New_Elements);
  //adds a matrix with new elements to the Hilbert basis
  void reduce_and_insert(const list< vector<Integer> >& New_Elements);
  //adds a list with new elements to the Hilbert basis
   void reduce_and_insert_extreme(const vector< Integer >& new_element);
  //select extreme rays  by reduction
  //used for the dual algoritm
  void find_new_face(); // to be used with a sheling in order to add to each simples the maximal new face
	void process_non_compressed(list< vector<int> > &non_compressed); //compute triangulations of the not compressed, not simplicial pieces and add them to Triangulation

  void only_hilbert_basis(); // computes only the Hilbert basis, support hyperplanes and triangulation must be computed in advance
  void global_reduction(set < vector<Integer> >& Candidates); //does the global reduction of the candidates
  vector<Integer> compute_degree_function() const;  //computes a degree function, s.t. every generator has value >0
//---------------------------------------------------------------------------
public:
//---------------------------------------------------------------------------
//				Construction and destruction
//---------------------------------------------------------------------------

  Full_Cone();
  Full_Cone(Matrix M);      //main constructor
  Full_Cone(const Full_Cone& C);  //copy constructor
  ~Full_Cone();            //destructor

//---------------------------------------------------------------------------
//						Data acces
//---------------------------------------------------------------------------

  void print() const;                   // to be modified, just for tests
  int read_dimension()const;         //returns dimension
  int read_nr_generators()const;    //returns the number of generators
  int read_hyp_size()const;         //returns hyp_size
  string read_status()const;       //returns status, may be:
								 // "non initialized" .......
  bool read_homogeneous() const; //returns homogeneous
  vector<Integer> read_linear_form() const; //returns the linear form
  Integer read_multiplicity() const; //returns multiplicity
  Matrix read_generators()const;                         //read the generators
  vector<bool>  read_extreme_rays()const;     //read the extreme rays
  Matrix read_support_hyperplanes()const; //read the support hyperplanes of the facets
  Matrix read_triangulation()const;      //read the triangulation
//  Matrix get_triangulation_list()const;      //read the triangulation
  Matrix read_triangulation_volume()const;      //read the triangulation and the volume
  //of each simplex , the volume is saved on the last column
  //the vectors coresponding to the generators of each simplex are sorted 
  Matrix read_hilbert_basis()const; //read the Hilbert basis
  Matrix read_homogeneous_elements()const; //read the homogeneous elements
  vector<Integer> read_h_vector() const; //returns the h-vector
  vector<Integer> read_hilbert_polynomial() const; //returns the Hilbert polynomial

//---------------------------------------------------------------------------
//					Algorithms
//---------------------------------------------------------------------------

  void extreme_rays(); //computes the extrem rays , the support hyperplanes must be known
  void support_hyperplanes(bool do_partial_triang=false); //computes the support hyperplanes
  void support_hyperplanes_pyramid(bool do_triang=false);  //computes the support hyperplanes via pyramids
  void support_hyperplanes_triangulation();//computes the support hyperplanes and triangulation
  void support_hyperplanes_triangulation_multiplicity();//calls support_hyperplanes_triangulation and computes the multiplicity
  void compute_multiplicity();//computes only the multiplicity
  void hilbert_basis(const bool compressed_test=false); // computes the Hilbert basis
  bool low_part_simplicial(); //computes the support hyperplanes and test if the lower part cone is simplicial
  void line_shelling(); //orders the lower part of the lifted cone after a line shelling
  void triangulation_lift(); //computed the triangulation by lifting
  vector<Integer> compute_e_vector(); //computes the e vector using the h vector
  void compute_polynomial(); // computes the Hilbert polynomial using the h-vector
  void hilbert_polynomial(); //computes the h-vector and the Hilbert polynomial
  void hilbert_basis_polynomial(); //computes the Hilbert basis and the Hilbert polynomial
  Integer primary_multiplicity() const; //computes the multiplicity of the ideal in case
  // of a Rees algebra (not the same as the multiplicity of the semigroup)
  void add_hyperplane(const int& hyp_counter, const bool& lifting, vector<Integer>& halfspace); //computes the Hilbert basis
  // after adding a support hyperplane with the dual algorithm
  Matrix add_hyperplane(const int& hyp_counter, const Matrix& Basis_Max_Subspace); //computes the Hilbert basis
  // after adding a support hyperplane with the dual algorithm , general case
  void extreme_rays_reduction(); //computes the extrem rays using reduction, used for the dual algorithm
  void extreme_rays_rank(); //computes the extrem rays using rank test, used for the dual algorithm
  void hilbert_basis_dual(); //computes the Hilbert basis with the dual algorithm

	bool check_compressed(); //checks if the cone is compressed, support hyperplanes must be computed
	void support_hyperplanes_dynamic();
//---------------------------------------------------------------------------
//						Error msg
//---------------------------------------------------------------------------
  void error(string s) const;
};
//class end *****************************************************************
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
 void lift(Full_Cone& Lifted,Matrix Extreme_Generators);//generates  a lifted cone with the lower part simplicial, needed for computing the  triangulation by lifting
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
