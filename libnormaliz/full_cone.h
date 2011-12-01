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

#ifndef FULL_CONE_H
#define FULL_CONE_H

#include <list>
#include <vector>
#include <set>
#include <boost/dynamic_bitset.hpp>

#include "libnormaliz.h"
#include "cone_property.h"
#include "integer.h"
#include "matrix.h"
#include "simplex.h"
#include "cone_dual_mode.h"
#include "HilbertSeries.h"

namespace libnormaliz {
using std::list;
using std::vector;
using std::set;
using std::pair;
using boost::dynamic_bitset;

template<typename Integer> class Cone;
template<typename Integer> class Full_Cone;

//generates a lifted cone with the lower part simplicial, needed for computing the triangulation by lifting
template<typename Integer>
void lift(Full_Cone<Integer>& Lifted, Matrix<Integer> Extreme_Generators);


template<typename Integer>
class Full_Cone {
    size_t dim;
    size_t nr_gen;
    size_t hyp_size;
    
    bool pointed;
    bool ht1_generated;
    bool ht1_extreme_rays;
    bool ht1_hilbert_basis;
    bool integrally_closed;
    
    bool do_triangulation;
    bool do_partial_triangulation;
    bool do_Hilbert_basis;
    bool do_ht1_elements;
    bool do_h_vector;
    bool keep_triangulation;
    bool is_pyramid;
    
    ConeProperties is_Computed;
    vector<Integer> Linear_Form;
    Integer multiplicity;
    Matrix<Integer> Generators;
    vector<bool> Extreme_Rays;
    list<vector<Integer> > Support_Hyperplanes;
        
    vector<bool> in_triang;
    
    list<vector<Integer> > Hilbert_Basis;
    list<vector<Integer> > Candidates;
    list<vector<Integer> > Ht1_Elements;
    vector<Integer> Hilbert_Polynomial;
    HilbertSeries Hilbert_Series;
    vector<long> gen_degrees;  // will contain the degrees of the generators

    friend class Cone<Integer>;
    friend class Simplex<Integer>;
    
    struct SHORTSIMPLEX{                   // type for simplex, short in contrast to class Simplex
        vector<size_t> key;                // full key of simplex
        Integer height;                    // height of last vertex over opposite facet
    };
    
    list <SHORTSIMPLEX> Triangulation;      // triangulation of cone
    size_t TriangulationSize;               // number of elements in Triangulation, for efficiency
    vector<typename list <SHORTSIMPLEX>::iterator> TriSectionFirst;   // first simplex with lead vertex i
    vector<typename list <SHORTSIMPLEX>::iterator> TriSectionLast;     // last simplex with lead vertex i
    vector<size_t> VertInTri;               // generators in the order in which they are inserted
       
    struct FACETDATA {
        vector<Integer> Hyp;               // linear form of the hyperplane
        boost::dynamic_bitset<> GenInHyp;  // incidence hyperplane/generators
        Integer ValNewGen;                 // value of linear form on the generator to be added
        Integer ValPrevGen;                 // value on last generator added
    };
    
    list<FACETDATA> Facets;  // contains the data for Fourier-Motzkin and extension of triangulation
    
    vector<Integer> Order_Vector;  // vector for inclusion-exclusion

    Full_Cone<Integer>* Top_Cone;     // Reference to cone on top level
    vector<size_t> Top_Key;  // Indices of generators w.r.t Top_Cone
    
    int pyr_level;  // 0 for top cone, increased by 1 for each level of pyramids
    
    size_t totalNrSimplices;   // total number of simplices evaluated

/* ---------------------------------------------------------------------------
 *              Private routines, used in the public routines
 * ---------------------------------------------------------------------------
 */
    void add_hyperplane(const size_t& ind_gen, const FACETDATA & positive,const FACETDATA & negative);
    void extend_triangulation(const size_t& new_generator);
    void find_new_facets(const size_t& ind_gen);
    void process_pyramids(const size_t ind_gen,const bool recursive);
    void process_pyramid(FACETDATA& l, const size_t ind_gen,const bool recursive);

    void find_and_evaluate_start_simplex();
    Simplex<Integer> find_start_simplex() const;
    void store_key(const vector<size_t>&, const Integer& height);
    
    void build_cone();
    
    // Returns true if new_element is reducible versus the elements in Irred
    bool is_reducible(list<vector<Integer> *> & Irred, const vector<Integer> & new_element);
    // reduce the Candidates against itself and stores the remaining elements in Hilbert_Basis */
    void global_reduction();
    /* computes a degree function, s.t. every generator has value >0 */
    vector<Integer> compute_degree_function() const;

    void extreme_rays_and_ht1_check();
    void set_degrees();
    void compute_support_hyperplanes();
    void evaluate_triangulation();
    void transfer_triangulation_to_top();
    void primal_algorithm(); 
     
    // void support_hyperplanes_partial_triangulation();
    // void compute_support_hyperplanes_pyramid(const bool do_triang = false);

    void compute_extreme_rays();
    void select_ht1_elements();
    void compute_hilbert_basis();
    void compute_ht1_elements();
    void compute_hilbert_polynomial();
    void compute_hilbert_basis_polynomial();

    void check_pointed();
    void check_ht1_generated();
    void check_ht1_extreme_rays();
    void check_ht1_hilbert_basis();
    void check_integrally_closed();

    void compute_multiplicity();
    bool low_part_simplicial();
    void line_shelling();
    void triangulation_lift();

    /* support hyperplanes computation for a dynamic lifting
     * adjusts the lifting if necessary, used in dual algorithm */
    void support_hyperplanes_dynamic();


    /* constructor for recursively generated subcones
     * int i is a dummy parameter to distinguish it from the standard constructor */
    Full_Cone(Matrix<Integer> M, int i);
    
    void reset_tasks();

public:
/*---------------------------------------------------------------------------
 *                      Constructors
 *---------------------------------------------------------------------------
 */
    Full_Cone();
    Full_Cone(Matrix<Integer> M);            //main constructor
    Full_Cone(const Cone_Dual_Mode<Integer> &C);
    Full_Cone(Full_Cone<Integer>& C, vector<size_t> Key); // for pyramids

/*---------------------------------------------------------------------------
 *                      Data access
 *---------------------------------------------------------------------------
 */
    void print() const;             //to be modified, just for tests
    size_t getDimension() const;       //returns dimension
    size_t getNrGenerators() const;    //returns the number of generators
    bool isPointed() const;
    bool isHt1ExtremeRays() const;
    bool isHt1HilbertBasis() const;
    bool isIntegrallyClosed() const;
    vector<Integer> getLinearForm() const; //returns the linear form
    Integer getMultiplicity() const; //returns multiplicity
    const Matrix<Integer>& getGenerators() const;
    vector<bool> getExtremeRays() const;
    Matrix<Integer> getSupportHyperplanes() const;
    /* saves the triangulation and the volume of each simplex
     * the vectors corresponding to the generators of each simplex are sorted and saved in Triang
     * the volumes are saved in TriangVol 
     * cleares the lists first */
    void getTriangulation(list< vector<size_t> >& Triang, list<Integer>& TriangVol) const;
    Matrix<Integer> getHilbertBasis() const;
    Matrix<Integer> getHt1Elements() const;
    vector<Integer> getHVector() const;
    vector<Integer> getHilbertPolynomial() const;
    
    bool isComputed(ConeProperty::Enum prop) const; 


/*---------------------------------------------------------------------------
 *              Computation Methods
 *---------------------------------------------------------------------------
 */
    void dualize_cone();
    void support_hyperplanes();
    void support_hyperplanes_triangulation();
    void support_hyperplanes_triangulation_pyramid();
    void triangulation_hilbert_basis();
    void hilbert_basis();
    void hilbert_polynomial();
    void hilbert_polynomial_pyramid();
    void hilbert_basis_polynomial();
    void hilbert_basis_polynomial_pyramid();
    void ht1_elements();

    /* computes the multiplicity of the ideal in case of a Rees algebra
     * (not the same as the multiplicity of the semigroup) */
    Integer primary_multiplicity() const;

    /* completes the computation when a Cone_Dual_Mode is given */
    void dual_mode();

    void error_msg(string s) const;
};
//class end *****************************************************************
//---------------------------------------------------------------------------

}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
