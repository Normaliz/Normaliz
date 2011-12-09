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

#include <algorithm>
#include <string>
#include <iostream>
#include <set>

#include "integer.h"
#include "vector_operations.h"
#include "matrix.h"
#include "simplex.h"
#include "list_operations.h"
#include "HilbertSeries.h"

//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//Private
//---------------------------------------------------------------------------
template<typename Integer>
void SimplexEvaluator<Integer>::reduce_and_insert_interior(const vector< Integer >& new_element){
    //implementing this function as a tree searching may speed up computations ...
    if (new_element[0]==0) {
        return; // new_element=0
    }
    else {
        size_t i,c=1,d=dim+1;
        typename list< vector<Integer> >::iterator j;
        for (j =Hilbert_Basis.begin(); j != Hilbert_Basis.end(); j++) {
            if (new_element[0]<2*(*j)[0]) {
                break; //new_element is not reducible;
            }
            else  {
                if ((*j)[c]<=new_element[c]){
                    for (i = 1; i < d; i++) {
                        if ((*j)[i]>new_element[i]){
                            c=i;
                            break;
                        }
                    }
                    if (i==d) {
                        Hilbert_Basis.push_front(*j);
                        Hilbert_Basis.erase(j);
                        return;
                    }
                    //new_element is not in the Hilbert Basis
                }
            }
        }
        Hilbert_Basis.push_back(new_element);
    }
}

//---------------------------------------------------------------------------
//Public
//---------------------------------------------------------------------------

template<typename Integer>
Simplex<Integer>::Simplex(const Matrix<Integer>& Map){
    dim=Map.nr_of_columns();
    key=Map.max_rank_submatrix_lex(dim);
    Generators=Map.submatrix(key);
    diagonal = vector< Integer >(dim);
    Support_Hyperplanes=Invert(Generators, diagonal, volume); //test for arithmetic
    //overflow performed
    v_abs(diagonal);
    Support_Hyperplanes = Support_Hyperplanes.transpose();
    multiplicators = Support_Hyperplanes.make_prime();
}

//---------------------------------------------------------------------------

template<typename Integer>
Simplex<Integer>::Simplex(const vector<size_t>& k, const Matrix<Integer>& Map){
    key=k;
    Generators=Map.submatrix(k);
    dim=k.size();
    diagonal = vector< Integer >(dim);
    Support_Hyperplanes=Invert(Generators, diagonal, volume);  //test for arithmetic
    //overflow performed
    v_abs(diagonal);
    Support_Hyperplanes=Support_Hyperplanes.transpose();
    multiplicators=Support_Hyperplanes.make_prime();
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t Simplex<Integer>::read_dimension() const{
    return dim;
}

//---------------------------------------------------------------------------

template<typename Integer>
void Simplex<Integer>::write_volume(const Integer& vol){
    volume=vol;
}

//---------------------------------------------------------------------------

template<typename Integer>
Integer Simplex<Integer>::read_volume() const{
    return volume;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<size_t> Simplex<Integer>::read_key() const{
    return key;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> Simplex<Integer>::read_generators() const{
    return Generators;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Simplex<Integer>::read_diagonal() const{
    return diagonal;
}

//---------------------------------------------------------------------------

template<typename Integer>
vector<Integer> Simplex<Integer>::read_multiplicators() const{
    return multiplicators;
}

//---------------------------------------------------------------------------


template<typename Integer>
Matrix<Integer> Simplex<Integer>::read_support_hyperplanes() const{
    return Support_Hyperplanes;
}

//---------------------------------------------------------------------------

template<typename Integer>
SimplexEvaluator<Integer>::SimplexEvaluator(Full_Cone<Integer>& fc)
: C(fc),
  dim(C.dim),
  Generators(dim,dim),
  diagonal(dim),
  RS(dim,1)
{
}

//---------------------------------------------------------------------------


/* evaluates a simplex in regard to all data, key must be initialized */
template<typename Integer>
Integer SimplexEvaluator<Integer>::evaluate(const vector<size_t>& key, const Integer& height) {

    bool do_only_multiplicity=!C.do_h_vector && !C.do_Hilbert_basis && !C.do_ht1_elements;

    bool unimodular=false;
    vector<Integer> Indicator;
    if(height <= 1 || do_only_multiplicity) {
        for(size_t i=0; i<dim; ++i) {
            Generators.write_column(i+1,C.Generators[key[i]-1]);
        } // already transposed Generators
        RS.write_column(1,C.Order_Vector);  // right hand side

        Matrix<Integer> Sol(Generators.solve_destructiv(RS,diagonal,volume));
        Indicator.resize(dim);
        for (size_t i=0; i<dim; i++){
            Indicator[i]=Sol[i][0];
        }
        if(volume==1)
            unimodular=true;
    }
            
    // in this case we have to add the volume and nothing else is to be done
    if ( do_only_multiplicity || (unimodular && !C.do_h_vector) ) {
        #pragma omp critical(MULTIPLICITY)
        C.multiplicity+=volume;
        return volume;
    }

    // compute degrees of the generators
    vector<long> gen_degrees(dim);
    HilbertSeries Hilbert_Series;

    if (C.do_h_vector || C.do_ht1_elements) {
        //degrees of the generators according to the Grading of C
        for (size_t i=0; i<dim; i++){
            gen_degrees[i] = C.gen_degrees[key[i]-1];
        }
        if (C.do_h_vector) {
            int max_degree = *max_element(gen_degrees.begin(),gen_degrees.end());
            vector<long64> denom(max_degree+1);
            for (size_t i=0; i<dim; i++) {
                denom[gen_degrees[i]]++;
            }

            Hilbert_Series = HilbertSeries(vector<long64>(dim+1),denom);
        }
    }

    bool decided=true; // true if order vector in no hyperplane of simplex
    size_t i,j;
    long64 Deg=0;    // Deg is the degree according to Grading in which the 0 vector is counted
    if(unimodular) {  // it remains to count the 0-vector in the h-vector 
        for(i=0;i<dim;i++){
            if(Indicator[i]<0) {       // facet opposite of vertex i excluded
                Deg += gen_degrees[i];
            }
            else if(Indicator[i]==0) { // Order_Vector in facet, to be decided later
                decided=false;
                break;
            }
        }
        if(decided){
            //only change in the H-vector, so done directly
            Hilbert_Series.add_to_num(Deg);
            #pragma omp critical(HSERIES) 
            C.Hilbert_Series += Hilbert_Series;
            #pragma omp critical(MULTIPLICITY)
            C.multiplicity += volume;
            return volume;               // if not we need lex decision, see below
        }
    } // We have tried to take care of the unimodular case WITHOUT the matrix inversion

    for(size_t i=0; i<dim; ++i) {
        Generators[i] = C.Generators[key[i]-1];
    }
    Matrix<Integer> InvGen=Invert(Generators, diagonal, volume);
    v_abs(diagonal);
    vector<bool> Excluded(dim,false);
    Integer Test; 
    
    if(C.do_h_vector){
        Deg=0;
        if (Indicator.size() != dim) { //it hasn't been computed yet
            Indicator = InvGen.VxM(C.Order_Vector);
        }
        for(i=0;i<dim;i++) // excluded facets and degree shift for 0-vector
        {
            Test=Indicator[i];
            if(Test<0)
            {
                Excluded[i]=true; // the facet opposite to vertex i is excluded
                Deg += gen_degrees[i];
            }
            if(Test==0){  // Order_Vector in facet, now lexicographic decision
                for(j=0;j<dim;j++){
                    if(InvGen[j][i]<0){ // COLUMNS of InvGen give supp hyps
                        Excluded[i]=true;
                        Deg += gen_degrees[i];
                        break;
                    }
                    if(InvGen[j][i]>0) // facet included
                        break;
                }
            }
        }
        Hilbert_Series.add_to_num(Deg);
        if(unimodular){     // and in the unimodular case nothing left to be done
            #pragma omp critical(HSERIES) 
            C.Hilbert_Series += Hilbert_Series;
            #pragma omp critical(MULTIPLICITY)
            C.multiplicity += volume;
            return volume;
        }
    }
    
    vector < Integer > norm(1);
    Integer normG;
    list < vector<Integer> > Candidates;
    typename list <vector <Integer> >::iterator c;
    size_t last;
    vector<Integer> point(dim,0);
 
    Matrix<Integer> elements(dim,dim); //all 0 matrix
    Matrix<Integer> V = InvGen; //Support_Hyperplanes.multiply_rows(multiplicators).transpose();
    V.reduction_modulo(volume); //makes reduction when adding V easier
    vector<Integer> help;

    //now we need to create the candidates
    while (true) {
        last = dim;
        for (int k = dim-1; k >= 0; k--) {
            if (point[k] < diagonal[k]-1) {
                last = k;
                break;
            }
        }
        if (last >= dim) {
            break;
        }

        point[last]++;
        v_add_to_mod(elements[last], V[last], volume);

        for (i = last+1; i <dim; i++) {
            point[i]=0;
            elements[i] = elements[last];
        }
        
        norm[0]=0; // norm[0] is just the sum of coefficients, = volume*degree for standard grading
        normG = 0;
        for (i = 0; i < dim; i++) {  // since generators have degree 1
            norm[0]+=elements[last][i];
            if(C.do_h_vector || C.do_ht1_elements) {
                normG += elements[last][i]*gen_degrees[i];
            }
        }

        if(C.do_h_vector){
            Deg = explicit_cast_to_long<Integer>(normG/volume);
            for(i=0;i<dim;i++) { // take care of excluded facets and increase degree when necessary
                if(elements[last][i]==0 && Excluded[i]) {
                    Deg += gen_degrees[i];
                }
            }
            
            Hilbert_Series.add_to_num(Deg);
        }
        
        if(C.do_ht1_elements && normG==volume && height >= 2) // found degree 1 element
        {                                                      // only added if height >=2
            help=Generators.VxM(elements[last]);
            v_scalar_division(help,volume);
            Ht1_Elements.push_back(help);
            continue;
        } 
        
        // now we are left with the case of Hilbert bases
        if(C.do_Hilbert_basis && height >= 2){                 // only added if height >=2
            Candidates.push_back(v_merge(norm,elements[last]));
        }
    }
    
    if(C.do_h_vector) {
        #pragma omp critical(HSERIES)
        C.Hilbert_Series += Hilbert_Series;
    }
    
    if(C.do_ht1_elements) {
        #pragma omp critical(HT1ELEMENTS)
        C.Ht1_Elements.splice(C.Ht1_Elements.begin(),Ht1_Elements);
    }
    

    #pragma omp critical(MULTIPLICITY)
    C.multiplicity+=volume;

    if(!C.do_Hilbert_basis)
        return volume;  // no local reduction in this case

    Candidates.sort();        
    typename list <vector <Integer> >::iterator cand=Candidates.begin();
    while(cand != Candidates.end()) {
        reduce_and_insert_interior((*cand));
        Candidates.pop_front();
        cand=Candidates.begin();
    }

    //inverse transformation
    //some test for arithmetic overflow may be implemented here

    l_cut_front(Hilbert_Basis,dim);
    typename list< vector<Integer> >::iterator jj;
    for (jj =Hilbert_Basis.begin(); jj != Hilbert_Basis.end(); jj++) {
        *jj=Generators.VxM(*jj);
        v_scalar_division(*jj,volume);
    } 
    
    #pragma omp critical(CANDIDATES)
    C.Candidates.splice(C.Candidates.begin(),Hilbert_Basis);
    
    return volume;
}

} /* end namespace */
