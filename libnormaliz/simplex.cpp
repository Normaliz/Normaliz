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

#include "my_omp.h"


//---------------------------------------------------------------------------

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------
//Private
//---------------------------------------------------------------------------
template<typename Integer>
bool SimplexEvaluator<Integer>::is_reducible_interior(const vector< Integer >& new_element){
    // the norm is at position dim
    if (new_element[dim]==0) {
        return true; // new_element=0
    }
    else {
        size_t i,c=0;
        typename list< vector<Integer> >::iterator j;
        for (j =Hilbert_Basis.begin(); j != Hilbert_Basis.end(); ++j) {
            if (new_element[dim]<2*(*j)[dim]) {
                break; //new_element is not reducible;
            }
            else  {
                if ((*j)[c]<=new_element[c]){
                    for (i = 0; i < dim; i++) {
                        if ((*j)[i]>new_element[i]){
                            c=i;
                            break;
                        }
                    }
                    if (i==dim) {
                        // move the reducer to the begin
                        Hilbert_Basis.splice(Hilbert_Basis.begin(), Hilbert_Basis, j);
                        return true;
                    }
                    //new_element is not in the Hilbert Basis
                }
            }
        }
        return false;
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
Simplex<Integer>::Simplex(const vector<key_t>& k, const Matrix<Integer>& Map){
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
vector<key_t> Simplex<Integer>::read_key() const{
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
  mult_sum(0),
  Generators(dim,dim),
  TGenerators(dim,dim),
  GenCopy(dim,dim),
  InvGenSelRows(dim,dim),
  InvGenSelCols(dim,dim),
  Sol(dim,dim+1),
  InvSol(dim,dim+1),
  GDiag(dim),
  TDiag(dim),
  Excluded(dim),
  Indicator(dim),
  gen_degrees(dim),
  RS(dim,1)
{
}

//---------------------------------------------------------------------------

size_t Unimod=0, Ht1NonUni=0, Gcd1NonUni=0, NonDecided=0, NonDecidedHyp=0;
    
/* evaluates a simplex in regard to all data, key must be initialized */
template<typename Integer>
Integer SimplexEvaluator<Integer>::evaluate(const vector<key_t>& key, const Integer& height) {
    

    bool do_only_multiplicity =
        (!C.do_h_vector && !C.do_Hilbert_basis && !C.do_ht1_elements)
        || (height==1 && C.do_partial_triangulation && !C.do_h_vector);
    
    size_t i,j;
    
    if(do_only_multiplicity){
        for(size_t i=0; i<dim; ++i)
            Generators[i] = C.Generators[key[i]];
        volume=Generators.vol_destructive();
        mult_sum += volume;
        return volume;         
    }  // done if only mult is asked for

    bool unimodular=false;
    bool GDiag_computed=false;
    bool potentially_unimodular=false;
    
    size_t g=0;
    if(height==1){
        for(i=0;i<dim;++i){
            g=gcd(g,gen_degrees[i]);
            if(g==1)
                break;        
        }
        potentially_unimodular=(g==1);
    }

    if(potentially_unimodular){ // very likely unimodular, Indicator computed first uses transpose of Gen
        for(i=0; i<dim; ++i)
            TGenerators.write_column(i,C.Generators[key[i]]); 
        RS.write_column(0,C.Order_Vector);  // right hand side
        TGenerators.solve_destructive_Sol(RS,TDiag,volume,Sol);
        for (i=0; i<dim; i++)
            Indicator[i]=Sol[i][0];
        if(volume==1){
            unimodular=true;
            #pragma omp atomic
            Unimod++;
            for(i=0;i<dim;i++)
                GDiag[i]=1;
        }
        else
            #pragma omp atomic
            Ht1NonUni++;
    }
    

    // we need the GDiag if not unimodular (to be computed from Gen)
    // if potentially unimodular, we combine its computation with that of the i-th support forms for Ind[i]==0
    // stored in InvSol (transferred to InvGenSelCols later)
    // if unimodular and all Ind[i] !=0, then nothing is done here
  
    vector<key_t> Ind0_key;  //contains the indices i as above 
    Ind0_key.reserve(dim-1);
    
    if(potentially_unimodular)
        for(i=0;i<dim;i++)
            if(Indicator[i]==0)
                Ind0_key.push_back(i);
    if(!unimodular || Ind0_key.size()>0){      
        for(i=0; i<dim; ++i)  // (uses Gen)
            Generators[i] = C.Generators[key[i]];
        if(!unimodular)
            GenCopy=Generators;
        if(Ind0_key.size()>0){
            Matrix<Integer> RSmult(dim,Ind0_key.size());
            for(i=0;i<Ind0_key.size();i++) // insert unit vectors
                RSmult[Ind0_key[i]][i]=1;
            Generators.solve_destructive_Sol(RSmult,GDiag,volume,InvSol);
            v_abs(GDiag);
            GDiag_computed=true;         
        }
        if(!GDiag_computed){
            Matrix<Integer> RSmult(dim,Ind0_key.size());
            Generators.solve_destructive_Sol(RSmult,GDiag,volume,InvSol);
            v_abs(GDiag);
            GDiag_computed=true;
        }
    }  
    
    mult_sum += volume;
    

    // now we must compute the matrix InvGenSelRows (selected rows of InvGen)
    // for those i for which Gdiag[i]>1 combined with computation
    // of Indicator in case of potentially_unimodular==false (uses transpose of Gen)
    

         
    vector<key_t> Last_key;
    Last_key.reserve(dim);       
    if (!unimodular) {
        for(i=0; i<dim; ++i) { 
            TGenerators.write_column(i,C.Generators[key[i]]);
            if(GDiag[i]>1)
                Last_key.push_back(i);
        }
        
        size_t RScol;
        if(potentially_unimodular)
            RScol=Last_key.size();
        else
            RScol=Last_key.size()+1;
        Matrix<Integer> RSmult(dim,RScol);
            
        for(i=0;i<Last_key.size();i++) // insert unit vectors
            RSmult[Last_key[i]][i]=1;
        if(!potentially_unimodular) // insert order vector if necessary
            RSmult.write_column(Last_key.size(),C.Order_Vector);
        TGenerators.solve_destructive_Sol(RSmult,TDiag,volume,Sol);
                // Sol.print(cout);
           
        for(i=0;i<Last_key.size();i++) // write solutions as selected rows of InvDen
            for(j=0;j<dim;j++){
                InvGenSelRows[Last_key[i]][j]=Sol[j][i]%volume; //makes reduction mod volume easier
                if(InvGenSelRows[Last_key[i]][j] <0)
                    InvGenSelRows[Last_key[i]][j]+=volume;
            }
        if(!potentially_unimodular) // extract Indicator
            for (i=0; i<dim; i++)
                Indicator[i]=Sol[i][Last_key.size()];
    }
    

        // InvGenSelRows.print(cout);exit(0);
    
    // if potentially_unimodular==false  it remains to compute support forms for i
    // with Ind[i]>0 (if there are any)

    
    if(!potentially_unimodular){
        for(i=0;i<dim;i++)
            if(Indicator[i]==0)
                Ind0_key.push_back(i);
        if(Ind0_key.size()>0){  
            Generators=GenCopy;
            Matrix<Integer> RSmult(dim,Ind0_key.size());
            for(i=0;i<Ind0_key.size();i++) // insert unit vectors
                    RSmult[Ind0_key[i]][i]=1;
            Generators.solve_destructive_Sol(RSmult,TDiag,volume,InvSol);  // kep GDiag from above     
        }
    }
    
    if(Ind0_key.size()>0){
        #pragma omp atomic
        NonDecided++;
        #pragma omp atomic
        NonDecidedHyp+=Ind0_key.size();
    }
    
    for(i=0;i<Ind0_key.size();i++) // insert selected columns of InvGen at right place
        for(j=0;j<dim;j++){
            InvGenSelCols[j][Ind0_key[i]]=InvSol[j][i];
        }
    
    // compute degrees of the generators and prepare Hilbert series if necessary
    vector<long> gen_degrees(dim);
    HilbertSeries Hilbert_Series;
    if (C.do_h_vector || C.do_ht1_elements) {
        //degrees of the generators according to the Grading of C
        for (size_t i=0; i<dim; i++){
            gen_degrees[i] = C.gen_degrees[key[i]];
        }
        if (C.do_h_vector) {
            int max_degree = *max_element(gen_degrees.begin(),gen_degrees.end());
            vector<denom_t> denom(max_degree+1);
            for (size_t i=0; i<dim; i++) {
                denom[gen_degrees[i]]++;
            }
            Hilbert_Series = HilbertSeries(vector<num_t>(dim+1),denom);
        }
    }
    
    Integer Test;
    size_t Deg=0;
    for(i=0;i<dim;i++)
        Excluded[i]=false;
    for(i=0;i<dim;i++){ // excluded facets and degree shift for 0-vector
        Test=Indicator[i];
        if(Test<0)
        {
            Excluded[i]=true; // the facet opposite to vertex i is excluded
            if(C.do_h_vector)
                Deg += gen_degrees[i];
        }
        if(Test==0){  // Order_Vector in facet, now lexicographic decision
            for(j=0;j<dim;j++){
                if(InvGenSelCols[j][i]<0){ // COLUMNS of InvGen give supp hyps
                    Excluded[i]=true;
                    if(C.do_h_vector)
                        Deg += gen_degrees[i];
                    break;
                }
                if(InvGenSelCols[j][i]>0) // facet included
                    break;
            }
        }
    }
    
    if(C.do_h_vector)
        Hilbert_Series.add_to_num(Deg); // count the 0-vector in k-vector with the right shift

    if(unimodular){  // do_h_vector==true automatically here
        C.HS[omp_get_thread_num()] += Hilbert_Series;
        return volume;
    } // the unimodular case has been taken care of
    

    // now we create and evaluate the points in par
    vector < Integer > norm(1);
    Integer normG;
    list < vector<Integer> > Candidates;
    typename list <vector <Integer> >::iterator c;
    size_t last;
    vector<Integer> point(dim,0);
 
    Matrix<Integer> elements(dim,dim); //all 0 matrix 
    vector<Integer> help;


    //now we need to create the candidates
    while (true) {
        last = dim;
        for (int k = dim-1; k >= 0; k--) {
            if (point[k] < GDiag[k]-1) {
                last = k;
                break;
            }
        }
        if (last >= dim) {
            break;
        }

        point[last]++;
        v_add_to_mod(elements[last], InvGenSelRows[last], volume);

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

        size_t Deg=0;
        if(C.do_h_vector){
            Deg = explicit_cast_to_long<Integer>(normG/volume);
            for(i=0;i<dim;i++) { // take care of excluded facets and increase degree when necessary
                if(elements[last][i]==0 && Excluded[i]) {
                    Deg += gen_degrees[i];
                }
            }
            
            Hilbert_Series.add_to_num(Deg);
        }
        
        // the case of Hilbert bases and height 1 elements, only added if height >=2
        if (!C.do_partial_triangulation || height >= 2) {
            if (C.do_Hilbert_basis) {
                Candidates.push_back(v_merge(elements[last],norm));
                continue;
            }
            if(C.do_ht1_elements && normG==volume && !isDuplicate(elements[last])) {
                help=GenCopy.VxM(elements[last]);
                v_scalar_division(help,volume);
                Ht1_Elements.push_back(help);
            }
        }
    }
        
    
    if(C.do_h_vector) {
        // #pragma omp critical(HSERIES)
        C.HS[omp_get_thread_num()] += Hilbert_Series;
    }
    
    if(C.do_ht1_elements) {
        #pragma omp critical(HT1ELEMENTS)
        C.Ht1_Elements.splice(C.Ht1_Elements.begin(),Ht1_Elements);
    }
 
    if(!C.do_Hilbert_basis)
        return volume;  // no local reduction in this case

    Candidates.sort(compare_last<Integer>);
    typename list <vector <Integer> >::iterator cand=Candidates.begin();
    while(cand != Candidates.end()) {
        if (is_reducible_interior(*cand)) // erase the candidate
            cand = Candidates.erase(cand);
        else // move it to the Hilbert basis
            Hilbert_Basis.splice(Hilbert_Basis.end(), Candidates, cand++);
    }

    //inverse transformation
    //some test for arithmetic overflow may be implemented here
    
    typename list< vector<Integer> >::iterator jj = Hilbert_Basis.begin();
    while (jj != Hilbert_Basis.end()) {
        if (isDuplicate(*jj)) { //delete the element
            jj = Hilbert_Basis.erase(jj);
        } else {
            jj->pop_back(); //remove the norm entry at the end
            *jj = GenCopy.VxM(*jj);
            v_scalar_division(*jj,volume);
            ++jj;
        }
    } 

    
    #pragma omp critical(CANDIDATES)
    C.Candidates.splice(C.Candidates.begin(),Hilbert_Basis);
        
    return volume;
}

template<typename Integer>
bool SimplexEvaluator<Integer>::isDuplicate(const vector<Integer>& cand) const {
    if (C.do_partial_triangulation) // we cannot use the criterion then
        return false;

    for (size_t i=0; i<dim; i++)
        if (cand[i]==0 && Excluded[i])
            return true;
    return false;
}

template<typename Integer>
Integer SimplexEvaluator<Integer>::getMultiplicitySum() const {
    return mult_sum;
}

} /* end namespace */
