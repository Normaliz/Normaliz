/*
 * Normaliz
 * Copyright (C) 2007-2019  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#include <iomanip>

#include "libnormaliz/cone.h"
#include "libnormaliz/descent.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/sublattice_representation.h"

namespace libnormaliz {

template <typename Integer>
DescentFace<Integer>::DescentFace() {
    // simplicial = false;
    coeff = 0;
    tree_size = 0;
    dead = false;
}

template <typename Integer>
DescentSystem<Integer>::DescentSystem() {
    descent_steps = 0;
    tree_size = 0;
    nr_simplicial = 0;
    system_size = 0;
    exploit_automorphisms = false;
}

template <typename Integer>
 void DescentSystem<Integer>::setExploitAutoms(bool on){
    exploit_automorphisms=on;
 }

template <typename Integer>
DescentSystem<Integer>::DescentSystem(const Matrix<Integer>& Gens_given,
                                      const Matrix<Integer>& SuppHyps_given,
                                      const vector<Integer>& Grading_given) {
    descent_steps = 0;
    tree_size = 0;
    nr_simplicial = 0;
    system_size = 0;
    exploit_automorphisms = false;

    Gens = Gens_given;
    SuppHyps = SuppHyps_given;
    Grading = Grading_given;

    nr_gens = Gens.nr_of_rows();
    nr_supphyps = SuppHyps.nr_of_rows();
    dim = Gens.nr_of_columns();

    GradGens.resize(nr_gens);
    GradGens_mpz.resize(nr_gens);
    for (size_t i = 0; i < nr_gens; ++i) {
        GradGens[i] = v_scalar_product(Grading, Gens[i]);
        convert(GradGens_mpz[i], GradGens[i]);
    }

    multiplicity = 0;

    SuppHypInd.resize(nr_supphyps);
    vector<size_t> NrFacetsContainingGen(nr_gens, 0);

    for (size_t i = 0; i < nr_supphyps; ++i) {
        INTERRUPT_COMPUTATION_BY_EXCEPTION

        SuppHypInd[i].resize(nr_gens);
        for (size_t j = 0; j < nr_gens; ++j)
            if (v_scalar_product(SuppHyps[i], Gens[j]) == 0) {
                SuppHypInd[i][j] = true;
                NrFacetsContainingGen[j]++;
            }
    }

    OldNrFacetsContainingGen.resize(nr_gens, 1);
    NewNrFacetsContainingGen.resize(nr_gens, 0);

    SimplePolytope = true;
    for (size_t j = 0; j < nr_gens; ++j) {
        if (NrFacetsContainingGen[j] > dim - 1) {
            SimplePolytope = false;
            break;
        }
    }

    OldNrFacetsContainingGen.resize(nr_gens, 1);
    NewNrFacetsContainingGen.resize(nr_gens, 0);
}


size_t nr_equalities = 0;

// The criteria for the optimal vertex in *this are discussed in the paper on the descent algorithm
template <typename Integer>
void DescentFace<Integer>::find_optimal_vertex(key_t& m_ind,    
                   const DescentSystem<Integer>& FF, const map<dynamic_bitset, dynamic_bitset>& FacetInds, const vector<key_t>& mother_key){     
    
    vector<size_t> count_in_facets(mother_key.size());
#pragma omp parallel for
    for (size_t i = 0; i < mother_key.size(); ++i) {
        size_t k = i;
        for (auto& FacetInd : FacetInds)
            if ((FacetInd.first)[k] == true)
                count_in_facets[k]++;
    }
    
    /* bool this_face_simple = true; // we test whether *this is a simple polytope
    for (size_t i = 0; i < mother_key.size(); ++i) { // could be used as a priori information for its children
        if(count_in_facets[i] > d){  // d-1){
            this_face_simple = false;
            break;
        }
    } */
   
    size_t m = count_in_facets[0];  // we must have at least one facet (actually 3, since dim 2 is simplicial)
    m_ind = 0;

    for (size_t i = 1; i < count_in_facets.size(); ++i) {
        if (count_in_facets[i] > m) {
            m = count_in_facets[i];
            m_ind = i;
            continue;
        }
        if (count_in_facets[i] == m &&
            FF.OldNrFacetsContainingGen[mother_key[i]] < FF.OldNrFacetsContainingGen[mother_key[m_ind]]) {
            m_ind = i;
        }
    }
}


// We want to find the sublattice defined by *this.
// mother_key is the vector of indices picking out the generators of *this
// The first goal is to find a small subset of the geberators spanning the sublattice up to saturation, i.e.
// a subset whose rank is equal to the dimension of *this. We try a random selection and increase it if necessary.
// The second goal, only computed if must_saturate = true, is the saturation, given back in Sublatt_this.
template <typename Integer>
void DescentFace<Integer>::find_sublattice(Matrix<Integer>& Gens_this, Sublattice_Representation<Integer>& Sublatt_this, 
                                           bool& sub_latt_computed, vector<key_t> mother_key, size_t dim,
                                           const Matrix<Integer>& FF_Gens){
    
    size_t FF_dim = FF_Gens[0].size();


    if (mother_key.size() > 3 * dim) { // 
        try {
            size_t nr_selected = 3 * dim;
            vector<key_t> selection;
            key_t j;
            size_t rk = 0;
            while (rk < dim && nr_selected <= mother_key.size()) {
                selection.resize(nr_selected);
                for (size_t i = 0; i < nr_selected; ++i) {
                    j = rand() % mother_key.size();
                    selection[i] = mother_key[j];
                }
                Gens_this = FF_Gens.submatrix(selection);
                rk = Gens_this.row_echelon();
                nr_selected *= 2;
            }
            if (rk < dim) {
                Gens_this = FF_Gens.submatrix(mother_key);
                Gens_this.row_echelon();
            }
        } catch (const ArithmeticException& e) {
            Gens_this = FF_Gens.submatrix(mother_key);
            Gens_this.row_echelon();
        }
    }
    else {
        Gens_this = FF_Gens.submatrix(mother_key);
        Gens_this.row_echelon();
    }

    bool must_saturate = false;

    for (size_t i = 0; i < Gens_this.nr_of_rows(); ++i) {
        for (size_t j = i; j < FF_dim; ++j) {
            if (Gens_this[i][j] == 0)
                continue;
            if (Gens_this[i][j] != 1 && Gens_this[i][j] != -1) {
                must_saturate = true;
            }
            break;
        }
        if (must_saturate)
            break;
    }

    sub_latt_computed = false;
    if (must_saturate){
        Sublatt_this = Sublattice_Representation<Integer>(Gens_this, true, false);  //  take saturation, no LLL
        sub_latt_computed = true;
    }
}

/*
template <typename Integer>
void DescentFace<Integer>::make_simplicial_facet(map<dynamic_bitset, vector<key_t> >& SimpKeys, map<dynamic_bitset, vector<bool> >& SimpInds,
                                                 map<dynamic_bitset, key_t>& CutOutBy, const bool ind_better_than_keys, 
                                                 const DescentSystem<Integer>& FF, const vector<key_t>& mother_key, 
                                                 map<dynamic_bitset, dynamic_bitset>& FacetInds, 
                                                 const dynamic_bitset& facet_ind, vector<key_t> facet_key){
    
    cout << "MMMMMM " << facet_ind << endl;
    bool set_FacetInds = false;

        if(facet_key.size() == 0){
            set_FacetInds = true;
            for(size_t i =0; i<mother_key.size(); ++i)
                if(facet_ind[i])
                    facet_key.push_back(i);            
        }

        // FacetInds[facet_ind] = dynamic_bitset(0);  // don't need support hyperplanes
        CutOutBy[facet_ind] = FF.nr_supphyps;  // signalizes "simplicial facet"
        if(set_FacetInds)
            FacetInds[facet_ind] = 
        if (ind_better_than_keys) {  // choose shorter representation
            vector<bool> gen_ind(FF.nr_gens);
            for (unsigned int k : facet_key)
                gen_ind[mother_key[k]] = 1;
            SimpInds[facet_ind] = gen_ind;
        }
        else {
            vector<key_t> trans_key;  // translate back to FF
            for (unsigned int k : facet_key)
                trans_key.push_back(mother_key[k]);
            SimpKeys[facet_ind] = trans_key;  // helps to pick the submatrix of its generators
        }
}
*/

/*
======= From signed dec mit merge übernommen.
    // Now we find the potential facets of *this.
    dynamic_bitset facet_ind(mother_key.size());    // lists Gens
    map<dynamic_bitset, dynamic_bitset> FacetInds;  // potential facets
    map<dynamic_bitset, key_t> CutOutBy;            // the facet citting it out
>>>>>>> e1b9050180645b2a330b743b353bffe66906d5cc*/

// We need to find the facets of *this.
// FacetInds are explained in compute, ditto CutOutBy.
// SimpInds and SimpKeys are special versions of FacetInds for the simplicial facets.
template <typename Integer>
void DescentFace<Integer>::find_facets(map<dynamic_bitset, dynamic_bitset>& FacetInds, 
                                       map<dynamic_bitset, key_t>& CutOutBy,
                                       map<dynamic_bitset, vector<key_t> >& SimpKeys, 
                                       map<dynamic_bitset, vector<bool> >& SimpInds,
                                       
                                       const bool ind_better_than_keys,                                       
                                       const DescentSystem<Integer>& FF, 
                                       const vector<key_t>& mother_key, 
                                       const dynamic_bitset& facets_cutting_mother_out, 
                                       size_t dim){ 
   
    size_t nr_supphyps = FF.nr_supphyps;
    size_t d = dim;
    
    dynamic_bitset facet_ind(mother_key.size());    // lists Gens, local variable for work

    for (size_t i = 0; i < nr_supphyps; ++i) {
        if (facets_cutting_mother_out[i] == true)  // contains *this
            continue;

        // we can identify the facet(*this) uniquely only via the Gens in it
        vector<libnormaliz::key_t> facet_key; // keys of extreme rays in current supphyp
        for (size_t k = 0; k < mother_key.size(); ++k) {
            if (FF.SuppHypInd[i][mother_key[k]] == true)
                facet_key.push_back(k);
        }
        if (facet_key.size() < d - 1)  // can't be a facet(*this)
            continue;

        // now we make facet_ind out of facet_key: key for gens in potatntial facet
        facet_ind.reset();
        for (unsigned int jj : facet_key)
            facet_ind[jj] = true;

        // next we check whether we have the intersection already
        // not necessary for simple polytopes and in top dimension
        // Note: if P is simple, F is a face of P and H a support hyperplave of P,
        // then F\cap H is either empty or a facet of F. Moreover H is eniquely determined
        // by F\cap H. This will again be used below.
        if (d < FF.dim && !FF.SimplePolytope) {
            if (FacetInds.find(facet_ind) != FacetInds.end()) {  // already found, we need it only once
                if (facet_key.size() > d - 1)
                    FacetInds[facet_ind][i] = true;
                // but in the nonsimplicial case we must add SuppHyps[i] to the facets(C) containing
                // the current facet(*this)
                continue;
            }
        }

        // now we have a new potential facet

        if (facet_key.size() == d - 1) {               // simplicial or not a facet
            FacetInds[facet_ind] = facets_cutting_mother_out;
            FacetInds[facet_ind][i] = true;
            CutOutBy[facet_ind] = i + FF.nr_supphyps;  // signalizes "simplicial facet"
            if (ind_better_than_keys) {  // choose shorter representation
                vector<bool> gen_ind(FF.nr_gens);
                for (unsigned int k : facet_key)
                    gen_ind[mother_key[k]] = 1;
                SimpInds[facet_ind] = gen_ind;
            }
            else {
                vector<key_t> trans_key;  // translate back to FF
                for (unsigned int k : facet_key)
                    trans_key.push_back(mother_key[k]);
                SimpKeys[facet_ind] = trans_key;  // helps to pick the submatrix of its generators
            }
        }
        else {
            FacetInds[facet_ind] = facets_cutting_mother_out;
            FacetInds[facet_ind][i] = true;  // plus the facet cutting out facet_ind
            CutOutBy[facet_ind] = i;         // memorize the facet that cuts it out
        }
    }

    if (d < FF.dim && !FF.SimplePolytope) {  // now we select the true facets of *this
        auto G = FacetInds.end();            // by taking those with a maximal set of gens
        for (--G; G != FacetInds.begin(); --G) {   // we use the lexicographic sorting by the "first" in the map
            for (auto F = FacetInds.begin(); F != G;) {   // superset ==> lex larger
                if (F->first.is_subset_of(G->first))
                    F = FacetInds.erase(F);
                else
                    ++F;
            }
        }
    }
}

// The same task as that of find_facets. But we can use extra information: For each facet F of *this we know a support hyperplane
// of the global cone that cuts out F from *this.
template <typename Integer>
void DescentFace<Integer>::find_facets_from_FacetsOfFace(map<dynamic_bitset, dynamic_bitset>& FacetInds, map<dynamic_bitset, key_t>& CutOutBy,
                                       map<dynamic_bitset, vector<key_t> >& SimpKeys, map<dynamic_bitset, vector<bool> >& SimpInds,
                                       vector<key_t>& Orbit,
                                       
                                       const bool ind_better_than_keys,                                       
                                       const DescentSystem<Integer>& FF, const vector<key_t>& mother_key, 
                                       const dynamic_bitset& facets_cutting_mother_out, size_t dim){


    // first we make the fecets
    for(size_t i=0; i< FF.nr_supphyps; ++i){
        if(!FacetsOfFace[i])
            continue;        
        dynamic_bitset facet_ind(mother_key.size());    // lists Gens, local variable for work 
        size_t nr_vertices = 0;
        for(size_t j=0; j< mother_key.size(); ++j){
            if(FF.SuppHypInd[i][mother_key[j]]){
                facet_ind[j] = true;
                nr_vertices++;
            }
        }
        FacetInds[facet_ind] = facets_cutting_mother_out; // furthers will be added
        CutOutBy[facet_ind] = i;
        if(nr_vertices == dim -1){  // simplicial facet
            CutOutBy[facet_ind] += FF.nr_supphyps;
            vector<key_t> facet_key = bitset_to_key(facet_ind);
            
            if (ind_better_than_keys) {  // choose shorter representation
            vector<bool> gen_ind(FF.nr_gens);
            for (unsigned int k : facet_key)
                gen_ind[mother_key[k]] = 1;
            SimpInds[facet_ind] = gen_ind;
            }
            else {
                vector<key_t> trans_key;  // translate back to FF
                for (unsigned int k : facet_key)
                    trans_key.push_back(mother_key[k]);
                SimpKeys[facet_ind] = trans_key;  // helps to pick the submatrix of its generators
            } 
        }
    }
    
    // Now wefind the global SuppHyps meeting the facets
    for(size_t i=0;i<FF.nr_supphyps; ++i){
        if(facets_cutting_mother_out[i]) // already covered
            continue;
        dynamic_bitset facet_ind(mother_key.size());
        for(size_t j=0; j< mother_key.size(); ++j){
            if(FF.SuppHypInd[i][mother_key[j]])
                facet_ind[j] = true;      
        }
        auto F = FacetInds.find(facet_ind);
        if(F == FacetInds.end()) // does not intresect in a facet
            continue;
        F->second[i] = true;
    }
    
    // cout << "FFFFFFFFFFFFFFFFFFFFFFFFFF " << FacetsOfFace << endl;
    

    // cout << FacetOrbits;

    if(FacetOrbits.empty())
        return;
    Orbit.resize(FF.nr_supphyps, FF.nr_supphyps);
    size_t k =0;
    for(auto& Orb: FacetOrbits){
        for(size_t i=0; i< FF.nr_supphyps; ++i){
            if(Orb[i])
                Orbit[i] = k;
        }
        k++;        
    }
    // cout << "GGGGGGGGGGGG " << Orbit;
}

    
template <typename Integer>
void DescentFace<Integer>::compute(DescentSystem<Integer>& FF, 
            size_t dim,  //  dim of *this
            const dynamic_bitset& facets_cutting_mother_out, // indicates the supphyps of which *this is the intersection
            size_t mother_tree_size) {
    long omp_start_level = omp_get_level();

    vector<key_t> mother_key;

    size_t nr_supphyps = FF.nr_supphyps;
    size_t nr_gens = FF.nr_gens;

    size_t d = dim;
    
    // cout << "FFFFFFFFFF " << coeff << endl;
 
    // -------------------------------------------------------------
    // find extreme_rays

    dynamic_bitset GensInd(nr_gens); // find the extreme rays of *this, indicated by GensInd,
    GensInd.set();                   // by intersecting the indicators of facets_cutting_mother_out

    for (size_t i = 0; i < nr_supphyps; ++i) {
        if (facets_cutting_mother_out[i] == true) {
            GensInd = GensInd & FF.SuppHypInd[i];
        }
    }

    for (size_t i = 0; i < nr_gens; ++i)
        if (GensInd[i])
            mother_key.push_back(i);
        
    // -------------------------------------------------------------
    // find sublattice (sometimes only Gens_this)

    Matrix<Integer> Gens_this;
    Sublattice_Representation<Integer> Sublatt_this;
    bool sub_latt_computed;
    
    find_sublattice(Gens_this, Sublatt_this,sub_latt_computed, mother_key,dim, FF.Gens);   

    // -------------------------------------------------------------
    // Now we find the potential facets of *this.

    map<dynamic_bitset, dynamic_bitset> FacetInds;  // potential facets, map: gens(potential facet)
                                                    // --> set of supphyps(C) containing these gens
                                                    // reference for gens(potential facet) is the selection via mother_key
    map<dynamic_bitset, key_t> CutOutBy;            // the global facet citting it out (not unique a priori)

    // We try to make the descent algorithm well-bahevd alsofor cones with many facets, provided not too many are non-simplicial
    // SimpKeys save memory in this case. Perhaps a luxury.
    map<dynamic_bitset, vector<key_t> > SimpKeys;  // generator keys for simplicial facets
    map<dynamic_bitset, vector<bool> > SimpInds;   // alternative: generator indices for simplicial facets (if less memory needed)
    bool ind_better_than_keys = (dim * 64 > FF.nr_gens); // decision between the alternatives
    
    vector<key_t> Orbit; // assigns global supphyp with index i its orbit in the facets of this face 
                        // (only for those in FacetsOfFace, of course)
    
    if(FacetsOfFace.size() == 0)
        find_facets(FacetInds, CutOutBy, SimpKeys,  SimpInds,
                ind_better_than_keys,FF, mother_key, facets_cutting_mother_out, dim);
    else
        find_facets_from_FacetsOfFace(FacetInds, CutOutBy, SimpKeys, SimpInds, Orbit,
                ind_better_than_keys,FF, mother_key, facets_cutting_mother_out, dim);
        
    
    // if we don't have the coordinate transformation and there is a simplicial facet, we must make the transformation
    if (!sub_latt_computed && (SimpKeys.size() > 0 || SimpInds.size() > 0)){
        Sublatt_this = Sublattice_Representation<Integer>(Gens_this, true, false);  //  take saturation, no LLL
        sub_latt_computed = true;
    }

    // At this point we know the facets of *this.
    // The map FacetInds assigns the set of containing SuppHyps to the facet_ind(Gens).
    // The set of containing SuppHyps is a unique signature as well.
    
    // -------------------------------------------------------------
    // Next we choose the "opposite" vertex
    
    key_t m_ind = 0; // will indicate the optimal vertex relative to mother_key
    find_optimal_vertex(m_ind, FF, FacetInds, mother_key); 

    key_t selected_gen = mother_key[m_ind];  // this is the global index of the selected generator 
    vector<Integer> embedded_selected_gen;
    if (sub_latt_computed)
        embedded_selected_gen = Sublatt_this.to_sublattice(FF.Gens[selected_gen]);
    
    // -------------------------------------------------------------
    // The last stp: process the facets opposite to the chosen vertex
    //
    // Thwe processing is divided into (I) nonsimplicial facets, (II) simplicial facets 
    //
    // Note: for (I) we search in this routune (i) an already computed ideantical copy of this face and
    // (ii) if automorphisms are exploired, an isomorphic copy. The copies are updated if found. Otherwise insertion
    // of this face.

    vector<vector<key_t> > OppositeOrbits; // opposite facets grouped in orbits

    
    if(FacetOrbits.size() > 0){
        OppositeOrbits.resize(FacetOrbits.size());
        auto G = FacetInds.begin();
        for (; G != FacetInds.end(); ++G) {
            if ((G->first)[m_ind] == true) 
                continue;
            size_t COB = CutOutBy[G->first];
            if(COB >= nr_supphyps){  // take care of simplicial facets
                COB -= nr_supphyps;
            }
            OppositeOrbits[Orbit[COB] ].push_back(COB);        
        }
    }
    /*
    cout << "------------------" << endl;
    cout << OppositeOrbits;
    cout << "==================" << endl;
    */

    vector<Integer> embedded_supphyp;
    Integer ht;    
    mpq_class divided_coeff = coeff / FF.GradGens_mpz[selected_gen]; // coeff = coeff(*this)
    vector<Integer> HeightSumOrbit(FacetOrbits.size()); // accumulates the heights in each orbit

    auto G = FacetInds.begin();
    for (; G != FacetInds.end(); ++G) { // loop over nonsimplicial facets
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        size_t COB = CutOutBy[G->first];
        

        if ((G->first)[m_ind] == true || COB >= FF.nr_supphyps)  // is not opposite or is simplicial
            continue; 
        
        // auto H = Children.insert(Children.begin(),make_pair(G->second,DescentFace<Integer>()) );         
        if (sub_latt_computed) {
            embedded_supphyp = Sublatt_this.to_sublattice_dual(FF.SuppHyps[COB]);
            ht = v_scalar_product(embedded_selected_gen, embedded_supphyp);
        }
        else {
            embedded_supphyp = Gens_this.MxV(FF.SuppHyps[COB]);
            Integer den = v_make_prime(embedded_supphyp);
            ht = v_scalar_product(FF.Gens[selected_gen], FF.SuppHyps[COB]) / den;
        }
        
        mpq_class new_coeff;
        
        if(FacetOrbits.size() > 0){
            size_t orbit = Orbit[COB];
#pragma omp critical(HEIGHTSUM)
            HeightSumOrbit[orbit] += ht;
            if(COB != OppositeOrbits[orbit].back()) // not the last one in its orbit
                continue;
            new_coeff = divided_coeff * convertTo<mpz_class>(HeightSumOrbit[orbit]);
        }
        else{
            new_coeff = divided_coeff * convertTo<mpz_class>(ht);
        }
        
        //---------------------------------------------------------------
        // give computation results to DescentSystem           
        
        auto H = FF.NewFaces.begin();
        bool inserted = false;
        bool found_identical = false;
        
#pragma omp critical(INSERT) // trying to find ideantical copy
        {
        H = FF.NewFaces.find(G->second);  // try to find identical face
        
        if(H != FF.NewFaces.end()){// identical face found
            (H->second).coeff += new_coeff;
            found_identical = true;
        }
        else{ // face itself has not yet appeared
            if(!FF.exploit_automorphisms){ // isomorphic copy can't be used
                H = FF.NewFaces.insert(FF.NewFaces.begin(), make_pair(G->second, DescentFace<Integer>()) );
                H->second.coeff = new_coeff;
                inserted = true;
            }
        }
//        } // critical(INSERT)  // trying to find ideantical copy
        
// #pragma omp critical(INSERT)  // trying to find isomorphic copy
//            {  
        
        if(!(inserted || found_identical)){ // preparations for using the automorphisms
                      
            // first we compute the facets of our face
            dynamic_bitset ExtRaysFacet(FF.nr_gens); // in the first step we must tranlate
            for(size_t kk=0; kk< G->first.size(); ++kk){  // G->first into an indictaor relative to the global 
                if( (G->first)[kk])                       // list of extreme rays
                    ExtRaysFacet[mother_key[kk]] = 1;
            }
            dynamic_bitset FacetCandidates = ~G->second; // indicates the gobal support bhyperplanes
                                                            // intersecting this facet in a proper subset
            vector<dynamic_bitset> Intersections(FF.nr_supphyps, dynamic_bitset(nr_gens));
            for(size_t i=0; i < FF.nr_supphyps; ++i){
                if(FacetCandidates[i] == 0)
                    continue;
                Intersections[i] = ExtRaysFacet & FF.SuppHypInd[i];
            }
            
            dynamic_bitset TheFacets = ~G->second; // some bits are preset. maximal_subsets takes care of it.
            maximal_subsets(Intersections, TheFacets); // indicator of global supphyps cutting out facets of this child
            
            auto N = FF.NewFaces.insert(FF.NewFaces.begin(), make_pair(G->second, DescentFace<Integer>()) );
            N->second.coeff = new_coeff;
            N->second.FacetsOfFace = TheFacets;
            inserted = true;
        }
        } // critical INSERT
        
        if (inserted) { // update statistics for optimization
            for (unsigned int& i : mother_key)
                if (FF.SuppHypInd[COB][i])
#pragma omp atomic
                    FF.NewNrFacetsContainingGen[i]++;
        }
        
#pragma omp atomic
        FF.descent_steps++;

#pragma omp atomic        
        (H->second).tree_size += mother_tree_size;
        
    }  // end loop over nonsimplicial facets
    
    //---------------------------------------------------------------
    
    if (SimpKeys.size() > 0 || SimpInds.size() > 0) {  // have to compute simplicial facets 
        G = FacetInds.begin();
        size_t loop_length = FacetInds.size();
        size_t fpos = 0;
        bool skip_remaining = false;
        vector<mpq_class> thread_mult(omp_get_max_threads(), 0);
        Matrix<Integer> Embedded_Gens(d, d);
        Matrix<Integer> Gens_this(d, FF.dim);

        std::exception_ptr tmp_exception;

#pragma omp parallel for firstprivate(G, fpos, Embedded_Gens, Gens_this, ht, embedded_supphyp)
        for (size_t ff = 0; ff < loop_length; ++ff) {
            if (skip_remaining)
                continue;
            for (; ff > fpos; ++fpos, ++G)
                ;
            for (; ff < fpos; --fpos, --G)
                ;

            int tn;
            if (omp_get_level() == omp_start_level)
                tn = 0;
            else
                tn = omp_get_ancestor_thread_num(omp_start_level + 1);

            try {
                INTERRUPT_COMPUTATION_BY_EXCEPTION
                
                size_t COB = CutOutBy[G->first];

                if ((G->first)[m_ind] == true || COB < nr_supphyps) // is not opposite or not simplicial
                    continue;
                
                COB -= nr_supphyps;
                
                if(FacetOrbits.size() > 0){
                        
                    if (sub_latt_computed) {  // ijn this case we need the height
                        embedded_supphyp = Sublatt_this.to_sublattice_dual(FF.SuppHyps[COB]);
                        ht = v_scalar_product(embedded_selected_gen, embedded_supphyp);
                    }
                    else {
                        embedded_supphyp = Gens_this.MxV(FF.SuppHyps[COB]);
                        Integer den = v_make_prime(embedded_supphyp);
                        ht = v_scalar_product(FF.Gens[selected_gen], FF.SuppHyps[COB]) / den;
                    }
                    
                    size_t orbit = Orbit[COB];
#pragma omp critical(HEIGHTSUM)
                    HeightSumOrbit[orbit] += ht;
                    if(COB != OppositeOrbits[orbit].back()) // not the last one in its orbit
                        continue;
                }

            
                if (ind_better_than_keys)
                    Gens_this = FF.Gens.submatrix(SimpInds[G->first]);
                else
                    Gens_this = FF.Gens.submatrix(SimpKeys[G->first]);
                Gens_this.append(FF.Gens[selected_gen]);
                Integer det;
                if (Sublatt_this.IsIdentity())
                    det = Gens_this.vol();
                else {
                    Embedded_Gens = Sublatt_this.to_sublattice(Gens_this);
                    det = Embedded_Gens.vol();
                }
                
                if(FacetOrbits.size() >0){
                    det /= ht; // must have multiplicity of opposite facet
                    det *= HeightSumOrbit[Orbit[COB]]; // and multiply it by the height sum of the orbit
                }                
                
                mpz_class mpz_det = convertTo<mpz_class>(det);
                mpq_class multiplicity = mpz_det;
                if (ind_better_than_keys) {
                    for (size_t i = 0; i < FF.nr_gens; ++i)
                        if (SimpInds[G->first][i] && FF.GradGens[i] > 1)
                            multiplicity /= FF.GradGens_mpz[i];
                }
                else {
                    for (size_t i = 0; i < Gens_this.nr_of_rows() - 1; ++i)
                        if (FF.GradGens[SimpKeys[G->first][i]] > 1)
                            multiplicity /= FF.GradGens_mpz[SimpKeys[G->first][i]];
                }
                if (FF.GradGens[selected_gen] > 1)
                    multiplicity /= FF.GradGens_mpz[selected_gen];
                thread_mult[tn] += multiplicity;
#pragma omp atomic
                FF.nr_simplicial++;
#pragma omp atomic
                FF.tree_size += tree_size;

            } catch (const std::exception&) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
        } // end loop over simplicial facets

        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);

        mpq_class local_multiplicity = 0;
        for (const auto& j : thread_mult)
            local_multiplicity += j;
#pragma omp critical(ADD_MULT)
        FF.multiplicity += local_multiplicity * coeff;

    } // end computation simplicial facets
    
    //---------------------------------------------------------------
}


template <typename Integer>
void DescentFace<Integer>::process_iso_class_of_face(const IsoType<Integer> IT,
             map< IsoType<Integer>, DescentFace<Integer>* , IsoType_compare<Integer> >& Isos){

#pragma omp critical(INSERT_ISO)
    {    
    auto G = Isos.find(IT); 
    if(G != Isos.end()){
        dead = true; // to be skipped in descent
        mpz_class index_source = convertTo<mpz_class>(IT.index);
        mpz_class index_traget = convertTo<mpz_class>(G->first.index);
        mpq_class corr = index_source;
        corr /= index_traget;
        // cout << "--------------------------" << endl;
        // cout << "Index Source " << index_source << " Index Target " << index_traget << " CORR " << corr << endl;
        G->second->coeff += corr*coeff;
        // cout << "coeff Source " << F->second.coeff << " Coeff Target " << G->second->coeff << endl;
        // cout << "--------------------------" << endl; */
        assert(corr == 1);
    }
    else{
        Isos[IT]= &(*this);
    }
    }    
}

template <typename Integer>
void DescentSystem<Integer>::collect_old_faces_in_iso_classes(){
    
    if(OldFaces.size() <= 1) // nothingt to do here
        return;
    
    if(verbose)
        verboseOutput() << "Collecting isomorphism classes" << endl;
    
    map< IsoType<Integer>, DescentFace<Integer>* , IsoType_compare<Integer> > Isos;
    
    size_t nr_F = OldFaces.size();
    auto F = OldFaces.begin();
    size_t kkpos = 0;
    std::exception_ptr tmp_exception;
    bool skip_remaining = false;
    
    // First we compute the isomorphism classes

#pragma omp parallel for firstprivate(kkpos, F) // schedule(dynamic)    
    for (size_t kk = 0; kk < nr_F; ++kk) {
        
        if(skip_remaining)
            continue;
        try {
            INTERRUPT_COMPUTATION_BY_EXCEPTION

            for (; kk > kkpos; kkpos++, F++)
                ;
            for (; kk < kkpos; kkpos--, F--)
                ;
            
            Matrix<Integer> Equations = SuppHyps.submatrix(bitset_to_key(F->first));
            Matrix<Integer> Inequalities = SuppHyps.submatrix(bitset_to_key(F->second.FacetsOfFace));
            IsoType<Integer> IT(Inequalities, Equations,Grading);
            
            vector<key_t> SupphypsCuttingOutFacets = bitset_to_key(F->second.FacetsOfFace);
            for(auto& Orb: IT.FacetOrbits){ // tranlate into indicator vectors using indices of global support hyperplanes
                dynamic_bitset orbit(nr_supphyps);
                for(size_t i=0; i < Orb.size(); ++i){
                    if(Orb[i])
                        orbit[SupphypsCuttingOutFacets[i]] = true;                    
                }
                F->second.FacetOrbits.push_back(orbit);
            
            }
            
            F->second.process_iso_class_of_face(IT, Isos);
            

        } catch (const std::exception&) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
#pragma omp flush(skip_remaining)
        }
    }  // parallel for kk
    
    if (!(tmp_exception == 0))
        std::rethrow_exception(tmp_exception);

    
    if(verbose)
        verboseOutput() << Isos.size() << " classes found" << endl;

    // cout << "Iso types " << Isos.size() << endl;
    /*for(auto& F: OldFaces){
        cout << "DDDD " << F.second.dead << " CCCC " << F.second.coeff << endl;
    }*/

}

template <typename Integer>
void DescentSystem<Integer>::compute() {
    
#ifdef NMZ_EXTENDED_TESTS
    if(!using_GMP<Integer>() && !using_renf<Integer>() && test_arith_overflow_descent)
        throw ArithmeticException(0);    
#endif
    
    if (verbose) {
        if (SimplePolytope)
            verboseOutput() << "Polytope is simple" << endl;
        else
            verboseOutput() << "Polytope is not simple" << endl;
    }

    const size_t ReportBound = 400;
    const size_t MaxBlocksize = 1000000;

    dynamic_bitset empty(nr_supphyps);
    DescentFace<Integer> top;
    OldFaces[empty] = top;
    OldFaces[empty].coeff = 1;
    OldFaces[empty].tree_size = 1;
    
    if(exploit_automorphisms){  // prepare orbits of global facets
        Matrix<Integer> Equations(0,dim);
        Matrix<Integer> Inequalities = SuppHyps;
        IsoType<Integer> IT(Inequalities, Equations,Grading);
        mpz_class index_source = convertTo<mpz_class>(IT.index);
        if(index_source == 1){
            OldFaces[empty].FacetsOfFace = dynamic_bitset(nr_supphyps);
            OldFaces[empty].FacetsOfFace.set();
            OldFaces[empty].FacetOrbits = IT.FacetOrbits;
        }        
    }
    
    long d = (long)dim;
    bool top_cone = true;

    while (!OldFaces.empty()) {

        size_t nr_F = OldFaces.size();
        system_size += nr_F;
        if (verbose)
            verboseOutput() << "Descent from dim " << d << ", size " << nr_F << endl;
        
        if(exploit_automorphisms && !top_cone)
            collect_old_faces_in_iso_classes();
    
        top_cone = false;

        bool in_blocks = false;
        if (nr_F > MaxBlocksize)
            in_blocks = true;
        if (in_blocks && verbose)
            verboseOutput() << "processing in blocks" << endl;

        size_t nr_remaining = nr_F;

        size_t nr_block = 0;

        while (nr_remaining > 0) {
            nr_block++;

            size_t block_size = min((long)MaxBlocksize, (long)nr_remaining);

            auto F = OldFaces.begin();

            size_t kkpos = 0;
            bool skip_remaining = false;

            const long VERBOSE_STEPS = 50;
            long step_x_size = block_size - VERBOSE_STEPS;
            size_t total = block_size;

            if (in_blocks && verbose)
                verboseOutput() << nr_block << ": " << flush;

            std::exception_ptr tmp_exception;
#pragma omp parallel for firstprivate(kkpos, F) schedule(dynamic) if (block_size > 1)
            for (size_t kk = 0; kk < block_size; ++kk) {
                if (skip_remaining)
                    continue;

                if (verbose && block_size >= ReportBound) {
#pragma omp critical(VERBOSE)
                    while ((long)(kk * VERBOSE_STEPS) >= step_x_size) {
                        step_x_size += total;
                        verboseOutput() << "." << flush;
                    }
                }

                try {
                    INTERRUPT_COMPUTATION_BY_EXCEPTION

                    for (; kk > kkpos; kkpos++, F++)
                        ;
                    for (; kk < kkpos; kkpos--, F--)
                        ;
                    
                    if(F->second.dead)
                        continue;
 
                    F->second.compute(*this, d, F->first, (F->second).tree_size);

                } catch (const std::exception&) {
                    tmp_exception = std::current_exception();
                    skip_remaining = true;
#pragma omp flush(skip_remaining)
                }
            }  // parallel for kk

            if (verbose && block_size >= ReportBound)
                verboseOutput() << endl;

            for (size_t i = 0; i < block_size; ++i)
                OldFaces.erase(OldFaces.begin());

            nr_remaining -= block_size;

        }  // while nr_remaining >0

        OldFaces.swap(NewFaces);
        NewFaces.clear();

        OldNrFacetsContainingGen.swap(NewNrFacetsContainingGen);
        for (size_t i = 0; i < nr_gens; ++i)
            NewNrFacetsContainingGen[i] = 0;

        d--;

    }  // while

    if (verbose) {
        verboseOutput() << "Mult (before NoGradingDenom correction) " << multiplicity << endl;
        verboseOutput() << "Mult (float) " << std::setprecision(12) << mpq_to_nmz_float(multiplicity) << endl;
        if(!exploit_automorphisms){
            verboseOutput() << "Full tree size (modulo 2^64)" << tree_size << endl;
            verboseOutput() << "Number of descent steps " << descent_steps << endl;
            verboseOutput() << "Determinants computed " << nr_simplicial << endl;
            verboseOutput() << "Number of faces in descent system " << system_size << endl;
        }
    }
}

template <typename Integer>
bool DescentSystem<Integer>::set_verbose(bool onoff) {
    bool old_verbose = verbose;
    verbose = onoff;
    return old_verbose;
}

template <typename Integer>
mpq_class DescentSystem<Integer>::getMultiplicity() {
    return multiplicity;
}

template class DescentFace<long>;
template class DescentFace<long long>;
template class DescentFace<mpz_class>;

template class DescentSystem<long>;
template class DescentSystem<long long>;
template class DescentSystem<mpz_class>;

}  // namespace libnormaliz
