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
    simplicial = false;
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
DescentSystem<Integer>::DescentSystem(const Matrix<Integer>& Gens_given,
                                      const Matrix<Integer>& SuppHyps_given,
                                      const vector<Integer>& Grading_given) {
    descent_steps = 0;
    tree_size = 0;
    nr_simplicial = 0;
    system_size = 0;
    exploit_automorphisms = true;

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

// One could return DecentFaces containing all the computed data that are not const in the parameter lists
// One could think about a return object, say class Pyramic<Integer>
template <typename Integer>
void DescentFace<Integer>::compute(DescentSystem<Integer>& FF,
                                   size_t dim,  //  dim of *this
                                   const dynamic_bitset& own_facets, // indicates the supphyps of which *this is the intersection
                                   vector<key_t>& mother_key,  // will indicate the extreme rays of *this
                                        // used after return from this function to count the number of faces containing an extreme ray
                                        // these data are used in this function for the choice of the optimal vertex
                                   vector<dynamic_bitset>& opposite_facets, // for each opposite facet ALL the supphyps defining it as intersection 
                                                                            // used as a signature in the descent system
                                   vector<key_t>& CuttingFacet, // the indices of facets opposite to selected extreme ray (not unique)
                                   vector<Integer>& heights,  // the heights of selected extreme ray over selected facets
                                   vector<dynamic_bitset>& FacetsOfFace, // for each facet of *this given back a selection of global support
                                                                        // hyperplanes cutting out its facets(one p0er facet, nit unique)
                                   key_t& selected_gen) { // index of selected extreme ray
    long omp_start_level = omp_get_level();

    mother_key.clear();
    opposite_facets.clear();
    heights.clear();
    CuttingFacet.clear();
    if(FF.exploit_automorphisms)
        FacetsOfFace.clear();

    size_t nr_supphyps = FF.nr_supphyps;
    size_t nr_gens = FF.nr_gens;

    size_t d = dim;
    
    // cout << "FFFFFFFFFF " << coeff << endl;

    dynamic_bitset GensInd(nr_gens); // find the extreme rays of *this, indicated by GensInd
    GensInd.set();                   // by intersecting the indicators of own_facets
    // vector<key_t> own_facets_key;
    for (size_t i = 0; i < nr_supphyps; ++i) {  // find Gens in this
        if (own_facets[i] == true) {
            GensInd = GensInd & FF.SuppHypInd[i];
        }
    }

    for (size_t i = 0; i < nr_gens; ++i)
        if (GensInd[i])
            mother_key.push_back(i);

    Matrix<Integer> Gens_this;

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
                Gens_this = FF.Gens.submatrix(selection);
                rk = Gens_this.row_echelon();
                nr_selected *= 2;
            }
            if (rk < dim) {
                Gens_this = FF.Gens.submatrix(mother_key);
                Gens_this.row_echelon();
            }
        } catch (const ArithmeticException& e) {
            Gens_this = FF.Gens.submatrix(mother_key);
            Gens_this.row_echelon();
        }
    }
    else {
        Gens_this = FF.Gens.submatrix(mother_key);
        Gens_this.row_echelon();
    }

    bool must_saturate = false;

    for (size_t i = 0; i < Gens_this.nr_of_rows(); ++i) {
        for (size_t j = i; j < FF.dim; ++j) {
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

    Sublattice_Representation<Integer> Sublatt_this;
    if (must_saturate)
        Sublatt_this = Sublattice_Representation<Integer>(Gens_this, true, false);  //  take saturation, no LLL

    // Now we find the potential facets of *this.

    dynamic_bitset facet_ind(mother_key.size());    // lists Gens, local variable for work
    map<dynamic_bitset, dynamic_bitset> FacetInds;  // potential facets, map from gens(potential facet)
                                                    //  to set of supphyps(C) containing these gens
                                                    // reference for gens(potential facet) is the selection via mother_key
    map<dynamic_bitset, key_t> CutOutBy;            // the facet citting it out

    map<dynamic_bitset, vector<key_t> > SimpKeys;  // generator keys for simplicial facets
    map<dynamic_bitset, vector<bool> > SimpInds;   // alternative: generator indices for simplicial facets (if less memory needed)
    bool ind_better_than_keys = (dim * 64 > FF.nr_gens); // decision between the alternatives

    for (size_t i = 0; i < nr_supphyps; ++i) {
        if (own_facets[i] == true)  // contains *this
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
            FacetInds[facet_ind] = dynamic_bitset(0);  // don't need support hyperplanes
            CutOutBy[facet_ind] = FF.nr_supphyps + 1;  // signalizes "simplicial facet"
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
            FacetInds[facet_ind] = own_facets;
            FacetInds[facet_ind][i] = true;  // plus the facet cutting out facet_ind
            CutOutBy[facet_ind] = i;         // memorize the facet that cuts it out
        }
    }

    // if we don't have the coordinate transformation and there is a simplicial facet, we must make it
    if (!must_saturate && (SimpKeys.size() > 0 || SimpInds.size() > 0))
        Sublatt_this = Sublattice_Representation<Integer>(Gens_this, true, false);  //  take saturation, no LLL

    if (d < FF.dim && !FF.SimplePolytope) {  // now we select the true facets of *this
        auto G = FacetInds.end();            // by taking those with a maximal set of gens
        for (--G; G != FacetInds.begin(); --G) {
            for (auto F = FacetInds.begin(); F != G;) {
                if (F->first.is_subset_of(G->first))
                    F = FacetInds.erase(F);
                else
                    ++F;
            }
        }
    }

    // At this point we know the facets of *this.
    // The map FacetInds assigns the set of containing SuppHyps to the facet_ind(Gens).
    // The set of containing SuppHyps is a unique signature as well.

    // Now we want to find the generator with the lrast number opf opposite facets(*this)
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
    libnormaliz::key_t m_ind = 0;

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

    selected_gen = mother_key[m_ind];  // this is the selected generator
    vector<Integer> embedded_selected_gen;
    if (must_saturate)
        embedded_selected_gen = Sublatt_this.to_sublattice(FF.Gens[selected_gen]);

    // now we must find the facets opposite to thge selected generator

    vector<Integer> embedded_supphyp;
    Integer ht;

    auto G = FacetInds.begin();
    for (; G != FacetInds.end(); ++G) {
        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if ((G->first)[m_ind] == false && CutOutBy[G->first] != FF.nr_supphyps + 1) {  // is opposite and not simplicial
            opposite_facets.push_back(G->second);
            if (must_saturate) {
                embedded_supphyp = Sublatt_this.to_sublattice_dual(FF.SuppHyps[CutOutBy[G->first]]);
                ht = v_scalar_product(embedded_selected_gen, embedded_supphyp);
            }
            else {
                embedded_supphyp = Gens_this.MxV(FF.SuppHyps[CutOutBy[G->first]]);
                Integer den = v_make_prime(embedded_supphyp);
                ht = v_scalar_product(FF.Gens[selected_gen], FF.SuppHyps[CutOutBy[G->first]]) / den;
            }
            // cout <<  ht << endl;
            heights.push_back(ht);
            // coeff *= convertTo<mpz_class>(ht);
            // coeff /= FF.GradGens_mpz[selected_gen];
            
            CuttingFacet.push_back(CutOutBy[G->first]);
            
            if(FF.exploit_automorphisms){
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
                
                vector<bool> TheFacets;
                maximal_subsets(Intersections, TheFacets);
                FacetsOfFace.push_back(bool_to_bitset(TheFacets));   // indicates exactly one support hyperplane cutting out our facet 
                                                                    // from *this
            }
        }
    }

    if (SimpKeys.size() > 0 || SimpInds.size() > 0) {
        G = FacetInds.begin();
        size_t loop_length = FacetInds.size();
        size_t fpos = 0;
        bool skip_remaining = false;
        vector<mpq_class> thread_mult(omp_get_max_threads(), 0);
        Matrix<Integer> Embedded_Gens(d, d);
        Matrix<Integer> Gens_this(d, FF.dim);

        std::exception_ptr tmp_exception;

#pragma omp parallel for firstprivate(G, fpos, Embedded_Gens, Gens_this)
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

                if ((G->first)[m_ind] == false && CutOutBy[G->first] == FF.nr_supphyps + 1) {  // is opposite and simplicial
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
                    // #pragma omp critical(ADD_MULT)
                    // FF.multiplicity+=multiplicity*coeff;
                    thread_mult[tn] += multiplicity;
#pragma omp atomic
                    FF.nr_simplicial++;
#pragma omp atomic
                    FF.tree_size += tree_size;
                }

            } catch (const std::exception&) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
        }

        if (!(tmp_exception == 0))
            std::rethrow_exception(tmp_exception);

        mpq_class local_multiplicity = 0;
        for (const auto& j : thread_mult)
            local_multiplicity += j;
#pragma omp critical(ADD_MULT)
        FF.multiplicity += local_multiplicity * coeff;
    }
}

template <typename Integer>
void DescentSystem<Integer>::collect_old_faces_in_iso_classes(){
    
    if(OldFaces.size() <= 1) // nothingt to do here
        return;
    
    // Isomorphism_Classes<Integer> Isos(AutomParam::rational_dual);
    map< IsoType<Integer>, DescentFace<Integer>* , IsoType_compare<Integer> > Isos;
    
    size_t nr_F = OldFaces.size();
    auto F = OldFaces.begin();
    size_t kkpos = 0;
    std::exception_ptr tmp_exception;
    bool skip_remaining = false;
    
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
            auto G = Isos.find(IT); 
            if(G != Isos.end()){
                F->second.dead = true; // to be skipped in descent
                mpz_class index_source = convertTo<mpz_class>(IT.index);
                mpz_class index_traget = convertTo<mpz_class>(G->first.index);
                mpq_class corr = index_source;
                corr /= index_traget;
                // cout << "--------------------------" << endl;
                // cout << "Index Source " << index_source << " Index Target " << index_traget << " CORR " << corr << endl;
                G->second->coeff += corr*F->second.coeff;
                // cout << "coeff Source " << F->second.coeff << " Coeff Target " << G->second->coeff << endl;
                // cout << "--------------------------" << endl; */
                if(corr!=1)
                    exit(0);
            }
            else{
                // cout << "========================" << endl;
                Isos[IT]= &(F->second);
                // cout << "New New New " << "Index " << IT.index << " Coeff " << F->second.coeff << endl;
                //IT.getCanType().pretty_print(cout);                
                // cout << "========================" << endl;
            }

        } catch (const std::exception&) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
#pragma omp flush(skip_remaining)
        }
    }  // parallel for kk
    
    if (!(tmp_exception == 0))
        std::rethrow_exception(tmp_exception);
    
    cout << "Iso types " << Isos.size() << endl;
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
    long d = (long)dim;

    while (!OldFaces.empty()) {
        size_t nr_F = OldFaces.size();
        system_size += nr_F;
        if (verbose)
            verboseOutput() << "Descent from dim " << d << ", size " << nr_F << endl;
        
        if(exploit_automorphisms)
            collect_old_faces_in_iso_classes();

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

            vector<key_t> mother_key;
            mother_key.reserve(nr_gens);
            vector<dynamic_bitset> opposite_facets;
            opposite_facets.reserve(nr_supphyps);
            vector<key_t> CuttingFacet;
            CuttingFacet.reserve(nr_supphyps);
            vector<Integer> heights;
            heights.reserve(nr_supphyps);
            key_t selected_gen = 0;
            
            vector<dynamic_bitset> FacetsOfFace;

            std::exception_ptr tmp_exception;
#pragma omp parallel for firstprivate(kkpos, F, mother_key, opposite_facets, CuttingFacet, heights, FacetsOfFace, selected_gen) \
    schedule(dynamic) if (block_size > 1)
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
                    // cout << "Rechne " << endl;
                    
                    F->second.compute(*this, d, F->first, mother_key, opposite_facets, CuttingFacet, heights, FacetsOfFace, selected_gen);
                    if (F->second.simplicial)
                        continue;

                    auto G = opposite_facets.begin();
                    // mpz_class deg_mpz=convertTo<mpz_class>(GradGens[selected_gen]);
                    // mpq_class divided_coeff=(F->second).coeff/deg_mpz;
                    mpq_class divided_coeff = (F->second).coeff / GradGens_mpz[selected_gen];
                    size_t j = 0;
                    for (; G != opposite_facets.end(); ++G) {
                        auto H = NewFaces.begin();
                        bool inserted = false;
#pragma omp critical(INSERT)
                        {
                            H = NewFaces.find(*G);
                            if (H == NewFaces.end()) {
                                H = NewFaces.insert(NewFaces.begin(), make_pair(*G, DescentFace<Integer>()));
                                inserted = true;
                            }
                        }
                        if (inserted) {
                            if(exploit_automorphisms)
                                H->second.FacetsOfFace = FacetsOfFace[j];
                            for (unsigned int& i : mother_key)
                                if (SuppHypInd[CuttingFacet[j]][i])
#pragma omp atomic
                                    NewNrFacetsContainingGen[i]++;
                        }
                        mpq_class dc = divided_coeff * convertTo<mpz_class>(heights[j]);
#pragma omp critical(ADD_COEFF)
                        {
                            (H->second).coeff += dc;
                            (H->second).tree_size += (F->second).tree_size;
                        }
                        ++j;
                        descent_steps++;
                    }

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
        verboseOutput() << "Full tree size (modulo 2^64)" << tree_size << endl;
        verboseOutput() << "Number of descent steps " << descent_steps << endl;
        verboseOutput() << "Determinants computed " << nr_simplicial << endl;
        verboseOutput() << "Number of faces in descent system " << system_size << endl;
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
