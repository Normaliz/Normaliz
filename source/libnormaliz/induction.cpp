/*
 * Normaliz
 * Copyright (C) 2007-2022  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#include <fstream>
#include <sstream>

#include "libnormaliz/cone.h"
#include "libnormaliz/induction.h"

namespace libnormaliz{
using std::vector;
using std::string;
using std::list;
using std::ifstream;
using std::ofstream;

template<typename Integer>
Induction<Integer>::Induction(){

}

template<typename Integer>
Induction<Integer>::Induction(const vector<Integer>& fus_type, const vector<key_t>& fus_duality ,
                              const vector<Integer>& FusRing, bool verb){

    verbose = verb;

    if(verbose)
        verboseOutput() << "Preparing induction matrices" << endl;

    mult_of_ev_ok = true;

    ImageRing = FusRing;

    FusBasic = FusionBasic();
    FusBasic.fusion_rank = fus_type.size();
    fusion_rank = FusBasic.fusion_rank;
    /* for(auto& f:fus_type){
        fusion_type_string += to_string(f) + " ";
    }
    FusBasic.fusion_type_string = fusion_type_string; */
    duality = fus_duality;
    fusion_type = fus_type;
    FusBasic.duality = duality;
    FusBasic.fusion_type = fusion_coincidence_pattern(fus_type);
    FusComp = FusionComp<Integer>(FusBasic);
    FusComp.make_CoordMap();
    Tables = FusComp.make_all_data_tables(FusRing);
    // cout << FusRing;
    /* for(auto& T: Tables)
        T.debug_print();*/

    FPdim = 0;
    for(auto& f: fus_type)
        FPdim += f*f;

    FPSquare = FPdim * FPdim;
    for(Integer t = 1; t <= FPdim; ++t){

        if(FPdim % t == 0)
            divisors.push_back(t);
    }

    // cout << FPdim << " " << divisors;

    /* for(auto& T: Tables)
        T.debug_print('$'); */

    EVMat.resize(fusion_rank, fusion_rank);
    for(size_t s =0; s < fusion_rank; ++s){
        for(size_t l = 0; l < fusion_rank; ++l){
            Integer S = 0;
            for(size_t t = 0; t < fusion_rank; ++t){
                for(size_t k = 0; k < fusion_rank; ++k){
                    // cout << t << " " << k << " " << N(t, duality[t], k) << " " << N(k,l,s)<< endl;
                    S += N(t, duality[t], k)*N(k,l,s);
                }
            }
            EVMat[s][l] = S;
        }
    }
    // EVMat.debug_print();

    MultEV.resize(divisors.size());

    size_t sum_mult = 0;
    if(verbose)
        verboseOutput() << "eigenvalues and their multiplicities" << endl;
    for(size_t i = 0; i< divisors.size(); ++i){
        MultEV[i] = EVMat.mult_of_eigenvalue(divisors[i]);
        if(MultEV[i] > 0 && verbose)
            verboseOutput() << divisors[i] << " mult " <<  MultEV[i] << endl;
        sum_mult += MultEV[i];
    }

    if(sum_mult < fusion_rank){
        if(verbose)
            verboseOutput() << "Sum of multiplicities of eigenvalues dividing FPdim < fusion_rank" << endl;
        mult_of_ev_ok = false;
        return;
    }

   Bounds.resize(fusion_rank, fusion_rank);

    for(size_t j = 0; j< fusion_rank; ++j){
        for(size_t k = 0; k < fusion_rank; ++k){
            Integer S =0;
            for (size_t s = 0; s < fusion_rank; ++s){
                for(size_t t = 0; t < fusion_rank; ++t )
                    S += N(t,j,s)*N(s, duality[t],k);
            }
            Bounds[j][k] = S;
        }
    }
    // Bounds.debug_print('&');

    if(verbose)
        verboseOutput() << "Computing reporesentations of divisors of FPdim in terms of type" << endl;

    HighRepresentations.resize(0, fusion_rank);

    for(auto& t: divisors){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if(t == FPdim)
            continue;
        LowRepresentations[t].resize(0, fusion_rank);
        Matrix<Integer> InhomEqu(1,fusion_rank +1);
        for(size_t i = 0; i < fusion_rank; ++i)
            InhomEqu[0][i] = fus_type[i];
        InhomEqu[0].back() = - t;
        Matrix<Integer> NonNeg(InhomEqu.nr_of_columns() - 1);
        Matrix<Integer> NeutralInEqu(1,fusion_rank +1);
        NeutralInEqu[0][0] = -1;
        NeutralInEqu[0].back() = 1;
        Cone<Integer> RepCone(Type::inhom_equations, InhomEqu,Type::inequalities, NonNeg,
                              Type::inhom_inequalities, NeutralInEqu);
        RepCone.setVerbose(false);
        Matrix<Integer> Reps = RepCone.getLatticePointsMatrix();
        size_t count_low = 0, count_high = 0;
        for(size_t i = 0; i < Reps.nr_of_rows(); ++i){
            bool too_large = false;
            for(size_t j = 0; j < fusion_rank; ++j){
                for(size_t k = j; k < fusion_rank; ++k){
                    if(Reps[i][j]*Reps[i][k] > Bounds[j][k]){
                        too_large = true;
                        break;
                    }
                }
                if(too_large)
                    break;
            }
            if(too_large)
                continue;
            vector<Integer> new_rep = Reps[i];
            new_rep.resize(new_rep.size() - 1);

            if(Reps[i][0] == 1){
                LowRepresentations[t].append(new_rep);
                count_low++;
            }
            if(Reps[i][0] == 0){
                HighRepresentations.append(new_rep);
                count_high++;
            }

        }
        if(verbose)
            verboseOutput() << "divisor " << t << " has " << count_low <<" low reps and " << count_high << " high" << endl;

       // LowRepresentations[t].debug_print('$');
    }
    // HighRepresentations.debug_print('%');


    /* cout << "Numbers of reporesentations low" << endl;
    for(auto& t: divisors){
        cout << t << " " << LowRepresentations[t].nr_of_rows() << endl;
    } */
    // cout << "high" << endl;


}

template<>
Induction<renf_elem_class>::Induction(const vector<renf_elem_class>& fus_type, const vector<key_t>& fus_duality , const vector<renf_elem_class>& FusRing, bool verbose){

    assert(false);
}

template<typename Integer>
Integer Induction<Integer>::N(const key_t i, const key_t j, const key_t k){
    return Tables[i][j][k];
}

template<typename Integer>
void Induction<Integer>::start_low_parts(){

    if(!mult_of_ev_ok)
        return;

    if(verbose)
        verboseOutput() << "Computing low parts (rows 1,...,r before sorting)" << endl;

    vector<Integer> EV;
    for(size_t i = 0; i < divisors.size(); ++i){
        for(size_t j = 0; j < MultEV[i]; ++j)
            EV.push_back(divisors[i]);
    }
    for(auto& t: EV)
        low_m.push_back(FPdim / t);
    sort(low_m.begin(), low_m.end());

    if(verbose){
        size_t upperbound = 1;
        for(size_t i = 0; i < fusion_rank; ++i)
            upperbound *= LowRepresentations[low_m[i]].nr_of_rows();
        verboseOutput() << "Number of low parts <= " << upperbound << endl;
    }


    Matrix<Integer> matrix_so_far(0, fusion_rank);
    Matrix<Integer> bounds_so_far(fusion_rank, fusion_rank);
    build_low_matrices(matrix_so_far, bounds_so_far);

    vector<Matrix<Integer> >  OrderedLowParts;
    for(auto& M: LowParts){
        /* for(size_t i = 0; i < fusion_rank; ++i)
            cout << v_scalar_product(M[i], fusion_type) << " ";
        cout << endl;*/
        // M.debug_print('C');
        bool ordered = true;
        for(size_t i = 0; i < M.nr_of_rows() - 1; ++i){
            if(low_m[i] == low_m[i+1]){
                if(M[i] > M[i+1]){
                    // cout << "ord " << i << " " << low_m[i] << " " << low_m[i+1] << endl;
                    ordered = false;
                    // break;
                }
            }
        }
        // cout << "OOOOOOO " << ordered << endl;
        if(ordered)
            OrderedLowParts.push_back(M);
    }

    // cout << "Old " << LowParts.size() << " New " << OrderedLowParts.size() << endl;
    swap(LowParts, OrderedLowParts);

    /* if(verbose)
        verboseOutput() << "Found " << LowParts.size() << " low parts"  << endl;*/

    /* for(auto& M: LowParts)
        M.debug_print('P'); */
    /* for(auto& M: LowPartsBounds)
        M.debug_print('&'); */

}


template<typename Integer>
void Induction<Integer>::build_low_matrices(Matrix<Integer> matrix_so_far, Matrix<Integer> bounds_so_far){

    key_t step = matrix_so_far.nr_of_rows();

    Integer m_step = low_m[step];
    for(size_t i = 0; i < LowRepresentations[m_step].nr_of_rows(); ++i){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        Matrix<Integer> NewMatrix = matrix_so_far;

        vector<Integer> cand_extension = LowRepresentations[m_step][i];
        bool potential_low_part = true;

        Matrix<Integer> check_bounds = bounds_so_far;

        for(size_t j = 0; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                check_bounds[j][k] +=  cand_extension[j] * cand_extension[k];
                if(check_bounds[j][k] > Bounds[j][k]){
                   potential_low_part = false;
                   break;
                }
                if(!potential_low_part)
                    break;
            }
        }

        if(!potential_low_part)
            continue;
        NewMatrix.append(cand_extension);

        if(NewMatrix.nr_of_rows() == fusion_rank){
            /* if(verbose && LowParts.size() % 10000 == 0 && LowParts.size() > 0)
                verboseOutput() << LowParts.size() << " low parts" << endl;
            LowParts.push_back(NewMatrix);*/

            from_low_to_full(NewMatrix);
            continue;
        }

        build_low_matrices(NewMatrix, check_bounds);
    }
}


template<typename Integer>
Matrix<Integer> Induction<Integer>::make_allowed_transpositions(Matrix<Integer> FusionMap){

    vector<Integer> type = FusionMap.MxV(fusion_type);

    Matrix<Integer> AllowedTranspositions(0,2);

    size_t rank_ZR = FusionMap.nr_of_rows();
    for(long i = 1; i < rank_ZR; ++i){
        for(long j = i; j < rank_ZR; ++j){
            if(type[i] != type[j])
                continue;
            bool allowed =true;
            for(long k = 0; k < fusion_rank; ++k){
                if(FusionMap[i][duality[k]] != FusionMap[j][k]){
                    allowed = false;
                    break;
                }
            }
            if(allowed){
                vector<long> Help = {i,j};
                vector<Integer> HelpInt;
                convert(HelpInt, Help);
                AllowedTranspositions.append(HelpInt);
            }

        }
    }
    return AllowedTranspositions;
}

template<typename Integer>
void Induction<Integer>::from_low_to_full(const Matrix<Integer>& ThisLowPart){


    if(!mult_of_ev_ok)
        return;

    if(verbose)
        verboseOutput() << "Extending low part to full induction matrices" << endl;


    Integer FPdim_so_far =0;
    for(size_t j = 0; j < fusion_rank; ++j)
        FPdim_so_far += low_m[j] * low_m[j];

    Integer MinusOne = -1;

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        // ThisLowPart.debug_print('#');
        // LowPartsBounds[iii].debug_print('B');

        Matrix<Integer> MinusLowBound(fusion_rank, fusion_rank);
        for(size_t i = 0; i< fusion_rank; ++i){
            vector<Integer> cand_extension = ThisLowPart[i];
            for(size_t j = 0; j < fusion_rank; ++j){
                for(size_t k = j; k < fusion_rank; ++k){
                    MinusLowBound[j][k] +=  cand_extension[j] * cand_extension[k];
                }
            }
        }
        MinusLowBound.scalar_multiplication(MinusOne);

        // Bounds.debug_print();

        Matrix<Integer> Remaining = Bounds.add(MinusLowBound);

        // Remaining.debug_print('R');

        Matrix<Integer> InhomEqu(0, HighRepresentations.nr_of_rows() + 1);

        // HighRepresentations.debug_print('H');

        for(size_t j = 1; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                vector<Integer> this_equ;
                for(size_t i = 0; i < HighRepresentations.nr_of_rows(); ++i){
                        this_equ.push_back(HighRepresentations[i][j] * HighRepresentations[i][k]);
                }
                this_equ.push_back(-Remaining[j][k]);
                // cout << j << " " << k << " " << "ttt " <<  this_equ;
                InhomEqu.append(this_equ);
            }
        }

        vector<Integer> this_equ;
        for(size_t i = 0; i < HighRepresentations.nr_of_rows(); ++i){
            Integer m_new = v_scalar_product(HighRepresentations[i], fusion_type);
            this_equ.push_back(m_new * m_new);
        }
        // cout << "FFFFFF " << FPdim_so_far << " " << -FPSquare << " " << (-FPSquare + FPdim_so_far) <<  endl;
        this_equ.push_back(-FPSquare + FPdim_so_far);
        InhomEqu.append(this_equ);

        // InhomEqu.debug_print('&');

        // we use inequalitiesqualities to avoid coordinate transformation
        Matrix<Integer> Copy = InhomEqu;
        Copy.scalar_multiplication(MinusOne);
        InhomEqu.append(Copy);

        Matrix<Integer> Unit(HighRepresentations.nr_of_rows());
        Cone<Integer> C(Type::inhom_inequalities, InhomEqu, Type::inequalities, Unit);
        C.setVerbose(false);
        C.compute(ConeProperty::HilbertBasis, ConeProperty::Projection);
        Matrix<Integer> LP = C.getLatticePointsMatrix();

        if(talkative && LP.nr_of_rows() == 0)
            verboseOutput() << "N" << endl;

        for(size_t k= 0; k< LP.nr_of_rows(); ++k){
            Matrix<Integer> IndMat = ThisLowPart;
            // cout << "CCC " << ThisLowPart.nr_of_columns() << endl;
            for(size_t j = 0; j < LP.nr_of_columns() - 1; ++j){
                Integer count = 0;
                while(count < LP[k][j]){
                    // cout << "HHH " << HighRepresentations[j].size() << endl;
                    IndMat.append(HighRepresentations[j]);
                    count = count +1;
                }

            }
            InductionMatrices.push_back(IndMat);
            if(verbose)
                IndMat.debug_print('I');
            if(talkative)
                verboseOutput() << "I" << endl;
        }
}

template<typename Integer>
void Induction<Integer>::augment_induction_matrices(){
    set<vector<vector<Integer> > > EquivHelp;
    vector<Matrix<Integer> > Representatives;
    for(auto& M: InductionMatrices){

        vector<vector<Integer> >  N = M.get_elements();
        sort(N.begin(), N. end());
        if(EquivHelp.find(N) == EquivHelp.end()){
            Representatives.push_back(M);
           EquivHelp.insert(N);
        }
    }
    swap(InductionMatrices, Representatives);

    for(auto& M: InductionMatrices){
        //we rteorder the rows by increasing m_i

        map<Integer, vector<vector<Integer> > > SortHelp;
        for(size_t i = 0; i < M.nr_of_rows(); ++i){
            Integer m = v_scalar_product(M[i], fusion_type);
            SortHelp[m].push_back(M[i]);
        }

        M.resize(0,fusion_rank);
        for(auto& T: SortHelp){
            for(auto& v: T.second)
                M.append(v);
        }
    }

    if(verbose)
        verboseOutput() << InductionMatrices.size() << " induction matrices found" << endl;

    vector<Matrix<Integer> > InductionMatWithType;
    // cout << "FFFFFFFFFFF " << ImageRing;
    InductionMatWithType.push_back(ImageRing);
    for(auto& M: InductionMatrices){
        InductionMatWithType.push_back(M);
        vector<Integer> type = M.MxV(fusion_type);
        InductionMatWithType.push_back(Matrix<Integer>(type));
        Matrix<Integer> AllowedTranspositions = make_allowed_transpositions(M);
        InductionMatWithType.push_back(AllowedTranspositions);
    }

    swap(InductionMatrices, InductionMatWithType);

    /* for(auto& M: InductionMatrices)
        M.debug_print('$');*/
}


template<typename Integer>
void Induction<Integer>::compute(){

    if(!mult_of_ev_ok)
        return;
    start_low_parts();
    augment_induction_matrices();
}

/* Not used anamore
template<typename Integer>
void Induction<Integer>::extend_matrix(Matrix<Integer> matrix_so_far, key_t rep_index,
                       Matrix<Integer> bounds_so_far, Integer FPdim_so_far){

    if(rep_index >= HighRepresentations.nr_of_rows())
        return;

    vector<Integer> cand_ext = HighRepresentations[rep_index];
    Integer new_m = v_scalar_product(cand_ext, fusion_type);

    while(true){

        extend_matrix(matrix_so_far, rep_index +1, bounds_so_far, FPdim_so_far);
        FPdim_so_far += new_m* new_m;
        if(FPdim_so_far > FPSquare){
            return;
        }
        bool potential_solution = true;
        if(FPdim_so_far < FPSquare)
            potential_solution = false;

        for(size_t j = 0; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                bounds_so_far[j][k] +=  cand_ext[j] * cand_ext[k];
                if(bounds_so_far[j][k] > Bounds[j][k]){
                    return;
                }
                if(bounds_so_far[j][k] < Bounds[j][k])
                    potential_solution = false;
            }
        }
        matrix_so_far.append(cand_ext);

        if(potential_solution){
            InductionMatrices.push_back(matrix_so_far);
            // matrix_so_far.debug_print('$');
            return;
        }

    }
}
*/

template class Induction<mpz_class>;
template class Induction<long long>;
template class Induction<long>;
#ifdef ENFNORMALIZ
template class Induction<renf_elem_class>;
#endif



} // namespace
