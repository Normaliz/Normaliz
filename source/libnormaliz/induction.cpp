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
#include "libnormaliz/normaliz_exception.h"
#include "libnormaliz/matrix.h"

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

    auto FusBasic = FusionBasic();
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
    /* for(auto& T: Tables);
        T.debug_print();*/

    test_commutativity();

    if(using_renf<Integer>()){
        size_t r = fusion_rank -1;
        kkk = N(r,r,r);
        d_plus = fusion_type.back();
        d_minus = kkk - d_plus;
    }

    FPdim = 0;
    for(size_t i = 0; i < fus_type.size(); ++i){
        FPdim += fusion_type[i]*fusion_type[i];
    }

    if(verbose){
        verboseOutput() << "Type " << fusion_type;
        verboseOutput() << "FPdim " << FPdim;
        verboseOutput() << endl;
    }

    FPSquare = FPdim * FPdim;

    // For the computation of eigenvalues and multiplicities
    size_t EVMat_size = fusion_rank;

    EVMat.resize(EVMat_size, EVMat_size);
    for(size_t s =0; s < EVMat_size; ++s){
        for(size_t l = 0; l < EVMat_size; ++l){
            Integer S = 0;
            for(size_t t = 0; t < EVMat_size; ++t){
                for(size_t k = 0; k < EVMat_size; ++k){
                    // cout << t << " " << k << " " << N(t, duality[t], k) << " " << N(k,l,s)<< endl;
                    S += N(t, duality[t], k)*N(k,l,s);
                }
            }
            EVMat[s][l] = S;
        }
    }

    EVMat.debug_print('E');

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

    Bounds.debug_print('B');


    make_divisors();

    if(commutative){
        codegrees_and_mult_commutative();
    }
    else{
        codegrees_and_mult_noncommutative();
    }

    if(verbose)
        verboseOutput() << "Computing reporesentations of divisors of FPdim in terms of type" << endl;

    HighRepresentations.resize(0, fusion_rank);

    for(auto& t: candidates_m_i){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if(t == FPdim)
            continue;
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
        // Reps.debug_print('R');
        size_t count_high = 0;
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

            if(Reps[i][0] == 0){
                HighRepresentations.append(new_rep);
                count_high++;
            }

        }
        if(verbose)
            verboseOutput() << "candidate " << t << " has " << count_high << " high" << " representations" << endl;

       // LowRepresentations[t].debug_print('$');
    }
}

template<typename Integer>
Integer Induction<Integer>::conjugate(const Integer& val){
    assert(false);
    return 0;
}

template<typename Integer>
bool Induction<Integer>::is_algebraic_integer(const Integer& val){
    assert(false);
    return true;
}

#ifdef ENFNORMALIZ
template<>
renf_elem_class Induction<renf_elem_class>::conjugate(const renf_elem_class& val){
    vector<mpz_class> num = val.num_vector();
    // cout << "vvvvvvvvv " << val << endl;
    assert(num.size() <= 2);
    if(num.size() <= 1)
        return val;
    renf_elem_class ret = (num[1]*d_minus + num[0])/val.den();
    return ret;
}


template<>
bool Induction<renf_elem_class>::is_algebraic_integer_old(const renf_elem_class& val){
    renf_elem_class conj = conjugate(val);
    if (!(val + conj).is_integer())
        return false;
    if(!(val*conj).is_integer())
        return false;
    return true;
}


template<>
bool Induction<renf_elem_class>::is_algebraic_integer(const renf_elem_class& val){
    if(val.is_rational()){
        if(val.is_integer())
            return true;
        return false;
    }
    // cout << "true renf" << endl;
    vector<renf_elem_class> powers;
    powers.push_back(1);
    vector<vector<mpz_class> > numerators;
    numerators.push_back(powers.back().num_vector());
    size_t max_deg = 0;
    while(true){
        powers.push_back(powers.back() * val);
        numerators.push_back(powers.back().num_vector());
        if(numerators.back().size() - 1 > max_deg)
            max_deg = numerators.back().size() - 1;
        for(auto& v: numerators)
            v.resize(max_deg + 1);
        Matrix<mpz_class> LinDep(numerators);
        // cout << powers;
        // cout << "rank " << LinDep.rank() << endl;
        if(LinDep.rank() == LinDep.nr_of_rows()) // powers still linearly independent
            continue;
        // now gedree found
        // get 5rid of denominators
        // LinDep.debug_print();
        mpz_class D = 1;
        for(auto& p: powers)
            D = libnormaliz::lcm(D, p.den());
        for(size_t i = 0; i < LinDep.nr_of_rows(); ++i){
            mpz_class fact = D/powers[i].den();
            v_scalar_multiplication(LinDep[i], fact);
        }
        // cout << "Den " << D << endl;
        // find linear elation on powers
        Matrix<mpz_class> T = LinDep.transpose();
        // T.debug_print('T');
        Matrix<mpz_class> Sol = T.kernel(false);
        assert(Sol.nr_of_rows() ==1); // since the rank pf the matrix is 1 nr_of_columns -1
        vector<mpz_class> SolVec = Sol[0];
        // cout << "Sol " << SolVec << endl;
        // now we must check whether we can divide by the coefficient
        // of the highesdt power.
        if(Iabs(SolVec.back()) == 1)
            return true;
        else
            return false;

    }
}
#endif

template<typename Integer>
void Induction<Integer>::make_divisors(){
    for(Integer t = 1; t <= FPdim; ++t){
        if(FPdim % t == 0){
            divisors.push_back(t);
            candidates_m_i.push_back(FPdim/t);
        }
    }
}
template<>
void Induction<renf_elem_class>::make_divisors(){

    vector<renf_elem_class> h = fusion_type;
    h.push_back(-FPdim);
    renf_elem_class MinusOne = -1;
    v_scalar_multiplication(h, MinusOne);
    Matrix<renf_elem_class> Hyp(h);
    Hyp.debug_print();
    Cone<renf_elem_class> CandCone(Type::inhom_inequalities, Hyp);
    CandCone.setNonnegative();
    Matrix<renf_elem_class> RawCands = CandCone.getLatticePointsMatrix();
    set<renf_elem_class> CandSet;

    for(size_t i = 0; i < RawCands.nr_of_rows(); ++i){
        auto c = RawCands[i];
        renf_elem_class cand = v_scalar_product_vectors_unequal_lungth(c, fusion_type);
        if(cand == 0)
            continue;
        if(CandSet.find(cand) != CandSet.end()) // already tested
            continue;
        CandSet.insert(cand);
        renf_elem_class inv = FPdim/cand;
        // cout << " inv " << inv << endl;
        // cout << "old " << is_algebraic_integer_old(inv) << endl;
        bool whow = is_algebraic_integer(inv);
        //cout << "new " <<whow << endl;
        // cout << "=============================================" << endl;
        if(whow){
            candidates_m_i.push_back(cand);
            divisors.push_back(inv);
        }
    }
}

template<typename Integer>
void Induction<Integer>::test_commutativity(){

    for(size_t i = 0; i < Tables.size(); ++i){
            for(size_t j = i + 1; j < Tables.size(); ++j){
                Matrix<Integer> Prod_1 = Tables[i].multiplication(Tables[j]);
                Matrix<Integer> Prod_2 = Tables[j].multiplication(Tables[i]);
                if(!Prod_1.equal(Prod_2)){
                    commutative = false;
                    return;
                }
            }
    }
    commutative =  true;
}

template<typename Integer>
void Induction<Integer>::codegrees_and_mult_commutative(){

    nr_rows_low_part = fusion_rank;

    size_t sum_mult = 0;
    if(verbose)
        verboseOutput() << "eigenvalues and their multiplicities in the commutative case" << endl;
    for(size_t i = 0; i< divisors.size(); ++i){
        size_t mult_ev = EVMat.mult_of_eigenvalue(divisors[i]);
        sum_mult += mult_ev;
        if(mult_ev > 0){
            EV_mult_n_i[divisors[i]] = make_pair(mult_ev,1);
            if(verbose){
                verboseOutput() << divisors[i] << " mult " <<  EV_mult_n_i[divisors[i]].first << endl;
            }
        }
    }
    if(sum_mult < fusion_rank){
        if(verbose)
            verboseOutput() << "Sum of multiplicities of eigenvalues dividing FPdim < fusion_rank" << endl;
        mult_of_ev_ok = false;
    }
}


template<typename Integer>
void Induction<Integer>::codegrees_and_mult_noncommutative(){

    if(fusion_rank > 8){
        throw BadInputException("Fusion rank must be <= 8 for induction matrices in the noncommutative case!");
    }

    nr_rows_low_part = fusion_rank -3;

    // First we must compute the eigenvalues and multiplicities .
    // Later on they can be distributed in pairs (d_i, n_i)
    // where the eigenvalue is n_id_i and the problem that we may havre
    // n_id_i = n_jd_j for i != j.
    map<Integer,size_t> MultEV_raw;
    set<Integer> already_tested;
    size_t sum_mult = 0;

    for(int i = 1; i <= fusion_rank; ++i){
        for(size_t j = 0; j < divisors.size(); ++j){
            Integer ev_cand = i*divisors[j];
            if(already_tested.find(ev_cand) != already_tested.end())
                continue;
            size_t mult_ev = EVMat.mult_of_eigenvalue(ev_cand);
            if(mult_ev > 0)
                sum_mult+= mult_ev;
            if(mult_ev > 0){
                 MultEV_raw[ev_cand] = mult_ev;
            }
            already_tested.insert(ev_cand);
        }
    }

    if(verbose){
        verboseOutput() <<  "Eigenvalues and their multiplicities" << endl;
        for(auto& mult_this: MultEV_raw){
            verboseOutput() <<  mult_this.first << " multiplicity " << mult_this.second << endl;
        }
    }

    if(sum_mult < fusion_rank){
        if(verbose)
            verboseOutput() << "Sum of multiplicities of eigenvalues dividing FPdim < fusion_rank" << endl;
        mult_of_ev_ok = false;
        return;
    }

    // now we want to identify the pairs (n_i, d_i) from their products n_id=i = n_jd_j
    // we do it only in case rank <= 8 where there is a unique solution: n_i = 2 for exactly
    // one i:
    // identify n_id_i of multiplicity >= 4
    for(auto& ev_x_mult: MultEV_raw){
        if(ev_x_mult.second < 4){
            EV_mult_n_i[ev_x_mult.first] = make_pair(ev_x_mult.second,1);
        }
        if(ev_x_mult.second >= 4){
            Integer codeg = ev_x_mult.first/2;
            EV_mult_n_i[codeg] =make_pair(1, 2);

            size_t k = ev_x_mult.second - 4;
            if(k > 0)
                EV_mult_n_i[ev_x_mult.first] = make_pair(k,1);
        }
    }

    if(verbose){
        verboseOutput() << "codegrees with dim of iorreducibles and multiplicities" << endl;
        for(auto& codeg_mult: EV_mult_n_i){
            verboseOutput() << codeg_mult.first << " mult " << codeg_mult.second.first << " dim irred " << codeg_mult.second.second << endl;
        }
    }

}

template<typename Integer>
Integer Induction<Integer>::N(const key_t i, const key_t j, const key_t k){
    return Tables[i][j][k];
}

template<typename Integer>
void Induction<Integer>::solve_system_low_parts(){

    // foirst row and first column are fixed
    Matrix<Integer> our_equs(0,(nr_rows_low_part -1)*(fusion_rank-1) + 1); // +1 for rhs

    // in the system we use that the entries in the first column are fixed
    // they are taken cae of in the right hand side

    // make row equations
    for(size_t i = 0; i < nr_rows_low_part -1; ++i){
        vector<Integer> this_equ((nr_rows_low_part -1)*(fusion_rank-1));
        for(size_t j = 0; j < fusion_rank -1; ++j){
            this_equ[i*(fusion_rank -1) +j] = fusion_type[j+1];
        }
        //bnowmthe nright hand side
        // note: entry in first column is n_i
        this_equ.push_back(-low_m[i+1].first +low_m[i+1].second);
        our_equs.append(this_equ);
    }

    // make column equations
    for(size_t j = 0; j < fusion_rank - 1; ++j){
        vector<Integer> this_equ((nr_rows_low_part -1)*(fusion_rank-1));
        for(size_t i = 0; i < nr_rows_low_part - 1; ++i){
            this_equ[j + i*(fusion_rank -1)] = low_m[i+1].second;  // ?????????
        }
        this_equ.push_back(- Bounds[0][j + 1]); // right hand side
        our_equs.append(this_equ);
    }

    // cout << "RRR " << our_equs.nr_of_rows() << endl;
    // cout << "CCC " << our_equs.nr_of_columns() << endl;
    // our_equs.debug_print('E');

    Cone<Integer> LP(Type::inhom_equations, our_equs);
    LP.setVerbose(false);
    Matrix<Integer> LowPartsRaw = LP.getLatticePointsMatrix();
    // LowPartsRaw.debug_print('#');

    for(size_t iii = 0; iii < LowPartsRaw.nr_of_rows(); iii++){
        Matrix<Integer> our_low_part(nr_rows_low_part, fusion_rank);
        for(size_t i = 0; i < our_low_part.nr_of_rows(); ++i){
            our_low_part[i][0] = low_m[i].second;
        }
        for(size_t i = 1; i < our_low_part.nr_of_rows(); ++i){
            for(size_t j = 1; j < our_low_part.nr_of_columns(); ++j){
                our_low_part[i][j] = LowPartsRaw[iii][(i-1)*(fusion_rank - 1)+ j-1];
            }
        }
        // our_low_part.debug_print('L');
        LowParts.push_back(our_low_part);
    }
}

template<typename Integer>
void Induction<Integer>::make_low_m_i(){

    // Very first we repeat the data accirding to their multiplicities
    // in second.first and don't need the multiplicities anymore

    vector<pair<Integer, size_t> > EV_n_i;

    for(auto& EV_mult: EV_mult_n_i){
        cout << "MMMMM " << EV_mult.first << "     "  << EV_mult.second.first << endl;
        for(size_t j = 0; j < EV_mult.second.first; ++j){
            EV_n_i.push_back(make_pair(EV_mult.first, EV_mult.second.second));
            // cout << EV_mult.first << " " << EV_mult.second.second << endl;
        }
    }

    cout << "EEEEEEEEEEE " << EV_n_i.size() << endl;

    // First we want to replace eigenvalues by codegrees
    // Note: eigenvalue n_if_i by f_i si9nce we want m_i = FPdim/f_i

    // convert the n_i to Integer because they are used as such in the following
    for(auto& t: EV_n_i){
        Integer dummy = convertTo<Integer>(static_cast<long long>(t.second));
        low_m.push_back(make_pair(FPdim  / t.first, dummy));

        // cout << "m " << FPdim  / t.first << " mult "  << dummy << endl;
    }
    sort(low_m.begin(), low_m.end());


      for(auto& mm: low_m){
        cout << mm.first << " ---- " << mm.second << endl;
    }

}

template<typename Integer>
void Induction<Integer>::build_low_parts(){

    if(!mult_of_ev_ok)
        return;

    if(verbose)
        verboseOutput() << "Computing low parts" << endl;

    make_low_m_i();

    solve_system_low_parts();

    // Low parts may come with permuted rows. We select those with
    // ordered rows
    vector<Matrix<Integer> >  OrderedLowParts;
    for(auto& M: LowParts){
        bool ordered = true;
        for(size_t i = 0; i < M.nr_of_rows() - 1; ++i){
            if(low_m[i].first == low_m[i+1].first){
                if(M[i] > M[i+1]){
                    // cout << "ord " << i << " " << low_m[i] << " " << low_m[i+1] << endl;
                    ordered = false;
                    // break;
                }
            }
        }
        if(ordered)
            OrderedLowParts.push_back(M);
    }

    swap(LowParts, OrderedLowParts);
    // LowParts now contains the otrdered ones

    if(verbose)
        verboseOutput() << "Found " << LowParts.size() << " low parts:"  << endl;

    if(verbose){
        for(auto& M: LowParts)
            M.debug_print('P');
    }
    /* for(auto& M: LowPartsBounds)
        M.debug_print('&'); */

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
void Induction<Integer>::from_low_to_full(){


    if(!mult_of_ev_ok)
        return;

    if(verbose)
        verboseOutput() << "Extending low parts to full induction matrices" << endl;

    //First we compute the contribution of the low part to the conditions for induction matrices

    size_t count_low_parts = 0;

    for(auto& ThisLowPart: LowParts){

        count_low_parts++;
        if(verbose){
            verboseOutput() << "Low part  " << count_low_parts << endl;
            ThisLowPart.debug_print('L');
        }

        Integer FPdim_so_far =0;
        for(size_t j = 0; j < nr_rows_low_part; ++j)
            FPdim_so_far += low_m[j].first * low_m[j].first;

        Integer MinusOne = -1;

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        // ThisLowPart.debug_print('#');
        // LowPartsBounds[iii].debug_print('B');

        Matrix<Integer> MinusLowBound(fusion_rank, fusion_rank);
        for(size_t i = 0; i< nr_rows_low_part; ++i){
            vector<Integer> cand_extension = ThisLowPart[i];
            // cout << "cccc " << cand_extension.size() << endl;
            for(size_t j = 0; j < fusion_rank; ++j){
                // cout << "j "  << j << endl;
                for(size_t k = j; k < fusion_rank; ++k){
                    MinusLowBound[j][k] +=  cand_extension[j] * cand_extension[k];
                }
            }
        }
        MinusLowBound.scalar_multiplication(MinusOne);

        // Now the conditions for the high part
        // Remaining: the difference between whe goal and the contribution of the low part

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

        // now the equation for the FPdim of the center

        vector<Integer> this_equ;
        for(size_t i = 0; i < HighRepresentations.nr_of_rows(); ++i){
            Integer m_new = v_scalar_product(HighRepresentations[i], fusion_type);
            this_equ.push_back(m_new * m_new);
        }
        // cout << "FFFFFF " << FPdim_so_far << " " << -FPSquare << " " << (-FPSquare + FPdim_so_far) <<  endl;
        this_equ.push_back(-FPSquare + FPdim_so_far);
        InhomEqu.append(this_equ);

        // InhomEqu.debug_print('&');

        // we use inequalities to avoid coordinate transformation
        Matrix<Integer> Copy = InhomEqu;
        Copy.scalar_multiplication(MinusOne);
        InhomEqu.append(Copy);

        Matrix<Integer> Unit(HighRepresentations.nr_of_rows());
        Cone<Integer> C(Type::inhom_inequalities, InhomEqu, Type::inequalities, Unit);
        C.setVerbose(false);
        C.compute(ConeProperty::LatticePoints);
        Matrix<Integer> LP = C.getLatticePointsMatrix();

        // LP.debug_print('P');

        if(verbose && LP.nr_of_rows() == 0)
            verboseOutput() << "No extension" << endl;

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
            if(verbose){
                IndMat.debug_print('I');
                verboseOutput() << endl;
            }
        }

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
    swap(InductionMatrices, Representatives); // matrices now oin InductionMatrices

    map<vector<Integer>, vector<Matrix<Integer> > > InductionMatricesByType;

    for(auto& M: InductionMatrices){
        vector<Integer> type = M.MxV(fusion_type);
        InductionMatricesByType[type].push_back(M);
    }
    if(verbose)
        verboseOutput() << InductionMatricesByType.size() << " fusion types defined by induction matrices" << endl;

    InductionMatrices.clear();
    for(auto& T: InductionMatricesByType){
        // cout << T.first;
        InductionMatrices.insert(InductionMatrices.end(), T.second.begin(), T.second.end());
    }

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
    build_low_parts();
    from_low_to_full();
    augment_induction_matrices();
}


template class Induction<mpz_class>;
template class Induction<long long>;
template class Induction<long>;
#ifdef ENFNORMALIZ
template class Induction<renf_elem_class>;
#endif



} // namespace
