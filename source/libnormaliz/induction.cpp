/*
 * Normaliz
 * Copyright (C) 2007-2022  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
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

// Helpers for algebraic numbers

template<typename Integer>
Integer power(const Integer& m, const size_t& k){
    if( k == 0)
        return 1;
    return m * power(m, k -1);
}

template<typename Integer>
long long our_ceil(const Integer& val ){
    return convertTo<long long>(val);
}

template<typename Integer>
long long our_floor(const Integer& val ){
    return convertTo<long long>(val);
}

template<typename Integer>
Integer our_scalar_product(const vector<long long>& av, const vector<Integer>& bv){
    Integer S = 0;
    for(size_t i = 0; i < av.size(); ++i)
        S += convertTo<Integer>(av[i]) *bv[i];
    return S;
}


template<typename Integer>
Integer Induction<Integer>::conjugate(const Integer& val){
    assert(false);
    return 0;
}


/* template<typename Integer>
bool is_algebraic_integer(const Integer& val){
    assert(false);
    return true;
}
*/

template<typename Integer>
Matrix<long long> SplitRepresentations(Integer val, vector<Integer> summands){
    size_t r = summands.size();
    Matrix<long long> InhomEqu(1,r +1);
    for(size_t i = 0; i < r; ++i)
        InhomEqu[0][i] = convertTo<long long>(summands[i]);
    InhomEqu[0].back() = convertTo<long long>(- val);
    return InhomEqu;
}

#ifdef ENFNORMALIZ

template<>
long long our_ceil(const renf_elem_class& val ){
    return convertTo<long long>(val.ceil());
}

template<>
long long our_floor(const renf_elem_class& val ){
    return convertTo<long long>(val.floor());
}

template<>
renf_elem_class our_scalar_product(const vector<long long>& av, const vector<renf_elem_class>& bv){
    renf_elem_class S = 0;
    for(size_t i = 0; i < av.size(); ++i)
        S += av[i] * bv[i];
    return S;
}

vector<mpz_class> minimal_polynomial(const renf_elem_class& val){

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
        if(LinDep.rank() == LinDep.nr_of_rows()) // powers still linearly independent
            continue;
        // now gedree found
        // get rid of denominators
        mpz_class D = 1;
        for(auto& p: powers)
            D = libnormaliz::lcm(D, p.den());
        for(size_t i = 0; i < LinDep.nr_of_rows(); ++i){
            mpz_class fact = D/powers[i].den();
            v_scalar_multiplication(LinDep[i], fact);
        }
        // find linear elation on powers
        Matrix<mpz_class> T = LinDep.transpose();
        Matrix<mpz_class> Sol = T.kernel(false);
        assert(Sol.nr_of_rows() ==1); // since the rank pf the matrix is 1 nr_of_columns -1
        vector<mpz_class> SolVec = Sol[0];
        if(SolVec.back() < 0){
            mpz_class MinusOne = -1;
            v_scalar_multiplication(SolVec, MinusOne);
        }
        return SolVec;
    }
}

/// template<>
bool is_algebraic_integer(const renf_elem_class& val){
    if(val.is_rational()){
        if(val.is_integer())
            return true;
        return false;
    }
    vector<mpz_class> MinPol = minimal_polynomial(val);
    if(MinPol.back() == 1)
        return true;
    else
        return false;
}

bool is_d_number(const renf_elem_class& val){
    if(val.is_rational()){
        if(val.is_integer())
            return true;
        return false;
    }
    vector<mpz_class> MinPol = minimal_polynomial(val);
    if(MinPol.back() != 1)
        return false;
    size_t n = MinPol.size() - 1;
    vector<mpz_class> Rev(n+1);
    for(size_t i = 0; i<= n; ++i)
        Rev[i] = MinPol[n - i];
    for(size_t i = 1; i < n; ++i){
        if( power(Rev[i], n)% power(Rev[n],i) != 0)
            return false;
    }
    return true;
}

template<>
Matrix<long long> SplitRepresentations(renf_elem_class val, vector<renf_elem_class> summands){

    // get rid of denominators
    mpz_class common_den = val.den();
    for(auto& s: summands){
        common_den = libnormaliz::lcm(common_den, s.den());
    }
    val *= common_den;
    for(auto& s: summands){
        s *= common_den;
    }
    vector<mpz_class> val_num = val.num_vector();
    mpz_class MinusOne = -1;
    v_scalar_multiplication(val_num, MinusOne);
    size_t comm_size = val_num.size();
    vector<vector<mpz_class> > Equ_mpz;
    for(auto& s: summands){
        Equ_mpz.push_back(s.num_vector());
        comm_size = max(comm_size, Equ_mpz.back().size());
    }
    for(auto& v: Equ_mpz){
        v.resize(comm_size);
    }
    val_num.resize(comm_size);
    Equ_mpz.push_back(val_num);
    Matrix<long long> Equ_ll(Equ_mpz.size(),comm_size);
    size_t i = 0;
    for(auto& v: Equ_mpz){
        convert(Equ_ll[i],v);
        i++;
    }
    return Equ_ll.transpose();
}


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

/*
template<>
bool Induction<renf_elem_class>::is_algebraic_integer_old(const renf_elem_class& val){
    renf_elem_class conj = conjugate(val);
    if (!(val + conj).is_integer())
        return false;
    if(!(val*conj).is_integer())
        return false;
    return true;
}
*/

#endif


template<typename Integer>
bool check_bounds(const vector<Integer> v, const Matrix<Integer> Bounds){
    assert(v.size() <= Bounds.nr_of_columns());
    assert(v.size() <= Bounds.nr_of_rows());
    for(size_t j = 0; j < v.size(); ++j){
        for(size_t k = j; k < v.size(); ++k){
            if(v[j]*v[k] > Bounds[j][k]){
                return false;
            }
        }
    }
    return true;
}

template<typename Integer>
vector<string> BoundsAsPolynomials(const Matrix<Integer> Bounds){
    size_t n = Bounds.nr_of_columns();
    assert(n == Bounds.nr_of_rows());
    vector<string> Polys;
    for(size_t j = 0; j < n; ++j){
        for(size_t k = j; k < n; ++k){
            // we count variables from 1 for polynomial bounds
            string B = "-x[" + to_string(j+1) + "]*x["+to_string(k+1)+"]+"+to_string(Bounds[j][k]);
            Polys.push_back(B);
        }
    }
    return Polys;
}
// constructors

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
    FusBasic.make_type_automs();
    type_automs.resize(FusBasic.type_automs.size());
    // transfer short_key_t to key_t
    for(size_t i = 0; i < FusBasic.type_automs.size(); ++i){
        type_automs[i].resize(FusBasic.type_automs[0].size());
        for(size_t j = 0; j < FusBasic.type_automs[0].size(); ++j)
            type_automs[i][j] = FusBasic.type_automs[i][j];
    }
    if(verbose)
        verboseOutput() << FusBasic.type_automs.size() << endl;
    // cout << FusRing;
    /* for(auto& T: Tables);
        T.debug_print();*/

    test_commutativity();

    near_integral = false;

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

    if(verbose)
        EVMat.debug_print('E');

    Bounds.resize(fusion_rank, fusion_rank);
    for(size_t j = 0; j< fusion_rank; ++j){
        for(size_t k = 0; k < fusion_rank; ++k){
            Integer S =0;
            for (size_t s = 0; s < fusion_rank; ++s){
                for(size_t t = 0; t < fusion_rank; ++t )
                    S += N(t,j,s)*N(s, duality[t],k);
            }
            Bounds[j][k] = convertTo<long long>(S);
        }
    }
    if(verbose)
        Bounds.debug_print('B');

    convert(Bounds_Int, Bounds);

    BoundsPolys =  BoundsAsPolynomials(Bounds);

    make_divisors();

    if(commutative){
        codegrees_and_mult_commutative();
    }
    else{
        codegrees_and_mult_noncommutative();
    }

    if(verbose)
        verboseOutput() << "Computing representations of divisors of FPdim in terms of type" << endl;

    HighRepresentations.resize(0, fusion_rank);

    for(auto& t: candidates_m_i){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        if(t == FPdim)
            continue;
        Matrix<long long> InhomEqu = SplitRepresentations(t, fusion_type);
        Matrix<long long> NeutralEqu(1,fusion_rank );
        NeutralEqu[0][0] = 1;

        Cone<long long> RepCone(Type::inhom_equations, InhomEqu, Type::equations, NeutralEqu);
        RepCone.setVerbose(false);
        RepCone.setPolynomialInequalities( BoundsPolys);
        RepCone.setNonnegative();
        Matrix<long long> Reps = RepCone.getLatticePointsMatrix();
        // cout << "Reps Reps Reps " << Reps.nr_of_rows() << endl;
        size_t count_high = 0;
        for(size_t i = 0; i < Reps.nr_of_rows(); ++i){
            vector<long long> new_rep = Reps[i];
            new_rep.resize(new_rep.size() - 1);
            if(!check_bounds(new_rep, Bounds))
                continue;

            //if(Reps[i][0] == 0){
                HighRepresentations.append(new_rep);
                count_high++;
            // }

        }

        if(verbose)
            verboseOutput() << "candidate " << t << " has " << count_high << " high" << " representations" << endl;

       // LowRepresentations[t].debug_print('$');
    }
}


template<typename Integer>
void Induction<Integer>::make_divisors(){
    for(Integer t = 1; t <= FPdim; ++t){
        if(FPdim % t == 0){
            divisors.push_back(t);
            candidates_m_i.push_back(FPdim/t);
        }
    }
}

template<typename Integer>
void Induction<Integer>::make_divisors_near_integral(){
    assert(false);
}

#ifdef ENFNORMALIZ

template<>
void Induction<renf_elem_class>::make_divisors_near_integral(){

    size_t r = fusion_type.size() - 1;
    kkk = N(r,r,r);
    d_plus = fusion_type.back();
    d_minus = kkk - d_plus;

    for(long long n = 0; n < FPdim; ++n){
        for(long long m = 0; m < (FPdim-n)/d_plus; ++m){

             INTERRUPT_COMPUTATION_BY_EXCEPTION

            renf_elem_class cand = n + m * d_plus;
            if(cand == 0)
                continue;
            //if(conjugate(cand) < 0)
             //   continue;
            if(is_d_number(cand)){
                if(is_algebraic_integer(FPdim/cand))
                    candidates_m_i.push_back(cand);
                if(is_d_number(FPdim/cand))
                    divisors.push_back(FPdim/cand);
            }
        }
    }
}

template<>
void Induction<renf_elem_class>::make_divisors(){

    near_integral = true;
    for(size_t i = 0; i< fusion_type.size() -1; ++i){
        if(!fusion_type[i].is_integer())
            near_integral = false;
    }
    if(near_integral){
        make_divisors_near_integral();
        return;
    }

    vector<renf_elem_class> h = fusion_type;
    vector<long long> floors;
    for(auto& f: fusion_type)
        floors.push_back(convertTo<long long>(f.floor()));
    floors.push_back(- convertTo<long long>(FPdim.ceil()));
    long long MinusOne = -1;
    v_scalar_multiplication(floors,MinusOne);
    Matrix<long long> Hyp(floors);
    // Hyp.debug_print();
    Cone<long long> CandCone(Type::inhom_inequalities, Hyp);
    CandCone.setNonnegative();
    CandCone.setPolynomialInequalities(BoundsPolys);
    Matrix<long long> RawCands = CandCone.getLatticePointsMatrix();
    set<renf_elem_class> CandSet;

    for(size_t i = 0; i < RawCands.nr_of_rows(); ++i){
        auto c = RawCands[i];
        renf_elem_class cand = v_scalar_product_vectors_unequal_lungth_mixed(c, fusion_type);
        if(cand == 0 || cand > FPdim)
            continue;
        if(CandSet.find(cand) != CandSet.end()) // already tested
            continue;
        CandSet.insert(cand);
        renf_elem_class inv = FPdim/cand;
        if(is_d_number(cand)){
            if(is_algebraic_integer(FPdim/cand))
                candidates_m_i.push_back(cand);
            if(is_d_number(FPdim/cand))
                divisors.push_back(FPdim/cand);
        }
    }
}

#endif

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
bool Induction<Integer>::column_normal(const Matrix<long long>& mat) const{
    for(auto& perm: type_automs){
        Matrix<long long> trans = mat.transpose();
        trans.order_rows_by_perm(perm);
        Matrix<long long> mat_perm = trans.transpose();
        if(mat_perm.get_elements() < mat.get_elements())
            return false;
    }
    return true;
}

template<typename Integer>
void Induction<Integer>::solve_system_low_parts(){

    // foirst row and first column are fixed
    Matrix<long long> our_equs(0,(nr_rows_low_part -1)*(fusion_rank-1) + 1); // +1 for rhs

    // in the system we use that the entries in the first column are fixed
    // they are taken cae of in the right hand side

    // make row equations
    for(size_t i = 0; i < nr_rows_low_part -1; ++i){

        Matrix<long long> RowSplit = SplitRepresentations(low_m[i+1].first, fusion_type);
        for(int k = 0; k < RowSplit.nr_of_rows(); ++k){
            vector<long long> this_equ((nr_rows_low_part -1)*(fusion_rank-1));
            for(size_t j = 0; j < fusion_rank -1; ++j){
                this_equ[i*(fusion_rank -1) +j] = RowSplit[k][j+1];
            }
            // now the nright hand side
            // note: entry in first column is n_i
            if(k == 0)
                this_equ.push_back(RowSplit[k].back() + low_m[i+1].second);
            else
                this_equ.push_back(RowSplit[k].back());
            our_equs.append(this_equ);
        }
    }

    // make column equations
    for(size_t j = 0; j < fusion_rank - 1; ++j){
        vector<long long> this_equ((nr_rows_low_part -1)*(fusion_rank-1));
        for(size_t i = 0; i < nr_rows_low_part - 1; ++i){
            this_equ[j + i*(fusion_rank -1)] = low_m[i+1].second; //  ?????????
        }
        this_equ.push_back(- Bounds[0][j + 1]); // right hand side
        our_equs.append(this_equ);
    }

    // cout << "RRR " << our_equs.nr_of_rows() << endl;
    // cout << "CCC " << our_equs.nr_of_columns() << endl;
    // our_equs.debug_print('E');

    Cone<long long> LP(Type::inhom_equations, our_equs);
    LP.setVerbose(false);
    Matrix<long long> LowPartsRaw = LP.getLatticePointsMatrix();
    // LowPartsRaw.debug_print('#');

    for(size_t iii = 0; iii < LowPartsRaw.nr_of_rows(); iii++){
        Matrix<long long> our_low_part(nr_rows_low_part, fusion_rank);
        for(size_t i = 0; i < our_low_part.nr_of_rows(); ++i){
            our_low_part[i][0] = low_m[i].second;
        }
        for(size_t i = 1; i < our_low_part.nr_of_rows(); ++i){
            for(size_t j = 1; j < our_low_part.nr_of_columns(); ++j){
                our_low_part[i][j] = LowPartsRaw[iii][(i-1)*(fusion_rank - 1)+ j-1];
            }
        }
        // our_low_part.debug_print('L');
        if(column_normal(our_low_part))
            LowParts.push_back(our_low_part);

    }
}

template<typename Integer>
void Induction<Integer>::make_low_m_i(){

    // Very first we repeat the data accirding to their multiplicities
    // in second.first and don't need the multiplicities anymore

    vector<pair<Integer, size_t> > EV_n_i;

    for(auto& EV_mult: EV_mult_n_i){
        // cout << "MMMMM " << EV_mult.first << "     "  << EV_mult.second.first << endl;
        for(size_t j = 0; j < EV_mult.second.first; ++j){
            EV_n_i.push_back(make_pair(EV_mult.first, EV_mult.second.second));
            // cout << EV_mult.first << " " << EV_mult.second.second << endl;
        }
    }

    // cout << "EEEEEEEEEEE " << EV_n_i.size() << endl;

    // First we want to replace eigenvalues by codegrees
    // Note: eigenvalue n_if_i by f_i si9nce we want m_i = FPdim/f_i

    for(auto& t: EV_n_i){
        // Integer dummy = convertTo<Integer>(static_cast<long long>(t.second));
        low_m.push_back(make_pair(FPdim  / t.first, t.second));

        // cout << "m " << FPdim  / t.first << " mult "  << dummy << endl;
    }
    sort(low_m.begin(), low_m.end());


    /* for(auto& mm: low_m){
        cout << mm.first << " ---- " << mm.second << endl;
    }*/

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
    vector<Matrix<long long> >  OrderedLowParts;
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
        verboseOutput() << "Found " << LowParts.size() << " low parts"  << endl;

    /* if(verbose){
        for(auto& M: LowParts)
            M.debug_print('P');
    }*/
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
void Induction<Integer>::high_parts_recursive(const Matrix<long long>& Remaining, size_t p, long start, const Matrix<long long>& Ind_so_far){

    bool not_zero = false;
    for(; p < fusion_rank; ++p){
        for( size_t k = 0; k <= p; ++k){
            if(Remaining[p][k]){
                not_zero = true;
                break;
            }
        }
        if(not_zero)
            break;
    }

    cout << "check " << not_zero << " ----- " << p << " --- " << start << endl;



    for(; start < HighRepsHere.nr_of_rows(); ++start){

        if(HighRepsHere[start][p] == 0)
            continue;

        bool good = true;
        for(size_t j = 0; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                if(Remaining[j][k] - HighRepsHere[start][j]*HighRepsHere[start][k] < 0){
                    good = false;
                    break;
                }
            }
            if(!good)
                break;
        }
        if(!good)
            continue;
        bool complete = true;
        for(size_t j = 0; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                if(Remaining[j][k] - HighRepsHere[start][j]*HighRepsHere[start][k] > 0){
                    complete = false;
                    break;
                }
            }
            if(!complete)
                break;
        }

        Matrix<long long> NewRemaining(Remaining);
        for(size_t j = 0; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                NewRemaining[j][k] -= HighRepsHere[start][j]*HighRepsHere[start][k];
            }
        }
        if(complete){
#pragma omp critical(INDUCTION)
{
            NewRemaining.debug_print('N');
            Matrix<long long> Final(Ind_so_far);
            Final.append(HighRepsHere[start]);
            Matrix<Integer> indmat;
            Final.debug_print('I');
            cout << "*********************************************" << endl;
            convert(indmat, Final);
            InductionMatrices.push_back(indmat);
}
        continue;
        }

        cout << "jump " << start << endl;
        Matrix<long long> NewInd_so_far(Ind_so_far);
        NewInd_so_far.append(HighRepsHere[start]);
        high_parts_recursive(NewRemaining, p, start, NewInd_so_far);
    }
}

template<typename Integer>
void Induction<Integer>::from_low_to_full(){


    if(!mult_of_ev_ok)
        return;

    if(verbose)
        verboseOutput() << "Extending low parts to full induction matrices" << endl;

    //First we compute the contribution of the low part to the conditions for induction matrices

    // size_t count_low_parts = 0;

// #pragma omp parallel for private(HighRepsHere) schedule(dynamic)
    for(size_t lll = 0; lll < LowParts.size(); ++lll){

       // if(lll + 1 != 1945)
       //   continue;

        Matrix<long long> ThisLowPart = LowParts[lll];

        // count_low_parts++;
        if(verbose){
            verboseOutput() << "Low part  " << lll +1 << endl;
            ThisLowPart.debug_print('L');
        }

        Integer FPdim_so_far =0;
        for(size_t j = 0; j < nr_rows_low_part; ++j)
            FPdim_so_far += low_m[j].first * low_m[j].first;

        long long MinusOne_ll = -1;

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        // ThisLowPart.debug_print('#');
        // LowPartsBounds[iii].debug_print('B');

        Matrix<long long> MinusLowBound(fusion_rank, fusion_rank);
        for(size_t i = 0; i< nr_rows_low_part; ++i){
            vector<long long> cand_extension = ThisLowPart[i];
            // cout << "cccc " << cand_extension.size() << endl;
            for(size_t j = 0; j < fusion_rank; ++j){
                // cout << "j "  << j << endl;
                for(size_t k = j; k < fusion_rank; ++k){
                    MinusLowBound[j][k] +=  cand_extension[j] * cand_extension[k];
                }
            }
        }
        MinusLowBound.scalar_multiplication(MinusOne_ll);

        // Now the conditions for the high part
        // Remaining: the difference between whe goal and the contribution of the low part

        Matrix<long long> Remaining = Bounds.add(MinusLowBound);
        // make it symmetric for nicer output
        for(size_t j = 0; j < fusion_rank; ++j){
            // cout << "j "  << j << endl;
            for(size_t k = j; k < fusion_rank; ++k){
                Remaining[k][j] = Remaining[j][k];
            }
        }
        // Remaining.debug_print('R');

        // HighRepresentations.debug_print('H');

        HighRepsHere.resize(0,HighRepresentations.nr_of_columns());
        for(size_t i = 0; i < HighRepresentations.nr_of_rows(); ++i){
            if(check_bounds(HighRepresentations[i], Remaining))
                HighRepsHere.append(HighRepresentations[i]);
        }
        if(verbose)
            verboseOutput() << "Old " << HighRepresentations.nr_of_rows() << " New " << HighRepsHere.nr_of_rows() << endl;

        if(HighRepsHere.nr_of_rows() == 0)
            continue;

        // if(HighRepsHere.nr_of_rows() != 343)
        //   continue;

        sort(HighRepsHere.access_elements().begin(), HighRepsHere.access_elements().end());
        // HighRepsHere.debug_print('H');
        for(size_t kk = 0; kk < HighRepsHere.nr_of_rows()/2; ++kk ){
            swap(HighRepsHere[kk], HighRepsHere[HighRepsHere.nr_of_rows()-1-kk]);
        }

        /*
         * size_t p = 0;
        long start = 0;
        Matrix<long long> Ind_so_far(0, fusion_rank);
        high_parts_recursive(Remaining, p, start, ThisLowPart);
        continue;
        */


        Matrix<long long> InhomEqu(0, HighRepsHere.nr_of_rows() + 1);


        for(size_t j = 1; j < fusion_rank; ++j){
            for(size_t k = j; k < fusion_rank; ++k){
                vector<long long> this_equ;
                for(size_t i = 0; i < HighRepsHere.nr_of_rows(); ++i){
                    this_equ.push_back(HighRepsHere[i][j] * HighRepsHere[i][k]);
                }
                this_equ.push_back(-Remaining[j][k]);
                // cout << j << " " << k << " " << "ttt " <<  this_equ;
                InhomEqu.append(this_equ);
            }
        }


        // now the equation for the FPdim of the center -- integer bound version

        vector<long long> this_inequ;
        for(size_t i = 0; i < HighRepsHere.nr_of_rows(); ++i){
            Integer m_new = our_scalar_product(HighRepsHere[i], fusion_type);
            long long m_int = our_floor(m_new);
            this_inequ.push_back(- m_int*m_int);
        }
        // cout << "FFFFFF " << FPdim_so_far << " " << -FPSquare << " " << (-FPSquare + FPdim_so_far) <<  endl;
        this_inequ.push_back(our_ceil(FPSquare - FPdim_so_far));

        // InhomEqu.debug_print('&');

        // we use inequalities to avoid coordinate transformation
        Matrix<long long > Copy = InhomEqu;
        Copy.scalar_multiplication(-1);
        InhomEqu.append(Copy);

        InhomEqu.append(this_inequ);  // InhimEqu not so nice name

        Cone<long long> C(Type::inhom_inequalities, InhomEqu);
        C.setVerbose(false);
        C.setNonnegative();
        // C.compute(ConeProperty::NumberLatticePoints, ConeProperty::NoPatching, ConeProperty::ShortInt);
        C.compute(ConeProperty::LatticePoints, ConeProperty::NoPatching, ConeProperty::ShortInt);
        Matrix<long long > LP = C.getLatticePointsMatrix();
        // long long NrLatt = C.getNumberLatticePoints(); // C.getLatticePointsMatrix();
        // cout << "NrExt " << NrLatt << endl;


        /* vector<long long> Statistics(HighRepresentations.nr_of_rows());
        for(long long i = 0; i < LP.nr_of_columns() -1; ++i){
            for(long long k = 0; k < LP.nr_of_rows(); ++k){
                Statistics[i] += LP[k][i];
            }
        }

        bool has_solution = (LP.nr_of_rows() > 0);
        cout << "Statistics" << endl;
        for(size_t i = 0; i < Statistics.size(); ++i){
            if(Statistics[i] > 0)
                cout << Statistics[i] << " -- " << our_scalar_product(HighRepsHere[i], fusion_type) << " --- " << HighRepsHere[i];
        }
        if(has_solution)
            exit(0);*/

        // LP.debug_print('P');

        if(verbose){
            if(LP.nr_of_rows() == 0)
                verboseOutput() << "No extension" << endl;
            else
                verboseOutput() << LP.nr_of_rows() << " extensions" << endl;
        }

        for(size_t k= 0; k< LP.nr_of_rows(); ++k){

            // C.setif(k >= 100)
            //  break;

            Matrix<Integer> IndMat;
            convert(IndMat, ThisLowPart);
            // cout << "CCC " << ThisLowPart.nr_of_columns() << endl;
            for(size_t j = 0; j < LP.nr_of_columns() - 1; ++j){
                long long count = 0;
                while(count < LP[k][j]){
                    // cout << "HHH " << HighRepresentations[j].size() << endl;
                    vector<Integer> HR_Int;
                    convert(HR_Int, HighRepsHere[j]);
                    IndMat.append(HR_Int);
                    count = count +1;
                }

            }
            InductionMatrices.push_back(IndMat);

            assert(IndMat.transpose().multiplication(IndMat).equal(Bounds_Int));
            /* if(verbose){
                IndMat.debug_print('I');
                verboseOutput() << endl;
            }*/
        }

    }
}

template<typename Integer>
void Induction<Integer>::augment_induction_matrices(){

    // Sort each induction matrix by increasing FPdim
    // and rows with equal FPdim lex
    if(verbose)
        verboseOutput() << "Sorting induction matrices individually" << endl;
    for(auto& M: InductionMatrices){
        vector<pair<Integer, vector<Integer> > > FPdim_row;
        for(size_t i = 0; i < M.nr_of_rows(); ++i){
            Integer FPdim = v_scalar_product(fusion_type, M[i]);
            FPdim_row.push_back(make_pair(FPdim, std::move(M[i])));
        }
        sort(FPdim_row.begin(), FPdim_row.end());
        for(size_t i = 0; i < M.nr_of_rows(); ++i){
            M[i] = std::move(FPdim_row[i].second);
        }
    }

    if(verbose)
        verboseOutput()<< "Sorting induction matrices lex" << endl;
    vector<vector<vector<Integer> > >  IndVecVecVec;

    for(auto& M: InductionMatrices){
        IndVecVecVec.push_back(std::move(M.get_elements()));
    }

    sort(IndVecVecVec.begin(), IndVecVecVec.end());

    InductionMatrices.clear();
    for(auto& VVV: IndVecVecVec){
        InductionMatrices.push_back(Matrix<Integer>(std::move(VVV)));
    }


    if(verbose)
        verboseOutput()<< "Grouping  induction matrices lby type" << endl;
    // collect induction matrices of equal type
    map<vector<Integer>, vector<Matrix<Integer> > > InductionMatricesByType;
    for(auto& M: InductionMatrices){
        vector<Integer> type = M.MxV(fusion_type);
        InductionMatricesByType[type].push_back(std::move(M));
    }
    if(verbose)
        verboseOutput() << InductionMatricesByType.size() << " fusion types defined by induction matrices" << endl;

    // reorder induction matrices by increasing type
    InductionMatrices.clear();
    for(auto& T: InductionMatricesByType){
        // cout << T.first;
        InductionMatrices.insert(InductionMatrices.end(), std::make_move_iterator(T.second.begin()), std::make_move_iterator(T.second.end()));
    }
    if(verbose)
        verboseOutput() << InductionMatrices.size() << " induction matrices found" << endl;
    InductionMatricesByType.clear();

    // Add additional info
    if(verbose)
        verboseOutput()<< "Augmenting induction matrices by fusion types and duality info" << endl;
    vector<Matrix<Integer> > InductionMatWithType;
    // cout << "FFFFFFFFFFF " << ImageRing;
    InductionMatWithType.push_back(ImageRing);
    for(auto& M: InductionMatrices){
        vector<Integer> type = M.MxV(fusion_type);
        Matrix<Integer> AllowedTranspositions = make_allowed_transpositions(M);
        InductionMatWithType.push_back(std::move(M));
        InductionMatWithType.push_back(Matrix<Integer>(type));
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
