#include <cstdlib>
#include <vector>
#include <fstream>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

#include "libnormaliz/libnormaliz.h"

using namespace libnormaliz;

template<typename Integer>
void print(const set<Integer>& M1){
    for(auto i: M1)
        cout << " " << i;
    cout << endl;
}

template<typename Integer>
set<Integer> set_union(const set<Integer>& M1, const set<Integer>& M2){
    set<Integer> U = M1;
    for(auto i:M2)
        U.insert(i);
    return U;
}

template<typename Integer>
set<Integer> intersection(const set<Integer>& M1, const set<Integer>& M2){
    set<Integer> I;
    for(auto i:M2)
        if(contains(M1,i))
            I.insert(i);
    return I;
}

template<typename Integer>
set<Integer> full_set(const size_t card){
    set<Integer> F;
    for(int i = 0; i < card; ++i)
        F.insert(i);
    return F;
}

// M1 \ M2
template<typename Integer>
set<Integer> set_difference(const set<Integer>& M1, const set<Integer>& M2){
    set<Integer> D = M1;
    for(auto i:M2){
        if(contains(M2,i))
            D.erase(i);
    }
    return D;
}


bool next_subset(vector<size_t>& key, const size_t& card)
{
    bool can_be_increased = false;
    size_t card_sub = key.size();
    size_t maxi = card-1;
    for(int i = card_sub-1; i >= 0; i--){
        if(key[i] == maxi){
            maxi --;
            continue;
        }
        can_be_increased = true;
        key[i]++;
        for(size_t j = i+1; j < card_sub; j++){
                key[j] = key[i] + j -i;
        }
        break;
    }
    return can_be_increased;
}


bool next_permutation(vector<size_t>& key){
    bool can_be_increased = false;
    size_t card = key.size();
    for(int i =  card -1; i > 0; i--){
        if(key[i] < key[i-1]){
            continue;
        }
        can_be_increased = true;
        for(int j = card - 1; j>i -1; j--){
            if(key[j] > key[i-1]){
                swap(key[j], key[i-1]);
                break;
            }
        }
        auto start = key.begin();
        for(int k = 0; k <= i-1; ++k){
            start++;
        }
        sort(start, key.end());
        break;
    }
    return can_be_increased;
}

void cancelGCD(vector< set<int> >& ThisDet){
    if(ThisDet.empty())
        return;

    set<int> GCD = ThisDet[0];
    for(auto& M: ThisDet)
        GCD = intersection(GCD, M);

    for(auto& M: ThisDet)
        M = set_difference(M, GCD);
}

void print_dets(const vector< vector<set <int> > >& Dets){
    for(auto& D: Dets){
        for(auto& M: D){
            cout << " [ ";
            for(auto i: M)
                cout << i << " ";
            cout << "]";
        }
        cout << endl;
    }
    cout << "---------------------------" << endl;
}

vector<set<int> > order;

void clean_up_dets(vector< vector<set <int> > >& Dets){

    for(auto D= Dets.begin(); D != Dets.end();){ // only marked term left?
        if(D->size() <= 1)
            D = Dets.erase(D);
        else
            D++;
    }

    for(auto& D: Dets)
        cancelGCD(D);


}

vector< vector<set <int> > > FullDets;

bool is_lex_compatible_inner(vector< vector<set <int> > > Dets){

    // print_dets(Dets);

    clean_up_dets(Dets);

    // print_dets(Dets);

    if(Dets.empty())
        return true;

    set<int> maximal_candidates_raw;

    for(auto& D: Dets){
        maximal_candidates_raw = set_union(maximal_candidates_raw, D[0]);

    }

    for(auto& D: Dets){
        set<int> non_max;
        for(size_t j = 1; j < D.size(); ++j){
            non_max = set_union(non_max, D[j]);
        }
        non_max =set_difference(non_max,D[0]);
        maximal_candidates_raw = set_difference(maximal_candidates_raw, non_max);
    }

    set<int> maximal_candidates = maximal_candidates_raw;

    // print_dets(Dets);
    // cout << "max candidates "; print(maximal_candidates);

    if(maximal_candidates.empty())
        return false;

    for(auto i: maximal_candidates){
        bool one_killed = false;
        auto TestDets = Dets;
        for(auto& D: TestDets){
            if(contains(D[0],i)){
                auto M = D.begin();
                M++;
                for(; M!= D.end();){
                    if(!contains(*M,i)){
                        M = D.erase(M);
                        one_killed = true;
                    }
                    else
                        M++;
                }
            }
        }

        if(!one_killed){
            cout << "Nothing killed" << endl;
            exit(0);
        }

        bool test = is_lex_compatible_inner(TestDets);
        if(test)
            return true;
    }
    return false;
}

bool is_lex_compatible_full(const vector< vector<set <int> > >& GivenDets, const size_t nr_vars,
                            const vector<size_t>& marking){

    auto Dets = GivenDets;
    Dets.resize(marking.size());
    for(size_t i = 0; i < marking.size(); ++i){
        swap(Dets[i][0], Dets[i][marking[i]]);
    }


    auto D = Dets.begin();
    auto F = FullDets.begin();
    for(; D != Dets.end(); D++, F++){
        for(size_t i = 0; i< F->size(); ++i){
            if((*F)[i] == (*D)[0]){
                swap((*F)[i], (*F)[0]);
                break;
            }
        }
    }


    bool result = is_lex_compatible_inner(FullDets);
    /* if(!result){
        print_dets(Dets);
        print_dets(FullDets);
        exit(0);
    }*/

    return result;
}

bool is_revlex_compatible_inner(vector< vector<set <int> > > Dets){

    clean_up_dets(Dets);

    // print_dets(Dets);

    if(Dets.empty())
        return true;

    set<int> vars_in_marked; // collect variables in marked monomials of living dets
    set<int> vars_in_tail;
    for(auto& D: Dets){
            vars_in_marked = set_union(vars_in_marked, D[0]);
            for(size_t j = 1; j < D.size(); ++j)
                vars_in_tail = set_union(vars_in_tail, D[j]);
    }

    set<int> minimal_candidates = set_difference( vars_in_tail, vars_in_marked);
    if(minimal_candidates.empty())
        return false;

    for(auto k: minimal_candidates){

        auto TestDets = Dets;
        for(auto& D: TestDets){
            for(auto M = D.begin(); M!= D.end(); ){
                if(contains(*M, k)){
                    M = D.erase(M);
                }
                else
                    M++;
            }
        }
        bool test = is_revlex_compatible_inner(TestDets);
        if(test)
            return true;
    }

    return false;

}

set<int> S1 = {1,2};
set<int> S2 = {3,};
vector<set<int> > M1 ={S1, S2};
vector<set<int> > M2 ={S2, S1};
 vector<vector<set<int> > > DDD = {M1,M2};

bool has_weight_vector(const vector< vector<set <int> > >& Dets,  const size_t nr_vars ){

    /*Dets = DDD;
    nr_vars = 6; */

    Matrix<mpz_class> StrictInequalities(0, nr_vars);
    for(auto& D: Dets){
        if(D.size() <= 1)
            continue;
        for(size_t i = 1; i < D.size(); ++i){
            vector<mpz_class> strict_ineq(nr_vars);
            for(auto j: D[0])
                strict_ineq[j] += 1;
            for(auto j: D[i])
                strict_ineq[j] += -1;
            StrictInequalities.append(strict_ineq);
        }
    }
    if(StrictInequalities.nr_of_rows() == 0)
        return true;
    // StrictInequalities.pretty_print(cout);
    Cone<mpz_class> TestCone(Type::strict_inequalities, StrictInequalities);
    // cout << "AAAAAAAAA " << TestCone.getAffineDim() << endl;

    if (TestCone.getAffineDim() >= 0){
        return true;
    }

    return false;

}

bool do_lex, do_revlex, do_all, do_Rees;

string monord_string, Rees_string;

bool is_compatible(const vector< vector<set <int> > >& GivenDets, const size_t nr_vars, const vector<size_t>& marking){

    // cout << "=============================" << endl;

    auto Dets = GivenDets;
    Dets.resize(marking.size());
    for(size_t i = 0; i < marking.size(); ++i){
        swap(Dets[i][0], Dets[i][marking[i]]);
    }
    bool test;
    if(do_revlex)
        test = is_revlex_compatible_inner(Dets);
    if(do_lex)
        test = is_lex_compatible_inner(Dets);
    if(do_all)
        test = has_weight_vector(Dets, nr_vars);

    // cout << "comp " << test << endl;

    return test;
}

set<vector<mpz_class> > HS_occurring;
set<vector<mpz_class> > Nonnormal_HS_occurring;
set<vector<mpz_class> > HS_Rees_occurring;
set<vector<mpz_class> > Nonnormal_HS_Rees_occurring;
bool ini_not_normal;
bool rees_not_normal;

vector<mpz_class> GoodHS;

// vector<mpz_class> TRR ={1,34,350,1505,3038,3039,1503,348,31};

void check_initial_algebra(const vector< vector<set <int> > >& Dets, const size_t& nr_vars,  vector<size_t> marking){
    Matrix<long long> A(0, nr_vars);
    size_t kk = 0;
    for(auto& D: Dets){
        vector<long long> gen(nr_vars,0);
        for(auto i: D[marking[kk]])
            gen[i] = 1;
        kk++;
        A.append(gen);
    }

    bool this_ini_is_good = true;
    vector<mpz_class> HS;
    Cone<long long> TestCone(Type::monoid, A);
    TestCone.compute(ConeProperty::HilbertSeries, ConeProperty::Multiplicity);
    HS = TestCone.getHilbertSeries().getNum();
    if(!TestCone.isIntegrallyClosed()){
        this_ini_is_good = false;
        if(!ini_not_normal){
            cout << "Not normal" << endl;
            cout << "$$$$$$$$$$" << endl;
            ini_not_normal = true;
        }
    }
    else{
        if(HS_occurring.empty())
            GoodHS = HS;

        if(HS != GoodHS){
            this_ini_is_good = false;
        }

        if(HS_occurring.find(HS) == HS_occurring.end()){
            cout << "*** New h-vector" << endl;
            cout << HS;
            cout << "KIrull dim " << TestCone.getRank() << " multiplicity " << TestCone.getMultiplicity() << endl;
            cout << "-----------------------------" << endl;
            HS_occurring.insert(HS);
        }
    }

    if(do_Rees && this_ini_is_good){  // this_ini_is_good &&
        Matrix<long long> A_Rees =A;
        Matrix<long long> UnitMat(A.nr_of_columns());
        A_Rees.append(UnitMat);
        vector<long long> ReesVar(A.nr_of_rows(),1);
        ReesVar.resize(A_Rees.nr_of_rows(),0);
        A_Rees.append_column(ReesVar);
        Cone<long long> ReesTestCone(Type::cone_and_lattice, A_Rees);
        ReesTestCone.compute(ConeProperty::HilbertSeries, ConeProperty::HilbertBasis);
        if(!ReesTestCone.isIntegrallyClosed() && !rees_not_normal){
            cout << "Rees not normal" << endl;
            cout << "&&&&&&&&&&&&&&& " << endl;
            rees_not_normal = true;
        }
        else{
            vector<mpz_class> HS_Rees = ReesTestCone.getHilbertSeries().getNum();

            if(HS_Rees_occurring.find(HS_Rees) == HS_Rees_occurring.end()){
                cout << "New Rees h-vector for " << HS;
                cout << HS_Rees;
                cout << "KIrull dim " << ReesTestCone.getRank() << " multiplicity " << ReesTestCone.getMultiplicity() << endl;
                cout << "+++++++++++++++++++++++++++++" << endl;
                HS_Rees_occurring.insert(HS_Rees);
            }
        }

    }
}

void build_maring(const vector< vector<set <int> > >& Dets, const size_t& nr_vars,  vector<size_t> marking, size_t& count_compatible){

    // cout << Dets.size() << " " << marking.size() << endl;

    size_t level = marking.size();

    // cout << "On level " << level << " -- " << marking << endl;

    if(level == Dets.size())
        return;

    marking.resize(level + 1);
    // cout << "LLL " << Dets[level].size() << endl;
    for(size_t i = 0; i< Dets[level].size(); ++i){

        // cout << "iii " << i << endl;

        // cout << marking;
        marking[level] = i;
        // cout << "mark " << marking;
        bool is_good = is_compatible(Dets, nr_vars, marking);
        if(! is_good)
            continue;
        if(level == Dets.size() -1){
            if(do_lex)
                is_good = is_lex_compatible_full(Dets, nr_vars, marking);
            if(!is_good)
                continue;
            count_compatible++;
            check_initial_algebra(Dets, nr_vars, marking);
        }
        build_maring(Dets, nr_vars, marking, count_compatible);

    }
    return;

}

int main(int argc, char* argv[]) {

    if(argc < 2){
        cout << "Not enough parameters" << endl;
        exit(1);
    }


    /* vector<size_t> subset_key = {0,1,2};
    size_t card = 5;
    while(true){
        cout << subset_key;
        if(!next_subset(subset_key,card))
            break;
    }
    cout << "-----------------" << endl;
    vector<size_t> perm = {0,1,2};
    while(true){
        cout << perm;
        if(!next_permutation(perm))
            break;
    }*/

    size_t m = 3;
    size_t n;

    string n_string(argv[1]);
    n = stoi(n_string);

    do_revlex = false;
    do_lex = false;
    do_all = false;
    do_Rees = false;

    for(int i = 2; i < argc; ++i){
        string para(argv[i]);
        if(para == "-r"){
            do_revlex = true;
            monord_string = " revlex ";
        }
        if(para == "-l"){
            do_lex = true;
            monord_string = " lex ";
        }
        if(para == "-a"){
            do_all = true;
            monord_string = " weight ";
        }
        if(para == "-R"){
            do_Rees = true;
            Rees_string = " with Rees";
        }
    }

    cout << " n = " << n  << " revlex " << do_revlex << " lex " << do_lex << " weight " << do_all << " Rees " << do_Rees << endl;;

    vector< set< pair<size_t, size_t> > > Patterns;
    Patterns.resize(4);
    Patterns[0] = { {1,0}, {2,0}, {0,1}, {2,1}, {0,2}, {1,2} };
    Patterns[1] = { {1,0}, {2,0}, {2,1}, {0,2}, {0,3}, {1,3} };
    Patterns[2] = { {1,0}, {2,0}, {2,1}, {1,2}, {0,3}, {0,4} };
    Patterns[3] = { {0,0}, {0,1}, {1,2}, {1,3}, {2,4}, {2,5}};

    ini_not_normal = false;
    rees_not_normal = false;


    for(size_t pat = 0; pat < Patterns.size() ; ++pat){

        if(do_lex && pat == Patterns.size()- 1)
            continue;

        FullDets.clear();
        HS_occurring.clear();
        HS_Rees_occurring.clear();
        ini_not_normal = false;
        rees_not_normal = false;

        vector< vector < long long  > > M (m, vector<long long> (n,0));

        size_t k =0; // we insert the variables into the matrix, -1 for "holes"
        int extra_vars = 0;
        for(size_t i = 0; i< m; i++){
            for(size_t j = 0; j < n; ++j){
                if(Patterns[pat].find(make_pair(i,j)) == Patterns[pat].end()){
                    M[i][j] = k;
                    k++;
                }
                else{
                    M[i][j] = -extra_vars  -1;
                    extra_vars++;
                }
            }
        }
        size_t nr_vars = k;

        for(size_t i = 0; i< m; i++){
            for(size_t j = 0; j < n; ++j){
                if(M[i][j]<0){
                    M[i][j] = nr_vars - M[i][j] -1;
                    k++;
                }
            }
        }
        Matrix<long long>(M).pretty_print(cout);

        // cout << "nr_vars " << nr_vars << endl;

        vector< vector<set <int> > > Dets;

        // make all mminors as vectors<set<int >> of monomials
        size_t count_dets = 0;
        vector<size_t> Cols(m);
        for(size_t i = 0; i< m; ++i)
            Cols[i] = i;
        do{ // go over column selections
            vector< set < int >> ThisDet; // using variables outside the pattern
            vector< set < int >> FullDet; // using all variables
            vector<size_t> Perm(m);
            for(size_t i = 0;  i < m; ++i)
                Perm[i] = i;
            do{ // go over permutations of rows
                set<int> ThisMon;
                // cout << "PPPP " << Perm;
                bool mon_is_zero = false;
                for(size_t i = 0; i < m; ++i){
                    // cout << Cols[i] << " "  << Perm[i] << endl;
                    ThisMon.insert(M[Perm[i]] [Cols[i]]);
                }
                FullDet.push_back(ThisMon);
                bool outside_pattern = false;
                for(auto u: ThisMon){
                    if(u >= nr_vars)
                        outside_pattern = true;
                }
                if(!outside_pattern)
                    ThisDet.push_back(ThisMon);
            } while (next_permutation(Perm));
            count_dets ++;
            assert(ThisDet.size() >= 1);

            Dets.push_back(ThisDet);
            FullDets.push_back(FullDet);
            /*for(auto& D: ThisDet){
                for(auto& I : D){
                    cout << I << " ";
                }
                cout << endl;
            }
            cout << "=================" << endl;*/
        } while(next_subset(Cols, n));

        // print_dets(Dets);
        // print_dets(FullDets);

        size_t count_compatible = 0;
        build_maring(Dets, nr_vars, vector<size_t>(0) , count_compatible);
        cout << "Type " << pat+1  << monord_string << "compatible " << count_compatible << endl;
        cout << "============================================" << endl;
    }

}
