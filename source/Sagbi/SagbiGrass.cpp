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

size_t iter = 0;

vector<set<int> > order;

bool is_revlex_compatible_old(const vector< vector<set <int> > >& GivenDets, const vector<size_t>& marking, size_t nr_vars){
    
    // cout << "marking " << marking.size() << " -- " << marking;
    // print_dets(GivenDets);
    
    /* cout << "In Revlex " << marking;
    
    for(auto& D: GivenDets){
        for(auto& M: D){
            print(M);
            
        }
        cout << endl;
        
    }
    cout << "********************" << endl;*/
    
    // so far marking assigned to Dets[i] with i < marking.size()
    
    set<size_t> all_dets = full_set<size_t>(marking.size());
    // cout << "Start Marking " << marking;
    // cout << "Start marked sets"; print(all_dets);
    
    vector< vector<set <int> > > Dets = GivenDets;
    Dets.resize(marking.size());
    
    /*for(auto& D: Dets){
        for(auto& M: D){
            print(M);
            
        }
        cout << endl;
        
    }
    cout << "********************" << endl;*/
    
    for(auto& D: Dets)
        cancelGCD(D);

    /* size_t kk = 0;
    for(auto& D: Dets){
        cout << "det " << kk << endl;
        kk++;
        for(auto& M: D){
            print(M);
            
        }
        cout << endl;
        
    }
    cout << "********************" << endl; */
    
    vector<vector< bool> > dead(marking.size()); // make dead lists and declare marked monomials dead
    for(auto& i: all_dets){
        dead[i].resize(Dets[i].size());
        dead[i][marking[i]] = true;        
    }
    
    set<int> all_vars = full_set<int>(nr_vars);
    // set<int> vars_used;
    
    order.clear();
    
    // set<int> new_vars;
    
    iter =0;
    
    while(true){
        
        // cout << "Hauptschleife " << endl;
        bool alive = false; // check whether all monomials are dead
        for(auto& i: all_dets){
            if(dead[i] != vector<bool>(Dets[i].size(),true)){
                    alive = true;
                    break;
            }                
        }
        
        // cout << "alive " << vars_used.size() << endl;
        if(!alive)
            return true;
         
       iter++;
        //if(iter > 8)
         //   exit(0);
        
        // cout << "************************" << endl;
        
        //cout << " *** all_vars"; print(all_vars);
        // cout << "*** used vars"; print(vars_used);
        // cout << "all dets dets"; print(all_dets);
        
        /* for(auto& i: all_dets){
            cout << "det " << i;
            if(dead[i] == vector<bool>(Dets[i].size(),true)){
                cout << " dead " << endl;
            }
            else{
                for(size_t j = 0; j < Dets[i].size(); ++j){
                    if(dead[i][j] == false)
                        print(Dets[i][j]);
                cout << endl;
                
                cout << "marked "; print(Dets[i][marking[i]]);
                }
            }
            
        }*/
        
        // if(vars_used.size() == nr_vars) // We have living dets, but
        //    return false;      // no variable by which we can kill them
        
        set<int> vars_in_marked; // collect variables in marked monomials of living dets
        for(auto& k: all_dets){
            if(dead[k] == vector<bool>(Dets[k].size(),true)){
                continue;
            }
            else{
                vars_in_marked = set_union(vars_in_marked, Dets[k][marking[k]]);                        
            }                
        }
        
        // cout << "vars_in_marked "; print(vars_in_marked);
        
        set<int> vars_in_tail;
        for(auto& k: all_dets){
            for(size_t j = 0; j< dead[k].size();++j){
                if(dead[k][j])
                    continue;
                vars_in_tail = set_union(vars_in_tail, Dets[k][j]);
                
            }                
        }
        // cout << "vars_in_tail "; print(vars_in_tail);

        set<int> minimal_vars = set_difference(vars_in_tail, vars_in_marked);
        
        // cout << "vars_in_minimal "; print(minimal_vars);
        
        if(minimal_vars.empty())
            return false;
        
        // cout << "minimal vars "; print(minimal_vars);

        order.push_back(minimal_vars);
        
        bool someone_killed = false;
        
        for(auto& i: all_dets){ // set minimal_vars 0
            for(size_t j = 0; j < Dets[i].size(); ++j){
                if(dead[i][j])
                    continue;
                if(!intersection(Dets[i][j], minimal_vars).empty()){
                    dead[i][j] =true;
                    someone_killed = true;
                }                
            }            
        }        
    }
}

bool is_revlex_compatible_inner(vector< vector<set <int> > >& GivenDets, size_t nr_vars){
    
    set<size_t> all_dets = full_set<size_t>(GivenDets.size());

    auto Dets = GivenDets;
    
    // print_dets(Dets);
    
    for(auto D= Dets.begin(); D != Dets.end();){
        if(D->size() <= 1)
            D = Dets.erase(D);
        else
            D++;
    }
    
    if(Dets.empty())
        return true;
    
    for(auto& D: Dets)
        cancelGCD(D);
    
    // cout << "cleaned" << endl;
    // print_dets(Dets);
    
    set<int> vars_in_marked; // collect variables in marked monomials of living dets
    set<int> vars_in_tail;
    for(auto& D: Dets){
            vars_in_marked = set_union(vars_in_marked, D[0]); 
            for(size_t j = 1; j < D.size(); ++j)
                vars_in_tail = set_union(vars_in_tail, D[j]);
    }
    
    set<int> minimal_candidates = set_difference( vars_in_tail, vars_in_marked);
    /* cout << "mared "; print(vars_in_marked);
    cout << "tail "; print(vars_in_tail);
    cout << "minimal "; print(minimal_candidates); */
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
        bool test = is_revlex_compatible_inner(TestDets, nr_vars);
        if(test)
            return true;
    }
    
    return false;
    
}

bool is_revlex_compatible(const vector< vector<set <int> > >& GivenDets, const vector<size_t>& marking, size_t nr_vars){
    
    // cout << "=============================" << endl;
    
    auto Dets = GivenDets;
    Dets.resize(marking.size());
    for(size_t i = 0; i < marking.size(); ++i){
        swap(Dets[i][0], Dets[i][marking[i]]);        
    }
    bool test = is_revlex_compatible_inner(Dets, nr_vars);
    
    return test;
}

set<vector<mpz_class> > HS_occurring;
set<vector<mpz_class> > HS_Rees_occurring;

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
    Cone<long long> TestCone(Type::cone_and_lattice, A);
    TestCone.compute(ConeProperty::HilbertSeries, ConeProperty::HilbertBasis);
    if(!TestCone.isIntegrallyClosed())
        cout << "Not normal" << endl;
    else{
        vector<mpz_class> HS = TestCone.getHilbertSeries().getNum();
        if(HS_occurring.find(HS) == HS_occurring.end()){
            cout << "New h-vector" << endl;
            cout << HS;
            cout << "KIrull dim " << TestCone.getRank() << " multiplicity " << TestCone.getMultiplicity() << endl;
            HS_occurring.insert(HS);
            
            /* cout << "Order " << endl;
            for(auto& S: order)
                print(S);
            cout << "--------" << endl;
            
            size_t kk = 0;
            for(auto& D: Dets){                  
                for(auto& M: D){
                    print(M);
                }
                cout << "Marking " << marking[kk] << endl;
                kk++;
                cout << "=================" << endl;
            }*/           
        }
    }
    
    Matrix<long long> A_Rees =A;
    Matrix<long long> UnitMat(A.nr_of_columns());
    A_Rees.append(UnitMat);
    vector<long long> ReesVar(A.nr_of_rows(),1);
    ReesVar.resize(A_Rees.nr_of_rows(),0);
    A_Rees.append_column(ReesVar);
    Cone<long long> ReesTestCone(Type::cone_and_lattice, A_Rees);
    ReesTestCone.compute(ConeProperty::HilbertSeries, ConeProperty::HilbertBasis);
    if(!ReesTestCone.isIntegrallyClosed())
        cout << "Rees not normal" << endl;
    else{
        vector<mpz_class> HS = ReesTestCone.getHilbertSeries().getNum();
        if(HS_Rees_occurring.find(HS) == HS_Rees_occurring.end()){
            cout << "New Rees h-vector" << endl;
            cout << HS;
            cout << "KIrull dim " << ReesTestCone.getRank() << " multiplicity " << ReesTestCone.getMultiplicity() << endl;
            HS_Rees_occurring.insert(HS);
        }
    }
}

void build_maring(const vector< vector<set <int> > >& Dets, const size_t& nr_vars,  vector<size_t> marking, size_t& count_compatible){
    
    size_t level = marking.size();
    
    // cout << "On level " << level << " -- " << marking << endl;
    
    if(level == Dets.size())
        return;
    
    marking.resize(level + 1);
    for(size_t i = 0; i< Dets[level].size(); ++i){
        marking[level] = i;
        bool is_good = is_revlex_compatible(Dets, marking, nr_vars);
        if(! is_good)
            continue;
        if(level == Dets.size() -1){
            count_compatible++;
            check_initial_algebra(Dets, nr_vars, marking);
        }
        build_maring(Dets, nr_vars, marking, count_compatible);
        
    }
    return;
    
}

/*
bool next_marking(const vector< vector<set <int> > >& Dets, vector<size_t>& marking){
    
    // fin d last non-maximal mark
    for(int i = Dets.size(); i >= 0, i--){
        if(marking[i] < Dets[i].size() - 1){
            next_exists = true;
            marking[i]++;
            for( int j = i+1; j < Dets[j],size(); j++)
                Marking[j] = 0;
            return true;
        }        
    }
    return false;
}
*/

int main(int argc, char* argv[]) {
    
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
    size_t n = 7;
    
    vector< set< pair<size_t, size_t> > > Patterns;
    Patterns.resize(4);
    Patterns[0] = { {1,0}, {2,0}, {0,1}, {2,1}, {0,2}, {1,2} };
    Patterns[1] = { {1,0}, {2,0}, {2,1}, {0,2}, {0,3}, {1,3} };
    Patterns[2] = { {1,0}, {2,0}, {2,1}, {1,2}, {0,3}, {0,4} };
    Patterns[3] = { {0,0}, {0,1}, {1,2}, {1,3}, {2,4}, {2,5}};

    
for(size_t pat = 0; pat < Patterns.size(); ++pat){
    
    vector< vector < long long  > > M (m, vector<long long> (n,0));
    
    size_t k =0; // we insert the variables into the matrix, -1 for "holes"
    for(size_t i = 0; i< m; i++){
        for(size_t j = 0; j < n; ++j){
            if(Patterns[pat].find(make_pair(i,j)) == Patterns[pat].end()){
                M[i][j] = k;
                k++;
            }
            else{
                M[i][j] = -1;
            }                
        }
    }
    size_t nr_vars = k;
    Matrix<long long>(M).pretty_print(cout);
    // cout << "nr_vars " << nr_vars << endl;
    
    vector< vector<set <int> > > Dets;
   
    // make all mminors as vectors<set<int >> of monomials
    size_t count_dets = 0;   
    vector<size_t> Cols(m);
    for(size_t i = 0; i< m; ++i)
        Cols[i] = i;    
    do{ // go over column selections
        vector< set < int >> ThisDet;
        vector<size_t> Perm(m);
        for(size_t i = 0;  i < m; ++i)
            Perm[i] = i;
        do{ // go over permutations of rows
            set<int> ThisMon;
            // cout << "PPPP " << Perm;
            bool mon_is_zero = false;
            for(size_t i = 0; i < m; ++i){
                // cout << Cols[i] << " "  << Perm[i] << endl;
                if(M[Perm[i]][Cols[i]] == -1){
                    mon_is_zero = true;
                    break;
                }
                else{
                    ThisMon.insert(M[Perm[i]] [Cols[i]]);
                }
            }
            // cout << "---------" << endl;
            if(mon_is_zero)
                continue;
            ThisDet.push_back(ThisMon);
        } while (next_permutation(Perm));
        count_dets ++;
        assert(ThisDet.size() >= 1);
 
        Dets.push_back(ThisDet);
        /*for(auto& D: ThisDet){
            for(auto& I : D){
                cout << I << " ";
            }
            cout << endl;
        }
        cout << "=================" << endl;*/
    } while(next_subset(Cols, n));
    
    size_t count_compatible = 0;
    build_maring(Dets, nr_vars, vector<size_t>(0) , count_compatible);
    cout << "pattern " << pat << " revlex compatible " << count_compatible << endl;
}

}

/*
void build_marking(Marking ){
    
    if(complete)
        return check_initial_ideal(Marking);
    
    for(Mons auf level)
        extend_marking(Marking, Mon);
        if(revlex_compatible(Marking)
            build_Maring)
        
    return;
}



        revlex_compatible (Dets, Marking);
        
        ziehe alle markierten auf Anfang.
        
        micjt benutzte  var = alle
        
        while( )
        
        wenn kein tail vhd. (d.h. alle size 1)
            return true
            
            
        Bilde Teilmenge von var die nicht in mark vorkommwen.
        
        Falls empty return false;
        
        haue alles weg, was davon teilbar
        
        nicht benutzte var -= Teilmenge
        
*/

