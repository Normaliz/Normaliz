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

#include "libnormaliz/fusion.h"
#include "libnormaliz/options.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/list_and_map_operations.h"

namespace libnormaliz{
using std::vector;
using std::string;
using std::list;
using std::ifstream;
using std::ofstream;


vector<dynamic_bitset> make_all_subsets(const size_t card){

    if(card == 0){
            return vector<dynamic_bitset>(1);
    }

    vector<dynamic_bitset> one_card_less = make_all_subsets(card -1);
    vector<dynamic_bitset> all_subsets;
    dynamic_bitset sub_one_less(card -1);
    dynamic_bitset sub(card);
    for(auto& s: one_card_less){
        for(size_t i = 0; i < card - 1; ++i)
            sub[i] = s[i];
        sub[card -1] = 0;
        all_subsets.push_back(sub);
        sub[card -1] = 1;
        all_subsets.push_back(sub);
    }
    return all_subsets;
}

vector<vector<key_t> > make_all_permutations(const size_t n){
    vector<vector <vector<key_t> > > Perms(1);
    Perms[0].resize(1);
    Perms[0][0].resize(1);
    Perms[0][0][0] = 0;
    for(key_t i=1; i<n; ++i){
        Perms.resize(i+1);
        for(key_t j=0;j<=i;++j){
            for(key_t k=0; k< (key_t) Perms[i-1].size(); ++k){
                vector<key_t> new_perm = Perms[i-1][k];
                new_perm.resize(i+1);
                new_perm[i] = i;
                swap(new_perm[j],new_perm[i]);
                Perms[i].push_back(new_perm);
            }
        }
    }
    sort(Perms[n-1].begin(),Perms[n-1].end());
    return Perms[n-1];
}

vector<vector<key_t> > make_all_permutations(const vector<key_t>& v){

    vector<vector<key_t> >Perms = make_all_permutations(v.size());
    for(auto& w: Perms){
        vector<key_t> w_new(v.size());
        for(size_t i = 0; i< w.size(); ++i)
            w_new[i] = v[w[i]];
        w = w_new;
    }
    return Perms;
}

vector<vector<key_t> > super_impose(const vector<vector<key_t> >& set_1, const vector<vector<key_t> >& set_2){

    vector<vector<key_t> > total;
    for(auto& v: set_1){
        for(auto& w: set_2)
            total.push_back(v_add(v,w));
    }
    return total;
}

vector<vector<key_t> > make_all_permutations(const vector<key_t>& type,const vector<key_t>& duality){

    auto type_1 = type;
    type_1[0] = 0; // to single out the unit

    auto coincidence_keys = collect_coincidence_subset_keys(type_1);
    vector<vector< vector<key_t> > > FullPermsByCoinc;

    for(auto& co: coincidence_keys){
        vector<vector<key_t> > ThisFullPerms;
        auto Perms = make_all_permutations(co);
        vector<key_t> FullPerm(type.size());
        for(auto& p: Perms){
            for(size_t i = 0; i< p.size(); ++i){
                FullPerm[co[i]] = p[i];
            }
            ThisFullPerms.push_back(FullPerm);
        }
        FullPermsByCoinc.push_back(ThisFullPerms);
    }

    vector<vector<key_t> > AllFullPerms;
    for(size_t i = 0; i < coincidence_keys.size(); ++i){
        if(i == 0){
            AllFullPerms=FullPermsByCoinc[0];
        }
        else{
            AllFullPerms = super_impose(AllFullPerms, FullPermsByCoinc[i]);
        }
    }

    if(duality == identity_key(type.size()))
        return AllFullPerms;

    vector<vector<key_t> > Compatible;
    for(auto& p: AllFullPerms){
        bool comp = true;
        for(size_t i = 0; i < p.size(); ++i){
            if(p[duality[i]] != duality[p[i]]){
                comp = false;
                break;
            }
        }
        if(comp)
            Compatible.push_back(p);
    }
    return Compatible;
}

vector<vector<key_t> > collect_coincidence_subset_keys(const vector<key_t>& type){

    vector<vector<key_t> > coincidence_keys;
    dynamic_bitset done(type.size());

    for(size_t i = 0; i < type.size(); ++i){
        if(done[i])
            continue;
        coincidence_keys.push_back(vector<key_t>(1,i));
        done[i] = 1;
        for(size_t j = i + 1; j < type.size(); ++j){
            if(done[j])
                continue;
            if(type[i] == type[j]){
                coincidence_keys.back().push_back(j);
                done[j] = 1;
            }
        }
    }
    return coincidence_keys;
}


// checks whether sol is contained in the "subring". If so, "false" is returned,
// meaning "not sinmple"

template <typename Integer>
void FusionData<Integer>::make_automorphisms(){

    make_CoordMap();

    auto type_automs = make_all_permutations(fusion_type,duality);

    for(auto& p: type_automs){
        vector<key_t> coord_perm(1); // must start with 0 !!!!
        for(auto& t: selected_ind_tuples){
            vector<key_t> image;
            for(auto& c: t)
                image.push_back(p[c]);
            coord_perm.push_back(coord(image));
        }
        Automorphisms.push_back(coord_perm);
    }
    if(verbose)
        verboseOutput() << "Automorphism group of order " << Automorphisms.size() << " computed" << endl;
}


template <typename Integer>
bool FusionData<Integer>::simplicity_check(const vector<key_t>& subring, const vector<Integer>& sol){

    for(auto& c: subring){
        if(sol[c] != 0)
            return true;
    }
    return false;
}

// checks whether sol is contained in one of the "subrings". If so "false" is returned.
template <typename Integer>
bool FusionData<Integer>::simplicity_check(const vector<vector<key_t> >& subrings, const vector<Integer>& sol){

    for(auto& sub: subrings){
        if(!simplicity_check(sub, sol)){
            return false;
        }
    }
    return true;
}

template <typename Integer>
FusionData<Integer>::FusionData(){
    initialize();
}

template <typename Integer>
void FusionData<Integer>::initialize(){
    check_simplicity = false;
    candidate_given = false;
    use_automorphisms = false;
    // select_iso_classes = false;
    verbose = false;
    activated = false;
    type_and_duality_set =false;
    commutative = false;
    nr_coordinates = 0;
}

template <typename Integer>
void FusionData<Integer>::set_options(const ConeProperties& ToCompute, const bool verb){

    verbose = verb;
    check_simplicity= ToCompute.test(ConeProperty::SimpleFusionRings);
    // select_simple = ToCompute.test(ConeProperty::SelectSimple);
    use_automorphisms = ToCompute.test(ConeProperty::FusionRings) || ToCompute.test(ConeProperty::SimpleFusionRings);
    // select_iso_classes = ToCompute.test(ConeProperty::FusionIsoClasses);
    if(check_simplicity || use_automorphisms)
        activated = true;
}

template <typename Integer>
void FusionData<Integer>::import_global_data(){

    fusion_type = fusion_type_coinc_from_input;
    fusion_rank = fusion_type.size();
    duality = fusion_duality_from_input;
    commutative = fusion_commutative_from_input;
    fusion_type_string = fusion_type_from_input;

    subring_base_key = candidate_subring_from_input;
    dynamic_bitset candidate_test = key_to_bitset(subring_base_key, fusion_rank);
    for(size_t i = 0; i < candidate_test.size(); ++i){
        if(candidate_test[i] && !candidate_test[duality[i]])
            throw BadInputException("Candidate subring not closed under duality");

    }
    candidate_given = (subring_base_key.size() >0);
    type_and_duality_set = true;
}


template <typename Integer>
pair<bool, bool>  FusionData<Integer>::read_data(const bool a_priori, const bool only_test) {

    bool dummy = false;

    if(!type_and_duality_set && fusion_type_coinc_from_input.size() > 0){
        import_global_data();
    }
    if(!type_and_duality_set){
            pair<bool, bool> standard_and_partition = data_from_string(global_project, only_test);
            if(only_test && (!standard_and_partition.first || standard_and_partition.second))
                return standard_and_partition;
    }

    set<key_t> subring_base_set;
    subring_base_set.insert(subring_base_key.begin(), subring_base_key.end());
    for(auto& kk: subring_base_key){
        if(subring_base_set.find(duality[kk]) == subring_base_set.end())
            throw BadInputException("Subring base not closed under duality");
    }

    if(fusion_type.size() == 0 || fusion_type.size() != duality.size() ||
                fusion_type[0] != 1 || duality[0] != 0){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("Fusion data corrupt");
    }

    fusion_rank = duality.size();
    dynamic_bitset dual_ind(fusion_rank);
    for(size_t i = 0; i < fusion_rank; ++i){
        if(duality[i] >= fusion_rank || fusion_type[i] != fusion_type[duality[i]]){
            if(only_test)
                return make_pair(false, dummy);
            throw BadInputException("Fusion data corrupt");
        }
         dual_ind[duality[i]] = 1;
    }
    if(dual_ind.count() != fusion_rank){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("Fusion data corrupt");
    }
    if(verbose){
        verboseOutput() << "rank " << fusion_rank << endl;
        verboseOutput() << "type " << fusion_type_string; // contains \n
        verboseOutput() << "duality " << duality;
        verboseOutput() << "commutative " << commutative << endl;
    }
    if(candidate_given && verbose)
        verboseOutput() << "candidate base of subring " << subring_base_key;

    if(!activated)
        return make_pair(true, dummy);

    if((use_automorphisms && a_priori) || (select_iso_classes && !a_priori) )
        make_automorphisms();

    if((check_simplicity && a_priori) || (select_simple && !a_priori) ) // after automorphisms !!
        prepare_simplicity_check();

    return make_pair(true, dummy);
}


template <typename Integer>
pair<bool, bool> FusionData<Integer>::data_from_string(const string& our_fusion, const bool only_test) {

    bool dummy = false;
    if(verbose)
        verboseOutput() << "Reading fusion data from string " << our_fusion << endl;
    string name = pureName(our_fusion);
    string clean_name;
    for(auto& c: name){
        if(c != ' ')
            clean_name.push_back(c);
    }
    size_t open_bracket = 0, close_bracket = 0;
    for(auto& c: clean_name){
        if(c == '[')
            open_bracket++;
        if(c == ']')
            close_bracket++;
    }
    if(clean_name.back() != ']'){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    if(open_bracket != close_bracket){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    if(open_bracket == 0 || open_bracket > 2){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    bool only_partition = false;
    if(open_bracket == 1)
        only_partition = true;

    stringstream data(clean_name);
    char c;
    data >> c;
    if(c !='['){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    vector<long> type;
    while(true){
        long nr;
        data >> nr;
        fusion_type.push_back(nr);
        data >> c;
        if(c == ']'){
            break;
        }
        else{
            if(c != ','){
                if(only_test)
                    return make_pair(false, dummy);
                throw BadInputException("String " + our_fusion +" not standard fusion");
            }
        }
    }
    if(only_partition){
        return make_pair(true, true);

    }
    data >> c;
    if(c !='['){
        if(only_test)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }

    while(true){
        long nr;
        data >> nr;
        if(nr == -1){
            commutative = true;
            nr = 0;
        }
        duality.push_back(nr);
        data >> c;
        if(c == ']'){
            break;
        }
        else{
            if(c != ','){
                if(only_test)
                    return make_pair(false, dummy);
                throw BadInputException("String " + our_fusion +" not standard fusion");
            }
        }
    }

    type_and_duality_set = true;
    return make_pair(true, false);
}

/*

template <typename Integer>
void FusionData<Integer>::read_data_from_file() {

    string file_name = global_project + ".fusion";
    ifstream in(file_name);

    fusion_rank = 0;

    string s;

    size_t duality_count = 0;

    while(true){
        in >> ws;
        int c = in.peek();
        if (c == EOF) {
            break;
        }
        in >> s;

        if(s == "rank"){
            in >> fusion_rank;
            if(fusion_rank == 0)
                throw BadInputException("Fusion rank must be > 0");
            duality = identity_key(fusion_rank);
            continue;
        }

        if(s == "type"){
            if(fusion_rank == 0)
                throw BadInputException("Need fusion rank before reading type");
            if(fusion_type.size() > 0)
                throw BadInputException("Only one fusion type allowed.");
            fusion_type.resize(fusion_rank);
            for(auto& b: fusion_type)
                in >> b;
            continue;
        }
        if(s == "duality"){
            if(fusion_rank == 0)
                throw BadInputException("Need fusion rank before reading duality");
            if(duality_count > 0)
                throw BadInputException("Only one duality input allowed.");
            duality_count++;
            size_t nr_transpositions;
            in >> nr_transpositions;
            for(size_t i = 0; i < nr_transpositions; ++i){
                key_t k,j;
                in >> k >> j;
                duality[k] = j;
                duality[j] = k;
            }
            continue;
        }
        if(s == "subring"){
        size_t rank_candidate;
            in >> rank_candidate;
            subring_base_key.resize(rank_candidate);
            for(size_t i = 0; i< rank_candidate; ++i)
                in >> subring_base_key[i];
            candidate_given = true;
            continue;
        }
        if(s == "commutative"){
            commutative = true;
            continue;
        }
        throw BadInputException("Illegal fusion keyword " + s);
    }

    if( (fusion_type.size() == 0) )
        throw BadInputException("No type in fusion file");

    type_and_duality_set = true;
}
*/

template <typename Integer>
set<vector<key_t> >  FusionData<Integer>::FrobRec(const vector<key_t>& ind_tuple){
    if(commutative)
        return FrobRec_12(ind_tuple);
    else
        return FrobRec_6(ind_tuple);
}


template <typename Integer>
set<vector<key_t> >  FusionData<Integer>::FrobRec_6(const vector<key_t>& ind_tuple){

    assert(ind_tuple.size() == 3);
    key_t i,j,k;
    i = ind_tuple[0];
    j = ind_tuple[1];
    k = ind_tuple[2];
    set< vector<key_t> > F;
    F = {
            {i,j,k},
            {duality[i],k,j},
            {j,duality[k],duality[i]},
            {duality[j],duality[i],duality[k]},
            {duality[k],i,duality[j]},
            {k,duality[j],i}
        };
    return F;
}

template <typename Integer>
set<vector<key_t> >  FusionData<Integer>::FrobRec_12(const vector<key_t>& ind_tuple){

    set< vector<key_t> > F = FrobRec_6(ind_tuple);
    vector<key_t> comm_tuple(3);
    comm_tuple[0] = ind_tuple[1];
    comm_tuple[1] = ind_tuple[0];
    comm_tuple[2] = ind_tuple[2];
    set< vector<key_t> > G = FrobRec_6(comm_tuple);
    for(auto& t: G)
        F.insert(t);
    return F;
}

template <typename Integer>
Integer FusionData<Integer>::value(const vector<Integer>& ring, vector<key_t>& ind_tuple){

    key_t i,j,k;
    i =ind_tuple[0];
    j = ind_tuple[1];
    k = ind_tuple[2];

    if(i == 0){
        if(j == k)
            return 1;
        else
            return 0;
    }
    if(j == 0){
        if(i == k)
            return 1;
        else
            return 0;
    }
    if(k == 0){
        if(duality[i] == j)
            return 1;
        else
            return 0;
    }

    return ring[coord_cone(ind_tuple)];
}


template <typename Integer>
key_t FusionData<Integer>::coord(vector<key_t>& ind_tuple){
    set<vector<key_t> > FR = FrobRec(ind_tuple);
    return coord(FR);
}

template <typename Integer>
key_t FusionData<Integer>::coord_cone(vector<key_t>& ind_tuple){
    key_t coord_compute = coord(ind_tuple);
    if(coord_compute == 0)
        return nr_coordinates;
    return coord_compute -1;
}

template <typename Integer>
key_t FusionData<Integer>::coord(set<vector<key_t> >& FR){
    return CoordMap[FR];
}


// makes the critical coordinates for the simplicity check
// bse_key is the vector of bases (by keys) of the potential subrings
template <typename Integer>
dynamic_bitset FusionData<Integer>::critical_coords(const vector<key_t>& base_key){
    set<key_t> cand_set;
    cand_set.insert(base_key.begin(), base_key.end());

    dynamic_bitset crit_coords(CoordMap.size() + 1); // coordinate 0 is omitted

    for(auto& ind_tuple: all_ind_tuples){
        if(cand_set.find(ind_tuple[0]) == cand_set.end() || cand_set.find(ind_tuple[1]) == cand_set.end()
                || cand_set.find(ind_tuple[2]) != cand_set.end() )
            continue;
        crit_coords[coord(ind_tuple)] = true;
    }
    return crit_coords;

}

template <typename Integer>
void FusionData<Integer>::make_all_ind_tuples(){
    for(key_t i = 1; i < fusion_rank; ++i){
        for(key_t j = 1; j < fusion_rank; ++j){
            for(key_t k = 1; k < fusion_rank; ++k){
                vector<key_t> ind_tuple = {i,j,k};
                all_ind_tuples.push_back(ind_tuple);
            }
        }
    }
}

template <typename Integer>
void FusionData<Integer>::make_CoordMap(){

    if(CoordMap.size() > 0)
        return;

    make_all_ind_tuples();

    key_t val = 1;  // coordinate 0 is the homogenizing one
    for(auto& ind_tuple: all_ind_tuples){
        set<vector<key_t> > F= FrobRec(ind_tuple);
        if(CoordMap.find(F) != CoordMap.end())
            continue;
        CoordMap[F] = val;
        val++;
    }

    nr_coordinates = CoordMap.size();

    // we also want the inverse i-th coordinate --> lex smallest index tuple
    for(auto m = CoordMap.begin(); m!= CoordMap.end(); ++m){
        selected_ind_tuples.push_back(*(m->first.begin()));
    }
}

template <typename Integer>
void  FusionData<Integer>::make_all_base_keys(){

    vector<dynamic_bitset> all_subsets = make_all_subsets(fusion_rank -1);
    for(auto& sub: all_subsets){
        if(sub.count() == 0 || sub.count() == fusion_rank -1) // must discard empty set and the full ring
            continue;
        vector<key_t> kk = bitset_to_key(sub);
        for(auto& c: kk)
            c++;
        bool duality_closed = true;
        for(auto& c: kk){
            if(!sub[duality[c] -1]){
                duality_closed = false;
                break;
            }
        }
        if(!duality_closed)
            continue;
        all_base_keys.push_back(kk);
    }
}

template <typename Integer>
bool FusionData<Integer>::automs_compatible(const vector<key_t>& cand ) const{

    for(auto& aa: Automorphisms){
        dynamic_bitset cand_ind = key_to_bitset(cand, Automorphisms.begin()->size());
        for(auto& c: cand){
            if(!cand_ind[aa[c]])
                return false;
        }
    }
    return true;
}

template <typename Integer>
void FusionData<Integer>::prepare_simplicity_check(){
    make_CoordMap();
    /* for(auto& t: all_ind_tuples){
        cout << coord(t) << " --- " <<t;
    }*/
    if(candidate_given){
        if(automs_compatible(subring_base_key)){
            coords_to_check_ind.push_back(critical_coords(subring_base_key));
            coords_to_check_key.push_back(bitset_to_key(coords_to_check_ind[0]));
        }
        else
            throw BadInputException("Candidate sunbring for non-simplicity not invarient under automorphisms.");
        return;
    }
    // now we must make all candidates
    make_all_base_keys();
    for(auto& bk: all_base_keys){
            coords_to_check_ind.push_back(critical_coords(bk));
            coords_to_check_key.push_back(bitset_to_key(coords_to_check_ind.back()));
    }
}

template <typename Integer>
Matrix<Integer> FusionData<Integer>::do_select_simple_inner(const Matrix<Integer>& LattPoints){
        prepare_simplicity_check();
    if(nr_coordinates != LattPoints.nr_of_columns() - 1)
        throw BadInputException("Wrong number of coordinates in fusion data. Mismatch of duality or commutativity.");

    for(auto& aa: coords_to_check_key){
        for(auto& c: aa) // homogenizing coordinate is no loonger at 0
            c--;
    }
    Matrix<Integer> SimplePoints;
    SimplePoints.resize(0,LattPoints.nr_of_columns());

    for(size_t i = 0; i < LattPoints.nr_of_rows(); ++i){
        if(simplicity_check(coords_to_check_key, LattPoints[i]))
            SimplePoints.append(LattPoints[i]);
    }

    string message = " simple fusion rings found";
    if(candidate_given)
        message = "fusion rings not containing candite subring";

    if(verbose)
        verboseOutput() << SimplePoints.nr_of_rows() << message << endl;

    return SimplePoints;
}


// We work with final format (last coordinate is homogenizing)
// This function protects *this
template <typename Integer>
Matrix<Integer> FusionData<Integer>::do_select_simple(const Matrix<Integer>& LattPoints) const {

    if(LattPoints.nr_of_rows() == 0 || !select_simple)
        return LattPoints;

    FusionData<Integer> work_fusion = *this;
    return work_fusion.do_select_simple_inner(LattPoints);
}

template <typename Integer>
Matrix<Integer> FusionData<Integer>::do_iso_classes_inner(const Matrix<Integer>& LattPoints){

   if(nr_coordinates != LattPoints.nr_of_columns() - 1)
        throw BadInputException("Wrong number of coordinates in fusion data. Mismatch of duality or commutativity.");

    Matrix<Integer> IsoClasses;
    IsoClasses.resize(0,LattPoints.nr_of_columns());

    for(auto& aa: Automorphisms){  // homogenizing coordinate is the last now
        vector<key_t> modified = aa;
        v_cyclic_shift_left(modified, modified.size() -1);
        modified.resize(modified.size()-1);
        for(auto& c: modified)
            c--;
        aa = modified;
    }

    set<vector<Integer> > Classes;

    for(size_t i = 0; i < LattPoints.nr_of_rows(); ++i){
        vector<Integer> max_conjugate = LattPoints[i];
        bool first = true;
        for(auto& aa: Automorphisms){
            vector<Integer> conjugate(LattPoints.nr_of_columns());
            for(size_t j = 0; j < aa.size(); ++ j){
                conjugate[j] = LattPoints[i][aa[j]];
            }
            conjugate.back() = 1;
            if(first || conjugate > max_conjugate){
                max_conjugate = conjugate;
                first = false;
            }
        }
        if(Classes.find(max_conjugate) == Classes.end())
            Classes.insert(max_conjugate);
    }

    for(auto& c: Classes)
        IsoClasses.append(c);

    if(verbose){
        verboseOutput() << IsoClasses.nr_of_rows() << " isomorphism classes computed" << endl;

    }

    return IsoClasses;
}

// We work with final format (last coordinate is homogenizing)
// This function protects *this
template <typename Integer>
Matrix<Integer> FusionData<Integer>::do_iso_classes(const Matrix<Integer>& LattPoints)const {

    if(LattPoints.nr_of_rows() == 0 || !select_iso_classes)
        return LattPoints;

    FusionData<Integer> work_fusion = *this;

    return work_fusion.do_iso_classes_inner(LattPoints);
}

template <typename Integer>
Matrix<Integer> FusionData<Integer>::make_linear_constraints(const vector<Integer>& d){

    make_CoordMap();

    Matrix<Integer> Equ(0, nr_coordinates + 1); // mudst accomodate right hand side in last coordinate

    vector<key_t> indices(3);
    for(key_t i = 1; i < fusion_rank; ++i){
        indices[0] = i;
        for(key_t j = 1; j < fusion_rank; ++j){
            indices[1] = j;
            vector<Integer> this_equ(nr_coordinates + 1);
            this_equ.back() = - d[i]*d[j];
            if(i == duality[j])
                this_equ.back() += 1;
            for(key_t k = 1; k < fusion_rank; ++k){
                indices[2] = k;
                this_equ[coord_cone(indices)] += d[k];
            }
            Equ.append(this_equ);
        }
    }

    Equ.remove_duplicate_and_zero_rows();

    // Equ.pretty_print(cout);
    return Equ;
}

template <typename Integer>
Matrix<Integer> FusionData<Integer>::make_linear_constraints_partition(const vector<Integer>& d,
                                                                       const vector<long>& card){
    make_CoordMap();

    /* cout << "DDDD " << d;
    cout << "CCCC " << card; */

    Matrix<Integer> Equ(0, nr_coordinates + 1); // mudst accomodate right hand side in last coordinate

    vector<key_t> indices(3);
    for(key_t i = 1; i < fusion_rank; ++i){
        indices[0] = i;
        for(key_t j = 1; j < fusion_rank; ++j){
            indices[1] = j;
            vector<Integer> this_equ(nr_coordinates + 1);
            this_equ.back() = - d[i]*d[j]*card[i]*card[j];
            if(i == j) // duality is trivial
                this_equ.back() += card[i];
            for(key_t k = 1; k < fusion_rank; ++k){
                indices[2] = k;
                this_equ[coord_cone(indices)] += d[k];
            }
            Equ.append(this_equ);
        }
    }

    Equ.remove_duplicate_and_zero_rows();

    // Equ.pretty_print(cout);
    return Equ;
}

// factor_1 = factor_1 * factor__2
template <typename Integer>
void prod(pair<Integer, vector<key_t> >& factor_1, const pair<Integer, vector<key_t> >& factor_2){

    if(factor_1.first == 0 || factor_2.first == 0){
        factor_1 = make_pair(0, vector<key_t>(0));
        return;
    }
    factor_1.first *= factor_2.first;
    factor_1.second.insert(factor_1.second.end(),factor_2.second.begin(), factor_2.second.end());
    sort(factor_1.second.begin(), factor_1.second.end());
}

template <typename Integer>
void add(map<vector<key_t>, Integer>& poly, const pair<Integer, vector<key_t> >& summand){

    // cout << summand.first << "+++++++++++++++" << summand.second;

    if(poly.find(summand.second) != poly.end()){
        poly[summand.second] += summand.first;
    }
    else{
        poly[summand.second] = summand.first;
    }
    // cout << "SSSSSSSSSSSSSS " << poly.size() << endl;
}

template <typename Integer>
void subtracct(map<vector<key_t>, Integer>& poly, const pair<Integer, vector<key_t> >& summand){
    pair<Integer, vector<key_t> > subtrahend = summand;
    subtrahend.first = -subtrahend.first;
    add(poly, subtrahend);
}

template <typename Integer>
pair<Integer, vector<key_t> >  FusionData<Integer>::term(const key_t& i, const key_t& j, const key_t& k){

    Integer coeff = -1;
    vector<key_t> exponent;
    if(k == 0){
        if(i == duality[j])
            coeff = 1;
        else
            coeff = 0;
    }
    if(coeff == -1 && i == 0 ){
        if(j == k)
            coeff = 1;
        else
            coeff = 0;
    }
    if(coeff == -1 && j == 0){
        if(i == k)
            coeff = 1;
        else
            coeff = 0;
    }
    if(coeff == -1){
        coeff = 1;
        vector<key_t> indices = {i,j,k};
        exponent.push_back(coord(indices));
    }

    return make_pair(coeff, exponent);
}

template <typename Integer>
set<map<vector<key_t>, Integer> > FusionData<Integer>::make_associativity_constraints(){

    make_CoordMap();

    // we produce associativity_equations b_i(b_j b_k) = (b_i b_j)b_k
    // parametrized by thge base vector b_t
    // if one of i,j,k is 0, no need for it since
    // the neutral element is automatically "associative"

    set<map<vector<key_t>, Integer> > Polys;

    for(key_t i = 1; i< fusion_rank; ++i){
        for(key_t j = 1; j < fusion_rank; ++j){
            for(key_t k = 1; k < fusion_rank; ++k){
                for(key_t t = 0; t < fusion_rank; ++t){
                    map<vector<key_t>, Integer> P;
                    for(key_t s = 0; s < fusion_rank; ++s){ // the poly is a sium over s
                        pair<Integer, vector<key_t> > t_1 = term(i,j,s);
                        pair<Integer, vector<key_t> > t_2 = term(s,k,t);
                        prod(t_1, t_2);
                        add(P, t_1);
                        pair<Integer, vector<key_t> > t_3 = term(j,k,s);
                        pair<Integer, vector<key_t> > t_4 = term(i,s,t);
                        prod(t_3, t_4);
                        subtracct(P, t_3);
                    }
                    bool nonzero = false;
                    for(auto& pp: P){
                        if(pp.second != 0){
                            // cout << pp.second << " -- " << pp.first;
                            nonzero = true;
                        }
                    }
                    if(nonzero)

                        Polys.insert(P);
                }
            }
        }
    }

    /*
    for(auto& q: Polys){
        cout << "****************" <<endl;
        for(auto& p: q){
            cout << p.second << " -- " << p.first;
        }
    }
    cout << "****************" <<endl;
    cout << "NR POLYS " << Polys.size() << endl;
    exit(0); */

    return Polys;
}

template <typename Integer>
void FusionData<Integer>::do_werite_input_file(InputMap<mpq_class>&  input){
    string name = global_project + ".in";
    ofstream out(name);
    if(!out.is_open())
        throw BadInputException("Cannot write input file");
    size_t rank;
    bool is_partition;
    if(contains(input, Type::fusion_type)){
        rank = input[Type::fusion_type].nr_of_columns();
        is_partition = false;
    }
    else{
        rank = input[Type::fusion_type_for_partition].nr_of_columns();
        is_partition = true;
    }
    out << "amb_space " << rank << endl << endl;
    if(is_partition){
        out << "fusion_type_for_partition" << endl;
        out << input[Type::fusion_type_for_partition][0];
    }
    else{
        out << "fusion_type" << endl;
        out << input[Type::fusion_type][0];
        out << endl;
        out << "fusion_duality" << endl;
        out << input[Type::fusion_duality][0];
    }
    out << endl;
    out.close();
    if(libnormaliz::verbose)
        verboseOutput() << "Wtote " << name << endl;
}

template <typename Integer>
void FusionData<Integer>::make_input_from_fusion_data(InputMap<mpq_class>&  input, const bool write_input_file){

    vector<long> bridge(fusion_type.size());
    Matrix<mpq_class> TypeInput(1, fusion_type.size());
    for(size_t i = 0; i< bridge.size(); ++i)
        bridge[i] = fusion_type[i];
    convert(TypeInput[0], bridge);
    for(size_t i = 0; i< bridge.size(); ++i)
        bridge[i] = duality[i];
    Matrix<mpq_class> DualityInput(1, fusion_type.size());
    convert(DualityInput[0], bridge);
    if(commutative)
        DualityInput[0][0] = -1;
    input[Type::fusion_type] = TypeInput;
    input[Type::fusion_duality] = DualityInput;
    if(write_input_file){
        do_werite_input_file(input);
    }
}

template <typename Integer>
void FusionData<Integer>::make_partition_input_from_fusion_data(InputMap<mpq_class>&  input,  const bool write_input_file){

    vector<long> bridge(fusion_type.size());
    Matrix<mpq_class> TypeInput(1, fusion_type.size());
    for(size_t i = 0; i< bridge.size(); ++i)
        bridge[i] = fusion_type[i];
    convert(TypeInput[0], bridge);
    input[Type::fusion_type_for_partition] = TypeInput;
    if(write_input_file){
        do_werite_input_file(input);
    }
}

template <typename Integer>
Matrix<Integer> FusionData<Integer>::data_table(const vector<Integer>& ring, const size_t i){

    Matrix<Integer> Table(fusion_rank, fusion_rank);

    for(key_t k = 0; k < fusion_rank; k++){
        for(key_t j= 0; j < fusion_rank; j++){
            key_t ii = i;
            vector<key_t>ind_tuple = {ii, j, k};
            Table[j][k] = value(ring, ind_tuple);
        }
    }
    // Table.debug_print();
    return Table;
}


template <typename Integer>
vector<Matrix<Integer> > FusionData<Integer>::make_all_data_tables(const vector<Integer>& ring){

    vector<Matrix<Integer> > Tables;

    for(size_t i = 0; i <fusion_rank; ++i){
        Tables.push_back(data_table(ring, i));
    }
    return Tables;
}

template <typename Integer>
void FusionData<Integer>::tables_for_all_rings(const Matrix<Integer>& rings){

    make_CoordMap();

    vector<vector<Matrix<Integer> > > AllTables;
    for(size_t i = 0; i < rings.nr_of_rows(); ++i)
        AllTables.push_back(make_all_data_tables(rings[i]));
}

template <typename Integer>
void FusionData<Integer>::write_all_data_tables(const Matrix<Integer>& rings, ostream& table_out){

    tables_for_all_rings(rings);

    table_out << "[" << endl;
    for(size_t kk = 0; kk < rings.nr_of_rows(); kk++){
        table_out << "  [" << endl;
        vector<Matrix<Integer> > Tables = AllTables[kk]; // for a fixed ring
        for(size_t nn = 0; nn < Tables.size(); ++nn){
            Matrix<Integer> table = Tables[nn];
            table_out << "    [" << endl;
            for(size_t jj = 0; jj < table.nr_of_rows(); ++jj){
                table_out << "      [";
                for(size_t mm = 0; mm < table.nr_of_columns(); ++mm){
                    table_out << table[jj][mm];
                    if(mm < table.nr_of_rows() - 1)
                        table_out << ",";
                    else
                        table_out << "]," << endl;
                }
            }
            if(nn == Tables.size() - 1)
                table_out << "    ]" << endl;
            else
                table_out << "    ]," << endl;
        }
        if(kk == rings.nr_of_rows() -1)
            table_out << "  ]" << endl;
        else
            table_out << "  ]," << endl;
    }
    table_out << "]" << endl;
}
//-------------------------------------------------------------------------------
// helper for fusion rings

// bridge to cone
template <typename Integer>
Matrix<Integer> select_simple(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb){

    FusionData<Integer> fusion;
    fusion.set_options(ToCompute, verb);
    fusion.read_data(false); // falsae = a posteriori
    return fusion.do_select_simple(LattPoints);
}

template <typename Integer>
void split_into_simple_and_nonsimple(Matrix<Integer>& SimpleFusionRings, Matrix<Integer>& NonsimpleFusionRings, const Matrix<Integer>& FusionRings, bool verb){

    if(verb)
        verboseOutput() << "Splitting fusion rings into simple and nonsimple" << endl;

    if(FusionRings.nr_of_rows() == 0){
        if(verb)
            verboseOutput() << "No fusion rings given" << endl;
        return;
    }

    FusionData<Integer> fusion;
    fusion.select_simple = true;
    fusion.activated = true;
    fusion.verbose = false;
    fusion.read_data(false); // falsae = a posteriori
    SimpleFusionRings = fusion.do_select_simple(FusionRings);
    string message = " simple fusion rings (or: not containing candidate subring)";
    /* if(candidate_given)
        message = "fusion rings not containing candite subring"; */

    if(verb)
        verboseOutput() << SimpleFusionRings.nr_of_rows() << message << endl;
    set<vector<Integer> > OurSimple;
    for(size_t i = 0; i < SimpleFusionRings.nr_of_rows(); ++i){
        OurSimple.insert(SimpleFusionRings[i]);
    }
    NonsimpleFusionRings.resize(0,FusionRings. nr_of_columns());
    for(size_t i = 0; i < FusionRings.nr_of_rows(); ++i){
        if(OurSimple.find(FusionRings[i]) == OurSimple.end())
            NonsimpleFusionRings.append(FusionRings[i]);
    }
    string message_1 = " nonsimple fusion rings (or: containing candidate subring)";
    /* if(candidate_given)
        message_1 = "fusion rings containing candite subring";*/

    if(verb)
        verboseOutput() << NonsimpleFusionRings.nr_of_rows() << message_1 << endl;
}

template <typename Integer>
Matrix<Integer> fusion_iso_classes(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb){

    FusionData<Integer> fusion;
    fusion.set_options(ToCompute, verb);
    fusion.read_data(false); // falsae = a posteriori
    return fusion.do_iso_classes(LattPoints);
}

Matrix<long long> extract_latt_points_from_out(ifstream& in_out){

    size_t nr_points;
    in_out >> nr_points;
    string s;
    in_out >> s;
    if(s != "lattice" && s != "fusion" && s!= "simple")
        throw BadInputException("out file not suitable for extraction of sim,ple fusion rtings");
    while(true){
        in_out >> s;
        if(s == "dimension")
            break;
    }
    in_out >> s; // skip = sign
    size_t emb_dim;
    in_out >> emb_dim;
    while(true){
        in_out >> s;
        if(s == "constraints:" || s == "isomorphism:" || s == "data:")
            break;
    }
    Matrix<long long> LattPoints(nr_points, emb_dim);
    for(size_t i = 0; i < nr_points; ++i)
        for(size_t j = 0; j < emb_dim; ++j)
            in_out >> LattPoints[i][j];

    if(in_out.fail())
        throw BadInputException("out file corrupt.");
    return LattPoints;
}


Matrix<long long> read_lat_points_from_file(bool our_verbose){

    string name = global_project + ".final.lat";
    Matrix<long long> LattPoints;
    ifstream in_final(name);
    if(in_final.is_open()){
        if(our_verbose)
            verboseOutput() << "Reading from " << name << endl;
        in_final.close();
        LattPoints = readMatrix<long long>(name);
    }
    else{
        name = global_project + ".out";
        ifstream in_out(name);
        if(!in_out.is_open())
            throw BadInputException("No file with lattice points found");
        if(our_verbose)
            verboseOutput() << "Reading from " << name << endl;
        LattPoints = extract_latt_points_from_out(in_out);
    }
    return LattPoints;
}

void post_process_fusion_file(const vector<string>& command_line_items,string our_project){

    bool non_simple_fusion_rings = true;
    bool verbose = false;
    for(auto& s: command_line_items){
        if(s == "--SimpleFusionRings")
            non_simple_fusion_rings = false;
        if(s == "-c" || s =="--verbose")
            verbose = true;
    }

    if(our_project.size() >= 11){
        if(our_project.substr(our_project.size()-10,10) == ".final.lat"){
            our_project = our_project.substr(0, our_project.size()-10);
        }
    }
    if(our_project.size() >= 5){
        if(our_project.substr(our_project.size()-4,4) == ".out"){
            our_project = our_project.substr(0, our_project.size()-4);
        }
    }
    if(our_project.size() >= 4){
        if(our_project.substr(our_project.size()-3,3) == ".in"){
            our_project = our_project.substr(0, our_project.size()-3);
        }
    }

    global_project = our_project;
    if(verbose)
        verboseOutput() << "Project " << global_project << endl;

    Matrix<long long> LattPoints = read_lat_points_from_file(verbose);
    // LatPoints.debug_print();
    LattPoints.sort_lex();
    size_t embdim = LattPoints.nr_of_columns();

    Matrix<long long> SimpleFusionRings, NonsimpleFusionRings;
    split_into_simple_and_nonsimple(SimpleFusionRings, NonsimpleFusionRings, LattPoints, verbose);

    string name = global_project + ".fusion";
    write_fusion_files(name, non_simple_fusion_rings, true, embdim,
                            SimpleFusionRings, NonsimpleFusionRings,false);
}

void post_process_fusion(const vector<string>& command_line_items){

    string our_project;
    bool list_processing = false;
    bool our_verbose = false;

    for(auto& s: command_line_items){
        if(s[0] != '-')
            our_project = s;
        if(s == "--List")
            list_processing = true;
       if(s == "-c" || s =="--verbose")
            our_verbose = true;
    }
    verbose = our_verbose;

    if(our_project.empty())
        throw BadInputException("No project derfined");
    if(verbose)
        verboseOutput() << "Given file " << our_project << endl;

    if(!list_processing){
        if(verbose)
            verboseOutput() << "Processing single file" << endl;
        post_process_fusion_file(command_line_items, our_project);
        return;
    }

    if(verbose)
        verboseOutput() << "Processing list of files" << endl;

     ifstream list(our_project);
     while(true){
        list >> ws;
        int c = list.peek();
        if (c == EOF) {
            break;
        }
        list >> our_project;
        post_process_fusion_file(command_line_items, our_project);
     }
}

template <typename Integer>
void make_full_input(InputMap<Integer>& input_data, set<map<vector<key_t>, Integer> >& Polys) {

    vector<Integer> full_type = input_data[Type::fusion_type][0];
    // cout << "FULL " << full_type;
    fusion_type_coinc_from_input = fusion_coincidence_pattern(full_type);
    // cout << "COINC " << fusion_type_coinc_from_input;
    size_t fusion_rank_from_input = fusion_type_coinc_from_input.size();
    stringstream for_type;
    for_type << full_type;
    fusion_type_from_input = for_type.str();


    fusion_commutative_from_input = false;

    if(contains(input_data, Type::fusion_duality)){
        vector<Integer> prel_duality = input_data[Type::fusion_duality][0];
        // cout << "PREL " << prel_duality  << " -- " << prel_duality.size() << " -- " <<  fusion_rank_from_input << endl;
        if(prel_duality.size() != fusion_rank_from_input || (prel_duality[0] != 0 && prel_duality[0] != -1))
            throw BadInputException("Fusion duality corrupt");
        if(prel_duality[0] == -1) {
            fusion_commutative_from_input = true;
            prel_duality[0] = 0;
        }
        fusion_duality_from_input.resize(fusion_rank_from_input);
        for(key_t i = 0; i < fusion_rank_from_input; ++i){
            bool in_range = false;
            for(long j = 0; j < fusion_rank_from_input; ++j){
                if(j == convertTo<long>(prel_duality[i])){
                    fusion_duality_from_input[i] = j;
                    in_range = true;
                    break;
                }
            }
            if(!in_range)
                throw BadInputException("Fusion duality corrupt");
        }
    }
    else{
        fusion_duality_from_input = identity_key(fusion_rank_from_input);
    }

    if(contains(input_data, Type::candidate_subring)){
        dynamic_bitset cand_indicator(input_data[Type::candidate_subring][0].size());
        if(cand_indicator.size() != full_type.size())
            throw BadInputException("Candidate subring corrupt");
        for(size_t i = 0; i < cand_indicator.size(); ++i){
            if(input_data[Type::candidate_subring][0][i] == 0){
                continue;
            }
            if(input_data[Type::candidate_subring][0][i] == 1){
                cand_indicator[i] = 1;
                continue;
            }
            throw BadInputException("Candidate subring corrupt");
        }
        if(!cand_indicator[0] || cand_indicator.count() <=1 || cand_indicator.count() == full_type.size())
            throw BadInputException("Candidate subring corrupt");
        candidate_subring_from_input = bitset_to_key(cand_indicator);
    }

    FusionData<Integer> OurFusion;
    OurFusion.import_global_data();
    OurFusion.read_data(true); // checks the duality

    if(verbose)
        verboseOutput() << "Making linear constraints for fusion rings" << endl;
    Matrix<Integer> Equ = OurFusion.make_linear_constraints(full_type);
    Matrix<Integer> InEqu = Equ;
    Integer MinusOne = -1;
    Equ.scalar_multiplication(MinusOne);
    InEqu.append(Equ);

    /* input_data.erase(Type::fusion_type);
    input_data.erase(Type::fusion_duality);
    input_data.erase(Type::candidate_subring);*/
    input_data.clear();
    input_data[Type::inhom_inequalities] = InEqu;
    input_data[Type::inequalities] = Matrix<Integer>(InEqu.nr_of_columns()-1);

    if(verbose)
    verboseOutput() << "Making accociativity constraints for fusion rings" << endl;
    Polys = OurFusion.make_associativity_constraints();

}

template <typename Integer>
void make_full_input_partition(InputMap<Integer>& input_data){

    vector<Integer> full_type = input_data[Type::fusion_type_for_partition][0];
    full_type[0] = 0; // tpo separate the neutral element
    map<Integer, long> blocks = count_in_map<Integer, long>(full_type);
    full_type[0] = 1; // restored;

    vector<Integer> d;
    vector<long> card;
    for(auto& b: blocks){
        card.push_back(b.second);
        d.push_back(b.first);
    }
    d[0] = 1; // restored

    FusionData<Integer> partition_fusion;
    partition_fusion.fusion_type = identity_key(blocks.size());
    partition_fusion.duality = identity_key(blocks.size());
    partition_fusion.commutative = true; // doesn't matter since duality = id
    partition_fusion.fusion_rank = d.size();

    if(verbose)
        verboseOutput() << "Making linear constraints for partition test of fusion rings" << endl;
    Matrix<Integer> Equ = partition_fusion.make_linear_constraints_partition(d, card);
    Matrix<Integer> InEqu = Equ;
    // Equ.pretty_print(cout);
    Integer MinusOne = -1;
    Equ.scalar_multiplication(MinusOne);
    InEqu.append(Equ);

    //input_data.erase(Type::fusion_type_for_partition);
    input_data.clear();
    input_data[Type::inhom_inequalities] = InEqu;
    input_data[Type::inequalities] = Matrix<Integer>(InEqu.nr_of_columns()-1);
}

template <typename Integer>
vector<key_t> fusion_coincidence_pattern(const vector<Integer>& v){

    vector<key_t> coinc;

    if(v.size() == 0)
        return coinc;

    coinc.resize(v.size());

    coinc[0] = 1;
    key_t last_new = 1;
    for(key_t i = 1; i < v.size(); ++i){
        for(key_t j = 1; j < i; ++j){
            if(v[i] == v[j]){
                coinc[i] = coinc[j];
                break;
            }
        }
        if(coinc[i] == 0){
            last_new++;
            coinc[i] = last_new;
        }
    }

    return coinc;
}

/*
template <typename Integer>
void FusionData<Integer>::set_global_fusion_data(){
    assert(false);
}

template <>
void FusionData<long long>::set_global_fusion_data(){
    fusion_type_coinc_from_input = fusion_type;
    fusion_type_from_input = fusion_type_string;
    fusion_duality_from_input = duality;
    fusion_commutative_from_input = commutative;
}
*/

template class FusionData<mpz_class>;
template class FusionData<long long>;
template class FusionData<long>;
#ifdef ENFNORMALIZ
template class FusionData<renf_elem_class>;
#endif

template Matrix<long> select_simple(const Matrix<long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<long long> select_simple(const Matrix<long long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<mpz_class> select_simple(const Matrix<mpz_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#ifdef ENFNORMALIZ
template Matrix<renf_elem_class> select_simple(const Matrix<renf_elem_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#endif

template void split_into_simple_and_nonsimple(Matrix<long long>& SimpleFusionRings, Matrix<long long>& NonsimpleFusionRings, const Matrix<long long>& FusionRings, bool verb);
template void split_into_simple_and_nonsimple(Matrix<long>& SimpleFusionRings, Matrix<long>& NonsimpleFusionRings, const Matrix<long>& FusionRings, bool verb);
template void split_into_simple_and_nonsimple(Matrix<mpz_class>& SimpleFusionRings, Matrix<mpz_class>& NonsimpleFusionRings, const Matrix<mpz_class>& FusionRings, bool verb);
#ifdef ENFNORMALIZ
template void split_into_simple_and_nonsimple(Matrix<renf_elem_class>& SimpleFusionRings, Matrix<renf_elem_class>& NonsimpleFusionRings, const Matrix<renf_elem_class>& FusionRings, bool verb);
#endif

template Matrix<long> fusion_iso_classes(const Matrix<long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<long long> fusion_iso_classes(const Matrix<long long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<mpz_class> fusion_iso_classes(const Matrix<mpz_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#ifdef ENFNORMALIZ
template Matrix<renf_elem_class> fusion_iso_classes(const Matrix<renf_elem_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#endif

template vector<key_t> fusion_coincidence_pattern(const vector<long>& v);
template vector<key_t> fusion_coincidence_pattern(const vector<long long>& v);
template vector<key_t> fusion_coincidence_pattern(const vector<mpz_class>& v);
#ifdef ENFNORMALIZ
template vector<key_t> fusion_coincidence_pattern(const vector<renf_elem_class>& v);
#endif

template void make_full_input(InputMap<long>& input_data, set<map<vector<key_t>, long> >& Polys);
template void make_full_input(InputMap<long long>& input_data, set<map<vector<key_t>, long long> >& Polys);
template void make_full_input(InputMap<mpz_class>& input_data, set<map<vector<key_t>, mpz_class> >& Polys);
#ifdef ENFNORMALIZ
template void make_full_input(InputMap<renf_elem_class>& input_data, set<map<vector<key_t>, renf_elem_class> >& Polys);
#endif

template void make_full_input_partition(InputMap<long>& input_data);
template void make_full_input_partition(InputMap<long long>& input_data);
template void make_full_input_partition(InputMap<mpz_class>& input_data);
#ifdef ENFNORMALIZ
template void make_full_input_partition(InputMap<renf_elem_class>& input_data);
#endif



} // namespace

