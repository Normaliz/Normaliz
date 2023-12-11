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
#include "libnormaliz/input.h"

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
void FusionData<Integer>::read_data(const bool a_priori) {

    string file_name = global_project + ".fusion";
    ifstream in(file_name);
    if(in.is_open()){
        in.close();
        read_data_from_file();
    }
    if(!type_and_duality_set)
        data_from_roject();

    set<key_t> subring_base_set;
    subring_base_set.insert(subring_base_key.begin(), subring_base_key.end());
    for(auto& kk: subring_base_key){
        if(subring_base_set.find(duality[kk]) == subring_base_set.end())
            throw BadInputException("Subring base not closed under duality");
    }

    if(fusion_type.size() == 0 || fusion_type.size() != duality.size() ||
                fusion_type[0] != 1 || duality[0] != 0)
        throw BadInputException("Filename not standard fusion");

    fusion_rank = duality.size();
    dynamic_bitset dual_ind(fusion_rank);
    for(size_t i = 0; i < fusion_rank; ++i){
        if(duality[i] >= fusion_rank || fusion_type[i] != fusion_type[duality[i]])
            throw BadInputException("Filename not standard fusion");
         dual_ind[duality[i]] = 1;
    }
    if(dual_ind.count() != fusion_rank)
        throw BadInputException("Filename not standard fusion");

    if(verbose){
        verboseOutput() << "rank " << fusion_rank << endl;
        verboseOutput() << "type " << fusion_type;
        verboseOutput() << "duality " << duality;
        verboseOutput() << "commutative " << commutative << endl;
    }

    if(!activated)
        return;

    if((use_automorphisms && a_priori) || (select_iso_classes && !a_priori) )
        make_automorphisms();

    if((check_simplicity && a_priori) || (select_simple && !a_priori) ) // after automorphisms !!
        prepare_simplicity_check();

    if(candidate_given && verbose)
        verboseOutput() << "candidate base of subring " << subring_base_key;
}


template <typename Integer>
void FusionData<Integer>::data_from_roject() {
    if(verbose)
        verboseOutput() << "Reading fusion data from project " << global_project << endl;
    string name = pureName(global_project);
    string clean_name;
    for(auto& c: name){
        if(c != ' ')
            clean_name.push_back(c);
    }
    stringstream data(clean_name);
    char c;
    data >> c;
    if(c !='[')
        throw BadInputException("Filename not standard fusion");
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
            if(c != ',')
                throw BadInputException("Filename not standard fusion");
        }
    }
    data >> c;
    if(c !='[')
        throw BadInputException("Filename not standard fusion");
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
            if(c != ',')
                throw BadInputException("Filename not standard fusion");
        }
    }


    type_and_duality_set = true;
}

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


/*
template <typename IntegerPL, typename IntegerRet>
key_t ProjectAndLift<IntegerPL,IntegerRet>::dual(const key_t i) const{
    return duality[i];
}
*/


template <typename Integer>
key_t FusionData<Integer>::coord(vector<key_t>& ind_tuple){
    set<vector<key_t> > FR = FrobRec(ind_tuple);
    return coord(FR);
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

    if(verbose)
        verboseOutput() << SimplePoints.nr_of_rows() << " simple fusion rings found" << endl;

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
    if(verb)
        verboseOutput() << SimpleFusionRings.nr_of_rows() << " simple fusion rings found" << endl;
    set<vector<Integer> > OurSimple;
    for(size_t i = 0; i < SimpleFusionRings.nr_of_rows(); ++i){
        OurSimple.insert(SimpleFusionRings[i]);
    }
    NonsimpleFusionRings.resize(0,FusionRings. nr_of_columns());
    for(size_t i = 0; i < FusionRings.nr_of_rows(); ++i){
        if(OurSimple.find(FusionRings[i]) == OurSimple.end())
            NonsimpleFusionRings.append(FusionRings[i]);
    }
    if(verb)
        verboseOutput() << NonsimpleFusionRings.nr_of_rows() << " nonsimple fusion rings found" << endl;
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

} // namespace

