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

// checks whether sol is contained in the "subring". If so, "false" is returned,
// meaning "not sinmple"
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
    select_iso_classes = false;
    verbose = false;
    activated = false;
    type_and_duality_set =false;
}

template <typename Integer>
void FusionData<Integer>::set_options(const ConeProperties& ToCompute, const bool verb){

    verbose = verb;
    check_simplicity= ToCompute.test(ConeProperty::OnlySimple);
    select_simple = ToCompute.test(ConeProperty::SelectSimple);
    use_automorphisms = ToCompute.test(ConeProperty::ExploitFusionAutoms);
    select_iso_classes = ToCompute.test(ConeProperty::FusionIsoClasses);
    if(check_simplicity || select_simple || use_automorphisms || select_iso_classes)
        activated = true;
}

template <typename Integer>
void FusionData<Integer>::read_data() {

    if(!activated)
        return;
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

    if(check_simplicity)
        prepare_simplicity_check();
    if(verbose){
        verboseOutput() << "rank " << fusion_rank << endl;
        verboseOutput() << "type " << fusion_type;
        verboseOutput() << "duality " << duality;
    }
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

    while(true){
        in >> ws;    int c = in.peek();
        if (c == EOF) {
            break;
        }
        in >> s;

        if(s == "rank"){
            in >> fusion_rank;
            if(fusion_rank == 0)
                throw BadInputException("Fusion rank must be > 0");
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
            if(verbose)
                verboseOutput() << "Fusion type " << fusion_type;
            continue;
        }
        if(s == "duality"){
            if(fusion_rank == 0)
                throw BadInputException("Need fusion rank before reading duality");
            if(duality.size() != 0)
                throw BadInputException("Only one duality input allowed.");
            duality = identity_key(fusion_rank);
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
        throw BadInputException("Illegal fusion keyword " + s);
    }

    if( (duality.size() == 0 && fusion_type.size() > 0) || (duality.size() > 0 && fusion_type.size() == 0) )
        throw BadInputException("Either both type and duality in fusion file or none");
    if(duality.size() > 0)
        type_and_duality_set = true;
}

template <typename Integer>
set<vector<key_t> >  FusionData<Integer>::FrobRec(const vector<key_t>& ind_tuple){

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
    key_t val = 1;  // coordinate 0 is the homogenizing one
    for(auto& ind_tuple: all_ind_tuples){
        set<vector<key_t> > F = FrobRec(ind_tuple);
        if(CoordMap.find(F) != CoordMap.end())
            continue;
        CoordMap[F] = val;
        val++;
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
void FusionData<Integer>::prepare_simplicity_check(){
    make_all_ind_tuples();
    make_CoordMap();
    /* for(auto& t: all_ind_tuples){
        cout << coord(t) << " --- " <<t;
    }*/
    if(candidate_given){
        coords_to_check_ind.push_back(critical_coords(subring_base_key));
        coords_to_check_key.push_back(bitset_to_key(coords_to_check_ind[0]));
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
void FusionData<Integer>::do_select_and_write_simple(const Matrix<Integer>& LattPoints){ // from out file or final.lat

    prepare_simplicity_check();
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

    string name = global_project + ".simple.lat";
    ofstream out(name);
    out << SimplePoints.nr_of_rows() << endl;
    out << SimplePoints.nr_of_columns() << endl;
    SimplePoints.pretty_print(out);
}

//-------------------------------------------------------------------------------
// helper for fusion rings

// bridge to cone
template <typename Integer>
void select_and_write_simple(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb){

    FusionData<Integer> fusion;
    fusion.set_options(ToCompute, verb);
    fusion.read_data();
    fusion.do_select_and_write_simple(LattPoints);
}

Matrix<long long> extract_latt_points_from_out(ifstream& in_out){

    size_t nr_points;
    in_out >> nr_points;
    string s;
    in_out >> s;
    if(s != "lattice")
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
        if(s == "constraints:")
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


/*
void select_simple_fusion_rings(){

    string name = global_project + ".final.lat";
    Matrix<long long> LattPoints;
    ifstream in_final(name);
    if(in_final.is_open()){
        in_final.close();
        LattPoints = readMatrix<long long>(name);
    }
    else{
        name = global_project + ".out";
        ifstream in_out(name);
        if(!in_out.is_open())
            throw BadInputException("No file with lattice points found");

        LattPoints = extract_latt_points_from_out(in_out);
    }
    FusionData<long long> fusion;
    fusion.read_data(); // true: allow using project name for fusion data
    fusion.select_and_write_simple(LattPoints);
}
*/

template class FusionData<mpz_class>;
template class FusionData<long long>;
template class FusionData<long>;
#ifdef ENFNORMALIZ
template class FusionData<renf_elem_class>;
#endif

template void select_and_write_simple(const Matrix<long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template void select_and_write_simple(const Matrix<long long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template void select_and_write_simple(const Matrix<mpz_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#ifdef ENFNORMALIZ
template void select_and_write_simple(const Matrix<renf_elem_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#endif

} // namespace
