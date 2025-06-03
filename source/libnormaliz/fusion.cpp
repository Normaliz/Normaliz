/*
 * Normaliz
 * Copyright (C) 2007-2025  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        for(size_t i = 0; i < card - 1; ++i)
            sub[i] = s[i];
        sub[card -1] = 0;
        all_subsets.push_back(sub);
        sub[card -1] = 1;
        all_subsets.push_back(sub);
    }
    return all_subsets;
}

vector<vector<shortkey_t> > make_all_permutations(const size_t n){
    vector<vector <vector<shortkey_t> > > Perms(1);
    Perms[0].resize(1);
    Perms[0][0].resize(1);
    Perms[0][0][0] = 0;
    for(shortkey_t i=1; i<n; ++i){
        Perms.resize(i+1);
        for(shortkey_t j=0;j<=i;++j){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            for(shortkey_t k=0; k< (shortkey_t) Perms[i-1].size(); ++k){
                vector<shortkey_t> new_perm = Perms[i-1][k];
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

template <typename Integer>
vector<vector<shortkey_t> > make_all_permutations(const vector<shortkey_t>& v, const vector<key_t>& duality,
                                                  const  Matrix<Integer>& fusion_ring_map){

    // the key contains the indices of a given FPdim (as indicated) in the coincidence pattern
    // The permutaions permute the key.

    vector<vector<shortkey_t> >Perms = make_all_permutations(v.size()); // prmutations of 0,...,v.size()-1
    vector<vector<shortkey_t> > KeyPerms;
    vector<shortkey_t> v_inv(duality.size()); // inv = inverse
    for(size_t i = 0; i < v.size(); ++i)
        v_inv[v[i]] = i;
    vector<shortkey_t> dual(v.size());
    for(size_t i = 0; i < v.size(); ++i)
        dual[i] = v_inv[duality[v[i]]];
    for(auto& w: Perms){

        INTERRUPT_COMPUTATION_BY_EXCEPTION

        bool comp = true;
        for(size_t i = 0; i < v.size(); ++i){
            if(dual[w[i]] != w[dual[i]]){
                comp = false;
                break;
            }
        }
        if(!comp)
            continue;

        vector<shortkey_t> w_new(v.size());
        for(size_t i = 0; i< w.size(); ++i)
            w_new[i] = v[w[i]];
        // we must check whether the paermutation leaves the images in fusion_ring_map invariant
        if(fusion_ring_map.nr_of_rows() > 0){
            for(size_t i = 0; i < v.size(); ++i){
                if(fusion_ring_map[v[i]] != fusion_ring_map[w_new[i]]){
                    comp = false;
                    break;
                }
            }
            if(!comp)
                continue;
            }
        KeyPerms.push_back(w_new);
    }
    return KeyPerms;
}


vector<vector<shortkey_t> > super_impose(const vector<vector<shortkey_t> >& set_1, const vector<vector<shortkey_t> >& set_2){

    size_t size_1 = set_1.size();
    size_t size_2 = set_2.size();
    size_t nr_coods = set_1[0].size();
    vector<vector<shortkey_t> > total(size_1*size_2, vector<shortkey_t>(nr_coods));

    bool skip_remaining = false;
    std::exception_ptr tmp_exception;

#pragma omp parallel for
    for(size_t i = 0; i < size_1; ++i){

        if(skip_remaining)
            continue;
        try{
            for(size_t j = 0; j < size_2; ++j){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

                total[i*size_2 +j] = v_add(set_1[i],set_2[j]);
            }
        }catch (const std::exception&) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
    }

    if (!(tmp_exception == 0))
        std::rethrow_exception(tmp_exception);

    return total;
}


/*
vector<vector<shortkey_t> > super_impose(const vector<vector<shortkey_t> >& set_1, const vector<vector<shortkey_t> >& set_2){

    vector<vector<shortkey_t> > total;
    for(auto& v: set_1){
        for(auto& w: set_2){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            total.push_back(v_add(v,w));
        }
    }
    return total;
}
*/

template <typename Integer>
vector<vector<shortkey_t> > make_all_full_permutations(const vector<key_t>& type,const vector<key_t>& duality,
                                                  const  Matrix<Integer>& fusion_ring_map){

    auto type_1 = type;
    type_1[0] = 0; // to single out the unit

    auto coincidence_keys = collect_coincidence_subset_keys(type_1);
    vector<vector< vector<shortkey_t> > > FullPermsByCoinc;

    for(auto& co: coincidence_keys){
        vector<vector<shortkey_t> > ThisFullPerms;
        auto Perms = make_all_permutations(co, duality, fusion_ring_map);
        vector<shortkey_t> FullPerm(type.size());
        for(auto& p: Perms){
            for(size_t i = 0; i< p.size(); ++i){
                FullPerm[co[i]] = p[i];
            }
            ThisFullPerms.push_back(FullPerm);
        }
        FullPermsByCoinc.push_back(ThisFullPerms);
    }

    vector<vector<shortkey_t> > AllFullPerms;
    for(size_t i = 0; i < coincidence_keys.size(); ++i){
        if(i == 0){
            AllFullPerms=FullPermsByCoinc[0];
        }
        else{
            AllFullPerms = super_impose(AllFullPerms, FullPermsByCoinc[i]);
        }
    }

    return AllFullPerms;
}

vector<vector<shortkey_t> > collect_coincidence_subset_keys(const vector<key_t>& type){

    vector<vector<shortkey_t> > coincidence_keys;
    dynamic_bitset done(type.size());

    for(size_t i = 0; i < type.size(); ++i){
        if(done[i])
            continue;
        coincidence_keys.push_back(vector<shortkey_t>(1,i));
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

//--------------------------------------------------------------
//
// FusionBasic
//
//--------------------------------------------------------------

FusionBasic::FusionBasic(){
    commutative = false;
    use_modular_grading = false;
    candidate_given = false;
    fusion_rank = 0;
    type_and_duality_set = false;
    type_automs_made = false;
    total_FPdim = 0;
}

template<typename Integer>
FusionBasic::FusionBasic(const FusionComp<Integer>& FC){
    commutative = FC.commutative;
    use_modular_grading = false;
    candidate_given = FC.candidate_given;
    fusion_rank = FC.fusion_rank;
    fusion_type = FC.fusion_type;
    duality = FC.duality;
    type_and_duality_set = true;
    total_FPdim = FC.total_FPdim;
    fusion_type_string = FC.fusion_type_string;
    subring_base_key = FC.subring_base_key;
    type_automs_made = FC.type_automs_made;
    type_automs = FC.type_automs;
    if(FC.fusion_ring_map.nr_of_rows() > 0){
        assert(!using_renf<Integer>());
        convert(fusion_image_type, FC.fusion_image_type);
        convert(fusion_image_ring, FC.fusion_image_ring);
        convert(fusion_ring_map, FC.fusion_ring_map);
        fusion_image_duality = FC.fusion_image_duality;
        fusion_image_commutative = FC.fusion_image_commutative;
        fusion_image_type_string = FC.fusion_image_type_string;
    }

    //Automorphisms = FC.Automorphisms; // not needed at pressent
    //type_automs = F:type_automs;
    // automorphisms_mde = FC.automorphisms_mde;
}

void  FusionBasic::data_from_mpq_input(ifstream& cone_in){
    InputMap<mpq_class> input;
    map<NumParam::Param, long> num_param_input;
    map<PolyParam::Param, vector<string> > poly_param_input;
    OptionsHandler options;
    renf_class_shared number_field;
    input = readNormalizInput<mpq_class>(cone_in, options, num_param_input, poly_param_input,  number_field);
    read_data_from_input<mpq_class>(input);
}

#ifdef ENFNORMALIZ
void  FusionBasic::data_from_renf_input(ifstream& cone_in){
    InputMap<renf_elem_class> input;
    map<NumParam::Param, long> num_param_input;
    map<PolyParam::Param, vector<string> > poly_param_input;
    OptionsHandler options;
    renf_class_shared number_field;
    input = readNormalizInput<renf_elem_class>(cone_in, options, num_param_input, poly_param_input,  number_field);
    read_data_from_input<renf_elem_class>(input);
}
#endif

void FusionBasic::restrict_type_automs_to_grading(){

    vector<key_t> in_comp(fusion_rank);
    for(size_t k = 0; k < group_order; ++k){
        for(size_t i = 0; i < chosen_modular_grading[k].size(); ++i){
            if(!chosen_modular_grading[k][i])
                continue;
            in_comp[i] = k;
        }
    }
    vector<vector<shortkey_t> >selection;
    for(auto& p: type_automs){
        // cout << p;
        vector<key_t> image;
        for(size_t i = 0; i< chosen_modular_grading.size(); ++i){
            size_t first_elem = chosen_modular_grading[i].find_first();
            image.push_back(in_comp[p[first_elem]]);
        }
        // cout << image;
        bool restrictable = true;
        for(size_t k = 0; k < fusion_rank; ++k){
            if(in_comp[p[k]] != image[in_comp[k]] ){
                restrictable = false;
                break;
            }
        }

        if(restrictable)
            selection.push_back(p);
    }
    swap(selection, type_automs);
    if(verbose)
        verboseOutput() << type_automs.size() << " type permutations selected" << endl;
}

void FusionBasic::make_type_automs(){
    if(type_automs_made)
        return;
    if(libnormaliz::verbose)
        verboseOutput() << "Making type automorphisms" << endl;
    type_automs = make_all_full_permutations(fusion_type,duality, fusion_ring_map);
    if(verbose)
        verboseOutput() << type_automs.size() << " type automorphisms made" << endl;
    type_automs_made = true;
}

bool  FusionBasic::data_from_file(const string& file_name){
    bool number_field_input = false;
    bool fusion_input = false;
    ifstream cone_in(file_name);
    string test;
    while (cone_in.good()) {
        cone_in >> test;
        // cout << test << endl;
        if (test == "number_field") {
            number_field_input = true;
        }
        if(test == "fusion_type"){
                fusion_input = true;
        }
    }
    cone_in.close();

    if(!fusion_input && number_field_input)
        throw BadInputException("Number filed input must be of fusion type tor fusion compoutation");

    if(!fusion_input)
        return false;

    cone_in.open(file_name.c_str(), ifstream::in);

#ifndef ENFNORMALIZ
    if(number_field_input)
        throw BadInputException("Number field input only possible with e-antic");
#else
    if(number_field_input){
        data_from_renf_input(cone_in);
        return true
        ;
    }
#endif
    data_from_mpq_input(cone_in);
    return true;
}


void  FusionBasic::data_from_file_or_string(const string& our_fusion){

    string file_name = our_fusion;
    if(file_name.size() < 3 || file_name.substr(file_name.size()-2,3) != ".in"){
        file_name += ".in";
    }
    ifstream cone_in(file_name);

    bool got_data_from_file = false;

    if(cone_in.is_open()){
        cone_in.close();
        got_data_from_file = data_from_file(file_name);
    }

    if(!got_data_from_file)
        data_from_string(our_fusion, false);
}

// Note: test_only = false produces data ready for FusionComp (coincidence, commutative)
// test_only = true copies the type and the duality 1-1
pair<bool, bool> FusionBasic::data_from_string(const string& our_fusion, const bool return_on_failure) {

    bool dummy = false;
    if(verbose)
        verboseOutput() << "Trying to read fusion data from string " << our_fusion << endl;
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
        if(return_on_failure)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    if(open_bracket != close_bracket){
        if(return_on_failure)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    if(open_bracket == 0 || open_bracket > 2){
        if(return_on_failure)
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
        if(return_on_failure)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    vector<long> type_input;
    while(true){
        long nr;
        data >> nr;
        if(data.fail()){
            if(return_on_failure)
                return make_pair(false, dummy);
            throw BadInputException("String " + our_fusion +" not standard fusion");
        }
        if(nr < 1){
            if(return_on_failure)
                return make_pair(false, dummy);
            throw BadInputException("String " + our_fusion +" not standard fusion");
        }
        type_input.push_back(nr);
        data >> c;
        if(c == ']'){
            break;
        }
        else{
            if(c != ','){
                if(return_on_failure)
                    return make_pair(false, dummy);
                throw BadInputException("String " + our_fusion +" not standard fusion");
            }
        }

    }
    if(type_input[0] != 1){
        if(return_on_failure)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }
    fusion_rank = type_input.size();
    total_FPdim = 0;
    for(size_t i = 0; i< fusion_rank; ++i){
        double this_FPdim = convertTo<double>(type_input[i]);
        total_FPdim += this_FPdim * this_FPdim;
    }
    stringstream for_type;
    for_type << type_input;
    fusion_type_from_command = type_input;
    // cout << "RRRRRRRRRRRRRR " << fusion_type_from_command;
    fusion_type_string = for_type.str();
    if(!return_on_failure)
        fusion_type = fusion_coincidence_pattern(type_input);
    else{ // we have checked positivity
        fusion_type.resize(fusion_rank);
        for(size_t i = 0; i < fusion_rank; ++i)
            fusion_type[i] = type_input[i];
    }

    if(only_partition){
        return make_pair(true, true);
    }

    vector<long> duality_input;
    data >> c;
    if(c !='['){
        if(return_on_failure)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +" not standard fusion");
    }

    while(true){
        long nr;
        data >> nr;
        if(data.fail()){
            if(return_on_failure)
                return make_pair(false, dummy);
            throw BadInputException("String " + our_fusion +" not standard fusion");
        }
        if(nr == -1 || nr == -3){
            commutative = true;
            if(nr == -3)
                use_modular_grading = true;
            nr = 0;
        }
        if(nr == -2){
            commutative = true;
           use_modular_grading = true;
           nr = 0;
        }
        duality_input.push_back(nr);
        data >> c;
        if(c == ']'){
            break;
        }
        else{
            if(c != ','){
                if(return_on_failure)
                    return make_pair(false, dummy);
                throw BadInputException("String " + our_fusion +" not standard fusion");
            }
        }
    }

    if(!check_duality(duality_input, type_input)){
        if(return_on_failure)
            return make_pair(false, dummy);
        throw BadInputException("String " + our_fusion +": duality does not fit type");
    };

    duality.resize(fusion_rank);
    for(key_t i = 0; i< fusion_rank; ++i){
        duality[i] = duality_input[i];
    }
    type_and_duality_set = true;
    return make_pair(true, false);
}


template <typename Integer>
vector<dynamic_bitset> make_subsets_FPdim(const vector<Integer>& d, const Integer&part_FPdim, const vector<dynamic_bitset>& AllSubsets){

    vector<dynamic_bitset> subsets_FPdim;
    for(auto& S: AllSubsets){
        Integer test_FPdim = 0;
        for(size_t i = 0; i < S.size(); ++i){
            if(!S[i])
                continue;
            test_FPdim += d[i]*d[i];
            if(test_FPdim > part_FPdim)
                break;
        }
        if(test_FPdim == part_FPdim)
            subsets_FPdim.push_back(S);
    }
    return subsets_FPdim;
}


void all_partitions(vector<dynamic_bitset>& Partitions, dynamic_bitset partition, const vector<dynamic_bitset>& subsets_FPdim,
               size_t index, dynamic_bitset already_covered, const size_t& group_order, size_t level){

    if(level == group_order){
        Partitions.push_back(partition);
        return;
    }

    for(size_t i = index; i < subsets_FPdim.size(); ++i){
        dynamic_bitset intersect = already_covered & subsets_FPdim[i];
        if(intersect.any())
            continue;
        dynamic_bitset covered = already_covered | subsets_FPdim[i];
        dynamic_bitset new_partition = partition;
        new_partition[i] = 1;
        all_partitions(Partitions, new_partition, subsets_FPdim, i + 1, covered, group_order, level +1);
    }
}

vector<dynamic_bitset> sort_part(const vector<dynamic_bitset>& to_be_sorted){

    vector<dynamic_bitset> sorted;
    map<int,int> sort_map;
    for(size_t i = 0; i < to_be_sorted.size(); ++i)
        sort_map[to_be_sorted[i].find_first()] = i;
    for(auto& q: sort_map)
        sorted.push_back(to_be_sorted[q.second]);
    return sorted;
}

template <typename Integer>
vector<vector<dynamic_bitset> > make_FPdim_partitions(const vector<Integer>& d, const Integer& part_FPdim, const size_t& group_order,
                                                      vector<dynamic_bitset>& AllSubsets){

    vector<dynamic_bitset> subsets_FPdim = make_subsets_FPdim(d, part_FPdim, AllSubsets);
    vector<dynamic_bitset> Partitions;
    dynamic_bitset partition(subsets_FPdim.size());
    dynamic_bitset already_covered(d.size());

    all_partitions(Partitions, partition, subsets_FPdim, 0, already_covered, group_order, 0);

    vector<vector<dynamic_bitset> > FPdimParts;
    for(auto& p :Partitions){
        vector<key_t> selection = bitset_to_key(p);
        vector<dynamic_bitset> FPdimPart;
        vector<dynamic_bitset> FPdimPart_sorted;
        for(auto& c: selection)
           FPdimPart.push_back(subsets_FPdim[c]);
        // sort(FPdimPart.begin(), FPdimPart.end());

        FPdimParts.push_back(sort_part(FPdimPart));
    }
    /*
    cout << "----------" << endl;
    for(auto& p: FPdimParts){
        for(auto& q:p)
            cout << bitset_to_key(q);
        cout << "------------" << endl;
    }*/
    return FPdimParts;
}

vector<vector<dynamic_bitset> > FusionBasic::make_part_classes(const vector<vector<dynamic_bitset> >& GradPartitions){

    vector<vector<dynamic_bitset> > PartNormalForms;
    for(auto& part: GradPartitions){
        bool is_normal_form = true;
        for(auto& aut: type_automs){
            vector<dynamic_bitset> cand;
            for(auto& p: part){
                dynamic_bitset image(part.front().size());
                for(size_t i = 0; i < p.size(); ++i){
                    if(p[i])
                        image[aut[i]] = 1;
                }
                cand.push_back(image);
            }
            cand = sort_part(cand);
            if(cand < part){
                is_normal_form = false;
                break;
            }
        }
        if(is_normal_form)
            PartNormalForms.push_back(part);
    }
    return PartNormalForms;

}

void FusionBasic::make_grad_mult_table(){

    vector<dynamic_bitset>& part = chosen_modular_grading;

    vector<key_t> in_comp(fusion_rank);
    for(size_t k = 0; k < group_order; ++k){
        for(size_t i = 0; i < part[k].size(); ++i){
            if(!part[k][i])
                continue;
            in_comp[i] = k;
        }
    }

     vector< vector<int> > MultTable(group_order, vector<int>(group_order,-1));
     for(size_t i = 0; i < group_order; ++i){
         MultTable[0][i] = i;
         MultTable[i][0] = i;
     }
     vector<key_t> inverse(group_order);
     for(size_t i = 0; i < group_order; ++i){
         size_t k = part[i].find_first();
         inverse[i] = in_comp[duality[k]];
     }
     for(size_t i = 0; i < group_order; ++i){
        MultTable[inverse[i]][i] = 0;
        MultTable[i][inverse[i]] = 0;
     }
    for(size_t i = 1; i < group_order; ++i){
        dynamic_bitset defined(group_order);
        dynamic_bitset used(group_order);
        for(size_t j = 1; j < group_order; ++j){
            if(MultTable[i][j] >= 0){
                defined[j] = 1;
                used[MultTable[i][j]] = 1;
            }
        }
        for(size_t j = 1; j < group_order; ++j){
            if(defined[j])
                continue;
            int third = -2;
            for(size_t k = 1; k < group_order; ++k){
                if(k == i || k == j)
                    continue;
                if(used[k])
                    continue;
                third = k;
            }
            MultTable[i][j] = third;
            defined[j] = 1;
            used[third] =1;
        }
    }

    GradMultTable = MultTable;
}




bool FusionBasic::compatible_duality(const vector<dynamic_bitset >& parts){
    if(group_type == "C2"){
        dynamic_bitset duals(fusion_rank);
        for(size_t i = 0; i < fusion_rank; ++i){
            if(parts[1][i] != 0)
                duals[duality[i]] = 1;
        }
        if(duals != parts[1]){
            return false;
        }
        return true;
    }
    if(group_type == "C3"){
        dynamic_bitset duals(fusion_rank);
        for(size_t i = 0; i < fusion_rank; ++i){
            if(parts[1][i] != 0)
                duals[duality[i]] = 1;
        }
        if(duals != parts[2])
            return false;
    }
    if(group_type == "C2xC2"){
        for(size_t k = 1; k < 4; ++k){
            dynamic_bitset duals(fusion_rank);
            for(size_t i = 0; i < fusion_rank; ++i){
                if(parts[k][i] != 0)
                    duals[duality[i]] = 1;
            }
            if(duals != parts[k])
                return false;
        }
    }
    if(group_type == "C4"){
        long order_2_at = -1;
        for(size_t k = 1; k < 4; ++k){
            dynamic_bitset duals(fusion_rank);
            for(size_t i = 0; i < fusion_rank; ++i){
                if(parts[1][i] != 0)
                    duals[duality[i]] = 1;
            }
            if(duals == parts[k]){
                order_2_at = k;
                break;
            }
        }
        if(order_2_at == -1)
            return false;

        vector<long> grad_dual;
        for(size_t k = 1; k < 4; ++k){
            if(k != order_2_at)
                grad_dual.push_back(k);
        }

        dynamic_bitset duals(fusion_rank);
        for(size_t i = 0; i < fusion_rank; ++i){
            if(parts[grad_dual[0]][i] != 0)
                duals[duality[i]] = 1;
        }
        if(duals != parts[grad_dual[1]])
            return false;
    }
    return true;
}


//--------------------------------------------------------------
// FusionComp
//--------------------------------------------------------------

template <typename Integer>
FusionComp<Integer>::FusionComp(){
    initialize();
}

template <typename Integer>
FusionComp<Integer>::FusionComp(const FusionBasic& basic){
    initialize();
    fusion_rank = basic.fusion_rank;
    commutative = basic.commutative;
    use_modular_grading = basic.use_modular_grading;
    candidate_given = basic.candidate_given;
    fusion_type = basic.fusion_type;
    fusion_type_string = basic.fusion_type_string;
    duality = basic.duality;
    subring_base_key = basic.subring_base_key;
    type_and_duality_set = basic.type_and_duality_set;
    total_FPdim = basic.total_FPdim;
    // Automorphisms = basic.Automorphisms;
    type_automs = basic.type_automs;
    type_automs_made = basic.type_automs_made;
    // automorphisms_mde = basic.automorphisms_mde;
    chosen_modular_grading = basic.chosen_modular_grading;
    if(basic.fusion_ring_map.nr_of_rows() > 0){
        assert(!using_renf<Integer>());
        convert(fusion_image_type, basic.fusion_image_type);
        convert(fusion_image_ring, basic.fusion_image_ring);
        convert(fusion_ring_map, basic.fusion_ring_map);
        fusion_image_duality = basic.fusion_image_duality;
        fusion_image_commutative = basic.fusion_image_commutative;
        fusion_image_type_string = basic.fusion_image_type_string;
    }
}

template <typename Integer>
void FusionComp<Integer>::initialize(){
    check_simplicity = false;
    candidate_given = false;
    use_automorphisms = false;
    type_automs_made = false;
    // select_iso_classes = false;
    verbose = false;
    // activated = false;
    automorphisms_mde = false;
    type_and_duality_set =false;
    commutative = false;
    use_modular_grading = false;
    nr_coordinates = 0;
    total_FPdim = 0;
}

template <typename Integer>
void FusionComp<Integer>::make_automorphisms(){

    if(automorphisms_mde)
        return;

    make_CoordMap();

    /* cout <<  "Coord " << CoordMap.size() << endl;
    cout << "Type " << fusion_type << endl;
    cout << "duality " << duality;*/

    if(!type_automs_made){
        if(libnormaliz::verbose)
            verboseOutput() << "Making type automorphisms" << endl;
        type_automs = make_all_full_permutations(fusion_type,duality, fusion_ring_map);
        if(libnormaliz::verbose)
            verboseOutput() << type_automs.size() << " type automorphisms" << endl;
        type_automs_made =true;
    }

    /*cout << "type automs " << type_automs.size() << endl;
     cout << "selected_ind_tuples " << selected_ind_tuples.size() << endl;*/

    if(libnormaliz::verbose)
        verboseOutput() << "Making coordinate automorphisms" << endl;

    Automorphisms.resize(type_automs.size());

    bool skip_remaining = false;
    std::exception_ptr tmp_exception;

#pragma omp parallel for
    for(size_t i = 0; i< type_automs.size(); ++i){

        if(skip_remaining)
            continue;
        try{
            vector<shortkey_t> coord_perm(1); // must start with 0 !!!!
            for(auto& t: selected_ind_tuples){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

                vector<key_t> image;
                for(auto& c: t)
                    image.push_back(type_automs[i][c]);
                coord_perm.push_back(coord(image));
            }
            Automorphisms[i] = coord_perm;
    }
    catch (const std::exception&) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
#pragma omp flush(skip_remaining)
            }
    }

    if (!(tmp_exception == 0))
        std::rethrow_exception(tmp_exception);

    if(libnormaliz::verbose)
        verboseOutput() << "Fusion data automorphism group of order " << Automorphisms.size() << " computed" << endl;
    automorphisms_mde = true;
}




// checks whether sol is contained in the "subring". If so, "false" is returned,
// meaning "not sinmple"
template <typename Integer>
bool FusionComp<Integer>::simplicity_check(const vector<key_t>& subring, const vector<Integer>& sol){

    for(auto& c: subring){
        if(sol[c] != 0)
            return true;
    }
    return false;
}

// checks whether sol is contained in one of the "subrings". If so "false" is returned.
template <typename Integer>
bool FusionComp<Integer>::simplicity_check(const vector<vector<key_t> >& subrings, const vector<Integer>& sol){

    for(auto& sub: subrings){
        if(!simplicity_check(sub, sol)){
            return false;
        }
    }
    return true;
}

template <typename Integer>
void FusionComp<Integer>::set_options(const ConeProperties& ToCompute, const bool verb){

    verbose = verb;
    check_simplicity= ToCompute.test(ConeProperty::SimpleFusionRings);
    use_automorphisms = ToCompute.test(ConeProperty::FusionRings) || ToCompute.test(ConeProperty::SimpleFusionRings);
    if(check_simplicity)
        prepare_simplicity_check();
    if(use_automorphisms)
        make_automorphisms();
    // cout << Automorphisms;
}

template <typename Integer>
set<vector<key_t> >  FusionComp<Integer>::FrobRec(const vector<key_t>& ind_tuple){
    if(commutative)
        return FrobRec_12(ind_tuple);
    else
        return FrobRec_6(ind_tuple);
}


template <typename Integer>
set<vector<key_t> >  FusionComp<Integer>::FrobRec_6(const vector<key_t>& ind_tuple){

    assert(ind_tuple.size() == 3);
    key_t i,j,k;
    i = ind_tuple[0];
    j = ind_tuple[1];
    k = ind_tuple[2];
    set< vector<key_t> > F;
    // cout << "FR "<< duality;
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
set<vector<key_t> >  FusionComp<Integer>::FrobRec_12(const vector<key_t>& ind_tuple){

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
Integer FusionComp<Integer>::value(const vector<Integer>& ring, vector<key_t>& ind_tuple){

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
key_t FusionComp<Integer>::coord(vector<key_t>& ind_tuple){
    set<vector<key_t> > FR = FrobRec(ind_tuple);
    return coord(FR);
}

template <typename Integer>
key_t FusionComp<Integer>::coord_cone(vector<key_t>& ind_tuple){
    key_t coord_compute = coord(ind_tuple);
    if(coord_compute == 0)
        return nr_coordinates;
    return coord_compute -1;
}

template <typename Integer>
key_t FusionComp<Integer>::coord(set<vector<key_t> >& FR){
    return CoordMap[FR];
}


// makes the critical coordinates for the simplicity check
// bse_key is the vector of bases (by keys) of the potential subrings
template <typename Integer>
dynamic_bitset FusionComp<Integer>::critical_coords(const vector<key_t>& base_key){
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
void FusionComp<Integer>::make_all_ind_tuples(){
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
void FusionComp<Integer>::make_CoordMap(){

    if(CoordMap.size() > 0)
        return;

    make_all_ind_tuples();
    // cout << "ind_tuples " << all_ind_tuples.size() << endl;

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
void  FusionComp<Integer>::make_all_base_keys(){

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
bool FusionComp<Integer>::automs_compatible(const vector<key_t>& cand ) const{

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
void FusionComp<Integer>::prepare_simplicity_check(){
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
            throw BadInputException("Candidate sunbring for non-simplicity not invariant under automorphisms.");
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
Matrix<Integer> FusionComp<Integer>::do_select_simple_inner(const Matrix<Integer>& LattPoints){
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
Matrix<Integer> FusionComp<Integer>::do_select_simple(const Matrix<Integer>& LattPoints) const {

    if(LattPoints.nr_of_rows() == 0 || !select_simple)
        return LattPoints;

    FusionComp<Integer> work_fusion = *this;
    return work_fusion.do_select_simple_inner(LattPoints);
}

template <typename Integer>
vector<Integer> FusionComp<Integer>::normal_form_of(const vector<Integer>& solution) const{

    vector<Integer> min_conjugate = solution;
    bool first = true;
    for(auto& aa: Automorphisms){
        vector<Integer> conjugate(solution.size());
        for(size_t j = 0; j < aa.size(); ++ j){
            conjugate[j] = solution[aa[j]];
        }
        conjugate[0] = 1;
        if(first || conjugate < min_conjugate){
            min_conjugate = conjugate;
            first = false;
        }
    }

    return min_conjugate;;
}

/*
template <typename Integer>
Matrix<Integer> FusionComp<Integer>::do_iso_classes_inner(const Matrix<Integer>& LattPoints){

   if(nr_coordinates != LattPoints.nr_of_columns() - 1)
        throw BadInputException("Wrong number of coordinates in fusion data. Mismatch of duality or commutativity.");

    Matrix<Integer> IsoClasses;
    IsoClasses.resize(0,LattPoints.nr_of_columns());

   // if(cone_coord){ // homogenizing coordinate is the last now
        for(auto& aa: Automorphisms){  // and we are allowed to destroy the original Automorphisms
            vector<shortkey_t> modified = aa;
            v_cyclic_shift_left(modified, modified.size() -1);
            modified.resize(modified.size()-1);
            for(auto& c: modified)
                c--;
            aa = modified;
        }
    // }
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

*/

// We work with final format (last coordinate is homogenizing)
// This function protects *this
/* template <typename Integer>
Matrix<Integer> FusionComp<Integer>::do_iso_classes(const Matrix<Integer>& LattPoints)const {

    if(LattPoints.nr_of_rows() == 0)
        return LattPoints;

    FusionComp<Integer> work_fusion = *this;

    return work_fusion.do_iso_classes_inner(LattPoints);
} */

template <typename Integer>
Matrix<Integer> FusionComp<Integer>::make_add_constraints_for_grading(){

    Matrix<Integer> GradEqu(0, nr_coordinates + 1);
    vector<key_t> indices(3);
    vector<key_t> in_comp(fusion_rank);
    for(size_t k = 0; k < chosen_modular_grading.size(); ++k){
        for(size_t i = 0; i < chosen_modular_grading[k].size(); ++i){
            if(!chosen_modular_grading[k][i])
                continue;
            in_comp[i] = k;
        }
    }

    // cout << chosen_modular_grading.size() << endl;

    for(key_t i = 1; i < fusion_rank; ++i){
        indices[0] = i;
        for(key_t j = 1; j < fusion_rank; ++j){
            indices[1] = j;
            key_t product_comp = GradMultTable[in_comp[i]][in_comp[j]];
            // cout << i << " " << j << " " << product_comp << endl;
            for(key_t k = 1; k < fusion_rank; k++){
                indices[2] =k;
                if(!chosen_modular_grading[product_comp][k]){
                    // cout << "zero" << endl;
                    ZeroCoords.insert(indices);
                    vector<Integer> this_equ(nr_coordinates + 1);
                    this_equ[coord_cone(indices)] = 1;
                    assert(coord_cone(indices) < nr_coordinates + 1);
                    GradEqu.append(this_equ);
                }
            }
        }
    }
    GradEqu.remove_duplicate_and_zero_rows();
    if(libnormaliz::verbose)
        verboseOutput() << GradEqu.nr_of_rows() << " coordinates set to 0" << endl;
    Matrix<Integer>  Neg = GradEqu;
    Integer MinusOne = -1;
    Neg.scalar_multiplication(MinusOne);
    GradEqu.append(Neg);
    return GradEqu;
}


            /*

                indices[2] = k;
                bool add_equ = false;
                // multiplication inside neutral component
                if((i <= half_at && j <= half_at)
                            && k > half_at){
                    add_equ = true;
                }
                // multiplication of non-neutral comp by neutral comp from left or right
                if( ( (i <= half_at && j> half_at) || (i > half_at && j<= half_at) )
                            && k <= half_at){
                    add_equ = true;
                }
                // multiplication of non-neutral component by itself
                if((i >  half_at && j > half_at) && k > half_at){
                    add_equ = true;
                }
                if(add_equ){
                    // cout << coord(indices) << " --- " << coord_cone(indices) << " ---- " << indices;
                    // counter++;
                    vector<Integer> this_equ(nr_coordinates + 1);
                    this_equ[coord_cone(indices)] = 1;
                    assert(coord_cone(indices) < nr_coordinates + 1);
                    GradEqu.append(this_equ);
                }
            }
        }
    }


    // cout << "CCCCC " << counter << endl;
    GradEqu.remove_duplicate_and_zero_rows();
    // cout << "Zero coords " << GradEqu.nr_of_rows() << " of " << GradEqu.nr_of_columns() << endl;

    vector<Integer> test_v(GradEqu.nr_of_columns());
    test_v.back() = 1;
    for(size_t kkn = 0; kkn < GradEqu.nr_of_rows(); ++kkn)
        if(test_v == GradEqu[kkn])
            assert(false);
    return GradEqu;
}*/

template <typename Integer>
void write_inhom_eq_as_lp(const Matrix<Integer>& Equ){

    if(verbose)
        verboseOutput() << "Writing inhomogeneous equations in lp format" << endl;

    string file_name = global_project+ ".lp";
    ofstream lp_out(file_name);
    size_t lhs_dim = Equ.nr_of_columns() -1;
    lp_out << "max:  ;" << endl;
    /*for(size_t i = 0; i < lhs_dim; ++i){
        lp_out << " + x" + to_string(i+1);
    }
    lp_out << ";" << endl;*/
    for(size_t i = 0; i < Equ.nr_of_rows(); ++i){
        for(size_t j = 0; j < lhs_dim; ++j){
            if(Equ[i][j] == 0)
                continue;
            if (Equ[i][j] > 0){
                if(Equ[i][j] == 1){
                    lp_out << " x" << to_string(j+1);
                    continue;
                }
                lp_out << " + " << Equ[i][j] << " x" << to_string(j+1);
            }
           if (Equ[i][j] < 0){
                if(Equ[i][j] == -1){
                    lp_out << " - x" << to_string(j+1);
                    continue;
                }
                lp_out << " - " << Iabs<Integer>(Equ[i][j]) << " x" << to_string(j+1);
            }
        }
        lp_out << " = " << -Equ[i][lhs_dim] << ";" << endl;
    }
    for(size_t j = 0; j < lhs_dim; ++j){
        lp_out << "x" << to_string(j+1) + " >= 0;" << endl;
    }
    lp_out << "int ";
    for(size_t j = 0; j < lhs_dim - 1; ++j){
        lp_out << "x" << to_string(j+1) << ",";
    }
    lp_out << "x" << to_string(lhs_dim) << ";" << endl;
}

template <typename Integer>
vector<Integer> FusionComp<Integer>::make_linear_equation(const map<vector<key_t>, Integer>& components, const Integer& rhs){

    vector<Integer> this_equ(nr_coordinates + 1);
    this_equ.back() = - rhs;
    for(auto& c: components){
        auto indices = c.first;
        auto coeff = c.second;
        if(indices[2] == 0){
            if(indices[0] == duality[indices[1]])
                this_equ.back() += coeff;
            continue;
        }
        this_equ[coord_cone(indices)] += coeff;
    }
    return this_equ;
}

template <typename Integer>
Matrix<Integer> FusionComp<Integer>::make_linear_constraints(const vector<Integer>& d){

    if(libnormaliz::verbose)
        verboseOutput() << "Making linear constraints for fusion rings" << endl;


    make_CoordMap();

    Matrix<Integer> Equ(0, nr_coordinates + 1); // mudst accomodate right hand side in last coordinate

    vector<key_t> indices(3);
    for(key_t i = 1; i < fusion_rank; ++i){
        indices[0] = i;
        for(key_t j = 1; j < fusion_rank; ++j){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

            indices[1] = j;
            map<vector<key_t>, Integer> components;
            Integer rhs = d[i] * d[j];
            for(key_t k = 0; k < fusion_rank; ++k){
                indices[2] = k;
                components[indices] = d[k];
            }
            auto this_equ = make_linear_equation(components, rhs);

            Equ.append(this_equ);
        }
    }

    if(write_lp_file)
        write_inhom_eq_as_lp(Equ);
    // Equ.print(global_project,"equ");

    Equ.remove_duplicate_and_zero_rows();
    if(libnormaliz::verbose)
        verboseOutput() << "Made " << Equ.nr_of_rows() << " inhom linear equations in " << Equ.nr_of_columns() -1 << " unknowns " << endl;

    // Equ.debug_print();
    return Equ;
}

template <typename Integer>
Matrix<Integer> FusionComp<Integer>::make_homomorphism_constraints(){

    if(libnormaliz::verbose)
        verboseOutput() << "Making homomorphism constraints for map of fusion rings" << endl;

    make_CoordMap();

    size_t fusion_image_rank = fusion_image_type.size();

    FusionBasic Image;
    Image.fusion_type = fusion_image_type;
    Image.duality = fusion_image_duality;
    Image.fusion_rank = fusion_image_rank;
    Image.commutative = commutative;

    FusionComp<Integer> ImageComp(Image);
    ImageComp.make_CoordMap();
    auto Tables = ImageComp.make_all_data_tables(fusion_image_ring);

    /* for(auto& M: Tables)
        M.debug_print();
    exit(0);*/

    vector<vector<vector<Integer> > > RHS;
    RHS.resize(fusion_rank);
    for(auto& r: RHS){
        r.resize(fusion_rank);
        for(auto& s: r)
            s.resize(fusion_image_rank);
    }

    for(size_t i = 0; i < fusion_rank; ++i){
        for(size_t j = 0; j < fusion_rank; ++j){
            for(size_t t = 0; t < fusion_image_rank; ++t){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

                Integer S = 0;
                for(size_t l = 0; l < fusion_image_rank; ++l){
                    for(size_t s = 0; s < fusion_image_rank; ++s)
                        S += fusion_ring_map[i][l] * fusion_ring_map[j][s] * Tables[l][s][t];
                }
                RHS[i][j][t] = S;
            }
        }
    }

    Matrix<Integer> Equ(0, nr_coordinates + 1); // mudst accomodate right hand side in last coordinate

    vector<key_t> indices(3);

    // the nneutral element is always mapped to the neutral element
    // So we can start the loops at 1.
    for(size_t i = 1; i < fusion_rank; ++i){
        indices[0] = i;
        for(size_t j = 1; j < fusion_rank; ++j){
            indices[1] = j;
            for(size_t t = 0; t < fusion_image_rank; ++t){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

                Integer rhs = RHS[i][j][t];
                map<vector<key_t>, Integer> components;
                for(size_t k = 0; k < fusion_rank; ++k){
                    indices[2] = k;
                    components[indices] = fusion_ring_map[k][t];
                }

                auto this_equ = make_linear_equation(components, rhs);
                Equ.append(this_equ);
            }
        }
    }

    // Equ.debug_print('+');

    Matrix<Integer> Help = Equ;
    Integer MinusOne = -1;
    Help.scalar_multiplication(MinusOne);
    Equ.append(Help);
    return Equ;
}

template <typename Integer>
Matrix<Integer> FusionComp<Integer>::make_linear_constraints_partition(const vector<Integer>& d,
                                                                       const vector<long>& card){
    make_CoordMap();

    /* cout << "DDDD " << d;
    cout << "CCCC " << card; */

    if(libnormaliz::verbose)
        verboseOutput() << "Making linear constraints for fusion rings partition" << endl;

    Matrix<Integer> Equ(0, nr_coordinates + 1); // mudst accomodate right hand side in last coordinate

    vector<key_t> indices(3);
    for(key_t i = 1; i < fusion_rank; ++i){
        indices[0] = i;
        for(key_t j = 1; j < fusion_rank; ++j){

            INTERRUPT_COMPUTATION_BY_EXCEPTION

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
    if(libnormaliz::verbose)
        verboseOutput() << "Made " << Equ.nr_of_rows() << " inhom linear equations in " << Equ.nr_of_columns() -1 << " unknowns " << endl;

    if(write_lp_file)
        write_inhom_eq_as_lp(Equ);

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

    if(summand.first == 0)
        return;

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
pair<Integer, vector<key_t> >  FusionComp<Integer>::term(const key_t& i, const key_t& j, const key_t& k){

    vector<key_t> indices = {i,j,k};
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
        exponent.push_back(coord(indices));
    }

    if(ZeroCoords.find(indices) != ZeroCoords.end())
        coeff = 0;

    return make_pair(coeff, exponent);
}

template <typename Integer>
set<map<vector<key_t>, Integer> > FusionComp<Integer>::make_associativity_constraints(){

    if(libnormaliz::verbose)
        verboseOutput() << "Making associativity constraints for fusion rings" << endl;

    make_CoordMap();

    // we produce associativity_equations b_i(b_j b_k) = (b_i b_j)b_k
    // parametrized by thge base vector b_t
    // if one of i,j,k is 0, no need for it since
    // the neutral element is automatically "associative"

    set<map<vector<key_t>, Integer> > Polys;

    for(key_t i = 1; i< fusion_rank; ++i){
        for(key_t j = 1; j < fusion_rank; ++j){
            for(key_t k = 1; k < fusion_rank; ++k){

                INTERRUPT_COMPUTATION_BY_EXCEPTION

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

    if(libnormaliz::verbose)
        verboseOutput() << "Made " << Polys.size() << " associativity constraints for fusion rings" << endl;

    return Polys;
}

void FusionBasic::do_write_input_file(InputMap<mpq_class>&  input) const{
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

void make_input_from_fusion_data(const FusionBasic& FusionInput, InputMap<mpq_class>&  input, const bool write_input_file){

    Matrix<mpq_class> TypeInput(1, FusionInput.fusion_rank);
    // cout << "TTTTTTT " << FusionInput.fusion_type_from_command;
    convert(TypeInput[0], FusionInput.fusion_type_from_command);
    vector<long> bridge(FusionInput.fusion_rank);
    for(size_t i = 0; i< bridge.size(); ++i)
        bridge[i] = FusionInput.duality[i];
    Matrix<mpq_class> DualityInput(1, FusionInput.fusion_rank);
    convert(DualityInput[0], bridge);
    if(FusionInput.commutative)
        DualityInput[0][0] = -1;
    if(FusionInput.use_modular_grading)
        DualityInput[0][0] -= 2;
    input[Type::fusion_type] = TypeInput;
    input[Type::fusion_duality] = DualityInput;
    if(write_input_file){
        FusionInput.do_write_input_file(input);
    }
}

void make_partition_input_from_fusion_data(const FusionBasic& FusionInput,InputMap<mpq_class>&  input,  const bool write_input_file){

    Matrix<mpq_class> TypeInput(1, FusionInput.fusion_rank);
    // cout << "TTTTTTT " << FusionInput.fusion_type_from_command;
    convert(TypeInput[0], FusionInput.fusion_type_from_command);
    input[Type::fusion_type_for_partition] = TypeInput;
    if(write_input_file){
        FusionInput.do_write_input_file(input);
    }
}

template <typename Integer>
Matrix<Integer> FusionComp<Integer>::data_table(const vector<Integer>& ring, const size_t i){

    Matrix<Integer> Table(fusion_rank, fusion_rank);

    for(key_t k = 0; k < fusion_rank; k++){
        for(key_t j= 0; j < fusion_rank; j++){
            key_t ii = i;
            vector<key_t>ind_tuple = {ii, j, k};
            Table[j][k] = value(ring, ind_tuple);
        }
    }
    // Table.debug_print('+');
    return Table;
}


template <typename Integer>
vector<Matrix<Integer> > FusionComp<Integer>::make_all_data_tables(const vector<Integer>& ring){

    vector<Matrix<Integer> > Tables;

    for(size_t i = 0; i <fusion_rank; ++i){
        Tables.push_back(data_table(ring, i));
        //Tables.back().debug_print('+');
    }
    return Tables;
}

template <typename Integer>
void FusionComp<Integer>::tables_for_all_rings(const Matrix<Integer>& rings){

    make_CoordMap();

    // vector<vector<Matrix<Integer> > > AllTables;
    for(size_t i = 0; i < rings.nr_of_rows(); ++i)
        AllTables.push_back(make_all_data_tables(rings[i]));
}

template <typename Integer>
void write_vec_vec_Mat(vector<vector<Matrix<Integer> > > AllTables, ostream& table_out){

    table_out << "[" << endl;
    for(size_t kk = 0; kk < AllTables.size(); kk++){
        table_out << "  [" << endl;
        vector<Matrix<Integer> > Tables = AllTables[kk]; // for a fixed ring
        for(size_t nn = 0; nn < Tables.size(); ++nn){
            Matrix<Integer> table = Tables[nn];
            table_out << "    [" << endl;
            for(size_t jj = 0; jj < table.nr_of_rows(); ++jj){
                table_out << "      [";
                for(size_t mm = 0; mm < table.nr_of_columns(); ++mm){
                    table_out << table[jj][mm];
                    if(mm < table.nr_of_columns() - 1)
                        table_out << ",";
                    else{
                        if(jj < table.nr_of_rows() -1)
                            table_out << "]," << endl;
                        else
                            table_out << "]" << endl;
                    }
                }
            }
            if(nn == Tables.size() - 1)
                table_out << "    ]" << endl;
            else
                table_out << "    ]," << endl;
        }
        if(kk == AllTables.size() -1)
            table_out << "  ]" << endl;
        else
            table_out << "  ]," << endl;
    }
    table_out << "]" << endl;


}

template <typename Integer>
void FusionComp<Integer>::write_all_data_tables(const Matrix<Integer>& rings, ostream& table_out){

    tables_for_all_rings(rings);
    write_vec_vec_Mat(AllTables, table_out);
}

/*
template <typename Integer>
Matrix<Integer> FusionComp<Integer>::data_table(const vector<Integer>& ring, const size_t i){

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
vector<Matrix<Integer> > FusionComp<Integer>::make_all_data_tables(const vector<Integer>& ring){

    vector<Matrix<Integer> > Tables;

    for(size_t i = 0; i <fusion_rank; ++i){
        Tables.push_back(data_table(ring, i));
    }
    return Tables;
}

template <typename Integer>
void FusionComp<Integer>::tables_for_all_rings(const Matrix<Integer>& rings){

    make_CoordMap();

    vector<vector<Matrix<Integer> > > AllTables;
    for(size_t i = 0; i < rings.nr_of_rows(); ++i)
        AllTables.push_back(make_all_data_tables(rings[i]));
}

template <typename Integer>
void FusionComp<Integer>::write_all_data_tables(const Matrix<Integer>& rings, ostream& table_out){

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
*/
//-------------------------------------------------------------------------------
// helper for fusion rings


// bridge to cone
template <typename Integer>
Matrix<Integer> select_simple(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb){

    FusionComp<Integer> fusion;
    fusion.set_options(ToCompute, verb);
    // fusion.read_data(false); // falsae = a posteriori
    return fusion.do_select_simple(LattPoints);
}


template <typename Integer>
void split_into_simple_and_nonsimple(const FusionBasic& basic, Matrix<Integer>& SimpleFusionRings, Matrix<Integer>& NonsimpleFusionRings, const Matrix<Integer>& FusionRings, bool verb){

    if(verb)
        verboseOutput() << "Splitting fusion rings into simple and nonsimple" << endl;

    if(FusionRings.nr_of_rows() == 0){
        if(verb)
            verboseOutput() << "No fusion rings given" << endl;
        return;
    }

    FusionComp<Integer> fusion(basic);
    // cout << fusion.fusion_rank << endl;
    fusion.select_simple = true;
    // fusion.activated = true;
    fusion.verbose = false;
    fusion.prepare_simplicity_check();
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

/*
template <typename Integer>
Matrix<Integer> fusion_iso_classes(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb){

    FusionComp<Integer> fusion;
    fusion.set_options(ToCompute, verb);
    fusion.read_data(false); // falsae = a posteriori
    return fusion.do_iso_classes(LattPoints);
}
*/

/*

Matrix<long long> extract_latt_points_from_out(ifstream& in_out){

    size_t nr_points;
    in_out >> nr_points;
    string s;
    in_out >> s;
    if(s != "lattice" && s != "fusion" && s!= "simple")
        throw BadInputException("out file not suitable for extraction of simple fusion rtings");
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
    FusionBasic blabla;
    split_into_simple_and_nonsimple(blabla, SimpleFusionRings, NonsimpleFusionRings, LattPoints, verbose);

    string name = global_project + ".fusion";
    write_fusion_files(blabla, name, non_simple_fusion_rings, true, embdim,
                            SimpleFusionRings, NonsimpleFusionRings,false,false);
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
        throw BadInputException("No project defined");
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
*/

/*
pair<bool, bool>  FusionBasic::read_data() {

    auto dummy = make_pair(true,true);

    return dummy;
}
*/




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

    FusionComp<Integer> partition_fusion;
    partition_fusion.fusion_type = identity_key(blocks.size());
    partition_fusion.duality = identity_key(blocks.size());
    partition_fusion.commutative = true; // doesn't matter since duality = id
    partition_fusion.fusion_rank = d.size();

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
void FusionComp<Integer>::set_global_fusion_data(){
    assert(false);
}

template <>
void FusionComp<long long>::set_global_fusion_data(){
    fusion_type_coinc_from_input = fusion_type;
    fusion_type_from_input = fusion_type_string;
    fusion_duality_from_input = duality;
    fusion_commutative_from_input = commutative;
}
*/

template class FusionComp<mpz_class>;
template class FusionComp<long long>;
template class FusionComp<long>;
#ifdef ENFNORMALIZ
template class FusionComp<renf_elem_class>;
#endif

/*
template class FusionBasic<mpz_class>;
template class FusionBasic<long long>;
template class FusionBasic<long>;
#ifdef ENFNORMALIZ
template class FusionBasic<renf_elem_class>;
#endif
*/

template Matrix<long> select_simple(const Matrix<long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<long long> select_simple(const Matrix<long long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<mpz_class> select_simple(const Matrix<mpz_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#ifdef ENFNORMALIZ
template Matrix<renf_elem_class> select_simple(const Matrix<renf_elem_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#endif

template void split_into_simple_and_nonsimple(const FusionBasic& basic, Matrix<long long>& SimpleFusionRings, Matrix<long long>& NonsimpleFusionRings, const Matrix<long long>& FusionRings, bool verb);
template void split_into_simple_and_nonsimple(const FusionBasic& basic, Matrix<long>& SimpleFusionRings, Matrix<long>& NonsimpleFusionRings, const Matrix<long>& FusionRings, bool verb);
template void split_into_simple_and_nonsimple(const FusionBasic& basic, Matrix<mpz_class>& SimpleFusionRings, Matrix<mpz_class>& NonsimpleFusionRings, const Matrix<mpz_class>& FusionRings, bool verb);
#ifdef ENFNORMALIZ
template void split_into_simple_and_nonsimple(const FusionBasic& basic, Matrix<renf_elem_class>& SimpleFusionRings, Matrix<renf_elem_class>& NonsimpleFusionRings, const Matrix<renf_elem_class>& FusionRings, bool verb);
#endif

/*
template Matrix<long> fusion_iso_classes(const Matrix<long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<long long> fusion_iso_classes(const Matrix<long long>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template Matrix<mpz_class> fusion_iso_classes(const Matrix<mpz_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#ifdef ENFNORMALIZ
template Matrix<renf_elem_class> fusion_iso_classes(const Matrix<renf_elem_class>& LattPoints, const ConeProperties& ToCompute, const bool verb);
#endif
*/

template vector<key_t> fusion_coincidence_pattern(const vector<long>& v);
template vector<key_t> fusion_coincidence_pattern(const vector<long long>& v);
template vector<key_t> fusion_coincidence_pattern(const vector<mpz_class>& v);
#ifdef ENFNORMALIZ
template vector<key_t> fusion_coincidence_pattern(const vector<renf_elem_class>& v);
#endif

/*
template bool FusionBasic::make_gradings(vector<long> d);
template bool FusionBasic::make_gradings(vector<long long> d);
template bool FusionBasic::make_gradings(vector<mpz_class> d);
#ifdef ENFNORMALIZ
template bool FusionBasic::make_gradings(vector<renf_elem_class> d);
#endif
*/

/*
template void make_full_input<long>(const FusionBasic& FusionInput, InputMap<long>& input_data, set<map<vector<key_t>, long> >& Polys);
template void make_full_input<long long>(const FusionBasic& FusionInput, InputMap<long long>& input_data, set<map<vector<key_t>, long long> >& Polys);
template void make_full_input<mpz_class>(const FusionBasic& FusionInput, InputMap<mpz_class>& input_data, set<map<vector<key_t>, mpz_class> >& Polys);
#ifdef ENFNORMALIZ
template void make_full_input<renf_elem_class>(const FusionBasic& FusionInput, InputMap<renf_elem_class>& input_data, set<map<vector<key_t>, renf_elem_class> >& Polys);
#endif
*/

template void make_full_input_partition(InputMap<long>& input_data);
template void make_full_input_partition(InputMap<long long>& input_data);
template void make_full_input_partition(InputMap<mpz_class>& input_data);
#ifdef ENFNORMALIZ
template void make_full_input_partition(InputMap<renf_elem_class>& input_data);
#endif

template void write_vec_vec_Mat(vector<vector<Matrix<long> > > AllTables, ostream& table_out);
template void write_vec_vec_Mat(vector<vector<Matrix<long long> > > AllTables, ostream& table_out);
template void write_vec_vec_Mat(vector<vector<Matrix<mpz_class> > > AllTables, ostream& table_out);
#ifdef ENFNORMALIZ
template void write_vec_vec_Mat(vector<vector<Matrix<renf_elem_class> > > AllTables, ostream& table_out);
#endif


template vector<vector<dynamic_bitset> > make_FPdim_partitions(const vector<long>& d, const long& part_FPdim, const size_t& group_order, vector<dynamic_bitset>& AllSubsets);
template vector<vector<dynamic_bitset> > make_FPdim_partitions(const vector<long long>& d, const long long& part_FPdim, const size_t& group_order, vector<dynamic_bitset>& AllSubsets);
template vector<vector<dynamic_bitset> > make_FPdim_partitions(const vector<mpz_class>& d, const mpz_class& part_FPdim, const size_t& group_order, vector<dynamic_bitset>& AllSubsets);
#ifdef ENFNORMALIZ
template vector<vector<dynamic_bitset> > make_FPdim_partitions(const vector<renf_elem_class>& d, const renf_elem_class& part_FPdim, const size_t& group_order, vector<dynamic_bitset>& AllSubsets);
#endif



} // namespace

/*
 *     if(check_fusion_grading){
        vector<vector<dynamic_bitset> > GradPartitions;
        bool has_grading = make_grading(d, GradPartitions, GradMultTable);
        if(!has_grading){
            if(running_input_file){
                FusionBasic aux_fusion(*this);
                Matrix<Integer> Empty(0,nr_coordinates);
                if(verbose)
                    verboseOutput() << "Writing output file with 0 fusion rings" << endl;
                write_fusion_files(aux_fusion,global_project,true,true,nr_coordinates, Empty, Empty,false);
                throw NoComputationException("No modular fusion grading found");
            }
        }
    }

    Matrix<Integer> GradEqu(0, nr_coordinates + 1);
    half_at = -1;
    if(use_modular_grading){
        find_grading(d);
        GradEqu = make_add_constraints_for_grading(d);
    }
*/
