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

#ifndef LIBNORMALIZ_FUSION_H_
#define LIBNORMALIZ_FUSION_H_

#include <vector>
#include <list>
#include <map>
#include <set>
#include <fstream>

#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/input_type.h"

namespace libnormaliz {
using std::vector;
using std::map;
using std::list;
using std::set;
using std::ifstream;
using std::pair;

class FusionBasic {

public:

    bool commutative;
    bool Z_2_graded;
    bool candidate_given;
    bool type_and_duality_set;

    size_t fusion_rank;
    vector<key_t> fusion_type;
    vector<long> fusion_type_from_command;
    string fusion_type_string;
    vector<key_t> duality;
    vector<key_t> subring_base_key;

    double total_FPdim;
    long half_at;

    // pair<bool, bool> read_data(const bool only_test);

    FusionBasic();

    template <typename Integer>
    void read_data_from_input(InputMap<Integer>& input_data);

    void  data_from_renf_input(ifstream& cone_in);
    void  data_from_mpq_input(ifstream& cone_in);
    void  data_from_file_or_string(const string& our_fusion);
    bool  data_from_file(const string& file_name);

    pair<bool, bool> data_from_string(const string& our_fusion, const bool return_on_failure);

    void do_write_input_file(InputMap<mpq_class>&  input) const;
};

template <typename Integer>
class FusionComp    {
    template <typename, typename>
    friend class ProjectAndLift;

public:

    bool activated;
    bool type_and_duality_set;

    bool verbose;

    bool commutative;
    bool Z_2_graded;

    bool check_simplicity;
    bool select_simple;
    bool candidate_given;

    bool use_automorphisms;
    bool select_iso_classes;
    bool write_mult_tables;

    size_t nr_coordinates;

    size_t fusion_rank;
    vector<key_t> fusion_type; // only coincidence pattern
    string fusion_type_string;
    vector<key_t> duality;

    double total_FPdim;

    long half_at; // temporarily used for ZZ_2-gradings

    void initialize();
    void import_global_data();
    vector<vector<vector<key_t> > > all_critical_coords_keys;
    vector<vector<key_t> > coords_to_check_key;
    vector<dynamic_bitset> coords_to_check_ind;
    vector<vector<key_t> > all_ind_tuples;
    vector<vector<key_t> > selected_ind_tuples; // the lex smallest in each FrobRec set
    map<set<vector<key_t> >, key_t> CoordMap;

    vector<vector<key_t> > Automorphisms;
    vector<dynamic_bitset> Orbits;

    vector<vector<Matrix<Integer> > > AllTables;

    FusionComp();
    FusionComp(const FusionBasic&);
    void set_options(const ConeProperties& ToCompute, const bool verb);
    //void read_data_from_file();

    // coordinates
    void make_CoordMap();
    set<vector<key_t> > FrobRec(const vector<key_t>& ind_tuple);
    set<vector<key_t> > FrobRec_6(const vector<key_t>& ind_tuple);
    set<vector<key_t> > FrobRec_12(const vector<key_t>& ind_tuple);
    key_t coord(set<vector<key_t> >& FR);
    key_t coord(vector<key_t>& ind_tuple);
    key_t coord_cone(vector<key_t>& ind_tuple);
    void make_all_ind_tuples();
    Integer value(const vector<Integer>& ring, vector<key_t>& ind_tuple);

    // for simplicity check
    vector<key_t> subring_base_key;
    dynamic_bitset critical_coords(const vector<key_t>& base_key);
    vector<vector<key_t> > all_base_keys;
    void make_all_base_keys();
    void prepare_simplicity_check();
    Matrix<Integer> do_select_simple_inner(const Matrix<Integer>& LattPoints);
    Matrix<Integer> do_iso_classes_inner(const Matrix<Integer>& LattPoints);
    Matrix<Integer> do_select_simple(const Matrix<Integer>& LattPoints) const;
    Matrix<Integer> do_iso_classes(const Matrix<Integer>& LattPoints) const;
    bool simplicity_check(const vector<key_t>& subring, const vector<Integer>& sol);
    bool simplicity_check(const vector<vector<key_t> >& subrings, const vector<Integer>& sol);
    bool automs_compatible(const vector<key_t>& cand) const;

    // for automosphisms
    void make_automorphisms();
    vector<Integer> norrmal_form(const vector<Integer> lattice_point);

    Matrix<Integer> make_linear_constraints(const vector<Integer>& d);
    Matrix<Integer> make_linear_constraints_partition(const vector<Integer>& d,
                                            const vector<long>& card);
    pair<Integer, vector<key_t> >  term(const key_t& i, const key_t& j, const key_t& k);
    set<map<vector<key_t>, Integer> > make_associativity_constraints();
    // void set_global_fusion_data();

    void find_grading(const vector<Integer>& d);
    Matrix<Integer> make_add_constraints_for_grading(const vector<Integer>& d);

    void write_all_data_tables(const Matrix<Integer>& rings, ostream& table_out);
    void tables_for_all_rings(const Matrix<Integer>& rings);
    vector<Matrix<Integer> > make_all_data_tables(const vector<Integer>& ring);
    Matrix<Integer> data_table(const vector<Integer>& ring, const size_t i);
};

// helpers
Matrix<long long> extract_latt_points_from_out(ifstream& in_out);
template <typename Integer>
Matrix<Integer> select_simple(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb);
template <typename Integer>
Matrix<Integer> fusion_iso_classes(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb);
//void select_simple_fusion_rings();
template <typename Integer>
void split_into_simple_and_nonsimple(const FusionBasic& basic, Matrix<Integer>& SimpleFusionRings, Matrix<Integer>& NonsimpleFusionRings, const Matrix<Integer>& FusionRings, bool verb);

template <typename Integer>
void make_full_input(const FusionBasic& FusionInput, InputMap<Integer>& input_data, set<map<vector<key_t>, Integer> >& Polys);
template <typename Integer>
void make_full_input_partition(InputMap<Integer>& input_data);

void make_input_from_fusion_data(const FusionBasic& FusionInput, InputMap<mpq_class>&  input, const bool write_input_file);
void make_partition_input_from_fusion_data(const FusionBasic& FusionInput, InputMap<mpq_class>&  input, const bool write_input_file);

vector<dynamic_bitset> make_all_subsets(const size_t card);
vector<vector<key_t> > make_all_permutations(size_t n);
vector<vector<key_t> > collect_coincidence_subset_keys(const vector<key_t>& type);
vector<vector<key_t> > make_all_permutations(const vector<key_t>& v);
vector<vector<key_t> > make_all_permutations(const vector<key_t>& type, const vector<key_t>& duality, const long& half_at);

// void remove_global_fusion_data();

// void post_process_fusion(const vector<string>& command_line_items);

template <typename Integer>
vector<key_t> fusion_coincidence_pattern(const vector<Integer>& v);

template<typename Integer>
bool check_duality(vector<Integer> test_duality, const vector<Integer>& test_type){
    if(test_duality[0] != 0 && test_duality[0] != -1)
        return false;
    test_duality[0] = 0;
    for(Integer i = 0; i< test_duality.size(); ++i){
        if(test_duality[i] < 0 || test_duality[i] >= test_duality.size())
            return false;
        if(test_duality[test_duality[i]] != i)
            return false;
        if(test_type[i] != test_type[test_duality[i]])
            return false;
    }
    return true;
}

// Note: the following routine must work for renf_elem_class
template <typename Integer>
void FusionBasic::read_data_from_input(InputMap<Integer>& input_data){

    vector<Integer> full_type = input_data[Type::fusion_type][0];
    total_FPdim = 0;
    for(size_t i = 0; i< full_type.size(); ++i){
        double this_FPdim = convertTo<double>(full_type[i]);
        total_FPdim += this_FPdim * this_FPdim;
    }
    // cout << "FULL " << full_type;
    fusion_type = fusion_coincidence_pattern(full_type);
    // cout << "COINC " << fusion_type_coinc_from_input;
    fusion_rank = full_type.size();
    stringstream for_type;
    for_type << full_type;
    fusion_type_string = for_type.str();

    commutative = false;
    if(contains(input_data, Type::fusion_duality)){
        vector<Integer> prel_duality = input_data[Type::fusion_duality][0];
        // std::cout << "PREL " << prel_duality; //" -- " << prel_duality.size() << " -- " <<  fusion_rank_from_input << endl;
        if(prel_duality.size() != fusion_rank || (prel_duality[0] != 0 && prel_duality[0] != -1 && prel_duality[0] != -2 && prel_duality[0] != -3))
            throw BadInputException("Fusion duality corrupt");
        if(prel_duality[0] == -1 || prel_duality[0] == -3) {
            commutative = true;
            if(prel_duality[0] == -3)
                Z_2_graded = true;
            prel_duality[0] = 0;
        }
        if(prel_duality[0] == -2) {
            Z_2_graded = true;
            prel_duality[0] = 0;
        }
        duality.resize(fusion_rank);
        for(key_t i = 0; i < fusion_rank; ++i){
            bool in_range = false;
            for(long j = 0; j < fusion_rank; ++j){  // conversion by counting
                if(convertTo<Integer>(j) == prel_duality[i]){
                    duality[i] = j;
                    in_range = true;
                    break;
                }
            }
            if(!in_range)
                throw BadInputException("Fusion duality out of range");
        }
        if(key_to_bitset(duality, fusion_rank).count() != fusion_rank)
            throw BadInputException("Fusion duality has repeated entries");
        if(!check_duality<key_t>(duality, fusion_type))
            throw BadInputException("Fusion duality does not fit type");
    }
    else{
        duality = identity_key(fusion_rank);
    }

    if(contains(input_data, Type::candidate_subring)){
        dynamic_bitset cand_indicator(input_data[Type::candidate_subring][0].size());
        if(cand_indicator.size() != fusion_rank)
            throw BadInputException("Candidate subring has wrong size");
        for(size_t i = 0; i < cand_indicator.size(); ++i){
            if(input_data[Type::candidate_subring][0][i] == 0){
                continue;
            }
            if(input_data[Type::candidate_subring][0][i] == 1){
                cand_indicator[i] = 1;
                continue;
            }
            throw BadInputException("Candidate subring not 0-1");
        }
        if(!cand_indicator[0] || cand_indicator.count() <=1 || cand_indicator.count() == full_type.size())
            throw BadInputException("Candidate subring corrupt");
        for(size_t i = 0; i < cand_indicator.size(); ++i){
            if(cand_indicator[i] && !cand_indicator[duality[i]])
                throw BadInputException("Candidate subring not closed iunder duality");
        }
        subring_base_key = bitset_to_key(cand_indicator);
    }

    type_and_duality_set = true;
}

template <typename Integer>
void make_full_input(FusionBasic& FusionInput, InputMap<Integer>& input_data, set<map<vector<key_t>, Integer> >& Polys) {

    FusionInput.read_data_from_input(input_data);
    FusionComp<Integer> OurFusion(FusionInput);
    vector<Integer> full_type = input_data[Type::fusion_type][0];
    Matrix<Integer> Equ = OurFusion.make_linear_constraints(full_type);
    FusionInput.half_at = OurFusion.half_at;
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

    Polys = OurFusion.make_associativity_constraints();

}

}  // end namespace libnormaliz



#endif /* FUSION_H_ */
