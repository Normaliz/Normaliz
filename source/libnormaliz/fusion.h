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
#include "libnormaliz/list_and_map_operations.h"

namespace libnormaliz {
using std::vector;
using std::map;
using std::list;
using std::set;
using std::ifstream;
using std::pair;

template <typename Integer>
class FusionComp;

class FusionBasic {

public:

    bool commutative;
    bool use_modular_grading;
    bool candidate_given;
    bool type_and_duality_set;

    size_t fusion_rank;
    vector<key_t> fusion_type;
    vector<long> fusion_type_from_command;
    string fusion_type_string;
    vector<key_t> duality;
    vector<key_t> subring_base_key;

    double total_FPdim;

    vector<vector<dynamic_bitset> > ModularGradings;
    size_t group_order;
    string group_type;
    vector<vector<int> > GradMultTable;
    vector<dynamic_bitset> chosen_modular_grading;

    vector<key_t> fusion_image_type;
    vector<key_t> fusion_image_duality;
    string fusion_image_type_string;
    vector<long> fusion_image_ring;
    Matrix<long> fusion_ring_map;
    bool fusion_image_commutative;

    vector<vector<shortkey_t> > type_automs; // permutations of the basis vectors
    bool type_automs_made;

    // pair<bool, bool> read_data(const bool only_test);

    FusionBasic();

    template<typename Integer>
    FusionBasic(const FusionComp<Integer>& FC);

    template <typename Integer>
    void read_data_from_input(InputMap<Integer>& input_data);

    void  data_from_renf_input(ifstream& cone_in);
    void  data_from_mpq_input(ifstream& cone_in);
    void  data_from_file_or_string(const string& our_fusion);
    bool  data_from_file(const string& file_name);

    pair<bool, bool> data_from_string(const string& our_fusion, const bool return_on_failure);

    void do_write_input_file(InputMap<mpq_class>&  input) const;

    template <typename Integer>
    void make_gradings(const vector<Integer>& d);
    vector<vector<dynamic_bitset> > make_part_classes(const vector<vector<dynamic_bitset> >& GradPartitions);
    bool compatible_duality(const vector<dynamic_bitset >& parts);
    void make_grad_mult_table();
    void restrict_type_automs_to_grading();

    void make_type_automs();
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
    bool use_modular_grading;

    bool check_simplicity;
    bool select_simple;
    bool candidate_given;
    bool automorphisms_mde;
    bool type_automs_made;

    bool use_automorphisms;
    bool write_mult_tables;

    size_t nr_coordinates;

    size_t fusion_rank;
    vector<key_t> fusion_type; // only coincidence pattern
    string fusion_type_string;
    vector<key_t> duality;

    vector<key_t> fusion_image_type;
    vector<key_t> fusion_image_duality;
    string fusion_image_type_string;
    vector<Integer> fusion_image_ring;
    Matrix<Integer> fusion_ring_map;
    bool fusion_image_commutative;

    vector<dynamic_bitset> chosen_modular_grading;
    vector<vector<int> > GradMultTable;
    set<vector<key_t> > ZeroCoords; // made 0 by grading

    double total_FPdim;

    void initialize();
    void import_global_data();
    vector<vector<vector<key_t> > > all_critical_coords_keys;
    vector<vector<key_t> > coords_to_check_key;
    vector<dynamic_bitset> coords_to_check_ind;
    vector<vector<key_t> > all_ind_tuples;
    vector<vector<key_t> > selected_ind_tuples; // the lex smallest in each FrobRec set
    map<set<vector<key_t> >, key_t> CoordMap;

    vector<vector<shortkey_t> > Automorphisms; // permutations of the coordinates
    vector<vector<shortkey_t> > type_automs; // permutations of the basis vectors
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
    // Matrix<Integer> do_iso_classes_inner(const Matrix<Integer>& LattPoints);
    Matrix<Integer> do_select_simple(const Matrix<Integer>& LattPoints) const;
    // Matrix<Integer> do_iso_classes(const Matrix<Integer>& LattPoints) const;
    vector<Integer> normal_form_of(const vector<Integer>& solution) const;
    bool simplicity_check(const vector<key_t>& subring, const vector<Integer>& sol);
    bool simplicity_check(const vector<vector<key_t> >& subrings, const vector<Integer>& sol);
    bool automs_compatible(const vector<key_t>& cand) const;

    // for automosphisms
    void make_automorphisms();
    vector<Integer> norrmal_form(const vector<Integer> lattice_point);

    Matrix<Integer> make_linear_constraints(const vector<Integer>& d);
    vector<Integer> make_linear_equation(const map<vector<key_t>, Integer>& components, const Integer& rhs);
    Matrix<Integer> make_linear_constraints_partition(const vector<Integer>& d,
                                            const vector<long>& card);
    pair<Integer, vector<key_t> >  term(const key_t& i, const key_t& j, const key_t& k);
    set<map<vector<key_t>, Integer> > make_associativity_constraints();
    // void set_global_fusion_data();

    Matrix<Integer> make_add_constraints_for_grading();

    void write_all_data_tables(const Matrix<Integer>& rings, ostream& table_out);
    void tables_for_all_rings(const Matrix<Integer>& rings);
    vector<Matrix<Integer> > make_all_data_tables(const vector<Integer>& ring);
    Matrix<Integer> data_table(const vector<Integer>& ring, const size_t i);

    Matrix<Integer> make_homomorphism_constraints();
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
void make_full_input(const FusionBasic& FusionInput, InputMap<Integer>& input_data);
template <typename Integer>
void make_full_input_partition(InputMap<Integer>& input_data);

void make_input_from_fusion_data(const FusionBasic& FusionInput, InputMap<mpq_class>&  input, const bool write_input_file);
void make_partition_input_from_fusion_data(const FusionBasic& FusionInput, InputMap<mpq_class>&  input, const bool write_input_file);

vector<dynamic_bitset> make_all_subsets(const size_t card);
vector<vector<shortkey_t> > make_all_permutations(size_t n);
vector<vector<shortkey_t> > collect_coincidence_subset_keys(const vector<key_t>& type);
template <typename Integer>
vector<vector<shortkey_t> > make_all_permutations(const vector<key_t>& v, const vector<key_t>& duality,
                                                  const  Matrix<Integer>& fusion_ring_map);
template <typename Integer>
vector<vector<shortkey_t> > make_all_permutations(const vector<key_t>& type, const vector<key_t>& duality,
                                                  const  Matrix<Integer>& fusion_ring_map);

template <typename Integer>
void write_vec_vec_Mat(vector<vector<Matrix<Integer> > > AllTables, ostream& table_out);

// void remove_global_fusion_data();

// void post_process_fusion(const vector<string>& command_line_items);

template <typename Integer>
void string_to_type(vector<Integer>& our_type, const string& our_type_string){

    istringstream type_stram(our_type_string);
    for(size_t i = 0; i < our_type.size(); ++i){
        type_stram >> our_type[i];
    }
}

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
        double this_FPdim;
        /* if(using_mpq_class<Integer>())
            this_FPdim = mpq_to_nmz_float(full_type[i]);
        else */
            this_FPdim = convertTo_nmz_float<Integer>(full_type[i]);
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
                use_modular_grading = true;
            prel_duality[0] = 0;
        }
        if(prel_duality[0] == -2) {
            commutative = true;
            use_modular_grading = true;
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

    bool has_fusion_image = false;
    if(contains(input_data, Type::fusion_image_ring)) {
        if(!contains(input_data, Type::fusion_image_type)
            || !contains(input_data, Type::fusion_ring_map) )
            throw BadInputException("Incomplete fusion image data");
        has_fusion_image = true;
    }

    if(!has_fusion_image)
        return;

    convert_vector_via_string(fusion_image_ring,input_data[Type::fusion_image_ring][0]);
    convert_matrix_via_string(fusion_ring_map, input_data[Type::fusion_ring_map]);

    if(!contains(input_data, Type::fusion_image_duality)){ // take the default
        input_data[Type::fusion_image_duality].resize(1);
        convert_vector_via_string(input_data[Type::fusion_image_duality][0],
                                  identity_key(fusion_ring_map.nr_of_columns()));
    }

    InputMap<Integer> Help;
    Help[Type::fusion_type]= input_data[Type::fusion_image_type][0];
    Help[Type::fusion_duality] = input_data[Type::fusion_image_duality][0];
    FusionBasic HB;
    HB.read_data_from_input(Help);
    fusion_image_type = HB.fusion_type;
    fusion_image_type_string = HB.fusion_type_string;
    fusion_image_duality = HB.duality;
    fusion_image_commutative = HB.commutative;

    if(fusion_image_type.size() != fusion_ring_map.nr_of_columns()
        || fusion_image_duality.size() != fusion_ring_map.nr_of_columns()
        || fusion_type.size() != fusion_ring_map.nr_of_rows())
        throw BadInputException("Formats of image data don't fit");

    for(size_t i = 0; i < fusion_ring_map.nr_of_rows(); ++i){
        for(size_t j =  0; j < fusion_ring_map.nr_of_columns(); ++j){
            if(fusion_ring_map[i][fusion_image_duality[j]] != fusion_ring_map[duality[i]][j])
                throw BadInputException("Fusion ring map not compatible with dualities:" + to_string(i) + "* = "
                        +to_string(duality[i]) + " not allowed");
        }
    }

    vector<long> full_image_type(fusion_image_duality.size());
    string our_type_string = fusion_image_type_string;
    string_to_type(full_image_type, our_type_string);

    vector<long> full_type_long(duality.size());
    convert_vector_via_string(full_type_long, full_type);

    vector<long> test_type = fusion_ring_map.MxV(full_image_type);
    if(test_type != full_type_long)
        throw BadInputException("Fusion type does not fit fusion ring map");
}

template <typename Integer>
void make_full_input(FusionBasic& FusionInput, InputMap<Integer>& input_data) {

    FusionInput.read_data_from_input(input_data);
    FusionComp<Integer> OurFusion(FusionInput);
    vector<Integer> full_type = input_data[Type::fusion_type][0];
    Matrix<Integer> Equ = OurFusion.make_linear_constraints(full_type);
    // FusionInput.type_aiutoms_mde = OurFusion.automorphisms_mde;
    // swap(FusionInput.Automorphisms,OurFusion.Automorphisms);
    // swap(FusionInput.type_automs, OurFusion.type_automs);
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

}

template <typename Integer>
inline bool is_divisible(Integer& part_FPdim,const Integer& full_FPdim, const size_t&group_order){
	long group_long = group_order; // to avoid problems with MS Windows gmpxx

    if(full_FPdim % group_long != 0)
        return false;
    part_FPdim = full_FPdim / group_long;
    return true;
}

#ifdef ENFNORMALIZ
template <>
inline bool is_divisible(renf_elem_class& part_FPdim,const renf_elem_class& full_FPdim, const size_t&group_order){

    part_FPdim = full_FPdim / group_order;
    return true;
}
#endif

template <typename Integer>
vector<vector<dynamic_bitset> > make_FPdim_partitions(const vector<Integer>& d, const Integer& part_FPdim, const size_t& group_order,
                                                      vector<dynamic_bitset>& AllSubsets);

template <typename Integer>
void FusionBasic::make_gradings(const vector<Integer>& d){

    group_order = 0;
    for(auto& t: d){
        if(t == 1)
            group_order++;
        if(t > 1)
            break;
    }
    for(size_t i = group_order; i < d.size(); ++i){
        if(d[i] == 1)
            throw BadInputException("Fusion type has 1 at wrong place");
    }
    if(group_order > 4)
        throw BadInputException("Group order > 4 not allowed for modular gradings");
    if(group_order == 1)
        throw BadInputException("Modular grading asked for perfect fusion rings");
    if(group_order == 2)
        group_type = "C2";
    if(group_order == 3){
        group_type = "C3";
        if(duality[1] != 2)
            throw BadInputException("Group " + group_type + " has wrong duality");
    }
    if(group_order == 4){
        if(duality[1] == 1 && duality[2] == 2)
            group_type = "C2xC2";
        else
            group_type = "C4";
    }
    if(verbose){
        verboseOutput() << "Modular grading group is " << group_type << endl;

    }

    Integer full_FPdim = 0;
    for(auto& t: d)
        full_FPdim += t*t;
    Integer part_FPdim;
    if(!is_divisible(part_FPdim,full_FPdim, group_order))
        throw BadInputException("Fusion type cannot be partitioned");

    vector<dynamic_bitset> AllSubsets = make_all_subsets(fusion_rank);
    vector<vector<dynamic_bitset> > FPdimParts = make_FPdim_partitions(d, part_FPdim, group_order, AllSubsets);
    vector< vector<dynamic_bitset> > GradPartitions;
    for(auto& P: FPdimParts){
        if(compatible_duality(P)){
            GradPartitions.push_back(P);
            continue;
        }
    }
    // identify partitions that are conjugate under automorphisms
    make_type_automs();
    ModularGradings = make_part_classes(GradPartitions);
    if(verbose){
        verboseOutput() << ModularGradings.size() << " grading partitions found:" << endl;
        size_t i = 0;
        for(auto& P: ModularGradings){
            i++;
            verboseOutput() << "Grading " << i << endl;
            for(auto& p: P)
                verboseOutput() << bitset_to_key(p);
        }
    }
}

}  // end namespace libnormaliz



#endif /* FUSION_H_ */
