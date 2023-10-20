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

namespace libnormaliz {
using std::vector;
using std::map;
using std::list;
using std::set;
using std::ifstream;

template <typename Integer>
class FusionData    {
    template <typename, typename>
    friend class ProjectAndLift;

public:

    bool activated;
    bool type_and_duality_set;

    bool verbose;

    bool check_simplicity;
    bool select_simple; // not used at present
    bool candidate_given;

    bool use_automorphisms;
    bool select_iso_classes;

    vector<key_t> duality;
    size_t fusion_rank;
    vector<size_t> fusion_type;

    vector<key_t> subring_base_key;
    vector<vector<key_t> > all_base_keys;
    vector<vector<vector<key_t> > > all_critical_coords_keys;
    vector<vector<key_t> > coords_to_check_key;
    vector<dynamic_bitset> coords_to_check_ind;
    vector<vector<key_t> > all_ind_tuples;

    map<set<vector<key_t> >, key_t> CoordMap;

    FusionData();
    void read_data();
    void read_data_from_file();
    void data_from_roject();

    void make_CoordMap();
    // key_t dual(key_t i);
    key_t dual(const key_t i) const;
    set<vector<key_t> > FrobRec(const vector<key_t>& ind_tuple);
    key_t coord(set<vector<key_t> >& FR);
    key_t coord(vector<key_t>& ind_tuple);
    void make_all_ind_tuples();
    dynamic_bitset critical_coords(const vector<key_t>& base_key);
    void make_all_base_keys();
    void prepare_simplicity_check();
    void do_select_and_write_simple(const Matrix<Integer>& LattPoints);
    bool simplicity_check(const vector<key_t>& subring, const vector<Integer>& sol);
    bool simplicity_check(const vector<vector<key_t> >& subrings, const vector<Integer>& sol);

    void set_options(const ConeProperties& ToCompute, const bool verb);

};

// helper
Matrix<long long> extract_latt_points_from_out(ifstream& in_out);
template <typename Integer>
void select_and_write_simple(const Matrix<Integer>& LattPoints, const ConeProperties& ToCompute, const bool verb);
//void select_simple_fusion_rings();

}  // end namespace libnormaliz

#endif /* FUSION_H_ */
