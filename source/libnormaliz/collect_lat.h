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

#ifndef NORMALIZ_COLLECTLAT_H
#define NORMALIZ_COLLECTLAT_H

#include "libnormaliz/general.h"

namespace libnormaliz {
using std::vector;
using std::string;


// This class contains two types of data:
// (i) data common to all
// (ii) data for "this" split determined this_split_index
class SplitData{

public:

    string project;

    size_t nr_split_levels;  // common to all splits
    vector<long> split_moduli; // they are common to all splits

    long this_refinement;
    vector<vector<long> > refinement_residues;  // for each split to record the history of residues leading to it
    vector<vector<long> > refinement_levels;  // ditto
    vector<vector<long> > refinement_total_indices;  // ditto
    vector<vector<long> > refinement_done_indices;  // ditto
    vector<vector<long> > refinement_predecessors;  // ditto

    // long max_nr_splits_per_round;
    long nr_splits_to_do;
    // long nr_rounds;

    long this_split_index; // for the index coming from -X=<...>
    vector<long> this_split_residues; // data depending on this_split_index
    vector<long> this_split_levels;  // TODO single out as a struct
    vector<long> this_split_total_indices;
    vector<long> this_split_done_indices;
    vector<long> this_split_min_returns;


    // long this_round;

    SplitData();
    SplitData(const string& this_project, const long& level, const size_t& nr_vectors);

    void read_data(const string& this_project);
    void set_this_split(const long& given_split);
    // void next_round() const;
    void set_default(const string& this_project);
    void write_data() const;
    long necessary_rounds() const;
};

void collect_lat(const string& project, const long given_nr_subsplits);

void next_round(const string& project);

}  // namespace libnormaliz

#endif  // NMZ_CHUNK_H
