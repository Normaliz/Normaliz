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
// (i) data common to all split_levels
// (ii) data for "this" split determined this_split_index
class SplitData{

public:

    vector<long> split_levels;
    size_t nr_split_levels;
    vector<long> split_moduli;

    long this_refinement;
    vector<vector<long> > refinement_residues;  // for each split to do the history of residues leading to it
    vector<long> refinement_predecessors; // for each split to do the last predecessor

    string project;

    long max_nr_splits_per_round;
    long nr_splits_to_do;

    long this_split_index; // for the index coming from -X=<...>
    vector<long> this_split_residues; // ditto

    long this_split_predecessor; // this_split_index of the "mother" split, needed to identify the
                      // the lat fdile produced by the mother
    vector<long> this_pred_min_return; // the min return level of the predecessors;
    vector< vector<long> > this_pred_done_indices; // indices of done elelemts of LatPopints on prwedecrssor levels
                                                   // to be read from predecessor lat file
    long this_round;
    long nr_rounds;

    SplitData();
    SplitData(const string& this_project, const long& level, const size_t& nr_vectors);

    void read_data(const string& this_project);
    void set_this_split(const long& given_split);
    // void next_round() const;
    void set_default(const string& this_project);
    void write_data() const;
    long necessary_rounds() const;
};

void collect_lat(const string& project);

void next_round(const string& project);

}  // namespace libnormaliz

#endif  // NMZ_CHUNK_H
