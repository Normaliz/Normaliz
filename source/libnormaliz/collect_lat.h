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

class SplitData{

public:

    vector<long> split_patches;
    size_t nr_split_patches;
    vector<long> split_moduli;

    long max_nr_splits_per_round;
    long nr_splits_to_do;

    long this_split;
    vector<long> this_split_residues;

    long this_round;
    long nr_rounds;
    long this_refinement;

    vector<vector<long> > refinement_residues;

    string project;

    SplitData();

    void read_data(const string& this_project);
    void set_this_split(const long& given_split);
    void next_round();
    void write_default(const string& this_project);
    void write_data();
};

void collect_lat();

}  // namespace libnormaliz

#endif  // NMZ_CHUNK_H
