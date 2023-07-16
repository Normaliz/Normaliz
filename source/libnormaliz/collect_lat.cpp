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

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "libnormaliz/collect_lat.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/input.h"

namespace libnormaliz {
using std::cout;
using std::endl;
using std::ifstream;

void collect_lat() {

    bool no_refinement = false;

    string name = global_project + ".split.data";
    ifstream split_control(name.c_str());
    if(!split_control.is_open())
        throw BadInputException(name + " does not exist");
    long nr_split_patches;
    long nr_splitPatches_all_rounds = 1;
    split_control >> nr_split_patches;
    split_patches.resize(nr_split_patches);
    split_moduli.resize(nr_split_patches);
    split_residues.resize(nr_split_patches);
    for(long i = 0; i < nr_split_patches; ++i){
            split_control >> split_patches[i] >> split_moduli[i];
            nr_splitPatches_all_rounds *= split_moduli[i];
    }
    size_t max_nr_splits, actual_nr_splits;
    split_control >> max_nr_splits >> actual_nr_splits;
    if(split_control.fail())
        throw BadInputException(name + " corrupt");
    split_control.close();

    name = global_project + ".rounds.data";
    ifstream rounds_control(name.c_str());
    if(rounds_control.is_open()){
        no_refinement = true;
        long total_rounds, this_round;
        rounds_control >> this_round >> total_rounds;
        if(total_rounds != this_round +1){
            throw BadInputException("Not all rounds done");
        }
        // cout << this_round << " " << total_rounds << " " << max_nr_splits << " " << nr_splitPatches_all_rounds << endl;
        actual_nr_splits = max_nr_splits * total_rounds;
        if( actual_nr_splits != nr_splitPatches_all_rounds)
            throw BadInputException("Numbers of splits for all rounds does not fit");
    }

    if(verbose)
        verboseOutput() << "Collecting lattice points from " << actual_nr_splits << " lat files" << endl;

    Matrix<long long> TotalLat;
    bool first = true;

    name = global_project + ".total.lat";
    ifstream lat_in(name.c_str());
    if(lat_in.is_open()){
        TotalLat = readMatrix<long long>(name);
        if(TotalLat.nr_of_rows() > 0){
                first = false;
        }
        lat_in.close();
    }

    vector<size_t> NotDone;

    for(size_t i = 0; i < actual_nr_splits; ++i){
        name = global_project + "." + to_string(i) + ".lat";
        if(verbose)
            verboseOutput()  << name << endl;

        const char* file_in = name.c_str();
        ifstream test_in;
        test_in.open(file_in, ifstream::in);
        if (!test_in.is_open()){
            if(verbose)
                verboseOutput() << name << "does not exist" << endl;
            NotDone.push_back(i);
            continue;
        }

        Matrix<long long> this_lat = readMatrix<long long>(name);
        if(this_lat.nr_of_rows() == 0)
            continue;
        if(first)
            TotalLat.resize(0, this_lat.nr_of_columns());
        first = false;
        TotalLat.append(this_lat);
    }

    if(NotDone.size() >0 && no_refinement){
        throw BadInputException("Incomplete computation after rounds");
    }

    if(NotDone.size() > 0){
        if(verbose)
            verboseOutput() << "Computation NOT complete" << endl;
        if(max_nr_splits >= 2*NotDone.size()){
            if(verbose)
                verboseOutput() << "Scheduling refinement" << endl;
            size_t nr_sub_splits = max_nr_splits/NotDone.size();
            if(verbose)
                verboseOutput() << nr_sub_splits << " subplits" << endl;

            name = global_project + ".split.dist";
            ofstream dist_out(name.c_str());

            for(auto& split_res:NotDone){
                long res = split_res;
                for(long i = 0; i < nr_split_patches; ++i){
                    split_residues[i] = res % split_moduli[i];
                    res /=  split_moduli[i];
                }
                vector<long> extended_res = split_residues;
                extended_res.resize(extended_res.size() +1);
                for(long i = 0; i < nr_sub_splits; ++i){
                    extended_res.back() = i;
                    dist_out << extended_res;
                }
            }
            long new_actual_nr_splits = NotDone.size()*nr_sub_splits;
            split_patches.push_back(split_patches.back()+1);
            split_moduli.push_back(nr_sub_splits);
            name = global_project + ".split.data";
            ofstream new_split_control(name.c_str());
            new_split_control << split_patches.size();
            for(size_t i = 0; i < split_patches.size(); ++i)
                new_split_control << " " << split_patches[i] << " " << split_moduli[i];
            new_split_control << endl;
            new_split_control << max_nr_splits << " " << new_actual_nr_splits << endl;
        }
        else{
            if(verbose)
                verboseOutput() << "Refinement not possible" << endl;
        }
    }
    else{
        if(verbose)
            verboseOutput() << "Computation complete" << endl;
    }

    name = global_project + ".total.lat";
    ofstream lat_out(name.c_str());
    TotalLat.print(lat_out);
}

}  // namespace libnormaliz
