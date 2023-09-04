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

SplitData::SplitData(){
}

// c instructor with default values
SplitData::SplitData(const string& this_project, const long& level){

    project = this_project;
    nr_split_levels = 1;
    split_levels.resize(nr_split_levels);
    split_levels[0] = level;
    split_moduli.resize(nr_split_levels);
    split_moduli[0] = 1000;
    max_nr_splits_per_round = 1000;
    nr_splits_to_do = 1000;
    this_round = 0;
    this_refinement = 0;
}

void SplitData::read_data(const string& this_project){

    project = this_project;
    string name = project + ".split.data";

    ifstream split_control(name.c_str());
    if(!split_control.is_open())
        throw BadInputException(name + " does not exist");

    long nr_splitPatches_all_rounds = 1;
    split_control >> nr_split_levels;
    split_levels.resize(nr_split_levels);
    split_moduli.resize(nr_split_levels);
    this_split_residues.resize(nr_split_levels);
    for(long i = 0; i < nr_split_levels; ++i){
            split_control >> split_levels[i] >> split_moduli[i];
            nr_splitPatches_all_rounds *= split_moduli[i];
    }
    split_control >> max_nr_splits_per_round >> nr_splits_to_do;
    split_control >> this_round >> this_refinement;

    if(split_control.fail())
        throw BadInputException(name + " corrupt");

    if(this_refinement == 0){
        if( nr_splits_to_do != nr_splitPatches_all_rounds)
            throw BadInputException("Numbers of splits for all rounds does not fit");
    }
    else{
        refinement_residues.resize(nr_splits_to_do);
        for(size_t i= 0; i < nr_splits_to_do; ++i){
            refinement_residues[i].resize(nr_split_levels);
            for(size_t j = 0; j < nr_split_levels; ++j)
                split_control >> refinement_residues[i][j];
        }
        // Matrix<long>(refinement_residues).debug_print();
    }

    assert(split_control.good());

    split_control.close();
}

void SplitData::write_data() const{

    string name = project + ".split.data";
    ofstream new_split_control(name.c_str());
    new_split_control << split_levels.size();
    for(size_t i = 0; i < split_levels.size(); ++i)
    new_split_control << " " << split_levels[i] << " " << split_moduli[i];
    new_split_control << endl;
    new_split_control << max_nr_splits_per_round << " " << nr_splits_to_do << endl;
    new_split_control << this_round << " " << this_refinement << endl;

    if(this_refinement == 0){
        new_split_control.close();
        return;
    }

    Matrix<long>(refinement_residues).pretty_print(new_split_control);
    new_split_control.close();
}

long SplitData::necessary_rounds() const{

    long requird_rounds = nr_splits_to_do / max_nr_splits_per_round;
    if(nr_splits_to_do % max_nr_splits_per_round != 0)
        requird_rounds++;
    long rounds_done = this_round +1;
    return requird_rounds - rounds_done;
}

/*
void SplitData::next_round() const{

    if(necessary_rounds() <= 0)
        throw BadInputException("All rounds done. No next round!");

    SplitData def_data = *this;
    def_data.this_round++;
    def_data.write_data();
}
*/

void SplitData::set_this_split(const long& given_index){
    if(given_index >= max_nr_splits_per_round)
        throw NoComputationException("Split index given by -X too large");
    this_split_index = given_index + this_round * max_nr_splits_per_round;
    split_index_rounds= this_split_index; // needed in cone for output
    if(this_split_index >= nr_splits_to_do)
        throw NoComputationException("Total split index too large");

    long res = this_split_index;
    if(this_refinement == 0){
        for(long i = 0; i < nr_split_levels; ++i){
            this_split_residues[i] = res % split_moduli[i];
            res /=  split_moduli[i];
        }
    }
    else
        this_split_residues = refinement_residues[this_split_index];
}

void next_round(const string& project) {

    SplitData our_split;
    our_split.read_data(project);
    if(our_split.necessary_rounds() <= 0)
        throw BadInputException("All rounds done. No next round!");
    // archive <project>.spli.data that was jiust read
    string command = "cp " + project + ".split.data " + project + "." + to_string(our_split.this_refinement) + "." + to_string(our_split.this_round) + ".split.data";
    int dummy = system(command.c_str());
    if(dummy > 0)
        throw NoComputationException("Problem in archiving <project.split.data");
    // now mwrite updated split data
    our_split.this_round++;
    our_split.write_data();
    if(verbose)
        verboseOutput() << "New round " << our_split.this_round << endl;
}

void collect_lat(const string& project) {

    string name;

    SplitData our_split;
    our_split.read_data(project);
    if(our_split.necessary_rounds() > 0)
        throw BadInputException("Last round not complete. Use --NextRound??");
    // archive <project>.split.data that was jiust read
    string command = "cp " + project + ".split.data " + project + "." + to_string(our_split.this_refinement) + "." + to_string(our_split.this_round) + ".split.data";
    int dummy = system(command.c_str());
    if(dummy > 0)
        throw NoComputationException("Problem in archiving project.split.data");

    if(verbose)
        verboseOutput() << "Collecting lattice points from " << our_split.nr_splits_to_do << " lat files" << endl;

    Matrix<long long> TotalLat;
    bool first = true;

    // First we must vread what has been compouted in previous refinements
    name = project + ".so_far.lat";
    ifstream lat_in(name.c_str());
    if(lat_in.is_open()){
        TotalLat = readMatrix<long long>(name);
        if(TotalLat.nr_of_rows() > 0){
                first = false;
        }
        lat_in.close();
    }

    vector<size_t> NotDone;


    // Now we read what has been compouted in the last refinement
    // and register the parts that have not been complete in NotDone
    for(size_t i = 0; i < our_split.nr_splits_to_do; ++i){
        name = project + "." + to_string(our_split.this_refinement) + "." +  to_string(i) + ".lat";
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

        // delete done file if lat exists
        name = global_project + "."  + v_to_point_list(our_split.this_split_residues) + "done";
        file_in = name.c_str();
        test_in.open(file_in, ifstream::in);
        if (test_in.is_open()){
            test_in.close();
            string command = "rm "+ name + " &> /dev/null";
            int dummy = system(command.c_str());
            if(dummy == 0){
                if(verbose)
                    verboseOutput() << name << " removed" << endl;
            }
            else{
                if(verbose)
                    verboseOutput() << name << " NOT removed" << endl;
            }
        }

        name = project + "." + to_string(our_split.this_refinement) + "." +  to_string(i) + ".lat";
        Matrix<long long> this_lat = readMatrix<long long>(name);
        if(this_lat.nr_of_rows() == 0)
            continue;
        if(first)
            TotalLat.resize(0, this_lat.nr_of_columns());
        first = false;
        TotalLat.append(this_lat);
    }

    name = global_project + ".so_far.lat";
    ofstream lat_out(name.c_str());
    TotalLat.print(lat_out);

    if(NotDone.size() == 0){
        string comp_name = global_project + ".final.lat";
        ofstream comp_out(comp_name.c_str());
        TotalLat.print(comp_out);
        if(verbose)
            verboseOutput() << "Computation of " + global_project + " complete " << endl;
        return;
    }

    if(verbose)
        verboseOutput() << "Computation of " + global_project + " NOT complete" << endl;
    SplitData new_split_data = our_split;
    size_t nr_sub_splits = 2;
    if(our_split.max_nr_splits_per_round/NotDone.size() > 2)
        nr_sub_splits = our_split.max_nr_splits_per_round/NotDone.size();
    new_split_data.nr_splits_to_do = NotDone.size() * nr_sub_splits;
    new_split_data.split_levels.push_back(our_split.split_levels.back() + 1);
    new_split_data.nr_split_levels++;
    new_split_data.split_moduli.push_back(nr_sub_splits);
    new_split_data.this_refinement++;
    new_split_data.this_round = 0;

    if(verbose)
        verboseOutput() << "Scheduling refinement" << endl;
    if(verbose)
        verboseOutput() << nr_sub_splits << " subplits" << endl;

    vector<long> split_residues(our_split.nr_split_levels);
    new_split_data.refinement_residues.clear();
    for(auto& split_index:NotDone){
        if(our_split.this_refinement == 0){ // no previous refinement
            long res = split_index;
            for(long i = 0; i < our_split.nr_split_levels; ++i){
                split_residues[i] = res % our_split.split_moduli[i];
                res /=  our_split.split_moduli[i];
            }
        }
        else{
            split_residues = our_split.refinement_residues[split_index];

        }
        vector<long> extended_res = split_residues;
        extended_res.resize(extended_res.size() +1);
        for(long i = 0; i < nr_sub_splits; ++i){
            extended_res.back() = i;
            new_split_data.refinement_residues.push_back(extended_res);
        }
    }
    new_split_data.write_data();
}

}  // namespace libnormaliz
