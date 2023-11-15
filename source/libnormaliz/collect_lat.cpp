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
SplitData::SplitData(const string& this_project, const long& level, const size_t& nr_vectors){

    project = this_project;
    nr_split_levels = 1;
    this_split_levels.resize(nr_split_levels);
    this_split_levels[0] = level;
    split_moduli.resize(nr_split_levels);
    split_moduli[0] = 1000;
    if(nr_vectors< 1000)
        split_moduli[0] = nr_vectors;
    max_nr_splits_per_round = 1000;
    nr_splits_to_do = split_moduli[0];
    this_round = 0;
    this_refinement = 0;
}

void SplitData::read_data(const string& this_project){

    project = this_project;
    string name = project + ".split.data";

    ifstream split_control(name.c_str());
    if(!split_control.is_open())
        throw BadInputException(name + " does not exist");

    char c;
    split_control >> ws;
    c = split_control.peek();
    bool start_split;
    if (c != 'r'){  // start split
        start_split = true;
    }
    else{
        start_split = false;
        string s;
        split_control >> s;
        if(s != "refinement")
            throw BadInputException("Split data for refinement corrupt");
    }

    split_control >> nr_split_levels;
    this_split_levels.resize(nr_split_levels);
    split_moduli.resize(nr_split_levels);
    this_split_residues.resize(nr_split_levels);
    split_control >> max_nr_splits_per_round >> nr_splits_to_do;
    split_control >> this_round >> this_refinement;

    if(start_split){
        long nr_splitPatches_all_rounds = 1;

        for(long i = 0; i < nr_split_levels; ++i){
                split_control >> this_split_levels[i] >> split_moduli[i];
                nr_splitPatches_all_rounds *= split_moduli[i];
        }

        if(split_control.fail())
            throw BadInputException(name + " corrupt");

        if( nr_splits_to_do != nr_splitPatches_all_rounds)
            throw BadInputException("Numbers of splits for all rounds does not fit");
    }
    else{
        refinement_residues.resize(nr_splits_to_do);
        refinement_levels.resize(nr_splits_to_do);
        refinement_done_indices.resize(nr_splits_to_do);
        for(long i = 0; i < nr_split_levels; ++i){
            split_control >> split_moduli[i];
        }
        for(size_t i= 0; i < nr_splits_to_do; ++i){
            refinement_levels[i].resize(nr_split_levels);
            refinement_residues[i].resize(nr_split_levels);
            for(size_t j = 0; j < nr_split_levels; ++j){
                split_control >> refinement_levels[i][j];
            }
            for(size_t j = 0; j < nr_split_levels; ++j){
                split_control >> refinement_residues[i][j];
            }
            refinement_done_indices[i].resize(this_refinement);
            for(size_t j = 0; j < this_refinement; ++j)
                split_control >> refinement_done_indices[i][j];
        }
        // Matrix<long>(refinement_residues).debug_print();
    }

    assert(!split_control.fail());

    split_control.close();
}

void SplitData::write_data() const{

    string name = project + ".split.data";
    ofstream new_split_control(name.c_str());
    if(this_refinement > 0){
        new_split_control << "refinement" << endl;
    }
    new_split_control << nr_split_levels << endl;
    new_split_control << max_nr_splits_per_round << " " << nr_splits_to_do << endl;
    new_split_control << this_round << " " << this_refinement << endl;
    if(this_refinement == 0){
        for(size_t i = 0; i < this_split_levels.size(); ++i)
            new_split_control << " " << this_split_levels[i] << " " << split_moduli[i];
        new_split_control << endl;
    }
    else{
            new_split_control << split_moduli;
    }

    if(this_refinement == 0){
        new_split_control.close();
        return;
    }

    new_split_control << endl;

    for(size_t i = 0; i < nr_splits_to_do; ++i){
        for(size_t j = 0; j < nr_split_levels; ++j)
            new_split_control << refinement_levels[i][j] << " ";
        for(size_t j = 0; j < nr_split_levels; ++j)
            new_split_control << refinement_residues[i][j] << " ";
       for(size_t j = 0; j < this_refinement; ++j)
           new_split_control << refinement_done_indices[i][j] << " ";
        new_split_control << endl;
    }

    new_split_control.close();
}

long SplitData::necessary_rounds() const{

    long requird_rounds = nr_splits_to_do / max_nr_splits_per_round;
    if(nr_splits_to_do % max_nr_splits_per_round != 0)
        requird_rounds++;
    long rounds_done = this_round +1;
    return requird_rounds - rounds_done;
}


//selects the split data for the given undex from the full data set
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
    else{
        this_split_residues = refinement_residues[this_split_index];
        this_split_levels = refinement_levels[this_split_index];
        this_split_done_indices = refinement_done_indices[this_split_index];
    }
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

void read_lat_file(ifstream& lat_in, const string& lat_name, size_t& min_return, size_t& done_indices,Matrix<long long>& TotalLat){

    string s1, s2;
    lat_in >> s1 >> s2;
    if(s1 != "preliminary_stage" || s2 != "min_return"){
        throw BadInputException("CollectLat failed because of corrupt file " + lat_name);
    }

    lat_in >> min_return;

    // now we simplify the preliminary lat file to make the work in project_and_lift easier
    Matrix<long long> solutions_so_far;

    while(true){
        lat_in >> ws;
        int c = lat_in.peek();
        if (c == EOF) {
            break;
        }
        string s1;
        lat_in >> s1;
        if(s1 != "done_indices")
            throw BadInputException(lat_name + " corrupt.");
        long prel_done_indices;
        lat_in >> prel_done_indices;
        lat_in >> s1;
        if(s1 != "found_solutions"){
            // cout << "Vergiss es" << endl;
            throw BadInputException(lat_name + " corrupt.");
        }
        size_t nr_rows, nr_cols;
        lat_in >> nr_rows >>  nr_cols;
        Matrix<long long> prel_solutions(nr_rows, nr_cols);
        for(size_t i = 0; i < nr_rows; ++i){
            for(size_t j = 0; j < nr_cols; ++j){
                lat_in >> prel_solutions[i][j];
            }
        }
        done_indices = prel_done_indices;
        solutions_so_far = prel_solutions;
    }

    lat_in.close();

    if(solutions_so_far.nr_of_rows() > 0){
        if(TotalLat.nr_of_rows() == 0){
            TotalLat.resize(0, solutions_so_far.nr_of_columns());
        }
        TotalLat.append(solutions_so_far);
        if(verbose)
            verboseOutput() << solutions_so_far.nr_of_rows() << " solutions_transferred" << endl;
    }
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

    if(verbose){
        verboseOutput() << "Collecting lattice points and preliminary data from " << our_split.nr_splits_to_do << " lat files" << endl;
    }

    Matrix<long long> TotalLat; // collects the solutions found so far

    vector<size_t> NotDone; // register the lat files of preliminary stage
    vector<size_t> MinReturnNotDone; // collects their min returns
    vector<size_t> DoneIndicesNotDone; // collects the indices already done

    // Now we read what has been compouted in the last refinement
    // and register the splits that have not been complete in NotDone
    // We read the data of the others and rewrite a simplified version.
    for(size_t i = 0; i < our_split.nr_splits_to_do; ++i){
        string lat_name = project + "." + to_string(our_split.this_refinement) + "." +  to_string(i) + ".lat";
        if(verbose)
            verboseOutput()  << lat_name << endl;

        ifstream lat_in(lat_name);
        if (!lat_in.is_open()){
            if(verbose)
                verboseOutput() << lat_name << " does not exist" << endl;
            throw BadInputException("Not all lat files computed. CollectLat not yet possible!");
        }

        char c;
        lat_in >> ws;
        c = lat_in.peek();
        if (c == 'p'){
            if(verbose)
                verboseOutput() << lat_name << " in preliminary stage" << endl;
            NotDone.push_back(i);
            size_t min_return;
            size_t done_indices;
            read_lat_file(lat_in, lat_name, min_return, done_indices, TotalLat);
            MinReturnNotDone.push_back(min_return);
            DoneIndicesNotDone.push_back(done_indices);

            continue; // next lat file
        } // done with preliminray stage

        // now the completed lat files
        lat_in.close();
        Matrix<long long> this_lat = readMatrix<long long>(lat_name);
        if(this_lat.nr_of_rows() == 0)
            continue;
        if(TotalLat.nr_of_rows() == 0)
            TotalLat.resize(0, this_lat.nr_of_columns());
        TotalLat.append(this_lat);

    } // loop over all lat files

    name = global_project + "." + to_string(our_split.this_refinement) + ".so_far.lat";
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

    // The computation is not complete and we must create the new split data

    if(verbose)
        verboseOutput() << "Computation of " + global_project + " NOT complete" << endl;
    SplitData new_split_data = our_split;
    size_t nr_sub_splits = 2;
    if(our_split.max_nr_splits_per_round/NotDone.size() > 2)
        nr_sub_splits = our_split.max_nr_splits_per_round/NotDone.size();
    new_split_data.nr_splits_to_do = NotDone.size() * nr_sub_splits;
    new_split_data.nr_split_levels++;
    new_split_data.split_moduli.push_back(nr_sub_splits);
    new_split_data.this_refinement++;
    new_split_data.this_round = 0;

    if(verbose)
        verboseOutput() << "Scheduling refinement" << endl;
    if(verbose)
        verboseOutput() << nr_sub_splits << " subplits" << endl;

    new_split_data.refinement_residues.clear();
    new_split_data.refinement_levels.clear();
    new_split_data.refinement_done_indices.clear();

    for(size_t spl = 0; spl < NotDone.size(); ++spl){

        vector<long> split_residues;
        vector<long> split_levels;
        vector<long> split_done_indices;

        size_t split_index = NotDone[spl];

        if(our_split.this_refinement == 0){ // no previous refinement
            long res = split_index;
            split_residues.resize(our_split.nr_split_levels);
            for(long i = 0; i < our_split.nr_split_levels; ++i){
                split_residues[i] = res % our_split.split_moduli[i];
                res /=  our_split.split_moduli[i];
            }
            split_levels = our_split.this_split_levels;
        }
        else{
            split_residues = our_split.refinement_residues[split_index];
            split_levels =  our_split.refinement_levels[split_index];
            split_done_indices = our_split.refinement_done_indices[split_index];
        }

        size_t next_split_level = split_levels.back();
        if(MinReturnNotDone[spl] > next_split_level)
            next_split_level = MinReturnNotDone[spl];

        // common to all subsplits
        split_levels.push_back(next_split_level);
        split_done_indices.push_back(DoneIndicesNotDone[spl]);

        for(long i = 0; i < nr_sub_splits; ++i){
            vector<long> extended_res = split_residues;
            extended_res.push_back(i); // the new subsplit
            new_split_data.refinement_residues.push_back(extended_res);
            new_split_data.refinement_levels.push_back(split_levels);
            new_split_data.refinement_done_indices.push_back(split_done_indices);
        }
    }
    new_split_data.write_data();
}

}  // namespace libnormaliz
