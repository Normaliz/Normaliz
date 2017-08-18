/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#include <vector>
#include <list>
#include <string>
#include <sstream>
using namespace std;

#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/cone.h"
using namespace libnormaliz;

//#include "Input.h"
#include "output.h"

#ifndef NMZ_OPTIONS_H
#define NMZ_OPTIONS_H

//---------------------------------------------------------------------------

class OptionsHandler {

public:
	OptionsHandler();

	// returns true if a help should be printed, false otherwise
    bool handle_commandline(int argc, char* argv[]);

    // returns true if default mode was activated, false otherwise
    bool activateDefaultMode();

    template<typename Integer>
    void applyOutputOptions(Output<Integer>& Out);

    bool isFilenameSet() const {
        return project_name_set;
    }

    bool isIgnoreInFileOpt() const {
        return ignoreInFileOpt;
    }

    int getNrThreads() const {
        return nr_threads;
    }

    void activateConeProperty(ConeProperty::Enum cp) {
        to_compute.set(cp, true);
    }

    void activateInputFileConeProperty(ConeProperty::Enum cp) {
        if (!ignoreInFileOpt) to_compute.set(cp, true);
    }
    /* void activateInputFileBigInt() {
        if (!ignoreInFileOpt) use_Big_Integer = true;
    }*/
    void activateInputFileLongLong() {
        if (!ignoreInFileOpt) use_long_long = true;
    }
    
    void activateNoExtRaysOutput() {
        if (!ignoreInFileOpt) no_ext_rays_output = true;
    }
    void activateExtRaysFloar() {
        if (!ignoreInFileOpt) ext_rays_float = true;
    }

    const ConeProperties& getToCompute() const {
        return to_compute;
    }

    /* bool isUseBigInteger() const {
        return use_Big_Integer;
    }*/
    bool isUseLongLong() const {
        return use_long_long;
    }
    
    bool isNoExtRaysOutput() const {
        return no_ext_rays_output;
    }
    bool isExtRaysFloat() const {
        return ext_rays_float;
    }

    const string& getProjectName() const {
        return project_name;
    }
    
    const string& getOutputDir() const {
        return output_dir;
    }

    void setProjectName(const string& s);
    void setOutputDirName(const string& s);

//---------------------------------------------------------------------------

private:
	bool project_name_set;
        bool output_dir_set;
	string project_name;
        string output_dir;
        string output_file;

	// bool use_Big_Integer; now in ConeProperty
	bool use_long_long;
        bool no_ext_rays_output;
        bool ext_rays_float;
        
    bool ignoreInFileOpt;

    int nr_threads;

    ConeProperties to_compute;

    bool write_extra_files, write_all_files;

    vector<string> OutFiles;

    //return true if help should be printed, false otherwise
    bool handle_options(vector<string>& LongOptions, string& ShortOptions);
};

//---------------------------------------------------------------------------

string pureName(const string& fullName); // extracts the pure filename from a path

#endif //NMZ_OPTIONS_H
