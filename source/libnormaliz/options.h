/*
 * Normaliz
 * Copyright (C) 2007-2025  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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

#ifndef NORMALIZ_OPTIONS_H
#define NORMALIZ_OPTIONS_H

#include "libnormaliz/general.h"
#include "libnormaliz/input_type.h"
#include "libnormaliz/output.h"

#include <vector>
#include <string>

#ifndef STRINGIFY
#define STRINGIFYx(Token) #Token
#define STRINGIFY(Token) STRINGIFYx(Token)
#endif

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------

class OptionsHandler {
   public:
    OptionsHandler();

    // returns true if a help should be printed, false otherwise
    // bool handle_commandline(int argc, char* argv[]);
    bool handle_commandline(vector<string> argv);

    // returns true if default mode was activated, false otherwise
    bool activateDefaultMode();

    template <typename Integer>
    void applyOutputOptions(Output<Integer>& Out);

    inline bool isFilenameSet() const {
        return project_name_set;
    }

    inline bool isIgnoreInFileOpt() const {
        return ignoreInFileOpt;
    }

    inline int getNrThreads() const {
        return nr_threads;
    }

    inline void activateConeProperty(ConeProperty::Enum cp) {
        to_compute.set(cp, true);
    }

    inline void activateInputFileConeProperty(ConeProperty::Enum cp) {
        if (!ignoreInFileOpt)
            to_compute.set(cp, true);
    }
    /* void activateInputFileBigInt() {
        if (!ignoreInFileOpt) use_Big_Integer = true;
    }*/
    inline void activateInputFileLongLong() {
        if (!ignoreInFileOpt)
            use_long_long = true;
    }

    inline void activateNoExtRaysOutput() {
        if (!ignoreInFileOpt)
            no_ext_rays_output = true;
    }

    inline void activateBinomialsPacked() {
        if (!ignoreInFileOpt)
            binomials_packed = true;
    }

    inline void activateNoMatricesOutput() {
        if (!ignoreInFileOpt)
            no_matrices_output = true;
    }
    inline void activateOutputOnInterrupt() {
        if (!ignoreInFileOpt)
            output_on_interrupt = true;
    }


    inline void activateNoSuppHypsOutput() {
        if (!ignoreInFileOpt)
            no_supp_hyps_output = true;
    }

    inline void activateNoHilbertBasisOutput() {
        if (!ignoreInFileOpt)
            no_hilbert_basis_output = true;
    }

    inline const ConeProperties& getToCompute() const {
        return to_compute;
    }

    /* bool isUseBigInteger() const {
        return use_Big_Integer;
    }*/
    inline bool isUseLongLong() const {
        return use_long_long;
    }

    inline bool isUseChunk() const {
        return use_chunk;
    }

    inline bool isUseCollectLat() const {
        return use_collect_lat;
    }
    inline bool isUseSaveLocalSolutions() const {
        return use_save_local_solutions;
    }

    inline bool isUseMakeFullInput() const {
        return use_make_full_input;
    }

    inline bool isUseAddChunks() const {
        return use_add_chunks;
    }

    inline bool isUseSplit() const {
        return use_split;
    }

    inline bool isNoExtRaysOutput() const {
        return no_ext_rays_output;
    }

    inline bool isBinomialsPacked() const {
        return binomials_packed;
    }
    inline bool isOutputOnInterrupt() const {
        return output_on_interrupt;
    }
    inline bool isNoMatricesOutput() const {
        return no_matrices_output;
    }


    inline bool isNoSuppHypsOutput() const {
        return no_supp_hyps_output;
    }

    inline bool isNoHilbertBasisOutput() const {
        return no_hilbert_basis_output;
    }

    inline const string& getProjectName() const {
        return project_name;
    }

    inline void setProjectName_from_list(const string& name){
        project_name = name;
        project_name_set = true;
    }

    inline const string& getOutputDir() const {
        return output_dir;
    }

    inline bool get_given_name_contains_in() const {
        return given_name_contains_in;
    }

    //---------------------------------------------------------------------------

   private:
    bool project_name_set;
    bool output_dir_set;
    string project_name;
    string output_dir;
    string output_file;

    // bool use_Big_Integer; now in ConeProperty
    bool use_long_long;
    bool use_chunk;
    bool use_collect_lat;
    bool use_save_local_solutions;
    bool use_make_full_input;
    bool use_add_chunks;
    bool no_ext_rays_output;
    bool no_supp_hyps_output;
    bool no_matrices_output;
    bool no_hilbert_basis_output;
    bool binomials_packed;
    bool use_split;
    bool given_name_contains_in;

    bool ignoreInFileOpt;

    int nr_threads;

    ConeProperties to_compute;

    bool write_extra_files, write_all_files;

    vector<string> OutFiles;

    // return true if help should be printed, false otherwise
    bool handle_options(vector<string>& LongOptions, string& ShortOptions);

    void setProjectName(const string& s);
    void setOutputDirName(const string& s);
};

//---------------------------------------------------------------------------

inline OptionsHandler::OptionsHandler() {
    project_name_set = false;
    output_dir_set = false;
    write_extra_files = false, write_all_files = false;
    // use_Big_Integer = false;
    use_long_long = false;
    use_chunk = false;
    use_collect_lat = false;
    use_save_local_solutions = false;
    use_add_chunks = false;
    ignoreInFileOpt = false;
    nr_threads = 0;
    no_ext_rays_output = false;
    no_supp_hyps_output = false;
    no_matrices_output = false;
    no_hilbert_basis_output = false;
    binomials_packed = false;
    use_split = false;
    use_make_full_input = false;
    given_name_contains_in = false;
}

template <typename Integer>
void OptionsHandler::applyOutputOptions(Output<Integer>& Out) {
    if (no_ext_rays_output)
        Out.set_no_ext_rays_output();
    if (no_supp_hyps_output)
        Out.set_no_supp_hyps_output();
    if (no_matrices_output)
        Out.set_no_matrices_output();
    if (no_hilbert_basis_output)
        Out.set_no_hilbert_basis_output();
    if (binomials_packed)
        Out.set_binomials_packed();
    if (write_all_files) {
        Out.set_write_all_files();
    }
    else if (write_extra_files) {
        Out.set_write_extra_files();
    }

    if (to_compute.test(ConeProperty::WritePreComp)) {
        Out.set_write_precomp(true);
    }

    for (const auto& OutFile : OutFiles) {
        if (OutFile == "gen") {
            Out.set_write_gen(true);
            continue;
        }
        if (OutFile == "cst") {
            Out.set_write_cst(true);
            continue;
        }
        if (OutFile == "inv") {
            Out.set_write_inv(true);
            continue;
        }
        if (OutFile == "ht1") {
            Out.set_write_ht1(true);
            continue;
        }
        if (OutFile == "ext") {
            Out.set_write_ext(true);
            continue;
        }
        if (OutFile == "egn") {
            Out.set_write_egn(true);
            continue;
        }
        if (OutFile == "esp") {
            Out.set_write_esp(true);
            continue;
        }
        if (OutFile == "typ") {
            Out.set_write_typ(true);
            continue;
        }
        if (OutFile == "lat") {
            Out.set_write_lat(true);
            continue;
        }
        if (OutFile == "msp") {
            Out.set_write_msp(true);
            continue;
        }
        if (OutFile == "precomp") {
            Out.set_write_precomp(true);
            continue;
        }
        if (OutFile == "mod") {
            Out.set_write_mod(true);
            continue;
        }
    }

    if (!project_name_set) {
        cerr << "ERROR: No project name set!" << endl;
        exit(1);
    }
    Out.set_name(output_file);
}

string pureName(const string& fullName);

inline string package_string() {
    string optional_packages;

#ifdef NMZ_COCOA
    optional_packages += " CoCoALib";
#endif
#ifdef NMZ_FLINT
#ifndef ENFNORMALIZ
    optional_packages += " Flint";
#endif
#endif
#ifdef ENFNORMALIZ
    optional_packages += " Flint e-antic";
#endif
#ifdef NMZ_NAUTY
    optional_packages += " nauty";
#endif
#ifdef NMZ_HASHLIBRARY
    optional_packages += " hash-library";
#endif
    return optional_packages;
}

}  // namespace libnormaliz

#endif  // NMZ_OPTIONS_H
