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

#include <stdlib.h>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <algorithm>
using namespace std;

#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/cone.h"
using namespace libnormaliz;

#include "normaliz.h"
#include "options.h"
#include "output.h"


OptionsHandler::OptionsHandler() {
    project_name_set = false;
    output_dir_set=false;
    write_extra_files = false, write_all_files = false;
    // use_Big_Integer = false;
    use_long_long = false;
    ignoreInFileOpt = false;
    nr_threads = 0;
    no_ext_rays_output=false;
}


bool OptionsHandler::handle_commandline(int argc, char* argv[]) {
        vector<string> LongOptions;
        string ShortOptions; //all options concatenated (including -)
        // read command line options
        for (int i = 1; i < argc; i++) {
                if (argv[i][0] == '-') {
                        if (argv[i][1] != '\0') {
                                if (argv[i][1] != 'x') {
                                        if (argv[i][1] == '-') {
                                                string LO = argv[i];
                                                LO.erase(0, 2);
                                                LongOptions.push_back(LO);
                                        } else
                                                ShortOptions = ShortOptions + argv[i];
                                } else if (argv[i][2] == '=') {
                        #ifdef _OPENMP
                        string Threads = argv[i];
                        Threads.erase(0,3);
                        if ( (istringstream(Threads) >> nr_threads) && nr_threads >= 0) {
                            set_thread_limit(nr_threads);
                            // omp_set_num_threads(nr_threads); -- now in cone.cpp
                        } else {
                            cerr<<"Error: Invalid option string "<<argv[i]<<endl;
                        exit(1);
                        }
                    #else
                                        cerr << "Warning: Compiled without OpenMP support, option "
                                                        << argv[i] << " ignored." << endl;
                                        #endif
                                } else {
                                        cerr << "Error: Invalid option string " << argv[i] << endl;
                                        exit(1);
                                }
                        }
                } else {
                    setProjectName(argv[i]);
                }
        }
        return handle_options(LongOptions, ShortOptions);
}

void OptionsHandler::setProjectName(const string& s) {
    if (project_name_set) {
        cerr << "Error: Second project name " << s << " in command line!" << endl;
        exit(1);
    }
    project_name = s;
    // check if we can read the .in file
    string name_in= project_name+".in";
    const char* file_in=name_in.c_str();
    ifstream in2;
    in2.open(file_in,ifstream::in);
    if (in2.is_open()==false) {
        //check if user added ".in" and ignore it in this case
        string suffix (".in");
        size_t found = project_name.rfind(suffix);
        if (found!=string::npos) {
            project_name.erase(found);
        }
    } else {
        in2.close();
    }
    project_name_set = true;
}

void OptionsHandler::setOutputDirName(const string& s) {
    output_dir=s;
    char slash='/';
    #ifdef _WIN32 //for 32 and 64 bit windows
        slash='\\';
    #endif
    if(output_dir[output_dir.size()-1]!=slash)
        output_dir+=slash;
    output_dir_set=true;
}

bool OptionsHandler::handle_options(vector<string>& LongOptions, string& ShortOptions) {
    //Analyzing short command line options
    for (size_t i = 1; i <ShortOptions.size(); i++) {
        switch (ShortOptions[i]) {
            case '-':
                break;
            case 'c':
                verbose=true;
                break;
            case 'f':
                write_extra_files = true;
                break;
            case 'a':
                write_all_files = true;
                break;
            case 'T':
                to_compute.set(ConeProperty::Triangulation);
                // to_compute.set(ConeProperty::Multiplicity);
                break;
            case 's':
                to_compute.set(ConeProperty::SupportHyperplanes);
                break;
            case 'S':
                to_compute.set(ConeProperty::Sublattice);
                break;
            case 't':
                to_compute.set(ConeProperty::TriangulationSize);
                break;
            case 'v':
                to_compute.set(ConeProperty::Multiplicity);
                break;
            case 'V':
                to_compute.set(ConeProperty::Volume);
                break;
            case 'n':
                to_compute.set(ConeProperty::HilbertBasis);
                to_compute.set(ConeProperty::Multiplicity);
                break;
            case 'N':
                to_compute.set(ConeProperty::HilbertBasis);
                break;
            case 'w':
                to_compute.set(ConeProperty::IsIntegrallyClosed);
                break;
            case '1':
                to_compute.set(ConeProperty::Deg1Elements);
                break;
            case 'q':
                to_compute.set(ConeProperty::HilbertSeries);
                break;
            case 'p':
                to_compute.set(ConeProperty::HilbertSeries);
                to_compute.set(ConeProperty::Deg1Elements);
                break;
            case 'h':
                to_compute.set(ConeProperty::HilbertBasis);
                to_compute.set(ConeProperty::HilbertSeries);
                break;
            case 'y':
                to_compute.set(ConeProperty::StanleyDec);
                break;
            case 'd':
                to_compute.set(ConeProperty::DualMode);
                break;
            case 'r':
                to_compute.set(ConeProperty::Approximate);
                break;
            case 'e':  //check for arithmetic overflow
                // test_arithmetic_overflow=true;
                cerr << "WARNING: deprecated option -e is ignored." << endl;
                break;
            case 'B':  //use Big Integer
                to_compute.set(ConeProperty::BigInt); // use_Big_Integer=true;
                break;
            case 'b':  //use the bottom decomposition for the triangulation
                to_compute.set(ConeProperty::BottomDecomposition);
                break;
            case 'C':  //compute the class group
                to_compute.set(ConeProperty::ClassGroup);
                break;
            case 'k':  //keep the order of the generators in Full_Cone
                to_compute.set(ConeProperty::KeepOrder);
                break;
            case 'o':  //suppress bottom decomposition in Full_Cone
                to_compute.set(ConeProperty::NoBottomDec);
                break;
            case 'M':  // compute minimal system of generators of integral closure
                       // as a module over original monoid
                to_compute.set(ConeProperty::ModuleGeneratorsOverOriginalMonoid);
                break;
            case '?':  //print help text and exit
                return true;
                break;
            case 'x': //should be separated from other options
                cerr<<"Error: Option -x=<T> has to be separated from other options"<<endl;
                exit(1);
                break;
            case 'I': 
                to_compute.set(ConeProperty::Integral);
                break;
            case 'L': 
                to_compute.set(ConeProperty::VirtualMultiplicity);
                break;
            case 'E': 
                to_compute.set(ConeProperty::WeightedEhrhartSeries);
                break;
            case 'i':
                ignoreInFileOpt=true;
                break;
            case 'H':
                to_compute.set(ConeProperty::IntegerHull);
                break;
            case 'D':
                to_compute.set(ConeProperty::ConeDecomposition);
                break;
            case 'P':
                to_compute.set(ConeProperty::PrimalMode);
                break;
            case 'Y':
                to_compute.set(ConeProperty::Symmetrize);
                break;
            case 'X':
                to_compute.set(ConeProperty::NoSymmetrization);
                break;
            case 'G':
                to_compute.set(ConeProperty::IsGorenstein);
                break;
            case 'j':
                to_compute.set(ConeProperty::Projection);
                break;
            case 'J':
                to_compute.set(ConeProperty::ProjectionFloat);
                break;
            default:
                cerr<<"Error: Unknown option -"<<ShortOptions[i]<<endl;
                exit(1);
                break;
        }
    }

    // Remember to update also the --help text and the documentation when changing this!
    vector<string> AdmissibleOut;
    string AdmissibleOutarray[]={"gen","cst","inv","ext","ht1","esp","egn","typ","lat","msp","mod"}; // "mod" must be last
    for(size_t i=0;i<11;++i)
        AdmissibleOut.push_back(AdmissibleOutarray[i]);
    assert(AdmissibleOut.back()=="mod");

    // analyzing long options
    for(size_t i=0; i<LongOptions.size();++i){ 
        size_t j;
        for(j=0;j<LongOptions[i].size();++j){
            if(LongOptions[i][j]=='=')
                break;            
        }
        if(j<LongOptions[i].size()){
            string OptName=LongOptions[i].substr(0,j);
            string OptValue=LongOptions[i].substr(j+1,LongOptions[i].size()-1);
            if(OptName=="OutputDir"){
                setOutputDirName(OptValue);
                continue;
            }
        }
        if(LongOptions[i]=="help"){
            return true; // indicate printing of help text
        }
        if(LongOptions[i]=="verbose"){
            verbose=true;
            continue;
        }
        if(LongOptions[i]=="version"){
            printVersion();
            exit(0);
        }
        /* if(LongOptions[i]=="BigInt"){
            use_Big_Integer=true;
            continue;
        }*/
        if(LongOptions[i]=="LongLong"){
            use_long_long=true;
            continue;
        }
        if(LongOptions[i]=="NoExtRaysOutput"){
            no_ext_rays_output=true;
            continue;
        }
        if(LongOptions[i]=="ignore"){
            ignoreInFileOpt=true;
            continue;
        }
        if(LongOptions[i]=="files"){
            write_extra_files = true;
            continue;
        }
        if(LongOptions[i]=="all-files"){
            write_all_files = true;
            continue;
        }
        if(find(AdmissibleOut.begin(),AdmissibleOut.end(),LongOptions[i])!=AdmissibleOut.end()){
            OutFiles.push_back(LongOptions[i]);
            continue;
        }
        try {
            to_compute.set(toConeProperty(LongOptions[i]));
            continue;
        } catch (const BadInputException& ) {};
        cerr << "Error: Unknown option --" << LongOptions[i] << endl;
        exit(1);
    }
    
    if(output_dir_set){
        output_file=output_dir+pureName(project_name);
    }
    else
        output_file=project_name;
    
    

    return false; //no need to print help text
}

template<typename Integer>
void OptionsHandler::applyOutputOptions(Output<Integer>& Out) {
    if(no_ext_rays_output)
        Out.set_no_ext_rays_output();
    if(write_all_files) {
        Out.set_write_all_files();
    } else if (write_extra_files) {
        Out.set_write_extra_files();
    }
    if (to_compute.test(ConeProperty::Triangulation) || to_compute.test(ConeProperty::ConeDecomposition)) {
        Out.set_write_tri(true);
        Out.set_write_tgn(true);
        Out.set_write_inv(true);
    }
    if (to_compute.test(ConeProperty::StanleyDec)) {
        Out.set_write_dec(true);
        Out.set_write_tgn(true);
        Out.set_write_inv(true);
    }
    for(size_t i=0;i<OutFiles.size();++i){
        if(OutFiles[i]=="gen"){
            Out.set_write_gen(true);
            continue;
        }
        if(OutFiles[i]=="cst"){
            Out.set_write_cst(true);
            continue;
        }
        if(OutFiles[i]=="inv"){
            Out.set_write_inv(true);
            continue;
        }
        if(OutFiles[i]=="ht1"){
            Out.set_write_ht1(true);
            continue;
        }
        if(OutFiles[i]=="ext"){
            Out.set_write_ext(true);
            continue;
        }
        if(OutFiles[i]=="egn"){
            Out.set_write_egn(true);
            continue;
        }
        if(OutFiles[i]=="esp"){
            Out.set_write_esp(true);
            continue;
        }
        if(OutFiles[i]=="typ"){
            Out.set_write_typ(true);
            continue;
        }
        if(OutFiles[i]=="lat"){
            Out.set_write_lat(true);
            continue;
        }
        if(OutFiles[i]=="msp"){
            Out.set_write_msp(true);
            continue;
        }
        if(OutFiles[i]=="mod"){
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

bool OptionsHandler::activateDefaultMode() {
    if (to_compute.goals().none() && !to_compute.test(ConeProperty::DualMode) 
                                  && !to_compute.test(ConeProperty::Approximate)
                                  && !to_compute.test(ConeProperty::ProjectionFloat)
                                  && !to_compute.test(ConeProperty::Projection)  ) {
        to_compute.set(ConeProperty::DefaultMode);
        return true;
    }
    return false;
}

string pureName(const string& fullName){
// extracts the pure filename

    string slash="/";
    #ifdef _WIN32 //for 32 and 64 bit windows
        slash="\\";
    #endif
    size_t found = fullName.rfind(slash);
    if(found==std::string::npos)
        return(fullName);
    found++;
    size_t length=fullName.size()-found;
    
    // cout << "**************************** " << fullName.substr(found,length) << endl;
    // exit(1);
    return(fullName.substr(found,length));  	

}
