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

#include "Options.h"
#include "output.h"


OptionsHandler::OptionsHandler() {
    filename_set = false;
    write_extra_files = false, write_all_files = false;
	do_bottom_dec=false;
	keep_order=false;
	use_Big_Integer = false;
	ignoreInFileOpt=false;
	nmzInt_E = false, nmzInt_I = false, nmzInt_L = false;
    nr_threads = 0;
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
                     if ( (istringstream(Threads) >> nr_threads) && nr_threads > 0) {
                         omp_set_num_threads(nr_threads);
                     } else {
                         cerr<<"Warning: Invalid option string "<<argv[i]<<endl;
                     }
                    #else
					cerr << "Error: Compiled without OpenMP support, option "
							<< argv[i] << " ignored." << endl;
					throw BadInputException();
					#endif
				} else {
					cerr << "Error: Invalid option string " << argv[i] << endl;
					throw BadInputException();
				}
			}
		} else if (!filename_set) {
			string s(argv[i]);
			output_name = s;
			filename_set = true;
		} else if (filename_set) {
			cerr << "Error: Second file name " << argv[i] << " in command line!" << endl;
			throw BadInputException();
		}
	}
	return handle_options(LongOptions, ShortOptions);
}

void OptionsHandler::setOutputName(const string& outputName) {
    output_name = outputName;
    // check if we can read the .in file
    string name_in= output_name+".in";
    const char* file_in=name_in.c_str();
    ifstream in2;
    in2.open(file_in,ifstream::in);
    if (in2.is_open()==false) {
        //check if user added ".in" and ignore it in this case
        string suffix (".in");
        size_t found = output_name.rfind(suffix);
        if (found!=string::npos) {
            output_name.erase(found);
        }
    } else {
        in2.close();
    }
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
            case 'V':
                to_compute.set(ConeProperty::Multiplicity);
                break;
            case 'n':
                to_compute.set(ConeProperty::HilbertBasis);
                to_compute.set(ConeProperty::Multiplicity);
                break;
            case 'N':
                to_compute.set(ConeProperty::HilbertBasis);
                break;
            case '1':
                to_compute.set(ConeProperty::Deg1Elements);
                break;
            case 'q':
                to_compute.set(ConeProperty::HilbertSeries);
                break;
            case 'p':
            case 'P':
                to_compute.set(ConeProperty::HilbertSeries);
                to_compute.set(ConeProperty::Deg1Elements);
                break;
            case 'h':
            case 'H':
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
                to_compute.set(ConeProperty::ApproximateRatPolytope);
                to_compute.set(ConeProperty::Deg1Elements);
                break;
            case 'e':  //check for arithmetic overflow
                // test_arithmetic_overflow=true;
                break;
            case 'B':  //use Big Integer
                use_Big_Integer=true;
                break;
			case 'b':  //use the bottom decomposition for the triangulation
				do_bottom_dec=true;
				break;
			case 'C':  //compute the class group
				to_compute.set(ConeProperty::ClassGroup);
				break;
			case 'k':  //keep the order of the generators in Full_Cone
				keep_order=true;
				break;
            case 'M':  // compute minimal system of generators of integral closure
                       // as a module over original monoid
                to_compute.set(ConeProperty::ModuleGeneratorsOfIntegralClosure);
                break;
            case '?':  //print help text and exit
                return true;
                break;
            case 'x': //should be separated from other options
                cerr<<"Error: Option -x=<T> has to be separated from other options"<<endl;
                throw BadInputException();
                break;
            case 'I':  //nmzIntegrate -I (integrate)
                nmzInt_I = true;
                to_compute.set(ConeProperty::Triangulation);
                to_compute.set(ConeProperty::Multiplicity);
                break;
            case 'L':  //nmzIntegrate -L (leading term)
                nmzInt_L = true;
                to_compute.set(ConeProperty::Triangulation);
                to_compute.set(ConeProperty::Multiplicity);
                break;
            case 'E':  //nmzIntegrate -E (Ehrhart series)
                nmzInt_E = true;
                to_compute.set(ConeProperty::StanleyDec);
                break;
            case 'i':
                ignoreInFileOpt=true;
                break;
            default:
                cerr<<"Error: Unknown option -"<<ShortOptions[i]<<endl;
                throw BadInputException();
                break;
        }
    }

    vector<string> ComputeLO;
    string ComputeLOarray[]={"SupportHyperplanes","HilbertBasis","Deg1Elements","ModuleGeneratorsOfIntegralClosure",
        "HilbertSeries","Multiplicity","ClassGroup","Triangulation","TriangulationSize","TriangulationDetSum",
        "StanleyDec","DualMode","ApproximateRatPolytope","BottomDecomposition","DefaultMode"}; // "DefaultMode" must be last
    for(size_t i=0;i<15;++i)
        ComputeLO.push_back(ComputeLOarray[i]);
    assert(ComputeLO.back()=="DefaultMode");

    vector<string> AdmissibleOut;
    string AdmissibleOutarray[]={"gen","cst","inv","ext","ht1","esp","egn","typ","lat","mod"}; // "mod" must be last
    for(size_t i=0;i<10;++i)
        AdmissibleOut.push_back(AdmissibleOutarray[i]);
    assert(AdmissibleOut.back()=="mod");

    // analyzing long options
    for(size_t i=0; i<LongOptions.size();++i){
        if(find(ComputeLO.begin(),ComputeLO.end(),LongOptions[i])!=ComputeLO.end()){
            to_compute.set(LongOptions[i]);
            continue;
        }
        if(find(AdmissibleOut.begin(),AdmissibleOut.end(),LongOptions[i])!=AdmissibleOut.end()){
            OutFiles.push_back(LongOptions[i]);
            continue;
        }
        if(LongOptions[i]=="Ignore"){
            ignoreInFileOpt=true;
            continue;
        }
        if(LongOptions[i]=="KeepOrder"){
            keep_order=true;
            continue;
        }
        if(LongOptions[i]=="BottomDec"){
            do_bottom_dec=true;
            continue;
        }
        if(LongOptions[i]=="Console"){
            verbose=true;
            continue;
        }
        if(LongOptions[i]=="Files"){
            write_extra_files = true;
            continue;
        }
        if(LongOptions[i]=="AllFiles"){
            write_all_files = true;
            continue;
        }
        if(LongOptions[i]=="BigInt"){
            use_Big_Integer=true;
            continue;
        }
        cerr<<"Warning: Unknown option --" << LongOptions[i]<<endl;
    }

    // activate default mode
    if (to_compute.none()) {
        to_compute.set(ConeProperty::DefaultMode);
    }

	if(keep_order)
		to_compute.set(ConeProperty::KeepOrder);

	if(do_bottom_dec)
		to_compute.set(ConeProperty::BottomDecomposition);

	return false; //no need to print help text
}

template<typename Integer>
void OptionsHandler::applyOutputOptions(Output<Integer>& Out) {
    if(write_all_files) {
        Out.set_write_all_files();
    } else if (write_extra_files) {
        Out.set_write_extra_files();
    }
    if (to_compute.test(ConeProperty::Triangulation)) {
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
        if(OutFiles[i]=="mod"){
            Out.set_write_mod(true);
            continue;
        }
    }

    Out.set_name(output_name);
}

bool OptionsHandler::anyNmzIntegrateOption() const {
    return nmzInt_E || nmzInt_I || nmzInt_L;
}

string OptionsHandler::getNmzOptions() const {
    string nmz_options;
    if (verbose) {
        nmz_options.append(" -c");
    }
    if (nr_threads > 0) {
        nmz_options.append(" -x=");
        ostringstream convert;
        convert << nr_threads;
        nmz_options.append(convert.str());
    }
    if (nmzInt_E) {
        nmz_options.append(" -E");
    }
    if (nmzInt_L) {
        nmz_options.append(" -L");
    }
    if (nmzInt_I) {
        nmz_options.append(" -I");
    }
    nmz_options.append(" \"");
    nmz_options.append(output_name);
    nmz_options.append("\"");
    return nmz_options;
}
