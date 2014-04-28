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

#include <sstream>
#include <algorithm>
using namespace std;

#include "Normaliz.h"
#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/cone.h"
//#include "libnormaliz/libnormaliz.cpp"
using namespace libnormaliz;
#include "Input.cpp"
#include "output.cpp"



void printHeader() {
    cout << "                                                    \\.....|"<<endl;
    cout << "                     Normaliz 2.11                   \\....|"<<endl;
    cout << "                                                      \\...|"<<endl;
    cout << "      (C) W. Bruns  B. Ichim  T. Roemer  C. Soeger     \\..|"<<endl;
    cout << "                       April 2014                       \\.|"<<endl;
    cout << "                                                         \\|"<<endl;
}
void printHelp(char* command) {
    cout << "usage: "<<command<<" [-sSvVnNpPhH1dBface?] [-x=<T>] [PROJECT]"<<endl;
    cout << "  runs normaliz on PROJECT.in"<<endl;
    cout << "options:"<<endl;
    cout << "  -?\tprint this help text and exit"<<endl;
    cout << "  -s\tcompute support hyperplanes"<<endl;
    cout << "  -t\tcompute triangulation"<<endl;
    cout << "  -v\tcompute volume"<<endl;
    cout << "  -n\tcompute Hilbert basis (with full triangulation)"<<endl;
    cout << "  -N\tcompute Hilbert basis (with partial triangulation)"<<endl;
    cout << "  -q\tcompute Hilbert (quasi-)polynomial"<<endl;
    cout << "  -p\tcompute Hilbert (quasi-)polynomial and degree 1 elements"<<endl;
    cout << "  -h\tcompute Hilbert basis and Hilbert polynomial (default)"<<endl;
    cout << "  -1\tcompute degree 1 elements"<<endl;
    cout << "  -y\tcompute Stanley decomposition"<<endl;
    cout << "  -d\tcomputation mode: dual"<<endl;
    cout << "  -r\tcomputation mode: approximate rational polytope"<<endl;
    cout << "  -f\tthe files .out .gen .inv .cst are written"<<endl;
    cout << "  -T\tthe file .tri is written (triangulation)"<<endl;
    cout << "  -a\tall output files are written (except .tri)"<<endl;
    cout << "  -e\tperform tests for arithmetic errors"<<endl;
    cout << "  -B\tuse indefinite precision arithmetic"<<endl;
    cout << "  -c\tverbose (prints control data)"<<endl;
    cout << "  -x=<T>\tlimit the number of threads to <T>"<<endl;
}

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{

    size_t i;       //used for iterations
    char c;
    string output_name;         //name of the output file(s) saved here
    long nr_threads = 0;

    // read command line options
    bool filename_set=false;
    string option;            //all options concatenated (including -)
    
    for (i = 1; i < (unsigned int)argc; i++) {
        if (argv[i][0]=='-') {
            if (argv[i][1]!='\0') {
                if (argv[i][1]!='x') {
                    option = option + argv[i];
                } else if (argv[i][2]=='=') {
                    #ifdef _OPENMP
                    string Threads = argv[i];
                    Threads.erase(0,3);
                    if ( (istringstream(Threads) >> nr_threads) && nr_threads > 0) {
                        omp_set_num_threads(nr_threads);
                    } else {
                        cerr<<"Warning: Invalid option string "<<argv[i]<<endl;
                    }
                    #else
                    cerr << "Warning: Compiled without OpenMP support, option "<<argv[i]<<" ignored."<<endl;
                    #endif
                } else {
                    cerr<<"Warning: Invalid option string "<<argv[i]<<endl;
                }
            }
        } else if (!filename_set) {
            string s(argv[i]);
            output_name=s;
            filename_set=true;
        }
    }



    //Analyzing the command line options
    bool write_extra_files = false, write_all_files = false;
    bool use_Big_Integer = false;
    ConeProperties to_compute;
    bool nmzInt_E = false, nmzInt_I = false, nmzInt_L = false;

    for (i = 1; i <option.size(); i++) {
        switch (option[i]) {
            case '-':
            case 'i':
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
            case 'S':
                to_compute.set(ConeProperty::SupportHyperplanes);
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
                test_arithmetic_overflow=true;
                break;
            case 'B':  //use Big Integer
                use_Big_Integer=true;
                break;
            case 'm':  //save memory / don't optimize for speed
            //    optimize_speed=false;
                break;
            case '?':  //print help text and exit
                printHeader();
                printHelp(argv[0]);
                exit(1);
                break;
            case 'x': //should be separated from other options
                cerr<<"Warning: Option -x=<T> has to be separated from other options"<<endl;
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
            case 'E':  //nmzIntegrate -E (ehrhart series)
                nmzInt_E = true;
                to_compute.set(ConeProperty::StanleyDec);
                break;
            default:
                cerr<<"Warning: Unknown option -"<<option[i]<<endl;
                break;
        }
    }
    // activate default mode
    if (to_compute.none()) {
        to_compute.set(ConeProperty::DefaultMode);
    }

    if (verbose || !filename_set) {
        printHeader();
    }
    if (!filename_set) {
        cout<<"Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger"<<endl
            <<"This program comes with ABSOLUTELY NO WARRANTY; This is free software,"<<endl
            <<"and you are welcome to redistribute it under certain conditions;"<<endl
            <<"See COPYING for details."
            <<endl<<endl;
        cout<<"Enter the input file name or -? for help: ";
        cin >>output_name;
        if (output_name == "-?") {
            printHelp(argv[0]);
            return 1;
        }
    }

    // check if we cand read the .in file
    string name_in=output_name+".in";
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

    int returnvalue;

    if(use_Big_Integer) {
        //if the program works with the indefinite precision arithmetic, no arithmetic tests are performed
        test_arithmetic_overflow=false;
        //Read and process Input
        returnvalue = process_data<mpz_class>(output_name, to_compute, write_extra_files, write_all_files);
    } else {
        //Read and process Input
        returnvalue = process_data<long long int>(output_name, to_compute, write_extra_files, write_all_files);
    }

    if (returnvalue == 0 && (nmzInt_E || nmzInt_I || nmzInt_L) ) {
        //cout << "argv[0] = "<< argv[0] << endl;
        string nmz_int_exec("\"");
        // the quoting requirements for windows are insane, one pair of "" around the whole command and one around each file
        #ifdef _WIN32 //for 32 and 64 bit windows
            nmz_int_exec.append("\"");
        #endif	
        nmz_int_exec.append(argv[0]);
        size_t found = nmz_int_exec.rfind("normaliz");
        if (found!=std::string::npos) {
            nmz_int_exec.replace (found,8,"nmzIntegrate");
        } else {
            cout << "ERROR: Could not start nmzIntegrate" << endl;
            return 2;
        }
        nmz_int_exec.append("\"");

        if (verbose) {
            nmz_int_exec.append(" -c");
        }
        if (nr_threads > 0) {
            nmz_int_exec.append(" -x=");
            ostringstream convert;
            convert << nr_threads;
            nmz_int_exec.append(convert.str());
        }
        if (nmzInt_E) {
            nmz_int_exec.append(" -E");
        }
        if (nmzInt_L) {
            nmz_int_exec.append(" -L");
        }
        if (nmzInt_I) {
            nmz_int_exec.append(" -I");
        }
        nmz_int_exec.append(" \"");
        nmz_int_exec.append(output_name);
        nmz_int_exec.append("\"");
        #ifdef _WIN32 //for 32 and 64 bit windows
            nmz_int_exec.append("\"");
        #endif

        cout << "executing: "<< nmz_int_exec << endl;
        returnvalue = system(nmz_int_exec.c_str());
    }
    //exit
    if (!filename_set) {
        cout<< "\nType something and press enter to exit.\n";
        cin >> c;
    }

    return returnvalue;
}

//---------------------------------------------------------------------------

template<typename Integer> int process_data(string& output_name, ConeProperties to_compute, bool write_extra_files, bool write_all_files ) {

    Output<Integer> Out;    //all the information relevant for output is collected in this object

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


    string name_in=output_name+".in";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if ( !in.is_open() ) {
        cerr<<"error: Failed to open file "<<name_in<<"."<<endl;
        return 1;
    }

    Out.set_name(output_name);

    //read the file
    map <Type::InputType, vector< vector<Integer> > > input = readNormalizInput<Integer>(in);

    in.close();

    if (verbose) {
        cout << "************************************************************" << endl;
        cout << "Compute: " << to_compute << endl;
    }

    Cone<Integer> MyCone = Cone<Integer>(input);
//    MyCone.compute(ConeProperties(ConeProperty::HilbertBasis,ConeProperty::HilbertSeries));
    MyCone.compute(to_compute);
    Out.setCone(MyCone);
    Out.write_files();

    return 0;
}
