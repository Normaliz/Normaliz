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

#include "Normaliz.h"
#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/cone.h"
//#include "libnormaliz/libnormaliz.cpp"
using namespace libnormaliz;
#include "Input.cpp"
#include "Options.cpp"
#include "output.cpp"

#ifndef STRINGIFY
#define STRINGIFYx(Token) #Token
#define STRINGIFY(Token) STRINGIFYx(Token)
#endif

void printHeader() {
    cout << "                                                    \\.....|"<<endl;
    cout << "                    Normaliz " << string( STRINGIFY(NMZ_VERSION) "           " ,11)
                                                 << "             \\....|"<<endl;
    cout << "                                                      \\...|"<<endl;
    cout << "      (C) W. Bruns  B. Ichim  T. Roemer  C. Soeger     \\..|"<<endl;
    cout << "                       May 2015                         \\.|"<<endl;
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

void printCopying() {
    cout<<"Copyright (C) 2007-2015  Winfried Bruns, Bogdan Ichim, Christof Soeger"<<endl
        <<"This program comes with ABSOLUTELY NO WARRANTY; This is free software,"<<endl
        <<"and you are welcome to redistribute it under certain conditions;"<<endl
        <<"See COPYING for details."<<endl;
}

void printVersion() {
    cout << "Normaliz " << string(STRINGIFY(NMZ_VERSION)) << endl;
    printCopying();
}

template<typename Integer> int process_data(OptionsHandler& options);

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{

    // read command line options

    OptionsHandler options;

    bool print_help = options.handle_commandline(argc, argv);

	if (print_help) {
        printHeader();
        printHelp(argv[0]);
        exit(0);
	}

	if (verbose) {
        printHeader();
    }

    process_data<mpz_class>(options);

    if (options.anyNmzIntegrateOption()) {
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
            cerr << "Error: Could not start nmzIntegrate" << endl;
            return 2;
        }
        nmz_int_exec.append("\"");

        nmz_int_exec.append(options.getNmzIntegrateOptions());

        #ifdef _WIN32 //for 32 and 64 bit windows
            nmz_int_exec.append("\"");
        #endif

        cout << "executing: "<< nmz_int_exec << endl;
        return system(nmz_int_exec.c_str());
    }
}

//---------------------------------------------------------------------------

template<typename Integer> int process_data(OptionsHandler& options) {

#ifndef NCATCH
    std::exception_ptr excep;
    try {
#endif

    Output<Integer> Out;    //all the information relevant for output is collected in this object

    options.applyOutputOptions(Out);

    string name_in=options.getOutputName()+".in";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if ( !in.is_open() ) {
        cerr<<"error: Failed to open file "<<name_in<<"."<<endl;
        return 1;
    }

    //read the file
    map <Type::InputType, vector< vector<Integer> > > input = readNormalizInput<Integer>(in, options);

    options.activateDefaultMode(); // only if no real cone property is given!

    Out.set_lattice_ideal_input(input.count(Type::lattice_ideal)>0);

    in.close();

    if (verbose) {
        cout << "************************************************************" << endl;
        cout << "Compute: " << options.getToCompute() << endl;
    }

    Cone<Integer> MyCone = Cone<Integer>(input);
    if (options.isUseBigInteger()) {
        MyCone.deactivateChangeOfPrecision();
    }
//    MyCone.compute(ConeProperty::HilbertBasis,ConeProperty::HilbertSeries));
    try {
        MyCone.compute(options.getToCompute());
    } catch(const NotComputableException& e) {
        std::cout << "Not all desired properties could be computed." << endl;
        std::cout << "Writing only available data." << endl;
    }
    Out.setCone(MyCone);
    Out.write_files();

#ifndef NCATCH
    } catch(const NormalizException& e) {
        excep = std::current_exception();
        std::cout << "NormalizException caught...";
    } catch(const std::exception& e) {
        excep = std::current_exception();
        std::cout << "std::exception caught...";
    }

    if (excep != std::exception_ptr()) {
        cout << "rethrowing the exception!" << endl;
        std::rethrow_exception (excep);
    }
#endif

    return 0;
}
