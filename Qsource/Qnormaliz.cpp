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

#include "Qnormaliz.h"
#include "libQnormaliz/Qinteger.h"
#include "libQnormaliz/libQnormaliz.h"
#include "libQnormaliz/Qcone.h"
#include "libQnormaliz/Qmy_omp.h"
//#include "libnormaliz/libnormaliz.cpp"
using namespace libQnormaliz;
#include "Qinput.cpp"
#include "Qoptions.cpp"
#include "Qoutput.cpp"

#ifndef STRINGIFY
#define STRINGIFYx(Token) #Token
#define STRINGIFY(Token) STRINGIFYx(Token)
#endif

void printHeader() {
    cout << "                                                    \\.....|"<<endl;
    cout << "                    Normaliz " << string( STRINGIFY(NMZ_VERSION) "           " ,11)
                                                 << "             \\....|"<<endl;
    cout << "                                                      \\...|"<<endl;
    cout << "     (C) The Normaliz Team, University of Osnabrueck   \\..|"<<endl;
    cout << "                    November  2016                      \\.|"<<endl;
    cout << "                                                         \\|"<<endl;
}
void printHelp(char* command) {
    cout << "Usage: "<<command<<" [options] PROJECT"<<endl;
    cout << "  runs normaliz on PROJECT.in"<<endl;
    cout << "Options:"<<endl;
    cout << "  -S\tcompute sublattice"<<endl;
    cout << "  -s\tcompute support hyperplanes"<<endl;
    cout << "  -t\tcompute triangulation"<<endl;
    cout << "  -v\tcompute volume"<<endl;
    cout << "  -n\tcompute Hilbert basis and volume (needs full triangulation)"<<endl;
    cout << "  -N\tcompute Hilbert basis (with partial triangulation)"<<endl;
    cout << "  -w\tcheck for integrally closed and compute witness if not"<<endl;
    cout << "  -q\tcompute Hilbert (quasi-)polynomial"<<endl;
    cout << "  -p\tcompute Hilbert (quasi-)polynomial and degree 1 elements"<<endl;
    cout << "  -h\tcompute Hilbert basis and Hilbert polynomial (default)"<<endl;
    cout << "  -1\tcompute degree 1 elements"<<endl;
    cout << "  -y\tcompute Stanley decomposition (output in file .dec)"<<endl;
    cout << "  -C\tcompute class group (default)"<<endl;
    cout << "  -T\tcompute triangulation  (output in file .tri)"<<endl;
    cout << "  -D\tcompute cone decomposition (includes -T)"<<endl;
    cout << "  -H\tcompute integer hull"<<endl;
    cout << "  -M\tcompute module generators over original monoid"<<endl;

    cout << endl;
    cout << "  -d\tcomputation mode: dual"<<endl;
    cout << "  -P\tcomputation mode: primal"<<endl;
    cout << "  -r\tcomputation mode: approximate"<<endl;
    cout << "  -b\tcomputation mode: bottom decomposition"<<endl;
    cout << "  -o\tcomputation mode: no bottom decomposition"<<endl;
    cout << "  -k\tcomputation mode: keep order"<<endl;
    cout << "  -Y\tcomputation mode: symmetrization"<<endl;
    cout << "  -X\tcomputation mode: no symmetrization"<<endl;


    cout << endl;
    cout << "      --<PROP>     compute the ConeProperty <PROP>"<<endl;

    cout << endl;
    cout << "  -f, --files      write the files .out .gen .inv .cst"<<endl;
    cout << "  -a, --all-files  write all output files (except .tri)"<<endl;
    cout << "      --<SUFFIX>   write the file .<SUFFIX> where <SUFFIX> can be one of"<<endl;
    cout << "                   cst, egn, esp, ext, gen, ht1, inv, lat, mod, typ"<<endl;

    cout << endl;
    cout << "  -B, --BigInt     directly use indefinite precision arithmetic"<<endl;
    cout << "      --LongLong   only use long long arithmetic, no conversion possible"<<endl;
    cout << "  -i, --ignore     ignore the compute options set in the input file"<<endl;
    cout << "  -x=<T>           limit the number of threads to <T>"<<endl;
    cout << "  --OutputDir=<path> set a path for the output files (relative to current directory)"<< endl;
    cout << "  -?, --help       print this help text and exit"<<endl;
    cout << "  -c, --verbose    verbose (prints control data)"<<endl;
    cout << "      --version    print version info and exit"<<endl;
    cout << endl;
    cout << "Please report bugs to <normaliz@uos.de> or directly to our issue tracker:" << endl;
    cout << "https://github.com/Normaliz/Normaliz/issues" << endl;
}

void printCopying() {
    cout<<"Copyright (C) 2007-2015  The Normaliz Team, University of Osnabrueck."<<endl
        <<"This program comes with ABSOLUTELY NO WARRANTY; This is free software,"<<endl
        <<"and you are welcome to redistribute it under certain conditions;"<<endl
        <<"See COPYING for details."<<endl;
}

void printVersion() {
    cout << "Normaliz " << string(STRINGIFY(NMZ_VERSION)) << endl;
    printCopying();
}

template<typename Number> int process_data(OptionsHandler& options, const string& command_line);

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{

    // read command line options

    OptionsHandler options;
    
    string command_line;
    for(int i=1; i< argc;++i)
        command_line=command_line+string(argv[i])+" ";

    bool print_help = options.handle_commandline(argc, argv);

    if (print_help) {
        //printHeader();
        printHelp(argv[0]);
        exit(0);
    }

    if (verbose) {
        printHeader();
    }

    if (!options.isUseLongLong()) {
        process_data<mpq_class>(options, command_line);
    }
    // the previous process_data might return unsuccessfully if the input file specifies to use long long

}

//---------------------------------------------------------------------------

template<typename Number> int process_data(OptionsHandler& options, const string& command_line) {

#ifndef NCATCH
    try {
#endif

    Output<Number> Out;    //all the information relevant for output is collected in this object

    options.applyOutputOptions(Out);

    string name_in=options.getOutputName()+".in";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if ( !in.is_open() ) {
        cerr << "error: Failed to open file "<<name_in<<"."<<endl;
        exit(1);
    }

    //read the file
    map <Type::InputType, vector< vector<Number> > > input = readNormalizInput<Number>(in, options);

    options.activateDefaultMode(); // only if no real cone property is given!

    Out.set_lattice_ideal_input(input.count(Type::lattice_ideal)>0);

    in.close();


    if (verbose) {
        cout << "************************************************************" << endl;
        cout << "Command line: " << command_line << endl;
        cout << "Compute: " << options.getToCompute() << endl;
    }

    Cone<Number> MyCone = Cone<Number>(input);
    long dim= (long) MyCone.getEmbeddingDim();
 #ifdef _OPENMP
    long max_threads=omp_get_max_threads();
    if(!options.nr_threads_explicitly_set && std::getenv("OMP_NUM_THREADS")==NULL){
        max_threads=min(max_threads,4*dim); // we limit the implicit number of threads
        omp_set_num_threads(max_threads);
    }
#endif
    /* if (options.isUseBigNumber()) {
        MyCone.deactivateChangeOfPrecision(); 
    } */
    try {
        MyCone.compute(options.getToCompute());
    } catch(const NotComputableException& e) {
        std::cout << "Not all desired properties could be computed." << endl;
        std::cout << e.what() << endl;
        std::cout << "Writing only available data." << endl;
    }
    Out.setCone(MyCone);
    Out.write_files();

#ifndef NCATCH
    } catch(const BadInputException& e) {
        cerr << e.what() << endl;
        cerr << "BadInputException caught... exiting." << endl;
        exit(1);
    } catch(const FatalException& e) {
        cerr << e.what() << endl;
        cerr << "FatalException caught... exiting." << endl;
        exit(2);
    } catch(const NormalizException& e) {
        cerr << e.what() << endl;
        cerr << "NormalizException caught... exiting." << endl;
        exit(3);
    } catch(const std::exception& e) {
        cerr << "std::exception caught... \""<< e.what()<<"\" ...  exiting." << endl;
        exit(4);
    }
#endif

    return 0;
}
