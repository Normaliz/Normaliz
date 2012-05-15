/*
 * Normaliz 2.7
 * Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger
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



void printHelp(char* command) {
    cout << "usage: "<<command<<" [-sSvVnNpPhH1dBface?] [-x=<T>] [PROJECT]"<<endl;
    cout << "  runs normaliz on PROJECT.in"<<endl;
    cout << "options:"<<endl;
    cout << "  -?\tprint this help text and exit"<<endl;
    cout << "  -s\tcomputation mode: support hyperplanes"<<endl;
    cout << "  -S\tcomputation mode: support hyperplanes (currently same as -s)"<<endl;
    cout << "  -t\tcomputation mode: triangulation size"<<endl;
    cout << "  -v\tcomputation mode: volume triangulation"<<endl;
    cout << "  -V\tcomputation mode: volume large"<<endl;
    cout << "  -n\tcomputation mode: Hilbert basis triangulation (previously normal)"<<endl;
    cout << "  -N\tcomputation mode: Hilbert basis (using a partial triangulation)"<<endl;
    cout << "  -p\tcomputation mode: Hilbert polynomial"<<endl;
    cout << "  -P\tcomputation mode: Hilbert polynomial large"<<endl;
    cout << "  -h\tcomputation mode: Hilbert basis polynomial (default)"<<endl;
    cout << "  -H\tcomputation mode: Hilbert basis polynomial large"<<endl;
    cout << "  -1\tcomputation mode: height 1 elements"<<endl;
    cout << "  -d\tcomputation mode: dual"<<endl;
    cout << "  -f\tthe files .out .gen .inv .typ .cst are written"<<endl;
    cout << "  -a\tall output files are written"<<endl;
    cout << "  -e\tperform tests for arithmetic errors"<<endl;
    cout << "  -B\tuse indefinite precision arithmetic"<<endl;
    cout << "  -c\tverbose (prints control data)"<<endl;
    cout << "  -m\tsave memory (currently has no effect)"<<endl;
    cout << "  -i\tobsolete option"<<endl;
    cout << "  -x=<T>\tlimit the number of threads to <T>"<<endl;
}

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{

    //libnormaliz::RecBoundFactor = 5000000;
    size_t i;       //used for iterations
    char c;
    ComputationMode computation_mode = Mode::hilbertBasisSeries;
    //for available modes see libnormaliz/libnormaliz.h
    //it is set by the option from the command line
    string output_name;         //name of the output file(s) saved here

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
                    size_t nr_threads;
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
            case 's':
                computation_mode = Mode::supportHyperplanes;
                break;
            case 'S':
                computation_mode = Mode::supportHyperplanes;
                break;
            case 't':
                computation_mode = Mode::triangulationSize;
                break;
            case 'v':
                computation_mode = Mode::volumeTriangulation;
                break;
            case 'V':
                computation_mode = Mode::volumeLarge;
                break;
            case 'n':
                computation_mode = Mode::hilbertBasisTriangulation;
                break;
            case 'N':
                computation_mode = Mode::hilbertBasisLarge;
                break;
            case '1':
                computation_mode = Mode::height1Elements;
                break;
            case 'p':
                computation_mode = Mode::hilbertSeries;
                break;
            case 'P':
                computation_mode = Mode::hilbertSeriesLarge;
                break;
            case 'h':
                computation_mode = Mode::hilbertBasisSeries;
                break;
            case 'H':
                computation_mode = Mode::hilbertBasisSeriesLarge;
                break;
            case 'd':
                computation_mode = Mode::dual;
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
                printHelp(argv[0]);
                exit(1);
                break;
            case 'x': //should be separated from other options
                cerr<<"Warning: Option -x=<T> has to be separated from other options"<<endl;
                break;
            default:
                cerr<<"Warning: Unknown option -"<<option[i]<<endl;
                break;
        }
    }

    if (!filename_set) {
        cout<<"Normaliz 2.7"<<endl
            <<"Copyright (C) 2007-2011  Winfried Bruns, Bogdan Ichim, Christof Soeger"<<endl
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

    int returnvalue;

    if(use_Big_Integer) {
        //if the program works with the indefinite precision arithmetic, no arithmetic tests are performed
        test_arithmetic_overflow=false;
        //Read and process Input
        returnvalue = process_data<mpz_class>(output_name, computation_mode, write_extra_files, write_all_files);
    } else {
        //Read and process Input
        returnvalue = process_data<long long int>(output_name, computation_mode, write_extra_files, write_all_files);
    }

    //exit
    if (!filename_set) {
        cout<< "\nType something and press enter to exit.\n";
        cin >> c;
    }
    return returnvalue;
}

//---------------------------------------------------------------------------

template<typename Integer> int process_data(string& output_name, ComputationMode computation_mode, bool write_extra_files, bool write_all_files ) {
    vector<Integer> Grading;

    Output<Integer> Out;    //all the information relevant for output is collected in this object

    if(write_all_files) {
        Out.set_write_all_files();
    } else if (write_extra_files) {
        Out.set_write_extra_files();
    }

    string name_in=output_name+".in";
    const char* file_in=name_in.c_str();
    ifstream in, in2;
    in2.open(file_in,ifstream::in);
    if (in2.is_open()==false) {
        //check if user added ".in" and ignore it in this case
        string suffix (".in");
        size_t found = output_name.rfind(suffix);
        if (found!=string::npos) {
            output_name.erase(found);
            name_in=output_name+".in";
            file_in=name_in.c_str();
        }
    } else {
        in2.close();
    }
    in.open(file_in,ifstream::in);
    if ( !in.is_open() ) {
        cerr<<"error: Failed to open file "<<name_in<<"."<<endl;
        return 1;
    }

    Out.set_name(output_name);

    //read the file
    map <Type::InputType, vector< vector<Integer> > > input = readNormalizInput (in, Out);

    in.close();

    //don't save the triangulation if the user doesn't want to see it
    //and we don't need it for the primary multiplicity later
    if (!write_all_files && input.count(Type::rees_algebra)==0) {
        if (computation_mode == Mode::volumeTriangulation)
            computation_mode  = Mode::volumeLarge;
        if (computation_mode == Mode::hilbertBasisTriangulation)
            computation_mode  = Mode::hilbertBasisMultiplicity;
        if (computation_mode == Mode::hilbertSeries)
            computation_mode  = Mode::hilbertSeriesLarge;
        if (computation_mode == Mode::hilbertBasisSeries)
            computation_mode  = Mode::hilbertBasisSeriesLarge;
    }

    if (verbose) {
        cout<<"\n************************************************************\n";
        cout<<"Running in computation mode "<<computation_mode<<"."<<endl;
    }

    Cone<Integer> MyCone = Cone<Integer>(input);
//    MyCone.compute(ConeProperties(ConeProperty::HilbertBasis,ConeProperty::HilbertSeries));
    MyCone.compute(computation_mode);
    Out.setCone(MyCone);
    Out.write_files();

    return 0;
}
