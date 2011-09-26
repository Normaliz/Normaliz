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
#include "output.cpp"



void printHelp(char* command) {
	cout << "usage: "<<command<<" [-sSvVnNpPhH1dBface?] [-x=<T>] [PROJECT]"<<endl;
	cout << "  runs normaliz on PROJECT.in"<<endl;
	cout << "options:"<<endl;
	cout << "  -?\tprint this help text and exit"<<endl;
	cout << "  -s\tcomputation mode: support hyperplanes"<<endl;
	cout << "  -S\tcomputation mode: support hyperplanes (currently same as -s)"<<endl;
	cout << "  -v\tcomputation mode: volume triangulation"<<endl;
	cout << "  -V\tcomputation mode: volume large"<<endl;
	cout << "  -n\tcomputation mode: Hilbert basis triangulation (previously normal)"<<endl;
	cout << "  -N\tcomputation mode: Hilbert basis (using a partial triangulation)"<<endl;
	cout << "  -p\tcomputation mode: Hilbert polynomial"<<endl;
	cout << "  -P\tcomputation mode: Hilbert polynomial large"<<endl;
	cout << "  -h\tcomputation mode: Hilbert basis polynomial"<<endl;
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
	ComputationMode computation_mode = Mode::hilbertBasisTriangulation;
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
				computation_mode = Mode::hilbertPolynomial;
				break;
			case 'P':
				computation_mode = Mode::hilbertPolynomialLarge;
				break;
			case 'h':
				computation_mode = Mode::hilbertBasisPolynomial;
				break;
			case 'H':
				computation_mode = Mode::hilbertBasisPolynomialLarge;
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
			//	optimize_speed=false;
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
	size_t i,j;

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
	if (in.is_open()==false) {
		cerr<<"error: Failed to open file "<<name_in<<"."<<endl;
		return 1;
	}




	string mode_string;
	size_t nr_rows,nr_columns;;
	int mode;
	InputType input_type = Type::integral_closure;
	Integer number;
	in >> nr_rows;
	in >> nr_columns;
	Matrix<Integer> M(nr_rows,nr_columns);
	for(i=1; i<=nr_rows; i++){
		for(j=1; j<=nr_columns; j++) {
			in >> number;
			M.write(i,j,number);
		}
	}

	in>>mode_string;
	if (mode_string=="0"||mode_string=="integral_closure") {
		mode=0;
		input_type = Type::integral_closure;
	} else
	if (mode_string=="1"||mode_string=="normalization") {
		mode=1;
		input_type = Type::normalization;
	} else
	if (mode_string=="2"||mode_string=="polytope") {
		mode=2;
		input_type = Type::polytope;
	} else
	if (mode_string=="3"||mode_string=="rees_algebra") {
		mode=3;
		input_type = Type::rees_algebra;
	} else
	if (mode_string=="4"||mode_string=="hyperplanes") {
		mode=4;
	} else
	if (mode_string=="5"||mode_string=="equations") {
		mode=5;
	} else
	if (mode_string=="6"||mode_string=="congruences") {
		mode=6;
	} else
	if (mode_string=="10"||mode_string=="lattice_ideal") {
		mode=10;
		input_type = Type::lattice_ideal;
	} else {
		cerr<<"Warning: Unknown mode "<<mode_string<<" and will be replaced with mode integral_closure."<<endl;
		mode=0;
		input_type = Type::integral_closure;
	}

	if ( in.fail() ) {
		cerr << "error: Failed to read file "<<name_in<<". May be a bad format of the input file."<<endl;
		return 1;
	}

	Out.set_name(output_name);
	
	if (mode >= 4 && mode <= 6) { //equations, inequalities, congruences
		size_t nc = nr_columns;
		if (mode == 6) {
			nc--;  //the congruence matrix has one extra column
		}
		multimap <Type::InputType, vector< vector<Integer> > > constraints;

		while (true) {
			switch(mode) {
				case 4:
					constraints.insert(pair<Type::InputType, vector< vector<Integer> > > (Type::hyperplanes, M.get_elements()));
					break;
				case 5:
					constraints.insert(pair<Type::InputType, vector< vector<Integer> > > (Type::equations, M.get_elements()));
					break;
				case 6:
					constraints.insert(pair<Type::InputType, vector< vector<Integer> > > (Type::congruences, M.get_elements()));
					break;
				default:
					cerr<<"Reached unreachable code in Normaliz.cpp. Please contact the developers"<<endl;
					return 1;
			}
			if (!in.good()) 
				break;
			in >> nr_rows;
			in >> nr_columns;
			if (in.eof())
				break;  //not enough data to read on
			M = Matrix<Integer>(nr_rows,nr_columns);
			for(i=1; i<=nr_rows; i++){
				for(j=1; j<=nr_columns; j++) {
					in >> number;
					M.write(i,j,number);
				}
			}

			in>>mode_string;
			if ( in.fail() ) {
				cerr << "error: Failed to read file "<<name_in<<". May be a bad format of the input file?"<<endl;
				return 1;
			}

			if	 (mode_string=="4"||mode_string=="hyperplanes") {
				mode=4;
				if(nr_columns!=nc) {
					cerr<<"Number of columns not matching. Check input file!";
					return 1;
				}
			} else
			if (mode_string=="5"||mode_string=="equations") {
				mode=5;
				if(nr_columns!=nc) {
					cerr<<"Number of columns not matching. Check input file!";
					return 1;
				}
			} else
			if (mode_string=="6"||mode_string=="congruences") {
				mode=6;
				if(nr_columns!=nc+1) {
					cerr<<"Number of columns not matching. Check input file!";
					return 1;
				}
			} else {
				cerr<<"Illegal input type \""<<mode_string<<"\" at this position."<<endl;
				return 1;
			}
		}
		
		in.close();
		if (verbose) {
			cout<<"\n************************************************************\n";
			cout<<"Running in computation mode "<<computation_mode<<" with input type constraints."<<endl;
		}
		Cone<Integer> MyCone = Cone<Integer>(constraints);
		MyCone.compute(computation_mode);
		Out.setCone(MyCone);
		Out.cone();
	} 
	else { // all other modes
		in.close();
		//main computations and output
		if (verbose) {
			cout<<"\n************************************************************\n";
			cout<<"Running in computation mode "<<computation_mode<<" with input type "<<mode<<"."<<endl;
		}
		Cone<Integer> MyCone = Cone<Integer>(M.get_elements(), input_type);
//		MyCone.compute(ConeProperties(ConeProperty::HilbertBasis,ConeProperty::HilbertPolynomial));
		MyCone.compute(computation_mode);
		Out.setCone(MyCone);
		if (mode == 2) {
			Out.polytop();
		} else if (mode == 3) {
			Out.rees();
		} else {
			Out.cone();
		}
	}
	return 0;
}
