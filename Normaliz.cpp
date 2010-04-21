/*
 * Normaliz 2.2
 * Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim
 * With contributions by Christof Soeger
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
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

#include "Normaliz.h"


// this function determinates if and how the program will be terminated in case of errors
void global_error_handling(){
	cout<<"Some error detected. The program will be terminated."<<endl;
	exit(1);
}


void printHelp(char* command) {
	cout << "usage: "<<command<<" [-acdefhimnpsv?] [PROJECT]"<<endl;
	cout << "  runs normaliz on PROJECT.in"<<endl;
	cout << "options:"<<endl;
	cout << "  -?\tprint this help text and exit"<<endl;
	cout << "  -s\tcomputation type: support_hyperplanes"<<endl;
	cout << "  -S\tcomputation type: support_hyperplanesvia pyramids"<<endl;
	cout << "  -v\tcomputation type: triangulation"<<endl;
	cout << "  -V\tcomputation type: triangulation via pyramids"<<endl;
	cout << "  -n\tcomputation type: triangulation_hilbert_basis (old type normal)"<<endl;
	cout << "  -N\tcomputation type: hilbert_basis (using a partial triangulation)"<<endl;
	cout << "  -1\tcomputation type: ht1_elements"<<endl;
	cout << "  -p\tcomputation type: hilbert_polynomial"<<endl;
	cout << "  -h\tcomputation type: hilbert_basis_polynomial"<<endl;
	cout << "  -d\tcomputation type: dual"<<endl;
	cout << "  -f\tthe files .out .gen .inv .typ .sup are written"<<endl;
	cout << "  -a\tall output files are written"<<endl;
	cout << "  -e\tperform tests for arithmetic errors"<<endl;
	cout << "  -c\tverbose (prints control data)"<<endl;
	cout << "  -m\tsave memory"<<endl;
	cout << "  -i\tthe config file normaliz.cfg will be ignored"<<endl;
}

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	int i,j;       //used for iterations
	char c;
	string computation_type="triangulation_hilbert_basis";
	//4  types available, "support_hyperplanes", "triangulation", "normal" and "hilbert_polynomial"
	//it is set by the setup file or by the options "s", "v", "n", "p", "h" and "d" in the command line
	//the type given in the command line overrides the type set by the setup file
	string output_name;         //name of the output file(s) saved here
	Output Out;                //all the information relevant for output is collected
	//in this object

	// read commandline options 
	bool read_setup=true;
	bool filename_set=false;
	string option;            //all options concatenated (including -)
	for (i = 1; i <argc; i++) {
		if (argv[i][0]=='-') {
			option = option + argv[i];
		} else if (!filename_set) {
			string s(argv[i]);
			output_name=s;
			filename_set=true;
		}
	}


	//determine whether the config file should be ignored
	for (i = 0; i <option.size(); i++) {
		if (option[i]=='i') {
			read_setup=false;
		}
	}


	//read config file

	if (read_setup==true) {
		ifstream setup("normaliz.cfg");
		if (setup.is_open()==false) {
			if (verbose) {
				cout << "normaliz.cfg not found. Using default values and command line parameters."<<endl;
			}
		}
		else{
			string buf;
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set test_arihmetic_overflow
				test_arithmetic_overflow=true;
			while (buf!="=")
				setup>>buf;
			setup>>overflow_test_modulus;             //set overflow_test_modulus
			setup>>buf;
			while (buf!="=")
				setup>>buf;
			setup>>lifting_bound;   //set lifting_bound 
			setup>>buf;
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set control data
				verbose=true;
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set save memory / don't optimization for speed
				optimize_speed=false;
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf!="support_hyperplanes" && buf!="triangulation" && buf!="normal" && buf!="triangulation_hilbert_basis" && buf!="hilbert_polynomial"&& buf!="hilbert_basis_polynomial"&& buf!="dual") {
				cerr<<"warning: Unknown \"Run mode type\" in file normaliz.cfg. May be a bad format of the file."<<endl;
				cerr<<"Running \"Run mode type\" = normal ..."<<endl;
			}
			else                     //set computation_type
				computation_type=buf;
			if (computation_type="normal") {
				computation_type="triangulation_hilbert_basis"
			}
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set out flag
				Out.set_write_out(true);
			else
				Out.set_write_out(false);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set inv flag
				Out.set_write_inv(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set ext flag
				Out.set_write_ext(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set esp flag
				Out.set_write_esp(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set typ flag
				Out.set_write_typ(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set egn flag
				Out.set_write_egn(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set gen flag
				Out.set_write_gen(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set sup flag
				Out.set_write_sup(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set tri flag
				Out.set_write_tri(true);
			while (buf!="=")
				setup>>buf;
			setup>>buf;
			if (buf=="YES")           //set ht1 flag
				Out.set_write_ht1(true);
			if (setup.fail()!= false ){
				cerr << "warning: Failed to read file setup.txt. May be a bad format of the file."<<endl;
				cerr << "The program will run in normal mode."<<endl;
			}
			setup.close();
		}
	}

	//Analyzing the command line options

	for (i = 1; i <option.size(); i++) {
		switch (option[i]) {
			case '-':
			case 'i':
				break;
			case 'c':
				verbose=true;
				break;
			case 'f':
				Out.set_write_extra_files();
				break;
			case 'a':
				Out.set_write_all_files();
				break;
			case 's':
				computation_type="support_hyperplanes";
				break;
			case 'S':
				computation_type="support_hyperplanes_pyramid";
				break;
			case 'v':
				computation_type="triangulation";
				break;
			case 'V':
				computation_type="triangulation_pyramid";
				break;
			case 'n':
				computation_type="triangulation_hilbert_basis";
				break;
			case 'N':
				computation_type="hilbert_basis";
				break;
			case '1':
				computation_type="ht1_elements";
				break;
			case 'p':
				computation_type="hilbert_polynomial";
				break;
			case 'h':
				computation_type="hilbert_basis_polynomial";
				break;
			case 'd':
				computation_type="dual";
				break;
			case 'e':  //check for arithmetic overflow
				test_arithmetic_overflow=true;
				break;
			case 'm':  //save memory / don't optimize for speed
				optimize_speed=false;
				break;
			case '?':  //print help text and exit
				printHelp(argv[0]);
				exit(1);
				break;
			default:
				cerr<<"Warning: Unknown option -"<<option[i]<<endl;
				break;
		}
	}


	//if the program works with the indefinite precision arithmetic, no arithmetic tests are performed
#ifdef normbig
	test_arithmetic_overflow=false;
#endif



	//Read Input

	if (!filename_set) {
		cout<<"Normaliz 2.2"<<endl
			<<"Copyright (C) 2007,2008,2009  Winfried Bruns, Bogdan Ichim"<<endl
			<<"With contributions by Christof Soeger"<<endl
			<<"This program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions; See COPYING for details."
			<<endl<<endl;
		cout<<"Enter the input file name or -? for help: ";
		cin >>output_name;
		if (output_name == "-?") {
			printHelp(argv[0]);
			exit(1);
		}
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
			in.open(file_in,ifstream::in);
		}
	} else {
		in2.close();
		in.open(file_in,ifstream::in);
	}
	if (in.is_open()==false) {
		cerr<<"error: Failed to open file "<<name_in<<"."<<endl;
		if (!filename_set) {
			cout<< "Type something and press enter to exit."<<endl;
			cin >> c;
		}
		return 1;
	}
	string mode_string;
	int nr_rows,nr_columns, mode;
	int nr_equations=0;
	Integer number;
	in >> nr_rows;
	in >> nr_columns;
	Matrix M(nr_rows,nr_columns);
	for(i=1; i<=nr_rows; i++){
		for(j=1; j<=nr_columns; j++) {
			in >> number;
			M.write(i,j,number);
		}
	}

	in>>mode_string;
	if (mode_string=="0"||mode_string=="integral_closure") {
		mode=0;
	} else
	if (mode_string=="1"||mode_string=="normalization") {
		mode=1;
	} else
	if (mode_string=="2"||mode_string=="polytope") {
		mode=2;
	} else
	if (mode_string=="3"||mode_string=="rees_algebra") {
		mode=3;
	} else
	if (mode_string=="4"||mode_string=="hyperplanes") {
		mode=4;
	} else
	if (mode_string=="5"||mode_string=="equations") {
		mode=5;
	} else
	if (mode_string=="6"||mode_string=="lattice_ideal") {
		mode=6;
	} else {
		cerr<<"Warning: Unknown mode "<<mode_string<<" and will be replaced with mode integral_closure."<<endl;
		mode=0;
	}

	if (in.fail()!= false ) {
		cerr << "error: Failed to read file "<<name_in<<". May be a bad format of the input file."<<endl;
		if (!filename_set) {
			cout<< "Type something and press enter to exit."<<endl;
			cin >> c;
		}
		return 1;
	}
	if (mode==4 && in.good()) {
		in >> nr_equations;
	}
	in.close();
	Out.set_name(output_name);
	//cout<<"test="<<test_arithmetic_overflow;

	//main computations and output
	if (verbose) {
		cout<<"\n************************************************************\n";
		cout<<"Running in mode "<<mode<<" and computation type "<<computation_type<<"."<<endl;
	}
	make_main_computation(mode, computation_type, M, nr_equations, Out);

	//exit
	if (!filename_set) {
		cout<< "\nProgram finished. Type something and press enter to exit.\n";
		cin >> c;
	}
	return 0;
}

//---------------------------------------------------------------------------
