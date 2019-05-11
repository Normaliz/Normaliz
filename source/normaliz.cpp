/*
* Normaliz
* Copyright (C) 2007-2019  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
#include <fstream>
#include <algorithm>
#include <csignal>
using namespace std;

#ifdef NMZ_GPERF
#include <gperftools/profiler.h>
#endif

#include "libnormaliz/integer.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/output.h"

using namespace libnormaliz;
#include "input.cpp"
#include "options.cpp"

#ifndef STRINGIFY
#define STRINGIFYx(Token) #Token
#define STRINGIFY(Token) STRINGIFYx(Token)
#endif

long CCCCCCC=0;

void printCopying() {
    cout<<"Copyright (C) 2007-2018  The Normaliz Team, University of Osnabrueck."<<endl
        <<"This program comes with ABSOLUTELY NO WARRANTY; This is free software,"<<endl
        <<"and you are welcome to redistribute it under certain conditions;"<<endl
        <<"See COPYING for details."<<endl;
}


void printVersion() {
    cout << "Normaliz " << string(STRINGIFY(NMZ_VERSION)) << endl;
    printCopying();
}

void printHeader() {
    cout << "                                                    \\.....|"<<endl;
    cout << "                    Normaliz " << string( STRINGIFY(NMZ_VERSION) "           " ,11)
                                                 << "             \\....|"<<endl;
    cout << "                                                      \\...|"<<endl;
    cout << "     (C) The Normaliz Team, University of Osnabrueck   \\..|"<<endl;
    cout << "                     May  2019                          \\.|"<<endl;
    cout << "                                                         \\|"<<endl;
    bool with_optional_packages=false;
    string optional_packages;
#ifdef NMZ_COCOA
    with_optional_packages=true;
    optional_packages+=" CoCoALib";
#endif
#ifdef NMZ_FLINT
#ifndef ENFNORMALIZ
    with_optional_packages=true;
    optional_packages+=" Flint";
#endif
#endif
#ifdef NMZ_SCIP
    with_optional_packages=true;
    optional_packages+=" SCIP";
#endif
#ifdef ENFNORMALIZ
    with_optional_packages=true;
    optional_packages+=" Flint antic arb e-antic";
#endif
    
    if(with_optional_packages){    
        cout << "------------------------------------------------------------" << endl;
        cout << "with package(s)" << optional_packages << endl; 
    }
}

void printHelp(char* command) {
    cout << "Usage: "<<command<<" [options] PROJECT"<<endl;
    cout << "  runs normaliz on PROJECT.in"<<endl;
    cout << "Options:"<<endl;
    cout << "  -S\tcompute sublattice"<<endl;
    cout << "  -s\tcompute support hyperplanes"<<endl;
    cout << "  -t\tcompute triangulation"<<endl;
    cout << "  -v\tcompute multiplicity"<<endl;
    cout << "  -v\tcompute volume"<<endl;
    cout << "  -F\tcompute volume by decent in the face lattice"<<endl;
    cout << "  -n\tcompute Hilbert basis and multiplicity (needs full triangulation)"<<endl;
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
    cout << "  -E\tcompute weighted Ehrhart series"<<endl;
    cout << "  -L\tcompute virtual multiplicity of weighted Ehrhart series"<<endl;
    cout << "  -I\tcompute integral"<<endl;
    cout << "  -G\tcheck Gorenstein"<<endl;

    cout << endl;
    cout << "  -d\tcomputation mode: dual"<<endl;
    cout << "  -P\tcomputation mode: primal"<<endl;
    cout << "  -j\tcomputation mode: project-and-lift"<<endl;
    cout << "  -J\tcomputation mode: project-and-lift with floating point arithmetic"<<endl;
    cout << "  -r\tcomputation mode: approximate"<<endl;
    cout << "  -b\tcomputation mode: bottom decomposition"<<endl;
    cout << "  -o\tcomputation mode: no bottom decomposition"<<endl;
    cout << "  -k\tcomputation mode: keep order"<<endl;
    cout << "  -Y\tcomputation mode: symmetrization"<<endl;

    cout << endl;
    cout << "      --<PROP>     compute the ConeProperty <PROP>"<<endl;

    cout << endl;
    cout << "  -f, --files      write the files .out .gen .inv .cst"<<endl;
    cout << "  -a, --all-files  write all output files (except  .dec .tri .typ)"<<endl;
    cout << "      --<SUFFIX>   write the file .<SUFFIX> where <SUFFIX> can be one of"<<endl;
    cout << "                   cst, egn, esp, ext, gen, ht1, inv, lat, mod, msp, typ"<<endl;

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


int process_data(OptionsHandler& options, const string& command_line, renf_class& number_field);

//---------------------------------------------------------------------------

int main(int argc, char* argv[])
{

#ifdef NMZ_GPERF
    ProfilerStart("normaliz.prof");
#endif
    
    /*cout << "Before AAA" << endl;
    renf_class K("a^2 - 5", "a", "2.0 +/- 1.0");
    cout << "Constructed" << endl;
   renf_elem_class b(K, "a+1");
   // std::cout << "number " << b << std::endl;
   cout << "Survived " << endl;*/
    
    // signal handler for interrupt
    signal(SIGINT, &interrupt_signal_handler);
    
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
    
    renf_class number_field; // is bool without e-antic
    
    process_data(options, command_line, number_field);
    
    if(nmz_interrupted)
        exit(10);
    
#ifdef NMZ_GPERF
    ProfilerStop();
#endif
    
    exit(0);
}

//---------------------------------------------------------------------------

template<typename ConeType, typename InputNumberType>
void compute_and_output(OptionsHandler& options, const map <Type::InputType, vector< vector<InputNumberType> > >& input, 
                        const map<NumParam::Param,long >& num_param_input, const string& polynomial, renf_class& number_field,
                        const map <Type::InputType, vector< vector<InputNumberType> > >& add_input)
{
    
    Output<ConeType> Out;    //all the information relevant for output is collected in this object

    options.applyOutputOptions(Out);
    
    options.activateDefaultMode();

    Out.set_lattice_ideal_input(input.count(Type::lattice_ideal)>0);

    Cone<ConeType> MyCone = Cone<ConeType>(input);
    MyCone.setPolynomial(polynomial);
    MyCone.setNumericalParams(num_param_input);
    /*MyCone.setNrCoeffQuasiPol(nr_coeff_quasipol);
    MyCone.setExpansionDegree(expansion_degree);
    MyCone.setFaceCodimBound(face_codim_bound);*/
    MyCone.setRenf(&number_field);
    MyCone.set_project(options.getProjectName());
    MyCone.set_output_dir(options.getOutputDir());
    try {
        MyCone.compute(options.getToCompute());
        if(add_input.size()>0){
            ConeProperties AddInputOptions;
            AddInputOptions.set(ConeProperty::SupportHyperplanes);
            MyCone.addInput(add_input);
            MyCone.compute(AddInputOptions);
        }
    } catch(const NotComputableException& e) {
        std::cout << "Not all desired properties could be computed." << endl;
        std::cout << e.what() << endl;
        std::cout << "Writing only available data." << endl;
    } catch(const InterruptException& e) {
        std::cout << endl;
        std::cout << "Computation was interrupted." << endl;
        std::cout << "Writing only available data." << endl;
    }
    Out.setCone(MyCone);
    Out.set_renf(&number_field);
    
    signal(SIGINT, SIG_DFL);
    
    Out.write_files();
    
    if(MyCone.isComputed(ConeProperty::IntegerHull)){
        Output<ConeType> IntHullOut;
        options.applyOutputOptions(IntHullOut);
        IntHullOut.set_name(options.getProjectName()+".IntHull");
        IntHullOut.setCone(MyCone.getIntegerHullCone());
        IntHullOut.set_renf(&number_field,true);
        IntHullOut.write_files();        
    }
    
    if(MyCone.isComputed(ConeProperty::ProjectCone)){
        Output<ConeType> ProjOut;
        options.applyOutputOptions(ProjOut);
        ProjOut.set_name(options.getProjectName()+".ProjectCone");
        ProjOut.setCone(MyCone.getProjectCone());
        ProjOut.set_renf(&number_field);
        ProjOut.write_files();        
    }

#ifdef NMZ_COCOA
    if(MyCone.isComputed(ConeProperty::Symmetrize)){
        Output<ConeType> SymmOut;
        options.applyOutputOptions(SymmOut);
        SymmOut.set_name(options.getProjectName()+".symm");
        SymmOut.setCone(MyCone.getSymmetrizedCone());
        SymmOut.write_files();        
    }
#endif    
}

//---------------------------------------------------------------------------

template<typename InputNumberType>    
map <Type::InputType, vector< vector<InputNumberType> > > 
    extract_additional_input(map <Type::InputType, vector< vector<InputNumberType> > >& input){
        map <Type::InputType, vector< vector<InputNumberType> > > add_input;
        size_t nr_add_input=0;
        auto M=input.find(Type::add_inequalities);
        if(M!=input.end()){
            add_input[Type::inequalities]=input[Type::add_inequalities];
            input.erase(Type::add_inequalities);
            nr_add_input++;
        }
        M=input.find(Type::add_inhom_inequalities);
        if(M!=input.end()){
            add_input[Type::inhom_inequalities]=input[Type::add_inhom_inequalities];
            input.erase(Type::add_inhom_inequalities);
            nr_add_input++;
        }
        M=input.find(Type::add_cone);
        if(M!=input.end()){
            add_input[Type::cone]=input[Type::add_cone];
            input.erase(Type::add_inequalities);
            nr_add_input++;
        }
        M=input.find(Type::add_vertices);
        if(M!=input.end()){
            add_input[Type::vertices]=input[Type::add_vertices];
            input.erase(Type::add_inequalities);
            nr_add_input++;
        }
        if(nr_add_input>1){
            throw BadInputException("Only one additional input tyoe allowed");
        }
        return add_input;        
}
//---------------------------------------------------------------------------

int process_data(OptionsHandler& options, const string& command_line, renf_class& number_field) {
    
    try {

    if(options.getProjectName()==""){
        cerr << "ERROR: No project name set!" << endl;
        exit(1);
    }

    string name_in=options.getProjectName()+".in";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if ( !in.is_open() ) {
        cerr << "error: Failed to open file "<<name_in<<"."<<endl;
        exit(1);
    }

    string polynomial="";  // these are default values
    /*long nr_coeff_quasipol=-1;
    long expansion_degree=-1;
    long face_codim_bound=-1;*/
    
    map <Type::InputType, vector< vector<mpq_class> > > input, add_input;
    map <Type::InputType, vector< vector<renf_elem_class> > > renf_input, renf_add_input;
    map <NumParam::Param, long > num_param_input;
    bool renf_read=false;
    
    try{
        
    input = readNormalizInput<mpq_class>(in, options,num_param_input,polynomial,number_field);
    if(nmz_interrupted)
        exit(10);
    }
#ifdef ENFNORMALIZ
    catch (const NumberFieldInputException& e) {
        if(verbose)
            verboseOutput() << "Input specifies a number field, trying again with number field implementation..." << endl;
        
        in.close();
        in.open(file_in,ifstream::in);
        renf_input = readNormalizInput<renf_elem_class>(in, options,num_param_input,polynomial,number_field);
        if(nmz_interrupted)
            exit(10);
        renf_read=true;
    }
#else
    catch (const NumberFieldInputException& e) {
        throw BadInputException("");
    }
#endif

    in.close();

    if (verbose) {
        cout << "------------------------------------------------------------" << endl;
        cout << "Command line: " << command_line << endl;
        cout << "Compute: " << options.getToCompute() << endl;
    }
    
    if(renf_read){
        if(options.isUseLongLong())
            throw BadInputException("LongLong not allowed for algebraic polyhedra");
        if(options.getToCompute().test(ConeProperty::Dynamic))
            renf_add_input=extract_additional_input<renf_elem_class>(renf_input);
            
        compute_and_output<renf_elem_class>(options, renf_input, num_param_input, polynomial, number_field, renf_add_input);
    }
    else{
        if(options.isUseLongLong()){
            if(options.getToCompute().test(ConeProperty::Dynamic))
                add_input=extract_additional_input<mpq_class>(input);
            compute_and_output<long long>(options, input, num_param_input,polynomial, number_field, add_input);
        }
        else{
            if(options.getToCompute().test(ConeProperty::Dynamic))
                add_input=extract_additional_input<mpq_class>(input);
            compute_and_output<mpz_class>(options, input, num_param_input,polynomial, number_field, add_input);
        }
    }

    } catch(const BadInputException& e) {
        cerr << e.what() << endl;
        cerr << "BadInputException caught... exiting." << endl;
        exit(1);
    } catch(const FatalException& e) {
        cerr << e.what() << endl;
        cerr << "FatalException caught... exiting." << endl;
        exit(2);
    } catch(const NmzCoCoAException& e) {
        cerr << e.what() << endl;
        cerr << "NmzCoCoAException caught... exiting." << endl;
        exit(3);
    } catch(const NormalizException& e) {
        cerr << e.what() << endl;
        cerr << "NormalizException caught... exiting." << endl;
        exit(4);
    } catch(const std::exception& e) {
        cerr << "std::exception caught... \""<< e.what()<<"\" ...  exiting." << endl;
        exit(5);
    }

    return 0;
}
