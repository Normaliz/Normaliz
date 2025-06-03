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

#include <cstdlib>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <csignal>
#include <memory>
using namespace std;

#ifdef NMZ_GPERF
#include <gperftools/profiler.h>
#endif

#include "libnormaliz/integer.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/output.h"
#include "libnormaliz/input.h"
#include "libnormaliz/options.h"
#include "libnormaliz/chunk.h"
#include "libnormaliz/collect_lat.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/project_and_lift.h"
#include "libnormaliz/fusion.h"
#include "libnormaliz/induction.h"

using namespace libnormaliz;

long CCCCCCC = 0;

void printHeader() {
    cout << "                                                     \\.....|" << endl;
#ifdef NMZ_DEVELOP
    cout << "                 Normaliz DEVELOPMENT "
         << "                \\....|" << endl;
#else
    cout << "                    Normaliz " << string(STRINGIFY(NMZ_VERSION) "           ", 11) << "              \\....|"
         << endl;
#endif
    cout << "                                                       \\...|" << endl;
    cout << "     (C) The Normaliz Team, University of Osnabrueck    \\..|" << endl;
    cout << "                      October  2024                      \\.|" << endl;
    cout << "                                                          \\|" << endl;
    string optional_packages = package_string();
    if (optional_packages.size() > 0) {
        cout << "-------------------------------------------------------------" << endl;
        cout << "with package(s)" << optional_packages << endl;
    }
}

void printHelp(char* command) {
    cout << "Usage: " << command << " [options] PROJECT" << endl;
    cout << "  runs normaliz on PROJECT.in" << endl;
    cout << "Computation goals with short options (selection):" << endl;
    cout << "  -s\tcompute support hyperplanes" << endl;
    cout << "  -v\tcompute multiplicity" << endl;
    cout << "  -V\tcompute volume" << endl;
    cout << "  -N\tcompute Hilbert basis (with partial triangulation)" << endl;
    cout << "  -w\tcheck for integrally closed and compute witness if not" << endl;
    cout << "  -q\tcompute Hilbert series" << endl;
    cout << "  -1\tcompute degree 1 elements" << endl;
    cout << "  -T\tcompute triangulation  (output in file .tri)" << endl;
    cout << "  -H\tcompute integer hull" << endl;
    cout << "  -M\tcompute module generators over original monoid" << endl;
    cout << "  -E\tcompute weighted Ehrhart series" << endl;
    cout << "  -L\tcompute virtual multiplicity of weighted Ehrhart series" << endl;
    cout << "  -I\tcompute integral" << endl;
    cout << "  -G\tcheck Gorenstein" << endl;

    cout << endl;
    cout << "Algorithmic variants with short options (selection):" << endl;
    cout << "  -d\t dual ode (includes Hilbert basis, unless combined with -1)" << endl;
    cout << "  -j\t project-and-lift" << endl;
    cout << "  -J\t project-and-lift with floating point arithmetic" << endl;
    cout << "  -k\t keep order" << endl;
    cout << "  -Y\t symmetrization" << endl;
    cout << "  -F\t multiplicity/volume by decent in the face lattice" << endl;

    cout << endl;
    cout << "For computation goals and variants not in the lists above use" << endl;
    // cout << endl;
    cout << "      --<PROP>     compute the ConeProperty <PROP>" << endl;
    // cout << endl;
    cout << "see doc/Normaliz.pdf or doc/NmzShortRef.pdf. Selection:" << endl;
    // cout << endl;
    cout << "      Automorphisms, EuclideanAutomorphisms, RationalA..., CombinatorialA..." << endl;
    cout << "      EhrhartSeries, LatticePoints, NumberLatticePoints" << endl;
    cout << "      FaceLattice, FVector (also with Orbits)" << endl;
    cout << "      DualFaceLattice, DualFVector (also with Orbits)" << endl;
    cout << "      Incidence, DualIncidence" << endl;
    cout << "      MasrkovBasis, GroebnerBasis, Lex, DegLex, RevLex" << endl;
    cout << "      FusionRings, SimpleFusionRings" << endl;

    cout << endl;
    cout << "Output and execution:" << endl;
    cout << "  -f, --files      write the files .out .gen .inv .cst" << endl;
    cout << "  -a, --all-files  write all optional output files" << endl;
    cout << "      --<SUFFIX>   write the file .<SUFFIX> where <SUFFIX> can be one of" << endl;
    cout << "                   cst, egn, esp, ext, gen, ht1, inv, lat, mod, msp, typ" << endl;

    cout << endl;
    cout << "  -B, --BigInt     directly use indefinite precision arithmetic" << endl;
    cout << "      --LongLong   only use long long arithmetic, no conversion possible" << endl;
    cout << "  -i, --ignore     ignore the compute options set in the input file" << endl;
    cout << "  -x=<T>           limit the number of threads to <T>, -x=0 switches the bound off" << endl;
    cout << "  --OutputDir=<path> set a path for the output files (relative to current directory)" << endl;
    cout << "  -?, --help       print this help text and exit" << endl;
    cout << "  -c, --verbose    verbose (prints log data on terminal)" << endl;
    cout << "      --version    print version info and exit" << endl;
    cout << endl;
    cout << "Please report bugs to <normaliz@uos.de> or directly to our issue tracker:" << endl;
    cout << "https://github.com/Normaliz/Normaliz/issues" << endl;
}

int process_data(OptionsHandler& options, const string& command_line);

//---------------------------------------------------------------------------

void set_normaliz_time(){

        GlobalTimeBound = -1.0;
        GlobalPredictionTimeBound = -1.0;
        string name_time = "normaliz.time";
        const char* file_time = name_time.c_str();
        ifstream in_time;
        in_time.open(file_time, ifstream::in);
        if (in_time.is_open()) {
            double test_input;
            in_time >> test_input;
            if(!in_time.fail()){
                GlobalTimeBound = test_input;
                if(verbose)
                    verboseOutput() << "TIME BOUND " << GlobalTimeBound << endl;
            }

            in_time >> test_input;
            if(!in_time.fail()){
                GlobalPredictionTimeBound = test_input;
                if(verbose)
                    verboseOutput() << "TIME PREDICTION BOUND " << GlobalPredictionTimeBound << endl;
            }
            in_time.close();
        }
}

//---------------------------------------------------------------------------

int main(int argc, char* argv[]){
#ifdef NMZ_GPERF
    ProfilerStart("normaliz.prof");
#endif

    running_input_file = true; // used to print output files directly

    verb_length = 0;

    // signal handler for interrupt
    signal(SIGINT, &interrupt_signal_handler);

    vector<string> command_line_items;
    string command_line;
    for (int i = 1; i < argc; ++i){
        command_line = command_line + string(argv[i]) + " ";
        command_line_items.push_back(string(argv[i]));
    }

/*
    for(auto& s: command_line_items){
        if(s == "--PostProcessFusion"){
            try{
                post_process_fusion(command_line_items);
            } catch (const BadInputException& e) {
                cerr << e.what() << endl;
                cerr << "BadInputException caught... exiting." << endl;
                exit(1);
            }
            exit(0);
        }
    }
*/

    string global_command_line = command_line;

    // read command line options
    OptionsHandler global_options;

    vector<string> arg_string_vector = to_string_vector(argc, argv);
    bool print_help = global_options.handle_commandline(arg_string_vector);

    if (print_help) {
        // printHeader();
        printHelp(argv[0]);
        exit(0);
    }

    if (verbose) {
        printHeader();
        verboseOutput() << "-------------------------------------------------------------" << endl;
        verboseOutput() << "Command line: " << command_line << endl;
    }

    /* vector<long long> our_type = {1,1,2,6};
    vector<unsigned int> our_dual = {0,1,2,3};
    vector<long long> our_ring = {0,0,0,1,0,1,1,0,2,5,1}; */

    /* vector<long long> our_type = {1,1,2};
    vector<unsigned int> our_dual = {0,1,2};
    vector<long long> our_ring ={0,0,1,1,1}; */

    /* vector<long long> our_type = {1,1,1,3};
    vector<unsigned int> our_dual = {0,2,1,3};
    vector<long long> our_ring ={0,1,0,0,0,1,2,1}; */

    /* vector<long long> our_type = {1,1,4,4,6};
    vector<unsigned int> our_dual = {0,1,3,2,4};
    vector<long long> our_ring ={0,0,0,1,0,0,1,0,1,1,0,2,1,1,2,3,1}; */

    /* vector<long long> our_type = {1,5,5,5,6,7,7};
    vector<unsigned int> our_dual = {0,1,2,3,4,5,6};
    vector<long long> our_ring ={1,1,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,1,1}; */


    /*

    Induction<long long> Indu(our_type, our_dual, our_ring);

    Indu.start_low_parts();
    Indu.from_low_to_full();

    exit(0);
    */



    /*Matrix<long long> FFT(1,4);
    FFT[0] = {1,1,2,2};
    Cone<long long> TT(Type::fusion_type, FFT);
    // TT.compute(ConeProperty::FusionRings);
    vector<vector<Matrix<long long> > > BB = TT.getFusionDataMatrix();
    BB[0][0].debug_print();
    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " << endl;
    exit(0);*/

    /*Matrix<long long> CC(3,3);
    CC[0] = {2,0,1};
    CC[1] = {0,2,1};
    CC[2] = {0,0,1};
    CC.debug_print('/');
    Cone<long long> LL(Type::vertices, CC);
    LL.compute(ConeProperty::LatticePoints);
    LL.getSupportHyperplanesMatrix().debug_print('$');
    LL.getLatticePointsMatrix().debug_print();
    LL.getModuleGeneratorsMatrix().debug_print('+');*/

    set_normaliz_time();

    vector<string> input_file_names;

    if(list_of_input_files){
        string name_of_file_with_input_file_names = global_options.getProjectName();
        ifstream list_file(name_of_file_with_input_file_names.c_str());
        if(!list_file.is_open())
            throw BadInputException("File with list of input files does not exist");
        while(true){
            string name_of_input_file;
            list_file >> name_of_input_file;
            if(!list_file.good()){
                break;
            }
            input_file_names.push_back(name_of_input_file);
        }
    }

    size_t total_nr_input_files = input_file_names.size();
    if(!list_of_input_files)
        total_nr_input_files = 1;

    size_t nr_input_files_this_instance;

    if(number_normaliz_instances > 0){ // -Z set, number_normaliz_instances = value of -Z
        if(input_file_option < 0) // -A not set
            throw BadInputException("-Z set, but no -A.");
        else{ // both -A and -Z set
            nr_input_files_this_instance = total_nr_input_files / number_normaliz_instances;
            if(input_file_option < total_nr_input_files % number_normaliz_instances)
                nr_input_files_this_instance ++;
        }
    }
    else{   // -Z not set
        if(input_file_option < 0){ // -A not set, we do the full list, otherwise value of -A
            number_normaliz_instances = 1;
            nr_input_files_this_instance = total_nr_input_files;
            input_file_option = 0;
        }
        else{  // -A set, we do a single file
            number_normaliz_instances = total_nr_input_files;
            nr_input_files_this_instance = 1;
        }
    }

    if(GlobalTimeBound != -1.0 && nr_input_files_this_instance > 1){ // all input_files  get the same time
        GlobalTimeBound /= nr_input_files_this_instance;
            if(verbose)
                verboseOutput() << "TIME BOUND PER FILE " << GlobalTimeBound << endl;
    }


    // main loop over input files
    for(size_t intput_file_index = 0; intput_file_index < total_nr_input_files; ++intput_file_index){
        OptionsHandler options;
        if(!list_of_input_files){
            options = global_options;
        }
        else{
            if(intput_file_index % number_normaliz_instances != input_file_option)
                continue;

            vector<string> local_arg_string_vector;
            for (auto& arg_string: arg_string_vector){
                if(arg_string[0] != '-'){
                    arg_string = input_file_names[intput_file_index];
                    local_arg_string_vector.push_back(arg_string);
                    continue;
                }

                if(arg_string.size() >= 2){
                    string test = arg_string.substr(0,2);
                    if(test == "-A" || test == "-Z" || arg_string == "--List")
                        continue;
                }
                local_arg_string_vector.push_back(arg_string);
            }
            print_help = options.handle_commandline(local_arg_string_vector);
            list_of_input_files = true; // must be restored !!
        }

        if(list_of_input_files && verbose)
            verboseOutput() << "*************************************************************" << endl;

        StartGlobalTime();

        process_data(options, command_line);

        if (verbose && GMP_hyp + GMP_scal_prod + GMP_mat > 0)
            verboseOutput() << "GMP transitions: matrices " << GMP_mat << " hyperplanes " << GMP_hyp << " vector operations "
                            << GMP_scal_prod << endl;


        MeasureGlobalTime(verbose);

        if (nmz_interrupted)
        exit(10);

    } // nr_input_files

#ifdef NMZ_GPERF
    ProfilerStop();
#endif

    exit(0);
}

//---------------------------------------------------------------------------

template <typename ConeType, typename InputNumberType>
void compute_and_output(OptionsHandler& options,
                        const InputMap<InputNumberType>& input,
                        const map<NumParam::Param, long>& num_param_input,
                        map<PolyParam::Param, vector<string> >& poly_param_input,
                        renf_class_shared number_field_ref,
                        InputMap<InputNumberType>& add_input) {
    Output<ConeType> Out;  // all the information relevant for output is collected in this object

    // const
    renf_class_shared number_field =
#ifdef ENFNORMALIZ
        number_field_ref.get();
#else
        number_field_ref;
#endif

    options.applyOutputOptions(Out);

    options.activateDefaultMode();

    // Out.set_lattice_ideal_input(input.count(Type::lattice_ideal) > 0);

    Cone<ConeType> MyCone = Cone<ConeType>(input);
    MyCone.setPolyParams(poly_param_input);
    MyCone.setNumericalParams(num_param_input);
    /*MyCone.setNrCoeffQuasiPol(nr_coeff_quasipol);
    MyCone.setExpansionDegree(expansion_degree);
    MyCone.setFaceCodimBound(face_codim_bound);*/
    MyCone.setRenf(number_field);
    MyCone.setProjectName(options.getProjectName());
    try {
        write_fusion_mult_tables_from_input = false;
        if(options.getToCompute().test(ConeProperty::FusionData))
            write_fusion_mult_tables_from_input = true;
        MyCone.compute(options.getToCompute());
        if (add_input.size() > 0) {
            ConeProperties AddInputOptions;
            AddInputOptions.set(ConeProperty::SupportHyperplanes);
            MyCone.modifyCone(add_input);
            MyCone.compute(AddInputOptions);
        }
    } catch (const NotComputableException& e) {
        cout << "Not all desired data could be computed." << endl;
        cout << e.what() << endl;
        cout << "Writing only available data." << endl;
    } catch (const InterruptException& e) {
        cout << endl;
        cout << "Computation was interrupted." << endl;
        cout << e.what() << endl;
        if(!output_on_interrupt){
            cout << "No output on inmterrupt" << endl;
            exit(10);
        }
        cout << "Output on interrupt activated. Writing available data." << endl;
    }

    if(is_split_patching){
        cout << "No file <project>.out for split computation" << endl;
        // MeasureGlobalTime(verbose);
        return;
    }

    Out.setCone(MyCone);
    Out.set_renf(number_field);

    signal(SIGINT, SIG_DFL);

    // Output may call extra computations. It does so for the Hilbert quasipolynomial
    // In order to throw the interrupt exception again, we disable it here.
    nmz_interrupted = 0;
    Out.write_files();

    if (MyCone.isComputed(ConeProperty::IntegerHull)) {
        Output<ConeType> IntHullOut;
        options.applyOutputOptions(IntHullOut);
        IntHullOut.set_name(options.getProjectName() + ".IntHull");
        IntHullOut.setCone(MyCone.getIntegerHullCone());
        IntHullOut.set_renf(number_field, true);
        IntHullOut.write_files();
    }

    if (MyCone.isComputed(ConeProperty::ProjectCone)) {
        Output<ConeType> ProjOut;
        options.applyOutputOptions(ProjOut);
        ProjOut.set_name(options.getProjectName() + ".ProjectCone");
        ProjOut.setCone(MyCone.getProjectCone());
        ProjOut.set_renf(number_field);
        ProjOut.write_files();
    }

#ifdef NMZ_COCOA
    if (MyCone.isComputed(ConeProperty::Symmetrize)) {
        Output<ConeType> SymmOut;
        options.applyOutputOptions(SymmOut);
        SymmOut.set_name(options.getProjectName() + ".symm");
        SymmOut.setCone(MyCone.getSymmetrizedCone());
        SymmOut.write_files();
    }
#endif
}

//---------------------------------------------------------------------------
// for testing only, not really useful in Normaliz
template <typename InputNumberType>
InputMap<InputNumberType> extract_additional_input(
    InputMap<InputNumberType>& input) {
    InputMap<InputNumberType> add_input;
    size_t nr_add_input = 0;
    auto M = input.find(Type::add_inequalities);
    if (M != input.end()) {
        add_input[Type::inequalities] = input[Type::add_inequalities];
        input.erase(Type::add_inequalities);
        nr_add_input++;
    }
    M = input.find(Type::add_equations);
    if (M != input.end()) {
        add_input[Type::equations] = input[Type::add_equations];
        input.erase(Type::add_equations);
        nr_add_input++;
    }
    M = input.find(Type::add_inhom_inequalities);
    if (M != input.end()) {
        add_input[Type::inhom_inequalities] = input[Type::add_inhom_inequalities];
        input.erase(Type::add_inhom_inequalities);
        nr_add_input++;
    }
    M = input.find(Type::add_inhom_equations);
    if (M != input.end()) {
        add_input[Type::inhom_equations] = input[Type::add_inhom_equations];
        input.erase(Type::add_inhom_equations);
        nr_add_input++;
    }
    M = input.find(Type::add_cone);
    if (M != input.end()) {
        add_input[Type::cone] = input[Type::add_cone];
        input.erase(Type::add_cone);
        nr_add_input++;
    }
    M = input.find(Type::add_subspace);
    if (M != input.end()) {
        add_input[Type::subspace] = input[Type::add_subspace];
        input.erase(Type::add_subspace);
        nr_add_input++;
    }
    M = input.find(Type::add_vertices);
    if (M != input.end()) {
        add_input[Type::vertices] = input[Type::add_vertices];
        input.erase(Type::add_vertices);
        nr_add_input++;
    }
    return add_input;
}

//---------------------------------------------------------------------------

int process_data(OptionsHandler& options, const string& command_line) {
    try {
        if (options.getProjectName() == "" && !options.isUseChunk()) {
            cerr << "ERROR: No project name set!" << endl;
            exit(1);
        }

        global_project = options.getProjectName();

        Check_Stop(); // check whether stop file has been set

        if(options.isUseCollectLat()){
            // cout << "PPPPPPPPPPPPPPPP " << split_index_option << endl;
            write_fusion_mult_tables_from_input = false;
            if(options.getToCompute().test(ConeProperty::FusionData))
                write_fusion_mult_tables_from_input = true;
            collect_lat(global_project, split_index_option);
            return 0;
        }

        if(options.isUseSaveLocalSolutions()){
            save_local_solutions = true;
            if(level_local_solutions == -1)
                throw BadInputException("SaveLocalSolutions requires level set by -Q");
        }

        /* if(options.isUseNextRound()){
            next_round(global_project);
            return 0;
        }*/

        if (options.isUseChunk()) {
            chunk();
            exit(0);
        }

        if (options.isUseAddChunks()) {
            add_chunks(options.getProjectName());
            exit(0);
        }

        if(options.isUseSplit()){
            is_split_patching = true;
        }

        bool standard_fusion_name = false;
        bool only_partition;
        // reset_global_fusion_data(); // because of list processing
        FusionBasic test_fusion;
        string name_in = options.getProjectName() + ".in";
        const char* file_in = name_in.c_str();
        ifstream in;
        in.open(file_in, ifstream::in);
        if (!in.is_open()) {
            if(options.get_given_name_contains_in() || is_split_patching){ // patching requires real input file
                cerr << "error: Failed to open file " << name_in << "." << endl;
                return 1;
            }
            pair<bool, bool> result = test_fusion.data_from_string(global_project, true);
            standard_fusion_name = result.first;
            only_partition = result.second;
            if(no_empty_output && standard_fusion_name && only_partition){
                set<unsigned long> test_dupl;
                for(auto& t: test_fusion.fusion_type)
                    test_dupl.insert(t);
                if(test_dupl.size() == test_fusion.fusion_type.size()){
                    only_partition = false;
                    test_fusion.duality = identity_key(test_fusion.fusion_type.size());
                }
            }

            if(!standard_fusion_name){
                cerr << "error: Failed to open file " << name_in << "." << endl;
                return 1;
            }
        }
        else{
            if(options.isUseMakeFullInput())
                throw BadInputException("MkaeFusionInput not allowed if input file exists");
        }

        string polynomial;
        vector<string> polynomial_equations;

        InputMap<mpq_class> input, add_input;
        InputMap<renf_elem_class> renf_input, renf_add_input;
        map<NumParam::Param, long> num_param_input;
        map<PolyParam::Param, vector<string> > poly_param_input;
        bool renf_read = false;
        renf_class_shared number_field;

        if(!standard_fusion_name){
            try {
                input = readNormalizInput<mpq_class>(in, options, num_param_input, poly_param_input,  number_field);
                if (nmz_interrupted)
                    exit(10);
            }
#ifdef ENFNORMALIZ
            catch (const NumberFieldInputException& e) {
                if (verbose)
                    verboseOutput() << "Input specifies a number field, trying again with number field implementation..." << endl;

                in.close();
                in.open(file_in, ifstream::in);
                renf_input = readNormalizInput<renf_elem_class>(in, options, num_param_input, poly_param_input, number_field);
                if (nmz_interrupted)
                    exit(10);
                renf_read = true;
            }
#else
            catch (const NumberFieldInputException& e) {
                throw BadInputException("");
            }
#endif
            in.close();
        }
        else{
            if(!only_partition)
                make_input_from_fusion_data(test_fusion, input, options.isUseMakeFullInput());
            else
                make_partition_input_from_fusion_data(test_fusion, input, options.isUseMakeFullInput());

            if(options.isUseMakeFullInput())
                return 0;
        }

        if (verbose) {
            // verboseOutput() << "-------------------------------------------------------------" << endl;
            if(list_of_input_files)
                verboseOutput() << "Input file: " << options.getProjectName() << endl;
            verboseOutput() << "Compute: ";
            if (options.getToCompute().none())
                verboseOutput() << "No computation goal/variant  set, using defaults for given input" << endl;
            else
                verboseOutput() << options.getToCompute() << endl;
            for(auto& P: num_param_input)
                verboseOutput() << numpar_to_string(P.first) << " = " << P.second << endl;
        }

        if (renf_read) {
            if (options.isUseLongLong())
                throw BadInputException("LongLong not allowed for algebraic polyhedra");
            // if(options.getToCompute().test(ConeProperty::Dynamic))
            renf_add_input = extract_additional_input<renf_elem_class>(renf_input);

            compute_and_output<renf_elem_class>(options, renf_input, num_param_input, poly_param_input, number_field, renf_add_input);
        }
        else {
            if (options.isUseLongLong()) {
                // if(options.getToCompute().test(ConeProperty::Dynamic))
                add_input = extract_additional_input<mpq_class>(input);
                compute_and_output<long long>(options, input, num_param_input, poly_param_input, number_field, add_input);
            }
            else {
                // if(options.getToCompute().test(ConeProperty::Dynamic))
                add_input = extract_additional_input<mpq_class>(input);
                compute_and_output<mpz_class>(options, input, num_param_input, poly_param_input, number_field, add_input);
            }
        }

    } catch (const BadInputException& e) {
        cerr << e.what() << endl;
        cerr << "BadInputException caught... exiting." << endl;
        if(!list_of_input_files)
            exit(1);
    } catch (const FatalException& e) {
        cerr << e.what() << endl;
        cerr << "FatalException caught... exiting." << endl;
        if(!list_of_input_files)
            exit(2);
    } catch (const NoComputationException& e) {
        cerr << e.what() << endl;
        cerr << "NoComputationException caught... exiting." << endl;
        if(!list_of_input_files)
            exit(7);
    } catch (const TimeBoundException& e) {
        cerr << e.what() << endl;
        cerr << "Time bound exceeded for " << global_project << " exiting." << endl;
        if(!is_split_patching && !no_empty_output){
            cerr << "Creating signal file with suffix exc" << endl;
            string exc_name = options.getProjectName() + ".exc";
            ofstream exc_file(exc_name.c_str());
        }
        if(!list_of_input_files)
            exit(6);
    } catch (const NmzCoCoAException& e) {
        cerr << e.what() << endl;
        cerr << "NmzCoCoAException caught... exiting." << endl;
        if(!list_of_input_files)
            exit(3);
    } catch (const NormalizException& e) {
        cerr << e.what() << endl;
        cerr << "NormalizException caught... exiting." << endl;
        if(!list_of_input_files)
            exit(4);
    } catch (const std::exception& e) {
        cerr << "std::exception caught... \"" << e.what() << "\" ...  exiting." << endl;
        if(!list_of_input_files)
            exit(5);
    }

    return 0;
}

    /*
    verbose = true;
    vector<mpq_class> ty = {1,1,2,2};
    vector<mpq_class> du = {-1,1,2,3};
    Matrix<mpq_class> ty_mat(0,4);
    ty_mat.append(ty);
    Matrix<mpq_class> du_mat(0,4);
    du_mat.append(du);
    Cone<long long> CT(Type::fusion_type, ty_mat, Type::fusion_duality, du_mat);
    // CT.compute(ConeProperty::FusionRings, ConeProperty::UseModularGrading);
    Matrix<long long> Res = CT.getFusionRingsMatrix();
    Res.debug_print();
    CT.compute(ConeProperty::FusionRings, ConeProperty::UseModularGrading);
    Res = CT.getFusionRingsMatrix();
    Res.debug_print();
    Res = CT.getSingleFusionRing();
    exit(0);
    Res.debug_print('%');
    CT.compute(ConeProperty::FusionRings, ConeProperty::UseModularGrading);
    Res = CT.getFusionRingsMatrix();
    Res.debug_print('@');
    CT.setModularGraing(0);
    cout << "-------------------------------------------------------------------------------------------------------" << endl;
    Res = CT.getFusionRingsMatrix();
    Res.debug_print();
    exit(0); */

