#include <cstdlib>
#include <vector>
#include <fstream>
#include <omp.h>
#include <iomanip>

using namespace std;

#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/full_cone.h"

using namespace libnormaliz;

void spoiled_basic_data(){
    cout << "File basic.data is spoiled" << endl;
    exit(1);
}

void set_parallelization(long thread_limit) {

    if (thread_limit <= 0){
        throw BadInputException("Invalid thread limit");
    }
    
    omp_set_max_active_levels(1);
    if (std::getenv("OMP_NUM_THREADS") == NULL) {
        omp_set_num_threads(thread_limit);
        cout << "SET SET SET" << endl;
    }
}

int main(int argc, char* argv[]) {
    
    if(argc < 2){
        cout << "Not enough parameters" << endl;
        exit(1);
    }    
    
    string type;
    
    size_t this_chunk;

    cin >> std::ws;
    int c = cin.peek();
    if (c == 'B') {
        cin >> type;
        if(type != "Block"){
            cout << "Hollow tri file spoiled" << endl;
            exit(1);
        }
        cin >> this_chunk;
    }
    else
        this_chunk = stoi(string(argv[1]));

    set_parallelization(8);    
    if(argc>2){
        long thread_limit = stoi(string(argv[2]));
        set_parallelization(thread_limit);
    }        

    string name_in = "basic.data";
    const char* file_in = name_in.c_str();    
    ifstream in;
    in.open(file_in, ifstream::in);
    if (in.is_open() == false){
        cout << "Cannot find basic.data" << endl;;
        exit(1);        
    }
    
    in >> type;
    if(type != "Dim")
        spoiled_basic_data();
    size_t dim;
    in >> dim;

    in >> type;
    if(type != "Gen")
        spoiled_basic_data();
    size_t nr_gen;
    in >> nr_gen;

    Matrix<mpz_class> Generators(nr_gen,dim);
    for(size_t i=0; i< nr_gen; ++i){
        for(size_t j=0; j< dim; ++j)
            in >> Generators[i][j];        
    }
    cout << "Generators" << endl;
    Generators.pretty_print(cout);
    
    in >> type;
    if(type != "Grad")
        spoiled_basic_data();
    vector<mpz_class> GradingOnPrimal(dim);
    for(size_t j=0; j< dim; ++j)
        in >> GradingOnPrimal[j];
    cout << "GradingOnPrimal" << endl;
    cout << GradingOnPrimal;

    in >> type;
    if(type != "Generic")
        spoiled_basic_data();
    vector<mpz_class> Generic(dim);
    for(size_t j=0; j< dim; ++j)
        in >> Generic[j];
    
    cout << "Generic" << endl;
    cout << Generic;

    in >> type;
    if(type != "Blocks")
        spoiled_basic_data();
    size_t nr_blocks,dummy;
    in >> nr_blocks;
    cout << "Blocks " << nr_blocks << endl;
    vector<size_t> block_start(nr_blocks), block_end(nr_blocks);
    for(size_t i = 0; i< nr_blocks; ++i){
        in >> dummy;
        if(dummy != i)
            spoiled_basic_data();
        in >> block_start[i] >> block_end[i];
        cout << i << "  " << block_start[i] << "  " << block_end[i] << endl;
    }
    if(this_chunk >= nr_blocks)
        spoiled_basic_data();
    cout << "This chunk " << this_chunk << endl;
    
    vector<pair<dynamic_bitset,dynamic_bitset> > Triangulation_ind(block_end[this_chunk] - block_start[this_chunk]);
 
    string input_string;
    for(size_t i= 0; i< Triangulation_ind.size(); ++i){
        cin >> input_string;
        Triangulation_ind[i].first.resize(nr_gen);
        for(size_t j=0; j < input_string.size(); ++j){
            if(input_string[j]=='1')
                Triangulation_ind[i].first[nr_gen-1-j] = 1;
        }
        cin >> input_string;
        Triangulation_ind[i].second.resize(nr_gen);
        for(size_t j=0; j < input_string.size(); ++j){
            if(input_string[j]=='1')
                Triangulation_ind[i].second[nr_gen-1-j] = 1;
        }
    }
    
    // cout << Triangulation_ind.back().first << endl;
    
    cin >> type;
    if(type != "End")
        spoiled_basic_data();
    
    int     omp_start_level = omp_get_level();
    
    SignedDec<mpz_class> SDMult(Triangulation_ind, Generators, GradingOnPrimal, omp_start_level);
    SDMult.verbose = true;
    SDMult.Generic = Generic;
    if(!SDMult.ComputeMultiplicity())
            assert(false);

    mpq_class multiplicity = SDMult.multiplicity;
    
    mpz_class corr_factor = v_gcd(GradingOnPrimal); // search in code for corr_factor to find an explanation
    multiplicity *= corr_factor; 
    
    cout << "Multiplicity" << endl;
    cout << multiplicity << endl;
    cout << "Mult (float) " << std::setprecision(12) << mpq_to_nmz_float(multiplicity) << endl;
    
    string file_name = "mult.";
    file_name += to_string(this_chunk);
    ofstream out(file_name.c_str());
    out << "multiplicity " << this_chunk << endl << endl;
    out << multiplicity<< endl;
    out.close();
    
    
    
}
