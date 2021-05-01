#include <cstdlib>
#include <vector>
#include <fstream>
#include <omp.h>
#include <iomanip>

#include "libnormaliz/libnormaliz.h"

using namespace std;

using namespace libnormaliz;

int main(int argc, char* argv[]) {

    if(argc < 3){
        cout << "Not enaough parameters" << endl;
        exit(1);        
    }
    
    string project_name(argv[1]);
    
    mpq_class total_mult = 0;
    
    size_t nr_blocks = stoi(string(argv[2]));
    cout << "Summing " << nr_blocks << " partial multiplicities" << endl;
    for(size_t i = 0; i< nr_blocks; ++i){
        cout << "Reading block " << i << endl;
        string name_in = project_name+".mult." + to_string(i);
        const char* file_in = name_in.c_str();    
        ifstream in;
        in.open(file_in, ifstream::in);
        string type;
        in >> type;
        if(type != "multiplicity"){
            cout << "spoiled mult " << i << endl;
            exit(1);
        }
        size_t this_chunk;
        in >> this_chunk;
        if(i != this_chunk){
            cout << "spoiled mult " << i << endl;
            exit(1);
        }
        mpq_class mult;
        in >> mult;
        total_mult += mult;        
    }
    cout << "Toatl miultiplicity" << endl;
    cout << total_mult << endl;
    cout << "Toatl miultiplicity (float) " << std::setprecision(12) << mpq_to_nmz_float(total_mult) << endl;
    
}
