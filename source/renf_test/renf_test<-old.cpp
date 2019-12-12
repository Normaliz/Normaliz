#include <sstream>
#include <cstdlib>
#include "e-antic/renfxx.h"

#include <gperftools/profiler.h>

using namespace std;

int main(void) {
    
    ProfilerStart("renf_time_new.prof");

    // istringstream is("min_poly (a2-2) embedding 1.4+/-0.1");
    istringstream is("min_poly (a^12 + a^6+a^5+a^2- 5) embedding [1.0 +/- 0.1]");
    renf_class NF;
    is >> NF;
    cout << NF << endl;

    renf_elem_class elem(NF.get_renf());    
    renf_elem_class squ(NF.get_renf());
    
    stringstream inout;
    inout >> set_renf(NF.get_renf());
    // inout << "( 37*a - 4)";
    inout << "( a^10-3a^9+a^6-a^5+a-22)";
    inout >> elem;
    
    for(long i=0;i< 100000000;++i){
        squ=elem*elem;        
    }
    
    cout << elem << endl;
    
    ProfilerStop();
}
