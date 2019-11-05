#include <sstream>
#include <cstdlib>
#include "e-antic/renfxx.h"

#include <gperftools/profiler.h>

using namespace std;

int main(void) {
    
    ProfilerStart("renf_time_old.prof");
    
    // string mp_string="a^2-2";
    // string emb_string="1.4+/-0.1";
    
    string mp_string="a^12 + a^6+a^5+a^2- 5";    
    string emb_string="2 +/- 1";
    
    renf_class Renf(mp_string, "a", emb_string);;

    char *res, *res1;
    res = fmpq_poly_get_str_pretty(Renf.get_renf()->nf->pol, "a");
    res1 = arb_get_str(Renf.get_renf()->emb, 64, 0);
    cout << "min_poly "
        << "(" << res << ")"
        << " embedding " << res1 << endl
        << endl;
    flint_free(res);
    flint_free(res1);

    
    stringstream inout;
    Renf.set_istream(inout);
    
    // inout << "( 37*a - 4)";
    
    inout << "( a^10-3a^9+a^6-a^5+a-22)";
    
    renf_elem_class elem, squ;
    inout >> elem;
    
    for(long i=0;i< 100000000;++i){
        squ=elem*elem;
        // cout << i << endl;
        
    }
    
    cout << elem << endl;
    
    ProfilerStop();
}
