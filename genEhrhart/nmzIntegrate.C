

#include "CoCoA/library.H"
using namespace CoCoA;
using namespace std;


#include <fstream> 
// #include <sstream>                

#include <omp.h>

#include "../libnormaliz/HilbertSeries.h"
// #include "../libnormaliz/general.h"
#include "../libnormaliz/matrix.h"
// #include "../libnormaliz/vector_operations.h"
#include "../libnormaliz/map_operations.h"
#include "../output.h"
 
    
#include "nmzIntInput.C"

#include "nmzIntPoly.C"

#include "cyclRatFunct.C"

#include "genEhrhart.C"

#include "nmzIntegral.C"                 


void printHelp(char* command) {
    cout << "usage: "<<command<<" [-cEIL?] [-x=<T>] [PROJECT]"<<endl;
    cout << "  runs nmzIntegrate on PROJECT.in"<<endl;
    cout << "options:"<<endl;
    cout << "  -?\tprint this help text and exit"<<endl;
    cout << "  -E\tcompute generalized Ehrhart series"<<endl;
    cout << "  -I\tcompute integral"<<endl;
    cout << "  -L\tcompute lead coefficient of quasipolynomial"<<endl;
    cout << "  -c\tverbose (prints control data)"<<endl;
    cout << "  -x=<T>\tlimit the number of threads to <T>"<<endl;
}

bool verbose_INT;

//----------------------------------------------------------------------
// Use main() to analyze options, start computations and 
// handle any uncaught exceptions.
int main(int argc, char* argv[])  
{
  try
  {
    size_t i;   
    string project;            

    bool filename_set=false;
    string option;            //all options concatenated (including -)
    for (i = 1; i < (unsigned int)argc; i++) {
        if (argv[i][0]=='-') {
            if (argv[i][1]!='\0') {
                if (argv[i][1]!='x') {
                    option = option + argv[i];
                } else if (argv[i][2]=='=') {
                    string Threads = argv[i];
                    Threads.erase(0,3);
                    size_t nr_threads;
                    if ( (istringstream(Threads) >> nr_threads) && nr_threads > 0) {
                        omp_set_num_threads(nr_threads);
                    } else {
                        cerr<<"Warning: Invalid option string "<<argv[i]<<endl;
                    }
                } else {
                    cerr<<"Warning: Invalid option string "<<argv[i]<<endl;
                }
            }
        } else if (!filename_set) {
            string s(argv[i]);
            project=s;
            filename_set=true;
        }
    }    
    

    bool do_genEhrhart=false, do_integral=false, do_leadCoeff=false;
    
    for (i = 1; i <option.size(); i++) {
        switch (option[i]) {
            case '-':
            case 'c':
                verbose_INT=true;
                break;
            case 'E':
                do_genEhrhart = true;
                break;
            case 'I':
                do_integral=true;
                break;
            case 'L': 
                do_leadCoeff=true;  
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
    
    if(!do_integral && !do_leadCoeff) // default is -E
        do_genEhrhart=true;
    
    cout << "+++ " << option << " " << do_genEhrhart << " " << do_integral << " " << do_leadCoeff << endl;

    if(do_genEhrhart)
        generalizedEhrhartSeries(project);
    else{
        if(do_integral)
            integrate(project,false);
        else{
            if(do_leadCoeff)
                integrate(project,true);
        }
    }
            
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }
  return 1;
}


