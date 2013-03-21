/*
 * nmzIntegrate
 * Copyright (C) 2012-2013  Winfried Bruns, Christof Soeger
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


#include "CoCoA/library.H"
using namespace CoCoA;

#include <fstream>
#include <sstream>

#include "../libnormaliz/my_omp.h"

#include "../libnormaliz/HilbertSeries.h"
#include "../libnormaliz/matrix.h"
// #include "../libnormaliz/vector_operations.h"
#include "../libnormaliz/map_operations.h"
 
using namespace std;
using namespace libnormaliz;

bool verbose_INT;

#include "nmzIntInput.C"
#include "nmzIntPoly.C"

#include "cyclRatFunct.C"
#include "genEhrhart.C"
#include "nmzIntegral.C"                 

void printHeader() {
    cout << "                                                    \\.....|"<<endl;
    cout << "                     nmzIntegrate 1.0                \\....|"<<endl;
    cout << "                                                      \\...|"<<endl;
    cout << "                 (C) W. Bruns  C. Soeger               \\..|"<<endl;
    cout << "                        March 2013                      \\.|"<<endl;
    cout << "                                                         \\|"<<endl;
}
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
                printHeader();
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

    if (verbose_INT) {
        printHeader();
    }

    if(!do_integral && !do_leadCoeff) // default is -E
        do_genEhrhart=true;
    
    // cout << "+++ " << option << " " << do_genEhrhart << " " << do_integral << " " << do_leadCoeff << endl;

    bool homogeneous=false;
    if(do_genEhrhart)
        generalizedEhrhartSeries(project,homogeneous);
    if(do_leadCoeff && !do_genEhrhart)
        integrate(project,true,homogeneous);
    if(do_integral && !homogeneous)
        integrate(project,false,homogeneous);
            
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


