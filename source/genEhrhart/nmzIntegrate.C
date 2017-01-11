/*
 * nmzIntegrate
 * Copyright (C) 2012-2014  Winfried Bruns, Christof Soeger
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


#include "CoCoA/library.H"
using namespace CoCoA;

#include <fstream>
#include <sstream>
#include<string>
#include <sys/stat.h>
#include <sys/types.h>

#include <boost/dynamic_bitset.hpp>

#include "../libnormaliz/my_omp.h"

#include "../libnormaliz/HilbertSeries.h"
#include "../libnormaliz/matrix.h"
// #include "../libnormaliz/vector_operations.h"
#include "../libnormaliz/map_operations.h"
 
using namespace std;
using namespace libnormaliz;

bool verbose_INT=false;

#include "nmzIntInput.C"
#include "cyclRatFunct.C"
#include "nmzIntPoly.C"
#include "genEhrhart.C"
#include "nmzIntegral.C"                 

void printHeader() {
    cout << "                                                    \\.....|"<<endl;
    cout << "                   nmzIntegrate 1.3.3                \\....|"<<endl;
    cout << "                                                      \\...|"<<endl;
    cout << "                (C) W. Bruns  C. Soeger                \\..|"<<endl;
    cout << "                      January 2017                      \\.|"<<endl;
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
    cout << "  -F=<pnm>\tspecifies file name <pnm> of polynomial" << endl;
    cout << "  --OutputDir=<path>\tnames the path to the output directory" << endl;
}

// the following two functions are copied from cone.cpp. The first has a sibbling in nmzIntInput.C
bool exists_file(string name_in){
//n check whether file name_in exists

    //b string name_in="nmzIntegrate";
    const char* file_in=name_in.c_str();
    
    struct stat fileStat;
    if(stat(file_in,&fileStat) < 0){
         return(false); 
    }
    return(true);
}

string command(const string& original_call, const string& to_replace, const string& by_this){
// in the original call we replace the program name to_replace by by_this
// wWe try variants with and without "lt-" preceding the names of executables 
// since libtools may have inserted "lt-" before the original name

    string copy=original_call;
    string search_lt="lt-"+to_replace;
    long length=to_replace.size();
    size_t found;
    found = copy.rfind(search_lt);
    if (found==std::string::npos) {
        found = copy.rfind(to_replace);
        if (found==std::string::npos){
            throw FatalException("Call "+ copy +" of "  +to_replace+" does not contain " +to_replace); 
        }
    }
    else{
            length+=3; //name includes lt-
    }
    string test_path=copy.replace (found,length,by_this);
    if(exists_file(test_path))
        return test_path;
    copy=original_call;
    string by_this_with_lt="lt-"+by_this;
    test_path=copy.replace (found,length,by_this_with_lt);
    cout << "TEST " << test_path << endl;
    if(exists_file(test_path))
        return test_path;
    return ""; // no executable found
}


//----------------------------------------------------------------------
// Use main() to analyze options, start computations and 
// handle any uncaught exceptions.
int main(int argc, char* argv[])  
{
  try
  {
    size_t i;   
    string project="",pnm="", output_dir="";
    
    string Threads;            

    string option;            //all options concatenated (including -)
    for (i = 1; i < (unsigned int)argc; i++) {
        
        string argument(argv[i]);
    
        if(argument=="-"){
            cerr << "Warning: empty option." << endl;
            continue;
        }

        if (argument[0]!='-') {
            if(project!=""){
                    cerr << "Fatal error: second project " << argument <<" specified." << endl;
                    exit(1);
            }
            project=argument;
            continue;
        }
        
        if(argument[1]=='-'){
            if(argument.substr(2,9)=="OutputDir" && argument[11]=='='){
                output_dir=argument.substr(12,argument.size()-1);
                if(output_dir.back()!='/')
                    output_dir+="/";
                continue;
            }            
        }
        
        if(argument[1]=='x'){
            if(argument.size()<=2 || argument[2]!='='){
                cerr << "Fatal error: -x not followed by =." << endl;
                exit(1);
            }
            #ifdef _OPENMP
            Threads = argument;
            Threads.erase(0,3);
            size_t nr_threads;
            if ( (istringstream(Threads) >> nr_threads) && nr_threads > 0) {
                omp_set_num_threads(nr_threads);
            } else {
                cerr<<"Fatal error: Invalid string following -x in "<< argument << endl;
                exit(1);
            }
            #else
            cerr << "Warning: Compiled without OpenMP support, option "<<argument<<" ignored."<<endl;
            #endif
            continue;
        }
        
        if(argument[1]=='F'){
            if(argument.size()<=2 || argument[2]!='='){
                cerr << "Fatal error: -F not followed by =." << endl;
                exit(1);
            }
            if(pnm!=""){
                    cerr << "Fatal error: second polynomial specified by " << argument << endl;
                    exit(1);
            }
            pnm=argument;
            pnm.erase(0,3);
            if(pnm==""){
                cerr << "Fatal error: option -F= without filename." << endl;
                exit(1);
            }
            if(pnm!=pureName(pnm)){
                cerr << "Fatal error: no path allowed in filename for polynomial." << endl;
                exit(1);
            }
            continue;        
        }

        option+=argument;
    } 
    
    // cout << "Options *****************************" << option << endl;   
    

    bool do_genEhrhart=false, do_integral=false, do_leadCoeff=false;
    
    for (i = 1; i <option.size(); i++) {
        switch (option[i]) {
            case '-': break;
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
                cerr<<"Fatal error : Option -x=<T> has to be separated from other options"<<endl;
                exit(1);
            case 'F': //should be separated from other options
                cerr<<"Fatal error: Option -F=<pnm> has to be separated from other options"<<endl;
                exit(1);
            default:
                cerr<<"Fatal error: Unknown option -"<<option[i]<<endl;
                exit(1);
        }
    }

    if (verbose_INT) {
        printHeader();
    }
    
    if(project==""){
        cerr << "Fatal error: No project specified" << endl;
        exit(1);
     }
     
     if(pnm==""){
        pnm=pureName(project);     
     }
     
     fullPnmName(project,pnm);

    if(!do_integral && !do_leadCoeff) // default is -E
        do_genEhrhart=true;
    bool do_int_or_lead=do_integral || do_leadCoeff;
    
    time_t pnmDate,mnz_inDate,invDate,tgnDate,decDate,triDate;
    
    existsFile(fullPnmName(project,pnm),"pnm",true,pnmDate); // true means: will exit if polynomial does not exist
                                                             // pnmDate irreleveant, only for completeness

    bool nmz_inExists=existsFile(project,"in",false,mnz_inDate); // checks if input file to Normaliz exists
    
    string nmz_output=project;
    if(output_dir!="")
        nmz_output=output_dir+pureName(project);
    /* if(verbose_INT)
        cout << "NMZ_OUT " << nmz_output << endl;*/
    bool invExists=existsFile(nmz_output,"inv",false,invDate);
    bool decExists=existsFile(nmz_output,"dec",false,decDate);
    bool triExists=existsFile(nmz_output,"tri",false,triDate); 
    bool tgnExists=existsFile(nmz_output,"tgn",false,tgnDate);
    
    bool makeInputFiles=!invExists || !tgnExists || (do_genEhrhart && !decExists)
                              || (do_int_or_lead && !triExists && ! decExists);
                              
    
    if(do_int_or_lead && !triExists)
        triDate=decDate;   
                              
    bool outdatedInputFiles= !makeInputFiles && nmz_inExists; // nevessary conditions
    if(outdatedInputFiles){
        outdatedInputFiles=invDate<mnz_inDate ||  tgnDate<mnz_inDate;
        if(!outdatedInputFiles && do_int_or_lead)
            outdatedInputFiles=triDate<mnz_inDate;
        if(!outdatedInputFiles && do_genEhrhart)
            outdatedInputFiles=decDate<mnz_inDate;
    }            
    
    
    // cout << "inv " << invExists << " dec " << decExists << " tri " << triExists << endl;

    if(makeInputFiles || outdatedInputFiles){
        if(verbose_INT)
            cout << "Input file(s) missing or  outdated. Trying to make them." << endl;
        existsFile(project,"in",true,mnz_inDate); // true ==> will exit if project.in does not exist
        
        string normalizExec("\""); // start with "
        
        // the quoting requirements for windows are insane, one pair of "" around the whole command and one around each file
        #ifdef _WIN32 //for 32 and 64 bit windows
            normalizExec.append("\"");
        #endif
        string normaliz_path=command(argv[0],"nmzIntegrate","normaliz");
        normalizExec.append(normaliz_path);
        normalizExec.append("\"");
        
        if(do_genEhrhart)
            normalizExec+=" -y";
        else
            normalizExec+=" -T";
                
        // normalizExec+="e";  // error check always activated // NO LONGER ncessary
        

        if(verbose_INT)
            normalizExec+="c";
        if(Threads!="")
            normalizExec+=(" -x="+Threads);
        if(output_dir!="")
            normalizExec+=(" --OutputDir="+output_dir);
            
            
        normalizExec+=(" \""); // enclose filename in by ""
        normalizExec+=project;
                    normalizExec+="\"";
                    
        #ifdef _WIN32 //for 32 and 64 bit windows
            normalizExec+="\"";
        #endif
        
        if(verbose_INT){
            cout << "executing: "<< normalizExec << endl;
            cout << "==========================================================" << endl;
        }        
        
        int returnvalue = system(normalizExec.c_str());
        if(returnvalue!=0){
            cerr <<  "Normaliz exited with error." << endl;
            exit(1);
        }
    }

    bool homogeneous=false;
    bool appendOutput=false;
    if (do_genEhrhart) {
        generalizedEhrhartSeries(project,output_dir,pnm,homogeneous);
        appendOutput=true;
        if(do_leadCoeff && verbose_INT){
            cout << "Suppressing computation LC since result contained in ES."  << endl;
        }
    }
    // cout << "hom " << homogeneous << endl;
    if (do_leadCoeff && !do_genEhrhart) {
        integrate(project,output_dir,pnm,true,homogeneous,appendOutput);
        appendOutput=true;
    }
    if(do_integral && homogeneous && verbose_INT){
            cout << "Suppressing computation Int since result contained in ES or LC."  << endl;
    }
    if (do_integral && !homogeneous) {
        integrate(project,output_dir,pnm,false,homogeneous,appendOutput);
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


