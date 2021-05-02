#include <cstdlib>
#include <vector>
#include <fstream>
#include <omp.h>
#include <iomanip>
#include <thread>
#include <chrono>
#include <iostream>
#include <set>
#include <vector>
#include <fstream>
#include <string>


using namespace std;

bool exists_file(string file_name){
    
   ifstream ifile;
   ifile.open(file_name);
   if(ifile) 
      return true;
   else 
      return false;   
}

// using namespace libnormaliz;

int main(int argc, char* argv[]) {
    
    if(argc < 6){
        cout << "Not enough parameters" << endl;
        exit(1);
    } 
    
    string project_name = string(argv[1]);
    
    long start_job = stoi(string(argv[2]));
    long end_job = stoi(string(argv[3]));
    long max_simultaneously = stoi(string(argv[4]));
    long nr_threads = stoi(string(argv[5]));
    
    vector<bool> done(end_job+1), running(end_job+1);
    long nr_done = 0;
    for(long i=0; i<start_job;++i){
        nr_done++;
        done[i] = true;
    }

    for(long j=start_job; j <=end_job;++j){
        string file_name = project_name+".hollow_tri."; // discard non-existing input files
        if(! exists_file(file_name + to_string(j)+".gz") ){
            done[j] = true;
            nr_done++;
            continue;
        }
        file_name = project_name+".mult."; // discard done jobs
        if(exists_file(file_name + to_string(j))){
            done[j] = true;
            nr_done++;
            continue;
        }
    }
    
    if(nr_done == end_job - start_job +1){
            cout << "Nothing to do" << endl;
            exit(0);
    }
    
    long nr_running = 0;
    
    while(true){
        for(long j=start_job; j <=end_job;++j){
            if(done[j])
                continue;
            string file_name = project_name+".mult.";
            if(exists_file(file_name + to_string(j))){
                done[j] = true;
                nr_done++;
                if(running[j]){
                    nr_running--;
                    running[j] = false;
                }
                continue;
            }
            std::this_thread::sleep_for (std::chrono::seconds(1));

            if(nr_running < max_simultaneously && !running[j]){
                long to_be_done = j;
                running[to_be_done] = true;
                nr_running++;
                string command = "nohup ../run_single.sh " + project_name + " ";
                command +=  to_string(to_be_done) + " " + to_string(nr_threads) + " &";
                int failure = system(command.c_str());
                if(failure > 0){
                    cout << "Command ended with exit code > 0";
                    exit(1);
                }
                break;
            }
            
        }
        if(nr_done == end_job - start_job +1 && nr_running == 0){
                cout << "All jobs done" << endl;
                exit(0);
        }
    
        std::this_thread::sleep_for (std::chrono::seconds(1));
        
    }    
}
