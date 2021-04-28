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
    
    if(argc < 5){
        cout << "Not enough parameters" << endl;
        exit(1);
    }    
    
    long start_job = stoi(string(argv[1]));
    long end_job = stoi(string(argv[2]));
    long max_simultaneously = stoi(string(argv[3]));
    long nr_threads = stoi(string(argv[4]));
    
    set<long> jobs;
    for(long i = start_job; i<=end_job;++i)
        jobs.insert(i);
    
    vector<bool> done(end_job+1), running(end_job+1);
    bool all_done = false;
    long nr_running = 0;
    long nr_done = 0;
    
    while(true){
        for(long j=start_job; j <=end_job;++j){
            if(done[j])
                continue;
            string file_name = "mult.";
            if(running[j] && exists_file(file_name + to_string(j))){
                running[j] = false;
                nr_running--;
                done[j] = true;;
                nr_done++;
            }
            std::this_thread::sleep_for (std::chrono::seconds(3));
            if(nr_done == end_job - start_job +1){
                all_done = true;
                break;
            }
            if(nr_running < max_simultaneously && jobs.size() > 0){
                long to_be_done = *(jobs.begin());
                jobs.erase(to_be_done);
                running[to_be_done] = true;
                string command = "nohup ../run_single.sh " + to_string(to_be_done) + " " + to_string(nr_threads) + " &>> " + "chunk.log."
                + to_string(to_be_done) +" &";
                nr_running++;
                int failure = system(command.c_str());
                if(failure > 0){
                    cout << "Command ended with exit code > 0";
                    exit(1);
                }
                break;
            }
            
        }
        if(all_done)
            break;
        std::this_thread::sleep_for (std::chrono::seconds(3));
        
    }    
}
