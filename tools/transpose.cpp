// reads a matrix in Normalizz format and prints the transpose to stdout
// again in Normaliz format

#include <stdlib.h>
#include <vector>
#include <list>
#include <map>

#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;


vector<vector<int> >  readMat(const string& project){
// reads one matrix from .in file

    string name_in=project;
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false){
        cerr << "Cannot find input file" << endl;
        exit(1);
    }

    int nrows,ncols;
    in >> nrows;
    in >> ncols;

    if(nrows==0 || ncols==0){
        cerr << "Matrix empty" << endl;
        exit(1);
    }


    int i,j,entry;
    vector<vector<int> > result(nrows);

    for(i=0;i<nrows;++i)
        for(j=0;j<ncols;++j){
            in >> entry;
            result[i].push_back(entry);
        }
    return(result);
}

int main(int argc, char* argv[])
{

    if(argc<2){
        cerr << "No input file given" << endl;
        exit(1);
    }
    string input_name=argv[1];
    vector< vector<int > > M;
    M=readMat(input_name);
    cout << M[0].size() << " " << M.size() << endl;
    for(size_t i=0;i<M[0].size();++i){
        for(size_t j=0;j<M.size();++j)
            cout << M[j][i] << " ";
       cout << endl;     
    }
    return(0);
}
