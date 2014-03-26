// Collects equal columns of a matrix and lists each with multiplicity

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
    vector< vector<int > > M, MT;
    M=readMat(input_name);
    MT.resize(M[1].size());
    for(size_t i=0;i<M.size();++i)
        for(size_t j=0;j<MT.size();++j)
            MT[j].push_back(M[i][j]);

    map< vector<int>, size_t > classes;
    map< vector<int>, size_t >::iterator C;

    for(size_t j=0;j<MT.size();++j){
        C=classes.find(MT[j]);
        if(C!=classes.end())
            C->second++;
        else
            classes.insert(pair<vector<int>, size_t>(MT[j],1));
    }
    
    for(size_t i=0;i<M.size();++i){
        for(C=classes.begin();C!=classes.end();++C)
            cout << C->first[i] << " ";
        cout << endl;
    }
    cout << "--------------" << endl;
    
    for(C=classes.begin();C!=classes.end();++C)
            cout << C->second << " ";
        
    cout << endl << endl;

    return(0);
}
