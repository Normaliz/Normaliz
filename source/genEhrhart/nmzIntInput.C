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

typedef int key_type;

struct STANLEYDATA_INT{
        vector<key_type> key;  // read from dec file
        vector<long> degrees;
        vector<vector<long> > offsets;  // ditto
        size_t classNr;  // number of class of this simplicial cone
};
    
struct TRIDATA{
        vector<key_type> key;
        long vol;  
};

struct SIMPLINEXDATA_INT{                        // local data of excluded faces
        boost::dynamic_bitset<> GenInFace;   // indicator for generators of simplex in face 
        long mult;                           // multiplicity of this face
        size_t card;                                // the cardinality of the face
        bool done;                           // indicates that this face has been done for a given offset
        vector<long> denom;
        vector<long> degrees;
        vector<long> key;
};

void fileMissing(const char* name_in){
    cerr << "Could not open file " << name_in << endl;
    exit(1);
}

void inputError(const char* name_in, const string message){
    cerr << "File " << name_in << ": " << message << endl;
    exit(1);
} 

long scalProd(const vector<long>& a, const vector<long>& b){
    long s=0;
    for(size_t i=0;i<a.size();++i)
        s+=a[i]*b[i];
    return(s);
}


bool existsFile(const string& project, const string& suffix, const bool& mustExist, time_t& fileDate){
//n check whether file project.suffix exists and retrieve last access time

    string name_in=project+"."+suffix;
    const char* file_in=name_in.c_str();
    
    struct stat fileStat;
    if(stat(file_in,&fileStat) < 0){
        if(mustExist)
            fileMissing(file_in);
         return(false); 
    }
    fileDate=fileStat.st_mtime;
    return(true);
}

bool existsFile(const string& project, const string& suffix, const bool& mustExist){
//n check whether file project.suffix exists

    time_t dummy;
    return(existsFile(project,suffix,mustExist,dummy));
}

string pureName(const string& fullName){
// extracts the pure filename

    string slash="/";
    #ifdef _WIN32 //for 32 and 64 bit windows
        slash="\\";
    #endif
    size_t found = fullName.rfind(slash);
    if(found==std::string::npos)
        return(fullName);
    found++;
    size_t length=fullName.size()-found;
    
    // cout << "**************************** " << fullName.substr(found,length) << endl;
    // exit(1);
    return(fullName.substr(found,length));  	

}

string fullPnmName(const string& project, const string& pnm){

    string name=project;
    string slash="/";
    #ifdef _WIN32 //for 32 and 64 bit windows
        slash="\\";
    #endif
    size_t found = name.rfind(slash);
    if(found==std::string::npos)
        return(pnm);
    found++;
    size_t length=project.size()-found;
    name.replace(found,length,pnm);
    // cout << "**************************** " << name << endl;
    // exit(1);
    return(name);  	
}


void getRankAndGrading(const string& project,long& rank,vector<long>& grading, long& gradingDenom){
// reads the grading data ffrom the inv file
    string name_in=project+".inv";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }
    string data;
    bool found=false;
    
    while(in.good()){
        in >> data;
        if(data=="rank"){
            found=true;
            in >> data; // skip = sign
            in >> rank;
            break;
        }
    }
    if(!found){
        inputError(file_in,"seems corrupted.");
    }
    
    found=false;
    long vectorSize;
    while(in.good()){
        in >> data;
        if(data=="vector")
            in >> vectorSize;
        if(data=="grading"){
            found=true;
            break;
        }
    }
    if(!found){
        inputError(file_in,"seems corrupted.");
    }

    in >> data; // skip = sign
    long g;
    for(long i=0;i<vectorSize;++i){
        in >> g;
        grading.push_back(g);
    }

    while(in.good()){
        in >> data;
        if(data=="grading_denom"){
            found=true;
            break;
        }
    }
    if(!found){
        inputError(file_in,"seems corrupted.");
    }
    in >> data; // skip = sign
    in >> gradingDenom;
    in.close();
}

void readGens(const string& project, vector<vector<long> >& gens, const vector<long>& grading){
// reads the generators from the tgn file
    string name_in=project+".tgn";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }
    long nrows, ncols;
    in >> nrows >> ncols;
    if(nrows==0 || ncols!=(long) grading.size()){
        inputError(file_in,"seems corrupted.");
    }


    long i,j;
    gens.resize(nrows);
    for(i=0;i<nrows;++i)
        gens[i].resize(ncols);

    long degree,prevDegree=1;
    for(i=0; i<nrows; i++){
        for(j=0; j<ncols; j++) {
            in >> gens[i][j];
        }
        degree=scalProd(gens[i],grading);
        if(degree<prevDegree){
               cerr << "Fatal error: degrees of generators not weakly ascending" << endl;
               cerr << "PLEASE CONTACT THE AUTHORS" << endl;
               exit(1);
        }
        prevDegree=degree;
    }
}



vector<RingElem> readFactorList(const string& project, const SparsePolyRing& R){
// reads factors of polynomial from file pnm
// This is now a single factor
    string name_in=project+".pnm";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }

    vector<RingElem> newPoly;
    RingElem InPoly(zero(R));
    InPoly=ReadExpr(in,R);
    newPoly.push_back(InPoly);
    return(newPoly);
}

RingElem readPolynomial(const string& project, const SparsePolyRing R){

    vector<RingElem> newPoly=readFactorList(project,R);
    RingElem p(one(R));
    for(size_t i=0;i<newPoly.size();++i)  // multiply the factors
        p*=newPoly[i];
    return(p);
}

void readInEx(ifstream& in, vector<pair<boost::dynamic_bitset<>, long> >& inExCollect, const size_t nrGen){

    size_t inExSize, keySize;
    in >> inExSize;
    key_type key;    
    long mult;
    boost::dynamic_bitset<> indicator(nrGen);
    for(size_t i=0;i<inExSize;++i){
        in >> keySize;
        indicator.reset();
        for(size_t j=0;j<keySize;++j){
            in >> key;
            indicator.set(key-1);
        }
        in >> mult;
        inExCollect.push_back(pair<boost::dynamic_bitset<>, long>(indicator,mult));       
    }
}

void readDecInEx(const string& project, const long& dim, list<STANLEYDATA_INT>& StanleyDec,
                vector<pair<boost::dynamic_bitset<>, long> >& inExCollect, const size_t nrGen){
// rads Stanley decomposition from file dec
    string name_in=project+".dec";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }
    
    string keyWord;
    in >> keyWord;
    if(keyWord!="in_ex_data" && keyWord!="Stanley_dec"){
        cerr << "Fatal error: dec file does not start with \"in_ex_data\" or \"Stanley_dec\"" << endl;
        cerr << "Potential reason: obsolete version of Normaliz" << endl;
        exit(1);
    }
    
    if(keyWord=="in_ex_data"){
        readInEx(in, inExCollect,nrGen);
        in >> keyWord;
        if(keyWord!="Stanley_dec"){
            cerr << "Fatal error: second keyword in dec file is not \"Stanley_dec\"" << endl;
            exit(1);
        }
    }
        
    size_t decSize;
    in >> decSize;

    STANLEYDATA_INT newSimpl;
    long i=0,j,det,dummy;
    newSimpl.key.resize(dim);
    
    long test;

    while(in.good()){
        in >> newSimpl.key[0];
        if(in.fail())
            break;

        for(i=1;i<dim;++i)
            in >> newSimpl.key[i];

        test=0;
        for(i=0;i<dim;++i){
            if(newSimpl.key[i]<=test){
                cerr << "Fatal error: Key of simplicial cone not ascending or out of range" << endl;
                cerr << "PLEASE CONTACT THE AUTHORS" << endl;
                exit(1);
            }
            test=newSimpl.key[i];
        }
        
        in >> det;
        in >> dummy;
        if(dummy!=dim){
            inputError(file_in,"wrong dimension in file.");
        }
        newSimpl.offsets.resize(det);
        for(i=0;i<det;++i)
            newSimpl.offsets[i].resize(dim);
        for(i=0;i<det;++i)
            for(j=0;j<dim;++j)
                in >> newSimpl.offsets[i][j];
        StanleyDec.push_back(newSimpl);
    }
    
    if(decSize!=StanleyDec.size()){
        cerr << "Fatal error: Actual size of Stanley decomposition does not match announced size" << endl;
        exit(1);
    }
    
}

void readTriFromDec(const string& project, const long& dim,list<TRIDATA>& triang){
// rads triangulation from Stanley decomposition file
    string name_in=project+".dec";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }

    TRIDATA newSimpl;
    long i=0,j,det,dummy;
    newSimpl.key.resize(dim);
    
    while(true){   // skip in_ex_data if present  
        string dummy;
        in >> dummy;
        if(dummy=="Stanley_dec")
            break;
    }
    
    size_t decSize;
    in >> decSize;

    while(in.good()){
        in >> newSimpl.key[0];
        if(in.fail())
            break;
        for(i=1;i<dim;++i)
            in >> newSimpl.key[i];
        sort(newSimpl.key.begin(),newSimpl.key.end()); // should come sorted, nevertheless for stability
        in >> det;
        newSimpl.vol=det;
        in >> dummy;
        if(dummy!=dim){
            inputError(file_in,"wrong dimension in file.");
        }

        for(i=0;i<det;++i)
            for(j=0;j<dim;++j)
                in >> dummy;
        triang.push_back(newSimpl);
    }
    if(decSize!=triang.size()){
        cerr << "Fatal error: Actual size of Stanley decomposition does not match announced size" << endl;
        exit(1);
    }
    
}

void readTri(const string& project, const long& dim, list <TRIDATA>& triang){
// rads triangulation from file tri
    string name_in=project+".tri";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false){
        if(verbose_INT) {
            cout << "Cannot find File " << name_in << ". Trying to read triangulation from dec." << endl;
        }
      readTriFromDec(project,dim,triang);
      return;   
    }
    
    long dummy;
    in >> dummy; // number of simpl in triang, not needed here
    in >> dummy;
    if(dim!=dummy-1){
            inputError(file_in,"wrong dimension in file.");
    }

    long i;
    TRIDATA newSimpl;
    newSimpl.key.resize(dim);
    
    while(in.good()){
        in >> newSimpl.key[0];
        if(in.fail())
            break;    
        for(i=1;i<dim;++i)
            in >> newSimpl.key[i];
        sort(newSimpl.key.begin(),newSimpl.key.end()); // should come sorted, nevertheless for stability
        in >> newSimpl.vol;
        triang.push_back(newSimpl);        
    }
}
