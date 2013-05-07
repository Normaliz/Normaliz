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

bool existsFile(const string& project, const string& suffix, const bool& mustExist){
//n check whether file project.suffix exists

    string name_in=project+"."+suffix;
    // cout << name_in << endl;
    const char* file_in=name_in.c_str();
    ifstream in2;
    in2.open(file_in,ifstream::in);
    if (in2.is_open()==false){
        if(mustExist)
            fileMissing(file_in); //and exit ...
        return(false);
    }
    in2.close();
    return(true);
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

void polyError(const string& token){
    cerr << "Error in polynomoial input at token " << token << endl;
    exit(1);
}

bool allDigits(const string& testString){

    for(size_t i=0;i<testString.size();++i)
        if(testString[i]<'0' || testString[i]>'9')
            return(false);
    return(true);

}

long analyzeToken(const string& token, BigRat& newCoeff, long& nrIndet, long& expo){
// analyzes a token in the polynomial input
// arguments 2-4 are used for returning values if type = 4 or 5
// types mof tokens:
// 1: *
// 2: +
// 3: -
// 4: (power of) indeterminate
// 5: rational number

    if(token=="*")
        return(1);
    if(token=="+")
        return(2);
    if(token=="-")
        return(3);
        
    if(token[0]=='x'){
        bool opening=false,closing=false;
        long openingAt=0,closingAt=0,expoAt=0;
        for(size_t i=1;i<token.size();++i){
            if(token[i]=='x')
                polyError(token);
            if(token[i]=='['){
                if(opening)
                    polyError(token);
                opening=true;
                openingAt=i;
                if(openingAt!=1)
                    polyError(token);
            }
            if(token[i]==']'){
                if(!opening || closing)
                    polyError(token);
                closing=true;
                closingAt=i;
            }
        }
        if(!opening || !closing)
            polyError(token);
        bool withexponent=false;
        for(size_t i=1;i<token.size();++i)
            if(token[i]=='^'){
                withexponent=true;
                expoAt=i;
                if(expoAt!=closingAt+1)  // ^ must directly follow ]
                    polyError(token);
        }
        string indetString=token.substr(openingAt+1,closingAt-openingAt-1);
        if(!allDigits(indetString))
            polyError(token);
        if(!withexponent && token[token.size()-1]!=']')
            polyError(token);
        if(withexponent){
            string expoString=token.substr(expoAt+1,token.size()-expoAt-1);
            if(!allDigits(expoString))
             polyError(token);
        }

        const char* token_C=token.c_str();
        expo=0;
        nrIndet=0;
        if(withexponent)
            sscanf(token_C,"x[%ld]^%ld",&nrIndet,&expo);
        else{
            expo=1;
            sscanf(token_C,"x[%ld]",&nrIndet);
        }
        if(expo<=0 || nrIndet<=0)
            polyError(token);
        return(4);
    }
    newCoeff=BigRat(token);
    return(5);
}

void tokenize(const string& inString, vector<string>& Tokens){
// cuts the inString into tokens
// 
// care must be taken since some characters that indicate the start of a new token
// are themselves tokens and others only (at present only 'x') start the new token

    bool isFullToken;
    string current;
      
    for(size_t i=0;i<inString.size();++i){
        isFullToken=false;
        if(inString[i]=='(' || inString[i]==')' || inString[i]=='*' || inString[i]=='+' || inString[i]=='-' || inString[i]=='x'){
            if(current!="")
                Tokens.push_back(current);
            current="";
            if(inString[i]!='x'){
                Tokens.push_back(inString.substr(i,1));
                isFullToken=true;    
            }
        }
        if(!isFullToken)
            current+=inString.substr(i,1);
    }
    if(current!="")
        Tokens.push_back(current);
}
           

vector<RingElem> readFactorList(const string& project, const SparsePolyRing& R){
// reads factors of polynomial from file pnm
    string name_in=project+".pnm";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }

    vector<RingElem> newPoly;
    newPoly.push_back(zero(R));
    long dim=NumIndets(R)-1;
    vector<long> newExpo(dim+1);
    BigRat newCoeff, readCoeff;
    string token, inString,readString;
    long typ, nrIndet,expo;
    vector<string> Tokens;
    size_t current;

    bool pm_preceding=false, first=true;
    newCoeff=1;

    while(in.good()){
        in >> readString;;
        if(in.fail())
            break;
        inString+=readString;
    }
    
    // cout <<  inString << endl;
    // exit(0);
    tokenize(inString,Tokens);
    
    if(Tokens.size()==0){
        cerr << "No polynomial given." << endl;
        exit(1);
    }
    
    
    /* cout <<"Tok " << Tokens.size() << endl;
    
    for(size_t i=0;i< Tokens.size();++i)
        cout << Tokens[i] << endl; */
    
    for(current=0;;++current){
    
        if(current==Tokens.size()){
            newPoly[newPoly.size()-1]+=monomial(R,newCoeff,newExpo); // add the very last token
            break;
        }

        token=Tokens[current];
        // cout << "*** " << token << endl;

        if(token.size()==0 || token=="(" || token==")") // ( and ) are skipped
            continue;
            
        if(current==Tokens.size()-1 && (token=="+" || token=="-" || token=="*")) // cannot end with such a token
            polyError(token);

        typ=analyzeToken(token,readCoeff,nrIndet,expo);

        if((typ<=3 && pm_preceding) || (typ==1 && first)) // pm_preceding indicates that
            polyError(token);               // the preceding token was + or -

        if(typ==5){                 // rational coefficient
            if(!first)              // must be the first token of the term
                polyError(token);
            newCoeff*=readCoeff;    // this takes care of the sign: the previous token has set it
            pm_preceding=false;
            first=false;
            continue;
        }
        if(typ==4){  // (power) indeterminate
            if(nrIndet>dim)
                polyError(token);
            newExpo[nrIndet]+=expo;
            pm_preceding=false;
            first=false;
            continue;
        }
        if(typ<=3 && !first){ // term finished, must be added to factor
            newPoly[newPoly.size()-1]+=monomial(R,newCoeff,newExpo);
            for(size_t i=0;i<newExpo.size();++i)  // reset exponent vector
                newExpo[i]=0;
        }
        if(typ==1){ // now even factor finished and we start a new one
            newPoly.push_back(zero(R));
        }

        newCoeff=1; // starting a new term
        first=true;
        pm_preceding=false;
        if(typ==2 || typ==3)   // the current token is + or -
            pm_preceding=true; 
        if(typ==3)     // the current token is -
            newCoeff=-1;

    } // current
    
    /*for(size_t i=0;i<newPoly.size();++i)
        cout << newPoly[i] << endl; */


    return(newPoly);
}

RingElem readPolynomial(const string& project, const SparsePolyRing R){

    vector<RingElem> newPoly=readFactorList(project,R);
    RingElem p(one(R));
    for(size_t i=0;i<newPoly.size();++i)  // multiply the factors
        p*=newPoly[i];
    return(p);
}

void readDec(const string& project, const long& dim, list<STANLEYDATA_INT>& StanleyDec){
// rads Stanley decomposition from file dec
    string name_in=project+".dec";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }

    STANLEYDATA_INT newSimpl;
    long i=0,j,det,dummy;
    newSimpl.key.resize(dim);
    
    long test;

    while(in.good()){
        in >> newSimpl.key[0];
        if(in.fail())
            break;
        test=-1;
        for(i=1;i<dim;++i){
            in >> newSimpl.key[i];
            if(newSimpl.key[i]<=test){
                cerr << "Fatal error: Key of simplicial cone not ascending" << endl;
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
}

void readTri(const string& project, const long& dim, list <TRIDATA>& triang){
// rads triangulation from file tri
    string name_in=project+".tri";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false && verbose_INT) {
      cout << "Cannot find File " << name_in << ". Trying to read triangulation from dec." << endl;
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
