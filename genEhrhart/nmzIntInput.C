typedef unsigned short key_type;

struct STANLEYDATA_INT{
        vector<key_type> key;
        vector<vector<long> > offsets;
        bool done;
    };
    
struct TRIDATA{
        vector<key_type> key;
        long vol;  
};

void fileMissing(const char* name_in){
    cerr << "Could not open file " << name_in << endl << flush;
    exit(1);
}

void inputError(const char* name_in, const string message){
    cerr << "File " << name_in << ": " << message << endl << flush;
    exit(1);
} 

void getGrading(const string& project, vector<long>& grading, long& gradingDenom){
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
    long vectorSize;
    while(1){
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

    while(1){
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

void readGens(const string& project, vector<vector<long> >& gens){
// reads the generaors from the tgn file
    string name_in=project+".tgn";
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false) {
        fileMissing(file_in);
    }
    long nrows, ncols;
    in >> nrows >> ncols;
    if(nrows==0 || ncols==0){
        inputError(file_in,"seems corrupted.");
    }


    long i,j;
    gens.resize(nrows);
    for(i=0;i<nrows;++i)
        gens[i].resize(ncols);

    for(i=0; i<nrows; i++){
        for(j=0; j<ncols; j++) {
            in >> gens[i][j];
        }
    }
}

void polyError(const string& token){
    cerr << "Error in polynomoial input at token " << token << endl << flush;
    exit(1);
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
        for(size_t i=1;i<token.size();++i){
            if(token[i]=='x')
                polyError(token);
            if(token[i]=='['){
                if(opening)
                    polyError(token);
                opening=true;
            }
            if(token[i]==']'){
                if(!opening || closing)
                    polyError(token);
                closing=true;
            }
        }
        if(!opening || !closing)
            polyError(token);
        bool withexponent=false;
        for(size_t i=1;i<token.size();++i)
            if(token[i]=='^')
                withexponent=true;
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
// cuts the instring into tokens
// 
// care must be taken since some characters that indicate the start of a new token
// are themselves tokens and others only start the new token
// oneback takes care of this problem
    bool first=true;
    size_t oneback=0;

    for(size_t i=0;i<inString.size();){
        size_t j;
        for(j=i;j<inString.size();++j)
            if(inString[j]=='(' || inString[j]==')' || inString[j]=='*' || inString[j]=='+' || inString[j]=='-' || (inString[j]=='x' && !first))
                break;
        first=false;
        if(j>i)
            Tokens.push_back(inString.substr(i-oneback,j-(i-oneback)));
        i=j;
        if(inString[j]!='x'){
            Tokens.push_back(inString.substr(i,1));
            oneback=0;
        }
        else
           oneback=1;
        i++;
    }
}


RingElem readPolynomial(const string& project, const SparsePolyRing R){
// reads polynomial from file pnm
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
    string token, inString;
    long typ, nrIndet,expo;
    vector<string> Tokens;
    size_t current=0;

    bool pm_preceding=false, first=true;
    newCoeff=1;

    while(in.good() || current < Tokens.size()){  // end of file not yet reached or remaining tokens to be processed

       if(current+1>Tokens.size()){ // indicates that the tokens already read
           Tokens.clear();          // have been completely processed
           in >> inString;          // and we must read a new string.
           if(in.fail()){           // eof reached: the last term must now be added to the current factor.
               if(first)            // first indicates that a new term must follow
                   polyError(token);
               newPoly[newPoly.size()-1]+=monomial(R,newCoeff,newExpo);
               break;
           }
           tokenize(inString,Tokens); // break the input string into tokens
           current=0;
        }

        token=Tokens[current];
        current++;

        if(token.size()==0 || token=="(" || token==")") // ( and ) are skipped
            continue;

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

    } // while in.good || current < Tokens.size()


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

    while(in.good()){
        in >> newSimpl.key[0];
        if(in.fail())
            break;
        for(i=1;i<dim;++i)
            in >> newSimpl.key[i];
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
        newSimpl.done=false;
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
      cout << "Cannot find File " << name_in << ". Trying to read triangulation from dec." << endl << flush;
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
        in >> newSimpl.vol;
        triang.push_back(newSimpl);        
    }
}
