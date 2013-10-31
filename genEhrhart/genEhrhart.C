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


vector<RingElem> power2ascFact(const SparsePolyRing& P, const long& k)
// computes the representation of the power x^n as the linear combination
// of (x+1)_n,...,(x+1)_0
// return value is the vector of coefficients (they belong to ZZ)
{
    RingElem t=indets(P)[0];
    const vector<long> ONE(NumIndets(P));
    RingElem f(P),g(P), h(P);
    f=power(t,k);
    long m;
    vector<RingElem> c(k+1,zero(P));
    while(f!=0)
    {
            m=deg(f);
            h=monomial(P,LC(f),ONE);
            c[m]=h;
            f-=h*ascFact(t,m);
    }
    return(c);
}



CyclRatFunct genFunctPower1(const SparsePolyRing& P, long k,long n)
// computes the generating function for
//  \sum_j j^n (t^k)^j
{
    vector<RingElem> a=power2ascFact(P,n);
    RingElem b(P);
    vector<long> u;
    CyclRatFunct g(zero(P)), h(zero(P));
    long i,s=a.size();
    for(i=0;i<s;++i)
    {
        u=makeDenom(k,i+1);
        b=a[i]*factorial(i);
        g.set2(b,u); 
        h.addCRF(g);
    }
    return(h);
}


CyclRatFunct genFunct(const vector<vector<CyclRatFunct> >& GFP, const RingElem& F, const vector<long>& degrees)
// writes \sum_{x\in\ZZ_+^n} f(x,t) T^x
// under the specialization T_i --> t^g_i
// as a rational function in t
{
    SparsePolyRing P=AsSparsePolyRing(owner(F));
    RingElem t=indets(P)[0];
    
    CyclRatFunct s(F); // F/1
    
    CyclRatFunct g(zero(P)),h(zero(P));
    
    long nd=degrees.size();  
    long i,k,mg;
    vector<RingElem> c;   

    for(k=1; k<=nd;k++)
    {
        c=ourCoeffs(s.num,k); // we split the numerator according 
                               // to powers of var k
        mg=c.size(); // max degree+1 in  var k

        h.set2(zero(P));
        for(i=0;i<mg;i++)     // now we replace the powers of var k
        {                      // by the corrseponding rational function,
                               // multiply, and sum the products

            h.num=(1-power(t,degrees[k-1]))*h.num+GFP[degrees[k-1]][i].num*c[i];
            h.denom=GFP[degrees[k-1]][i].denom;
        }
        s.num=h.num;
        s.denom=prodDenom(s.denom,h.denom);
    }
    return(s);   
}

CyclRatFunct evaluateFaceClasses(const vector<vector<CyclRatFunct> >& GFP,
                                    map<vector<long>,RingElem>& faceClasses){
// computes the generating rational functions
// for the denominator classes collected from proper faces and returns the sum

    SparsePolyRing R=AsSparsePolyRing(owner(faceClasses.begin()->second));
    CyclRatFunct H(zero(R));
    // vector<CyclRatFunct> h(omp_get_max_threads(),CyclRatFunct(zero(R)));
    // vector<CyclRatFunct> h(1,CyclRatFunct(zero(R)));
    
    long mapsize=faceClasses.size();
    if(verbose_INT){    
        cout << "--------------------------------------------" << endl;
        cout << "Evaluating " << mapsize <<" face classes" << endl;
        cout << "--------------------------------------------" << endl;
    }
    #pragma omp parallel
    {
    
    map<vector<long>,RingElem>::iterator den=faceClasses.begin();
    long mpos=0;
    CyclRatFunct h(zero(R));
   
    #pragma omp for schedule(dynamic)
    for(long dc=0;dc<mapsize;++dc){
        for(;mpos<dc;++mpos,++den);
        for(;mpos>dc;--mpos,--den);
        // cout << "mpos " << mpos << endl;
        
        h = genFunct(GFP,den->second,den->first);
        h.simplifyCRF();
        if(verbose_INT){
            #pragma omp critical(VERBOSE)
            {
            cout << "Class ";
            for(size_t i=0;i<den->first.size();++i)
                cout << den->first[i] << " ";
            cout  << "NumTerms " << NumTerms(den->second) << endl;
        
            // cout << "input " << den->second << endl;
            }
        }
        
        // h.showCoprimeCRF();
        #pragma omp critical(ADDCLASSES)
        H.addCRF(h);
    }
    
    } // parallel 
    faceClasses.clear();
    H.simplifyCRF();
    return(H);        
}

struct denomClassData{
    vector<long> degrees;
    size_t simplDue;
    size_t simplDone;
  };

CyclRatFunct evaluateDenomClass(const vector<vector<CyclRatFunct> >& GFP,
                                    pair<denomClassData,vector<RingElem> >& denomClass){
// computes the generating rational function
// for a denominator class and returns it

    SparsePolyRing R=AsSparsePolyRing(owner(denomClass.second[0]));
    
    if(verbose_INT){
    #pragma omp critical(PROGRESS)
    {
        cout << "--------------------------------------------" << endl;
        cout << "Evaluating denom class ";
        for(size_t i=0;i<denomClass.first.degrees.size();++i)
            cout << denomClass.first.degrees[i] << " ";
        cout  << "NumTerms " << NumTerms(denomClass.second[0]) << endl;
        // cout << denomClass.second << endl;
        cout << "--------------------------------------------" << endl;
    }
    }

    CyclRatFunct h(zero(R));
    h = genFunct(GFP,denomClass.second[0],denomClass.first.degrees);

    denomClass.second[0]=0;  // to save memory
    h.simplifyCRF();
    return(h);
}

mpz_class ourFactorial(const long& n){
    mpz_class fact=1;
    for(long i=1;i<=n;++i)
        fact*=i;
    return(fact);
}

void transferFacePolys(deque<pair<vector<long>,RingElem> >& facePolysThread, 
                            map<vector<long>,RingElem>& faceClasses){


    // cout << "In Transfer " << facePolysThread.size() << endl;
    map<vector<long>,RingElem>::iterator den_found;                            
    for(size_t i=0;i<facePolysThread.size();++i){
        den_found=faceClasses.find(facePolysThread[i].first);
        if(den_found!=faceClasses.end()){
                den_found->second+=facePolysThread[i].second;    
        }
        else{
            faceClasses.insert(facePolysThread[i]);
            if(verbose_INT){
                #pragma omp critical(VERBOSE)
                {
                    cout << "New face class " << faceClasses.size() <<    " degrees ";
                    for(size_t j=0;j<facePolysThread[i].first.size();++j)
                        cout << facePolysThread[i].first[j] << " ";
                    cout << endl << flush;
                    }
            }
        } // else
    }
    facePolysThread.clear();
} 

void writeGenEhrhartSeries(const string& project, const factorization<RingElem>& FF,
                   const libnormaliz::HilbertSeries& HS, const long& virtDeg, const mpz_class & commonDen){
                            
    string name_open=project+".intOut";                              
    const char* file=name_open.c_str();
    ofstream out(file);
    
    out <<"Factorization of polynomial:" << endl;  // we show the factorization so that the user can check
    for(size_t i=0;i<FF.myFactors.size();++i)
        out << FF.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
    out << "Remaining factor " << FF.myRemainingFactor << endl << endl;
    
    out << "Generalized Ehrhart series:" << endl;
    vector<mpz_class> num( HS.getNum());
    for(size_t i=0;i<num.size();++i)
        out << num[i] << " ";
    out << endl << "Common denominator of coefficients: ";
    out << commonDen << endl;
    map<long, long> HS_Denom = HS.getDenom();
    long nr_factors = 0;
    for (map<long, long>::iterator it = HS_Denom.begin(); it!=HS_Denom.end(); ++it) {
        nr_factors += it->second;
    }
    out << "Series denominator with " << nr_factors << " factors:" << endl;
    out << HS.getDenom();
    out << endl;
    long period = HS.getPeriod();
    if (period == 1) {
        out << "Generalized Ehrhart polynomial:" << endl;
        for(size_t i=0; i<HS.getHilbertQuasiPolynomial()[0].size();++i)
            out << HS.getHilbertQuasiPolynomial()[0][i] << " ";
        out << endl;
        out << "with common denominator: ";
        out << HS.getHilbertQuasiPolynomialDenom();
        out << endl<< endl; 
    } else {
        // output cyclonomic representation
        out << "Generalized Ehrhart series with cyclotomic denominator:" << endl;
        num=HS.getCyclotomicNum();
        for(size_t i=0;i<num.size();++i)
            out << num[i] << " ";
        out << endl << "Common denominator of coefficients: ";
        out << commonDen << endl;
        out << "Series cyclotomic denominator:" << endl;
        out << HS.getCyclotomicDenom();
        out << endl;
        // Generalized Ehrhart quasi-polynomial
        vector< vector<mpz_class> > hilbert_quasi_poly = HS.getHilbertQuasiPolynomial();
        if (hilbert_quasi_poly.size() > 0) { // == 0 means not computed
            out<<"Generalized Ehrhart quasi-polynomial of period " << period << ":" << endl;
            Matrix<mpz_class> HQP(hilbert_quasi_poly);
            HQP.pretty_print(out,true);
            out<<"with common denominator: "<<HS.getHilbertQuasiPolynomialDenom();
        }
    }
    
   long deg=HS.getHilbertQuasiPolynomial()[0].size()-1;
   out  << endl << endl << "Degree of (quasi)polynomial: " << deg << endl;

   out << endl << "Expected degree: " << virtDeg << endl;
        
   out << endl << "Virtual multiplicity: ";

   mpq_class genMultQ;
   if(deg==virtDeg)   
       genMultQ=HS.getHilbertQuasiPolynomial()[0][virtDeg];
   genMultQ=genMultQ/(HS.getHilbertQuasiPolynomialDenom()*commonDen)*ourFactorial(virtDeg);
   out << genMultQ << endl;
}

mpz_class mpz(BigInt B) {
    return(mpz_class(mpzref(B)));  
}

libnormaliz::HilbertSeries nmzHilbertSeries(const CyclRatFunct& H, mpz_class& commonDen)
{ 

  size_t i;
  vector<RingElem> HCoeff0=ourCoeffs(H.num,0); // we must convert the coefficients
  BigInt commonDenBI(1);                         // and find the common denominator 
  vector<BigRat> HCoeff1(HCoeff0.size());
  for(i=0;i<HCoeff0.size();++i){
    IsRational(HCoeff1[i],HCoeff0[i]);          // to BigRat
    commonDenBI=lcm(den(HCoeff1[i]),commonDenBI);
  }
  
  commonDen=mpz(commonDenBI);   // convert it to mpz_class
  
  BigInt HC2;
  vector<mpz_class> HCoeff3(HCoeff0.size());
  for(i=0;i<HCoeff1.size();++i){
    HC2=num(HCoeff1[i]*commonDenBI);        // to BigInt
    HCoeff3[i]=mpz(HC2);      // to mpz_class 
  }

  vector<long> denomDeg=denom2degrees(H.denom);
  libnormaliz::HilbertSeries HS(HCoeff3,count_in_map<long, long>(denomDeg)); 
  HS.simplify();
  return(HS);
}

bool compareDegrees(const STANLEYDATA_INT& A, const STANLEYDATA_INT& B){

    return(A.degrees < B.degrees);
}

bool compareFaces(const SIMPLINEXDATA_INT& A, const SIMPLINEXDATA_INT& B){

    return(A.card > B.card);
}

void prepare_inclusion_exclusion_simpl(const STANLEYDATA_INT& S,
      const vector<pair<boost::dynamic_bitset<>, long> >& inExCollect, 
      vector<SIMPLINEXDATA_INT>& inExSimplData) {

    size_t dim=S.key.size();
    vector<key_type> key=S.key;
    for(size_t i=0;i<dim;++i)  // BECAUSE OF INPUT
        key[i]--;
    
    boost::dynamic_bitset<> intersection(dim), Excluded(dim);
    
    Excluded.set();
    for(size_t j=0;j<dim;++j)  // enough to test the first offset (coming from the zero vector)
        if(S.offsets[0][j]==0)
            Excluded.reset(j); 

    vector<pair<boost::dynamic_bitset<>, long> >::const_iterator F;    
    map<boost::dynamic_bitset<>, long> inExSimpl;      // local version of nExCollect   
    map<boost::dynamic_bitset<>, long>::iterator G;

    for(F=inExCollect.begin();F!=inExCollect.end();++F){
        // cout << "F " << F->first << endl;
       bool still_active=true;
       for(size_t i=0;i<dim;++i)
           if(Excluded[i] && !F->first.test(key[i])){
               still_active=false;
               break;
           }
       if(!still_active)
           continue;
       intersection.reset();
       for(size_t i=0;i<dim;++i){
           if(F->first.test(key[i]))
               intersection.set(i);
       }    
       G=inExSimpl.find(intersection);
       if(G!=inExSimpl.end())
           G->second+=F->second;
       else
           inExSimpl.insert(pair<boost::dynamic_bitset<> , long>(intersection,F->second)); 
    } 
    
    SIMPLINEXDATA_INT HilbData;
    inExSimplData.clear();
    vector<long> degrees;
    
    for(G=inExSimpl.begin();G!=inExSimpl.end();++G){
       if(G->second!=0){
           HilbData.GenInFace=G->first;
           HilbData.mult=G->second;
           HilbData.card=G->first.count();
           degrees.clear();
           for(size_t j=0;j<dim;++j)
             if(G->first.test(j))
                degrees.push_back(S.degrees[j]);
           HilbData.degrees=degrees;
           HilbData.denom=degrees2denom(degrees);
           inExSimplData.push_back(HilbData);
       }
    }
    
    sort(inExSimplData.begin(),inExSimplData.end(),compareFaces);
    
    /* for(size_t i=0;i<inExSimplData.size();++i)
        cout << inExSimplData[i].GenInFace << " ** " << inExSimplData[i].card << " || " << inExSimplData[i].mult << " ++ "<< inExSimplData[i].denom <<  endl;
    cout << "InEx prepared" << endl; */
        
}

void generalizedEhrhartSeries(const string& project, bool& homogeneous){
  GlobalManager CoCoAFoundations;
  
  long i,j;
  
  if(verbose_INT){
    cout << "==========================================================" << endl;
    cout << "Generalized Ehrhart series " << project << endl;
    cout << "==========================================================" << endl << endl;
  }

  vector<long> grading;
  long gradingDenom;
  long rank;
  getRankAndGrading(project,rank,grading,gradingDenom);
  
  vector<vector<long> > gens;
  readGens(project,gens,grading);
  if(verbose_INT)
    cout << "Generators read" << endl;
  long dim=gens[0].size();
  long maxDegGen=scalProd(gens[gens.size()-1],grading)/gradingDenom; 
  
  list<STANLEYDATA_INT> StanleyDec;
  vector<pair<boost::dynamic_bitset<>, long> > inExCollect;
  readDecInEx(project,rank,StanleyDec,inExCollect,gens.size());
  if(verbose_INT)
    cout << "Stanley decomposition (and in/ex data) read" << endl;
    
  size_t dec_size=StanleyDec.size();
    
  // Now we sort the Stanley decomposition by denominator class (= degree class)

  list<STANLEYDATA_INT>::iterator S = StanleyDec.begin();

  vector<long> degrees(rank);
  vector<vector<long> > A(rank);
  
  // prepare sorting by computing degrees of generators

  BigInt lcmDets(1); // to become the lcm of all dets of simplicial cones
  
  for(;S!=StanleyDec.end();++S){
      for(i=0;i<rank;++i)    // select submatrix defined by key
        A[i]=gens[S->key[i]-1];
          degrees=MxV(A,grading);
      for(i=0;i<rank;++i)
        degrees[i]/=gradingDenom; // must be divisible
      S->degrees=degrees;
      lcmDets=lcm(lcmDets,S->offsets.size());
  }
  
  if(verbose_INT)
    cout << "lcm(dets)=" << lcmDets << endl;
  
  StanleyDec.sort(compareDegrees);
  
  SparsePolyRing R=NewPolyRing(RingQQ(),dim+1,lex);
  SparsePolyRing RZZ=NewPolyRing(RingZZ(),dim+1,lex);
  const RingElem& t=indets(RZZ)[0];

  if(verbose_INT)
    cout << "Stanley decomposition sorted" << endl; 

  vector<pair<denomClassData, vector<RingElem> > > denomClasses;
  denomClassData denomClass;
  vector<RingElem> ZeroVectRingElem;
  for(int j=0;j<omp_get_max_threads();++j)
    ZeroVectRingElem.push_back(zero(RZZ));
  
  map<vector<long>,RingElem> faceClasses; // denominator classes for the faces
                 // contrary to denomClasses these cannot be sorted beforehand
                 
  vector<deque<pair<vector<long>,RingElem> > > facePolys(omp_get_max_threads()); // intermediate storage
  bool evaluationActive=false;

  // we now make class 0 to get started
  S=StanleyDec.begin();
  denomClass.degrees=S->degrees;  // put degrees in class
  denomClass.simplDone=0;
  denomClass.simplDue=1;           // already one simplex to be done 
  denomClasses.push_back(pair<denomClassData,vector<RingElem> >(denomClass,ZeroVectRingElem));
  size_t dc=0;
  S->classNr=dc; // assignment of class 0 to first simpl in sorted order

  list<STANLEYDATA_INT>::iterator prevS = StanleyDec.begin();

  for(++S;S!=StanleyDec.end();++S,++prevS){
    if(S->degrees==prevS->degrees){                     // compare to predecessor
        S->classNr=dc;              // assign class to simplex
        denomClasses[dc].first.simplDue++;         // number of simplices in class ++
    }
    else{
        denomClass.degrees=S->degrees;  // make new class
        denomClass.simplDone=0;
        denomClass.simplDue=1;
        denomClasses.push_back(pair<denomClassData,vector<RingElem> >(denomClass,ZeroVectRingElem));
        dc++;
        S->classNr=dc;
    }
  }

  if(verbose_INT)
    cout << denomClasses.size() << " denominator classes built" << endl;

  vector<RingElem> primeFactors;
  vector<RingElem> primeFactorsNonhom;
  vector<long> multiplicities;
  RingElem remainingFactor(one(R));
  
  RingElem F=processInputPolynomial(project,R,RZZ,primeFactors, primeFactorsNonhom,
                multiplicities,remainingFactor,homogeneous,false);
                
  vector<BigInt> Factorial(deg(F)+dim); // precomputed values
  for(i=0;i<deg(F)+dim;++i)
      Factorial[i]=factorial(i);
  
  factorization<RingElem> FF(primeFactors,multiplicities,remainingFactor); // assembeles the data
  factorization<RingElem> FFNonhom(primeFactorsNonhom,multiplicities,remainingFactor); // for output

  long nf=FF.myFactors.size();
  if(verbose_INT){
    cout <<"Factorization" << endl;  // we show the factorization so that the user can check
    for(i=0;i<nf;++i)
        cout << FFNonhom.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
    cout << "Remaining factor " << FF.myRemainingFactor << endl << endl;
  }

  vector<vector<CyclRatFunct> > GFP; // we calculate the table of generating functions
  vector<CyclRatFunct> DummyCRFVect; // for\sum i^n t^ki vor various values of k and n
  CyclRatFunct DummyCRF(zero(RZZ));
  for(j=0;j<=deg(F);++j)
    DummyCRFVect.push_back(DummyCRF);
  for(i=0;i<=maxDegGen;++i){
    GFP.push_back(DummyCRFVect);
    for(j=0;j<=deg(F);++j)
        GFP[i][j]=genFunctPower1(RZZ,i,j);
  }

  CyclRatFunct H(zero(RZZ)); // accumulates the series
  
  if(verbose_INT){
    cout << "********************************************" << endl;
    cout << dec_size <<" simplicial cones to be evaluated" << endl;
    cout << "********************************************" <<  endl;
  }
 
  size_t nrSimplDone=0;

  #pragma omp parallel private(i)
  {

  long degree_b;
  long det;
  bool evaluateClass;
  vector<long> degrees;
  vector<vector<long> > A(rank);
  list<STANLEYDATA_INT>::iterator S=StanleyDec.begin();

  RingElem h(zero(RZZ));     // for use in a simplex
  CyclRatFunct HClass(zero(RZZ)); // for single class
  

  size_t s,spos=0;  
  #pragma omp for schedule(dynamic) 
  for(s=0;s<dec_size;++s){
    for(;spos<s;++spos,++S);
    for(;spos>s;--spos,--S);

    det=S->offsets.size();
    degrees=S->degrees;
    
    for(i=0;i<rank;++i)    // select submatrix defined by key
        A[i]=gens[S->key[i]-1];
        
    vector<SIMPLINEXDATA_INT> inExSimplData;
    if(inExCollect.size()!=0)    
        prepare_inclusion_exclusion_simpl(*S,inExCollect,inExSimplData);

    h=0;
    long iS=S->offsets.size();    // compute numerator for simplex being processed   
    for(i=0;i<iS;++i){
        degree_b=scalProd(degrees,S->offsets[i]);
        degree_b/=det;
        h+=power(t,degree_b)*affineLinearSubstitutionFL(FF,A,S->offsets[i],det,RZZ,degrees,lcmDets,inExSimplData, facePolys);
    }
    
    evaluateClass=false; // necessary to evaluate class only once
    
    int tn;
    if(omp_get_level()==0)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(1);
        
    // #pragma omp critical (ADDTOCLASS) 
    { 
        denomClasses[S->classNr].second[tn]+=h;
        #pragma omp critical (ADDTOCLASS)
        {
        denomClasses[S->classNr].first.simplDone++;
        
        if(denomClasses[S->classNr].first.simplDone==denomClasses[S->classNr].first.simplDue)
            evaluateClass=true;
        }
    }
    if(evaluateClass)
    {
    
        for(int j=1;j<omp_get_max_threads();++j){
            denomClasses[S->classNr].second[0]+=denomClasses[S->classNr].second[j];
            denomClasses[S->classNr].second[j]=0;
        }        
            
        // denomClasses[S->classNr].second=0;  // <------------------------------------- 
        HClass=evaluateDenomClass(GFP,denomClasses[S->classNr]);
        #pragma omp critical(ACCUMULATE)
        {
            H.addCRF(HClass);
        }
        
    }
    
    if(!evaluationActive && facePolys[tn].size() >= 20){
        #pragma omp critical(FACEPOLYS)
        {
            evaluationActive=true;
            transferFacePolys(facePolys[tn],faceClasses);
            evaluationActive=false;
        }
     }
    
    #pragma omp critical(PROGRESS) // a little bit of progress report
    {
    if((++nrSimplDone)%10==0 && verbose_INT)
        cout << nrSimplDone << " simplicial cones done  " << endl; // nrActiveFaces-nrActiveFacesOld << " faces done" << endl;
        // nrActiveFacesOld=nrActiveFaces;
    }
 
  }  // Stanley dec
    
  } // parallel
  
  // collect the contribution of proper fases from inclusion/exclusion as far as not done yet
  
    for(int i=0;i<omp_get_max_threads();++i)
        transferFacePolys(facePolys[i],faceClasses);
  
  if(!faceClasses.empty())
    H.addCRF(evaluateFaceClasses(GFP,faceClasses));
    
    // now we must return to rational coefficients 
 
  CyclRatFunct HRat(zero(R));
  HRat.denom=H.denom;
  HRat.num=makeQQCoeff(H.num,R); 
   
  HRat.num*=FF.myRemainingFactor;
  HRat.num/=power(lcmDets,deg(F));
  
  HRat.showCoprimeCRF();
  
  mpz_class commonDen; // common denominator of coefficients of numerator of H  
  libnormaliz::HilbertSeries HS(nmzHilbertSeries(HRat,commonDen));
  writeGenEhrhartSeries(project, FFNonhom,HS,deg(F)+rank-1,commonDen);
}



  /* if(project=="ele2symm6")
        F= binomial(x[1]+5,5)*binomial(x[2]+5,5)*binomial(x[3]+2,2)*binomial(x[4]+2,2)
                  *binomial(x[5]+2,2)*binomial(x[6]+2,2);
  if(project=="ele1symm8")
        F=binomial(x[1]+5,5)*(x[2]+1)*(x[3]+1)*(x[4]+1)*(x[5]+1)*(x[6]+1)*(x[7]+1)
        *binomial(x[8]+5,5);
  if(project=="ele3symm13")
        F=binomial(x[1]+5,5)*(x[2]+1)*(x[4]+1)*(x[6]+1)*(x[8]+1)*(x[10]+1)*(x[12]+1);
  if(project=="ele3pr18sy9")
          F=binomial(x[1]+5,5)*(x[2]+1)*(x[4]+1)*(x[6]+1)*(x[8]+1);
  if(project=="ele3pr20sy10")
        F=binomial(x[1]+5,5)*(x[2]+1)*(x[4]+1)*(x[6]+1)*(x[8]+1)*(x[10]+1);*/
