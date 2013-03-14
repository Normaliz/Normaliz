

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

RingElem orderExpos(const RingElem& F, const vector<long>& degrees){
 // orders the exponent vectors v of the terms of F
 // the exponents v[i] and v[j], i < j,  are swapped if
 // (1) degrees[i]==degrees[j] and (2) v[i] < v[j]
 // so that the exponents are descending in each degree block
 // the ordered exponent vectors are inserted into a map
 // and their coefficients are added
 // at the end the polynomial is rebuilt from the map

    SparsePolyRing P=AsSparsePolyRing(owner(F));
    vector<long> v(NumIndets(P));

    vector<long> denom=degrees2denom(degrees); // first we must find the blocks of equal degree
    vector<long> St,End;
    St.push_back(1);
    for(size_t i=0;i<denom.size();++i)
        if(denom[i]!=0){
            End.push_back(St[St.size()-1]+denom[i]-1);
            if(i<denom.size()-1)
                St.push_back(End[End.size()-1]+1);
        }

    // now the main job

    map<vector<long>,RingElem> orderedMons;  // will take the ordered exponent vectors
    map<vector<long>,RingElem>::iterator ord_mon;
    long p,s,pend,pst;
    bool ordered;
    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    for (; !IsEnded(mon); ++mon){
      exponents(v,PP(mon)); // this function gives the exponent vector back as v
      for(size_t j=0;j<St.size();++j){  // now we go over the blocks
        pst=St[j];
        pend=End[j];
        while(1){
                ordered=true;
                for(p=pst;p<pend;++p){
                    if(v[p]<v[p+1]){
                        ordered=false;
                        s=v[p];
                        v[p]=v[p+1];
                        v[p+1]=s;
                    }
                }
            if(ordered)
                break;
            pend--;
        }
      }
        ord_mon=orderedMons.find(v); // insert into map or add coefficient
        if(ord_mon!=orderedMons.end()){
            ord_mon->second+=coeff(mon);
        }
        else{
            orderedMons.insert(pair<vector<long>,RingElem>(v,coeff(mon)));
        }
    }

    // now we must reconstruct the polynomial
    // we use that the natural order of vectors in C++ STL is inverse
    // to lex. Therefore push_front

    RingElem r(zero(P));
    for(ord_mon=orderedMons.begin();ord_mon!=orderedMons.end();++ord_mon){
        PushFront(r,ord_mon->second,ord_mon->first);
    }
    return(r);
}

CyclRatFunct genFunctPower1(const SparsePolyRing& P, long k,long n)
// computes the generating function for
//  \sum_j j^n (t^k)^j
{
    vector<RingElem> a=power2ascFact(P,n);
    // cout << a << endl;
    RingElem b(P);
    vector<long> u;
    CyclRatFunct g(zero(P)), h(zero(P));
    long i,s=a.size();
    for(i=0;i<s;++i)
    {
        u=makeDenom(k,i+1);
        /* cout << "k " << k << " u ";
        for(size_t j=0;j<u.size();++j)
            cout << u[j]<< " "; 
        cout << endl;*/
        b=a[i]*factorial(i);
        g.set2(b,u); 
        // cout << "g inner "; g.showCRF(); cout << endl << flush;  
        h.addCRF(g);
    }
    return(h);
}


CyclRatFunct genFunct(const vector<vector<CyclRatFunct> >& GFP, const RingElem& F, const vector<long>& degrees)
// writes \sum_{x\in\ZZ_+^n} f(x,t) T^x
// under the specialization T_i --> t^g_i
// as a rational function in t
{
    // cout << "Start genFunct" << endl;
    SparsePolyRing P=AsSparsePolyRing(owner(F));
    RingElem t=indets(P)[0];
    
    CyclRatFunct s(F); // F/1
    
    CyclRatFunct g(zero(P)),h(zero(P));
    
    long nd=degrees.size();  
    long i,k,mg;
    vector<RingElem> c;   

    for(k=1; k<=nd;k++)
    {
        // cout << "k " << k<< endl;
        c=ourCoeffs(s.num,k); // we split the numerator according 
                               // to powers of var k
        mg=c.size(); // max degree+1 in  var k

        h.set2(zero(P));
        for(i=0;i<mg;i++)     // now we replace the powers of var k
        {                      // by the corrseponding rational function,
                               // multiply, and sum the products
            // g=GFP[degrees[k-1]][i];
            // g=genFunctPower1(P,degrees[k-1],i);
            // cout << "g "; g.showCRF(); cout << endl;
            // g.num*=c[i];
            h.num=(1-power(t,degrees[k-1]))*h.num+GFP[degrees[k-1]][i].num*c[i];
            h.denom=GFP[degrees[k-1]][i].denom;
        }
        s.num=h.num;
        // cout << "genFunct NumTerms " << NumTerms(s.num) << endl;
        s.denom=prodDenom(s.denom,h.denom);
    }
    return(s);   
}

struct denomClassData{
    vector<long> degrees;
    size_t simplDue;
    size_t simplDone;
  };

CyclRatFunct evaluateDenomClass(const vector<vector<CyclRatFunct> >& GFP,
                                    pair<denomClassData,RingElem>& denomClass){
// computes the generating rational function
// for a denominator class and return it

    SparsePolyRing R=AsSparsePolyRing(owner(denomClass.second));
    
    if(verbose_INT){
    #pragma omp critical(PROGRESS)
    {
        // cout << "--------------------------------------------" << endl;
        cout << "Evaluating denom class ";
        for(size_t i=0;i<denomClass.first.degrees.size();++i)
            cout << denomClass.first.degrees[i] << " ";
        cout  << "NumTerms " << NumTerms(denomClass.second) << endl;
        // cout << "--------------------------------------------" << endl;
    }
    }

    CyclRatFunct h(zero(R));
    h = genFunct(GFP,denomClass.second,denomClass.first.degrees);

    denomClass.second=0;  // to save memory
    h.simplifyCRF();
    return(h);
}

mpz_class ourFactorial(const long& n){
    mpz_class fact=1;
    for(long i=1;i<=n;++i)
        fact*=i;
    return(fact);
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
    HCoeff3[i]=mpz_class(mpzref(HC2));      // to mpz_class 
  }

  vector<long> denomDeg=denom2degrees(H.denom);
  libnormaliz::HilbertSeries HS(HCoeff3,count_in_map<long, long>(denomDeg)); 
  HS.simplify();
  return(HS);
}

bool compareDegrees(STANLEYDATA_INT& A, STANLEYDATA_INT& B){

    return(A.degrees < B.degrees);
}

void generalizedEhrhartSeries(const string& project, bool& homogeneous){
  GlobalManager CoCoAFoundations;
  
  long i,j;
  
  if(verbose_INT){
    cout << "Generalized Ehrhart series " << project << endl;
    cout << "=======================================" << endl << endl;
  }

  vector<long> grading;
  long gradingDenom;
  getGrading(project,grading,gradingDenom);
  
  vector<vector<long> > gens;
  readGens(project,gens);
  if(verbose_INT)
    cout << "Generators read" << endl;
  long dim=gens[0].size();
  long maxDegGen=scalProd(gens[gens.size()-1],grading)/gradingDenom; 
  
  SparsePolyRing R=NewPolyRing(RingQQ(),dim+1,lex);
  
  list<STANLEYDATA_INT> StanleyDec;
  readDec(project,dim,StanleyDec);
  if(verbose_INT)
    cout << "Stanley decomposition read" << endl;
    
  size_t dec_size=StanleyDec.size();
    
  // Now we sort the Stanley decomposition by denominator class (= degree class)

  list<STANLEYDATA_INT>::iterator S = StanleyDec.begin();
  
  long rank=S->key.size();
  vector<long> degrees(rank);
  vector<vector<long> > A(rank);
  
  // prepare sorting by computing degrees of generators

  for(;S!=StanleyDec.end();++S){
      for(i=0;i<rank;++i)    // select submatrix defined by key
        A[i]=gens[S->key[i]-1];
          degrees=MxV(A,grading);
      for(i=0;i<rank;++i)
        degrees[i]/=gradingDenom; // must be divisible
      S->degrees=degrees;
  }
  
  StanleyDec.sort(compareDegrees);


  if(verbose_INT)
    cout << "Stanley decomposition sorted" << endl; 

  vector<pair<denomClassData,RingElem> > denomClasses;
  denomClassData denomClass;

  // we now make class 0 to get started
  S=StanleyDec.begin();
  denomClass.degrees=S->degrees;  // put degrees in class
  denomClass.simplDone=0;
  denomClass.simplDue=1;           // already one simplex to be done 
  denomClasses.push_back(pair<denomClassData,RingElem>(denomClass,zero(R)));
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
        denomClasses.push_back(pair<denomClassData,RingElem>(denomClass,zero(R)));
        dc++;
        S->classNr=dc;
    }
  }

  if(verbose_INT)
    cout << denomClasses.size() << " denominator classes built" << endl;

  const RingElem& t=indets(R)[0];
  // const vector<RingElem>& x = indets(R);


  RingElem F(one(R)); // to have something

  F=readPolynomial(project,R);        
  if(verbose_INT)
    cout << "Polynomial read" << endl;

  vector<RingElem> compsF= homogComps(F);
  homogeneous=true;
  if(F!=compsF[compsF.size()-1])
    homogeneous=false;
  for(int i=0;i< (long) compsF.size();++i) // no longer needed
    compsF[i]=0;

  factorization<RingElem> FF=factor(F);
  long nf=FF.myFactors.size();
  if(verbose_INT){
    cout <<"Factorization" << endl;  // we show the factorization so that the user can check
    for(i=0;i<nf;++i)
        cout << FF.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
    cout << "Remaining factor " << FF.myRemainingFactor << endl << endl;
  }

  vector<vector<CyclRatFunct> > GFP; // we calculate the table of generating functions
  vector<CyclRatFunct> DummyCRFVect; // for\sum i^n t^ki vor various values of k and n
  CyclRatFunct DummyCRF(zero(R));
  for(j=0;j<=deg(F);++j)
    DummyCRFVect.push_back(DummyCRF);
  for(i=0;i<=maxDegGen;++i){
    GFP.push_back(DummyCRFVect);
    for(j=0;j<=deg(F);++j)
        GFP[i][j]=genFunctPower1(R,i,j);
  }

  CyclRatFunct H(zero(R)); // accumulates the series
  
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

  RingElem h(zero(R));     // for use in a simplex
  CyclRatFunct HClass(zero(R)); // for single class
  

  size_t s,spos=0;  
  #pragma omp for schedule(dynamic) 
  for(s=0;s<dec_size;++s){
    for(;spos<s;++spos,++S);
    for(;spos>s;--spos,--S);

    det=S->offsets.size();
    degrees=S->degrees;
    
    for(i=0;i<rank;++i)    // select submatrix defined by key
        A[i]=gens[S->key[i]-1];

    h=0;
    long iS=S->offsets.size();    // compute numerator for simplex being processed   
    for(i=0;i<iS;++i){
        degree_b=scalProd(degrees,S->offsets[i]);
        degree_b/=det;
        h+=power(t,degree_b)*affineLinearSubstitutionFL(FF,A,S->offsets[i],det,F,degrees);
    }
    // h=orderExpos(h,degrees);
    // cout << "Simpl " << s+1 << " NumTerms " << NumTerms(h) << endl;
    
    evaluateClass=false; // necessary to evaluate class only once
    #pragma omp critical (ADDTOCLASS) 
    { 
        denomClasses[S->classNr].second+=h;
        denomClasses[S->classNr].first.simplDone++;
        
        if(denomClasses[S->classNr].first.simplDone==denomClasses[S->classNr].first.simplDue)
            evaluateClass=true;
    }
    if(evaluateClass)
    {
        HClass=evaluateDenomClass(GFP,denomClasses[S->classNr]);
        #pragma omp critical(ACCUMULATE)
        {
            H.addCRF(HClass);
        }
        
    }
    
    #pragma omp critical(PROGRESS) // a little bit of progress report
    {
    if((++nrSimplDone)%10==0 && verbose_INT)
        cout << nrSimplDone << " simplicial cones done" << endl;
    }
 
  }  // Stanley dec
    
  } // parallel
  
  
  
  H.showCoprimeCRF();
  
  mpz_class commonDen; // common denominator of coefficients of numerator of H  
  libnormaliz::HilbertSeries HS(nmzHilbertSeries(H,commonDen));
  writeGenEhrhartSeries(project, FF,HS,deg(F)+dim-1,commonDen);
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
