RingElem IntegralUnitSimpl(RingElem F, vector<BigInt> Factorial){


    SparsePolyRing P=AsSparsePolyRing(owner(F));
    vector<long> v(NumIndets(P));
    long dim=v.size()-1;

    RingElem I(zero(RingQQ()));
    long deg;
    // RingElem facProd(one(RingQQ()));
    BigInt facProd;
    // RingElem dummy(zero(RingQQ()));
    // vector<long> ONE(NumIndets(P),0);

    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    for (; !IsEnded(mon); ++mon){
      exponents(v,PP(mon)); // this function gives the exponent vector back as v
      deg=0;
      facProd=1;
      for(long i=1;i<=dim;++i){
          deg+=v[i];
          facProd*=Factorial[v[i]];
       }
       // dummy=facProd*monomial(P,coeff(mon),ONE);
       I+=facProd*coeff(mon)/Factorial[deg+dim-1];
    }

    return(I);
}

void  writeIntegral(const string&  project,const factorization<RingElem>& FF,
                             const RingElem& I,const bool& do_leadCoeff, const long& virtDeg){
                             
    string name_open=project+".intOut";                              
    const char* file=name_open.c_str();
    ofstream out(file);
    
    out <<"Factorization of polynomial:" << endl;  // we show the factorization so that the user can check
    for(size_t i=0;i<FF.myFactors.size();++i)
        out << FF.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
    out << "Remaining factor " << FF.myRemainingFactor << endl << endl;
    
    if(do_leadCoeff){
        out << "Virtual leading coefficient: " << I << endl;
        out << endl << "Virtual multiplicity: " << I*factorial(virtDeg) << endl;
    }
    else
       out << "Integral: " << I << endl;
}

void integrate(const string& project, const bool do_leadCoeff){
  GlobalManager CoCoAFoundations;

  cout << "Integration for " << project << endl;
  cout << "=======================================" << endl << endl;;

  vector<long> grading;
  long gradingDenom;
  getGrading(project,grading,gradingDenom);

  vector<vector<long> > gens;
  readGens(project,gens);
  cout << "Generators read" << endl;
  long dim=gens[0].size();

  list<TRIDATA> triang;
  readTri(project,dim,triang);
  cout << "Triangulation read" << endl;

  SparsePolyRing R=NewPolyRing(RingQQ(),dim+1,lex);
  // const RingElem& t=indets(R)[0];
  // const vector<RingElem>& x = indets(R);

  map<vector<long>,RingElem> denomClasses;

  RingElem F(one(R)); // to have something

  F=readPolynomial(project,R);
  cout << "Polynomial read" << endl;
  
  if(do_leadCoeff){
    vector<RingElem> compsF= homogComps(F);
    // cout << "comps " << compsF << endl;
    if(F!=compsF[compsF.size()-1]){
        cout << "Polynomial is inhomogeneous. Replacing it by highest hom comp." << flush << endl;
        F=compsF[compsF.size()-1];
    }
  }
  

  long i;

  vector<BigInt> Factorial(deg(F)+dim);
  for(i=0;i<deg(F)+dim;++i)
      Factorial[i]=factorial(i);

  factorization<RingElem> FF=factor(F);
  long nf=FF.myFactors.size();
  cout <<"Factorization" << endl;  // we show the factorization so that the user can check
  for(i=0;i<nf;++i)
        cout << FF.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
  cout << "Remaining factor " << FF.myRemainingFactor << endl << endl;

  RingElem I(zero(RingQQ())); // accumulates the integral

  size_t tri_size=triang.size();

  cout << "********************************************" << endl;
  cout << tri_size <<" simplicial cones to be evaluated" << endl;
  cout << "********************************************" << endl;

  #pragma omp parallel private(i)
  {

  list<TRIDATA>::iterator S = triang.begin();
  long det, rank=S->key.size();
  vector<long> degrees(rank);
  vector<vector<long> > A(rank);
  RingElem ISimpl(RingQQ()); // integral over a simplex
  BigInt prodDeg; // product of the degrees of the generators
  RingElem h(zero(R));

  size_t spos=0,s;
  #pragma omp for schedule(dynamic) // parallelization will be interrupted as soon as
  for(s=0;s<tri_size;++s){          // 50 denominator classes have run up
      for(;spos<s;++spos,++S);
      for(;spos>s;--spos,--S);

    det=S->vol;
    for(i=0;i<rank;++i)    // select submatrix defined by key
        A[i]=gens[S->key[i]-1]; // will not be changed

    degrees=MxV(A,grading);
    prodDeg=1;
    for(i=0;i<rank;++i){
        degrees[i]/=gradingDenom;
        prodDeg*=degrees[i];
    }

    h=homogeneousLinearSubstitutionFL(FF,A,degrees,F);
    ISimpl=(det*IntegralUnitSimpl(h,Factorial))/prodDeg;

    #pragma omp critical(INTEGRAL)
    I+=ISimpl;

    #pragma omp critical(PROGRESS) // a little bit of progress report
    {
    if((s+1)%10==0)
        cout << "Simpl " << s+1 << " done" << endl << flush;
    }

  }  // triang

  } // parallel
  
  string result="Integral";
  if(do_leadCoeff)
    result="(Virtual) leading coefficient of quasipol";


   cout << "********************************************" << endl;
   cout << result << " is " << endl << I << endl;
   cout << "********************************************" << endl;
   
   if(do_leadCoeff){
    cout << "Virtual multiplicity  is " << endl << I*factorial(deg(F)+dim-1) << endl;
    cout << "********************************************" << endl;
   }
   
   writeIntegral(project,FF,I,do_leadCoeff,deg(F)+dim-1);
}
