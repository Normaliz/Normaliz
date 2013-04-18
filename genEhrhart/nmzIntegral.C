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

void writeIntegral(const string& project, const factorization<RingElem>& FF,
                   const RingElem& I, const bool& do_leadCoeff,
                   const long& virtDeg, const bool& appendOutput) {

    string name_open=project+".intOut";                              
    const char* file=name_open.c_str();
    ofstream out;
    if(appendOutput){
        out.open(file,ios_base::app);
        out << endl
            << "============================================================"
            << endl << endl;
    }
    else {
        out.open(file);
    }
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

void integrate(const string& project, const bool do_leadCoeff, bool& homogeneous, const bool& appendOutput) {
  GlobalManager CoCoAFoundations;

  if (verbose_INT) {
    cout << "============================================================" << endl;
    cout << "Integration for " << project << endl;
    cout << "============================================================" << endl << endl;
  }
  vector<long> grading;
  long gradingDenom;
  getGrading(project,grading,gradingDenom);

  vector<vector<long> > gens;
  readGens(project,gens);
  if(verbose_INT) 
    cout << "Generators read" << endl;
  long dim=gens[0].size();

  list<TRIDATA> triang;
  readTri(project,dim,triang);
  if(verbose_INT)
     cout << "Triangulation read" << endl;

  SparsePolyRing R=NewPolyRing(RingQQ(),dim+1,lex);
  // const RingElem& t=indets(R)[0];
  // const vector<RingElem>& x = indets(R);

  map<vector<long>,RingElem> denomClasses;

  RingElem F(one(R)); // to have something

  F=readPolynomial(project,R);
  if(verbose_INT)
    cout << "Polynomial read" << endl;
  
  if(do_leadCoeff){
    vector<RingElem> compsF= homogComps(F);
    // cout << "comps " << compsF << endl;
    if(F!=compsF[compsF.size()-1]){
        homogeneous=false;
        if(verbose_INT) 
            cout << "Polynomial is inhomogeneous. Replacing it by highest hom comp." << endl;
        F=compsF[compsF.size()-1];
        for(size_t i=0;i<compsF.size();++i) // no longer needed
            compsF[i]=0;     
    }
  }
  

  long i;

  vector<BigInt> Factorial(deg(F)+dim);
  for(i=0;i<deg(F)+dim;++i)
      Factorial[i]=factorial(i);

  factorization<RingElem> FF=factor(F);
  long nf=FF.myFactors.size();
  if(verbose_INT) 
    cout <<"Factorization" << endl;  // we show the factorization so that the user can check
  if(verbose_INT)
    for(i=0;i<nf;++i)
        cout << FF.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
  if(verbose_INT)
    cout << "Remaining factor " << FF.myRemainingFactor << endl << endl;

  RingElem I(zero(RingQQ())); // accumulates the integral

  size_t tri_size=triang.size();

  if(verbose_INT){
    cout << "********************************************" << endl;
    cout << tri_size <<" simplicial cones to be evaluated" << endl;
    cout << "********************************************" << endl;
  }

  size_t nrSimplDone=0;

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

    // h=homogeneousLinearSubstitutionFL(FF,A,degrees,F);
    ISimpl=(det*substituteAndIntegrate(FF,A,degrees,F,Factorial))/prodDeg;

    #pragma omp critical(INTEGRAL)
    I+=ISimpl;

    #pragma omp critical(PROGRESS) // a little bit of progress report
    {
    if ((++nrSimplDone)%10==0 && verbose_INT)
        cout << nrSimplDone << " simplicial cones done" << endl;
    }

  }  // triang

  } // parallel
  
  string result="Integral";
  if(do_leadCoeff)
    result="(Virtual) leading coefficient of quasipol";


   if(verbose_INT){
    cout << "********************************************" << endl;
    cout << result << " is " << endl << I << endl;
    cout << "********************************************" << endl;
   }
   
   if(do_leadCoeff && verbose_INT){
    cout << "Virtual multiplicity  is " << endl << I*factorial(deg(F)+dim-1) << endl;
    cout << "********************************************" << endl;
   }
   
   writeIntegral(project,FF,I,do_leadCoeff,deg(F)+dim-1,appendOutput);
}
