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

BigRat IntegralUnitSimpl(const RingElem& F, const vector<BigInt>& Factorial,const long& rank){

    SparsePolyRing P=AsSparsePolyRing(owner(F));
    vector<long> v(NumIndets(P));
    
    BigRat Irat;
    long deg;
    BigInt facProd,I,maxFact;
    I=0;
    maxFact=Factorial[Factorial.size()-1];

    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    for (; !IsEnded(mon); ++mon){
      exponents(v,PP(mon)); // this function gives the exponent vector back as v
      deg=0;
      IsInteger(facProd,coeff(mon)); // start with coefficient and multipliy by Factorials
      for(long i=1;i<=rank;++i){
          deg+=v[i];
          facProd*=Factorial[v[i]];
       }
       I+=facProd*maxFact/Factorial[deg+rank-1];
    }
    Irat=I;
    return(Irat/maxFact);
}

void writeIntegral(const string& project, const factorization<RingElem>& FF,
                   const BigRat& I, const bool& do_leadCoeff,
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

void integrate(const string& project, const bool& do_leadCoeff, bool& homogeneous, const bool& appendOutput) {
  GlobalManager CoCoAFoundations;
  
  long i;

  if (verbose_INT) {
    cout << "==========================================================" << endl;
    cout << "Integration for " << project << endl;
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

  list<TRIDATA> triang;
  readTri(project,rank,triang);
  if(verbose_INT)
     cout << "Triangulation read" << endl;
     
  list<TRIDATA>::iterator S = triang.begin(); // first we compute the lcm of all generator degrees
  vector<long> degrees(rank); 
  BigInt lcmDegs(1);
  vector<vector<long> > A(rank);
  for(;S!=triang.end();++S){          
    for(i=0;i<rank;++i)    // select submatrix defined by key
        A[i]=gens[S->key[i]-1]; 
    degrees=MxV(A,grading);
    for(i=0;i<rank;++i){
        degrees[i]/=gradingDenom;
        lcmDegs=lcm(lcmDegs,degrees[i]);
    }
  }

  SparsePolyRing R=NewPolyRing(RingQQ(),dim+1,lex);
  SparsePolyRing RZZ=NewPolyRing(RingZZ(),dim+1,lex);
  vector<RingElem> primeFactors;
  vector<RingElem> primeFactorsNonhom;
  vector<long> multiplicities;
  RingElem remainingFactor(one(R));
   
  RingElem F=processInputPolynomial(project,R,RZZ,primeFactors, primeFactorsNonhom,
                multiplicities,remainingFactor,homogeneous,do_leadCoeff);
                
  vector<BigInt> Factorial(deg(F)+dim); // precomputed values
  for(i=0;i<deg(F)+dim;++i)
      Factorial[i]=factorial(i);
  
  factorization<RingElem> FF(primeFactors,multiplicities,remainingFactor); // assembels the data
  factorization<RingElem> FFNonhom(primeFactorsNonhom,multiplicities,remainingFactor); // for output

  long nf=FF.myFactors.size();
  if(verbose_INT){
    cout <<"Factorization" << endl;  // we show the factorization so that the user can check
    for(i=0;i<nf;++i)
        cout << FFNonhom.myFactors[i] << "  mult " << FF.myExponents[i] << endl;
    cout << "Remaining factor " << FF.myRemainingFactor << endl << endl;
  }


  if(verbose_INT)
    cout << "Polynomial read" << endl;

  BigRat I; // accumulates the integral
  I=0;

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
  BigRat ISimpl; // integral over a simplex
  BigInt prodDeg; // product of the degrees of the generators
  RingElem h(zero(R));

  size_t spos=0,s;
  #pragma omp for schedule(dynamic) 
  for(s=0;s<tri_size;++s){         
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
    ISimpl=(det*substituteAndIntegrate(FF,A,degrees,RZZ,Factorial,lcmDegs))/prodDeg;

    #pragma omp critical(INTEGRAL)
    I+=ISimpl;

    #pragma omp critical(PROGRESS) // a little bit of progress report
    {
    if ((++nrSimplDone)%10==0 && verbose_INT)
        cout << nrSimplDone << " simplicial cones done" << endl;
    }

  }  // triang

  } // parallel
  
  I/=power(lcmDegs,deg(F));
  BigRat RFrat;
  IsRational(RFrat,remainingFactor); // from RingQQ to BigRat
  I*=RFrat;
  
  string result="Integral";
  if(do_leadCoeff)
    result="(Virtual) leading coefficient of quasipol";


   if(verbose_INT){
    cout << "********************************************" << endl;
    cout << result << " is " << endl << I << endl;
    cout << "********************************************" << endl;
   }
   
   if(do_leadCoeff && verbose_INT){
    cout << "Virtual multiplicity  is " << endl << I*factorial(deg(F)+rank-1) << endl;
    cout << "********************************************" << endl;
   }
   
   writeIntegral(project,FFNonhom,I,do_leadCoeff,deg(F)+rank-1,appendOutput);
}
