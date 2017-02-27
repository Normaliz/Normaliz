/*
 * nmzIntegrate
 * Copyright (C) 2012-2014  Winfried Bruns, Christof Soeger
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
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#include <fstream>
#include <sstream>
#include<string>

#include "libnormaliz/nmz_integrate.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/map_operations.h"

using namespace CoCoA;

#include <boost/dynamic_bitset.hpp>

#include "../libnormaliz/my_omp.h"

namespace libnormaliz {


BigRat IntegralUnitSimpl(const RingElem& F, const vector<BigInt>& Factorial,
                const vector<BigInt>& factQuot, const long& rank){
                
    SparsePolyRing P=owner(F);
    long dim=NumIndets(P);
    vector<long> v(dim);
    
    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    map<vector<long>,RingElem> orderedMons;  // will take the ordered exponent vectors
    map<vector<long>,RingElem>::iterator ord_mon;

    for (; !IsEnded(mon); ++mon){
      exponents(v,PP(mon)); // this function gives the exponent vector back as v
      sort(v.begin()+1,v.begin()+rank+1);
      ord_mon=orderedMons.find(v); // insert into map or add coefficient
      if(ord_mon!=orderedMons.end()){
          ord_mon->second+=coeff(mon);
      }
      else{
          orderedMons.insert(pair<vector<long>,RingElem>(v,coeff(mon)));
      }
    }


    long deg;
    BigInt facProd,I;
    I=0;
    for(ord_mon=orderedMons.begin();ord_mon!=orderedMons.end();++ord_mon){
      deg=0;
      v=ord_mon->first;
      IsInteger(facProd,ord_mon->second); // start with coefficient and multipliy by Factorials
      for(long i=1;i<=rank;++i){
          deg+=v[i];
          facProd*=Factorial[v[i]];
       }
       I+=facProd*factQuot[deg+rank-1];// maxFact/Factorial[deg+rank-1];
    }   
    
    BigRat Irat;
    Irat=I;
    return(Irat/Factorial[Factorial.size()-1]);            
}

BigRat substituteAndIntegrate(const ourFactorization& FF,const vector<vector<long> >& A,
                     const vector<long>& degrees, const SparsePolyRing& R, const vector<BigInt>& Factorial, 
                     const vector<BigInt>& factQuot,const BigInt& lcmDegs){
// we need F to define the ring
// applies linar substitution y --> y*(lcmDegs*A/degrees) to all factors in FF 
// where row A[i] is divided by degrees[i]
// After substitution the polynomial is integrated over the unit simplex
// and the integral is returned


    size_t i;
    size_t m=A.size();
    long rank=(long) m; // we prefer rank to be of type long
    vector<RingElem> v(m,zero(R));
    
    BigInt quot;
    for(i=0;i<m;i++){
        quot=lcmDegs/degrees[i];
        v[i]=indets(R)[i+1]*quot;
    }
    vector<RingElem> w=VxM(v,A);
    vector<RingElem> w1(w.size()+1,zero(R));
    w1[0]=RingElem(R,lcmDegs);
    for(i=1;i<w1.size();++i) // we have to shift w since the (i+1)st variable
        w1[i]=w[i-1];        // corresponds to coordinate i (counted from 0)
        
    
    RingHom phi=PolyAlgebraHom(R,R,w1);
    
    RingElem G1(zero(R));    
    list<RingElem> sortedFactors;
    for(i=0;i<FF.myFactors.size();++i){
        G1=phi(FF.myFactors[i]);
        for(int nn=0;nn<FF.myMultiplicities[i];++nn)         
                sortedFactors.push_back(G1);
    }
    
    list<RingElem>::iterator sf;
    sortedFactors.sort(compareLength);
    
    RingElem G(one(R));
    
    for(sf=sortedFactors.begin();sf!=sortedFactors.end();++sf)
        G*=*sf;

    // cout << "Evaluating integral over unit simplex" << endl;
    // boost::dynamic_bitset<> dummyInd;
    // vector<long> dummyDeg(degrees.size(),1);
    return(IntegralUnitSimpl(G,Factorial,factQuot,rank));  // orderExpos(G,dummyDeg,dummyInd,false)
}

template<typename Integer>
void readGens(Cone<Integer>& C, vector<vector<long> >& gens, const vector<long>& grading, bool check_ascending){
// get  from C for nmz_integrate functions

    size_t i,j;
    size_t nrows, ncols;
    nrows=C.getNrGenerators();
    ncols=C.getEmbeddingDim();
    gens.resize(nrows);
    for(i=0;i<nrows;++i)
        gens[i].resize(ncols);

    for(i=0; i<nrows; i++){
        for(j=0; j<ncols; j++) {
            convert(gens[i],C.getGenerators()[i]);
        }
        if(check_ascending){
            long degree,prevDegree=1;
            degree=v_scalar_product(gens[i],grading);
            if(degree<prevDegree){
                throw FatalException( " Degrees of generators not weakly ascending!");
            }
            prevDegree=degree;
        }
    }
}

template<typename Integer>
void readTri(Cone<Integer>& C, list<TRIDATA>& triang){
// get triangulation from C for nmz_integrate functions

    size_t rank=C.getRank();
    TRIDATA new_entry;
    new_entry.key.resize(rank);
    for(size_t i=0;i<C.getTriangulation().size();++i){
        for(size_t j=0;j<rank;++j)
            new_entry.key[j]=C.getTriangulation()[i].first[j]+1;
        convert(new_entry.vol,C.getTriangulation()[i].second);
        triang.push_back(new_entry);
    }
}


template<typename Integer>
void integrate(Cone<Integer>& C, const bool do_leadCoeff) {
  GlobalManager CoCoAFoundations;
 
  bool verbose_INTsave=verbose_INT;
  verbose_INT=C.get_verbose();
  
  long i;

  if (verbose_INT) {
    verboseOutput() << "==========================================================" << endl;
    verboseOutput() << "Integration" << endl;
    verboseOutput() << "==========================================================" << endl << endl;
  }
  vector<long> grading;
  convert(grading,C.getGrading());
  long gradingDenom;
  convert(gradingDenom,C.getGradingDenom());
  long rank=C.getRank();

  vector<vector<long> > gens;
  readGens(C,gens,grading,false);
  if(verbose_INT) 
    verboseOutput() << "Generators read" << endl;
  long dim=C.getEmbeddingDim();

  list<TRIDATA> triang;
  readTri(C,triang);
  if(verbose_INT)
     verboseOutput() << "Triangulation read" << endl;
     
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

  SparsePolyRing R=NewPolyRing_DMPI(RingQQ(),dim+1,lex);
  SparsePolyRing RZZ=NewPolyRing_DMPI(RingZZ(),PPM(R)); // same indets and ordering as R
  vector<RingElem> primeFactors;
  vector<RingElem> primeFactorsNonhom;
  vector<long> multiplicities;
  RingElem remainingFactor(one(R));

  bool homogeneous;
  RingElem F=processInputPolynomial(C.getPolynomial(),R,RZZ,primeFactors, primeFactorsNonhom,
                multiplicities,remainingFactor,homogeneous,do_leadCoeff);
                
  vector<BigInt> Factorial(deg(F)+dim); // precomputed values
  for(i=0;i<deg(F)+dim;++i)
      Factorial[i]=factorial(i);
      
  vector<BigInt> factQuot(deg(F)+dim); // precomputed values
  for(i=0;i<deg(F)+dim;++i)
      factQuot[i]=Factorial[Factorial.size()-1]/Factorial[i];
  
  ourFactorization FF(primeFactors,multiplicities,remainingFactor); // assembels the data
  ourFactorization FFNonhom(primeFactorsNonhom,multiplicities,remainingFactor); // for output

  long nf=FF.myFactors.size();
  if(verbose_INT){
    verboseOutput() <<"Factorization" << endl;  // we show the factorization so that the user can check
    for(i=0;i<nf;++i)
        verboseOutput() << FFNonhom.myFactors[i] << "  mult " << FF.myMultiplicities[i] << endl;
    verboseOutput() << "Remaining factor " << FF.myRemainingFactor << endl << endl;
  }


  if(verbose_INT)
    verboseOutput() << "Polynomial read" << endl;

  BigRat I; // accumulates the integral
  I=0;

  size_t tri_size=triang.size();

  if(verbose_INT){
    verboseOutput() << "********************************************" << endl;
    verboseOutput() << tri_size <<" simplicial cones to be evaluated" << endl;
    verboseOutput() << "********************************************" << endl;
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
    ISimpl=(det*substituteAndIntegrate(FF,A,degrees,RZZ,Factorial,factQuot,lcmDegs))/prodDeg;

#pragma omp critical(INTEGRAL)
    I+=ISimpl;

#pragma omp critical(PROGRESS) // a little bit of progress report
    {
    if ((++nrSimplDone)%10==0 && verbose_INT)
        verboseOutput() << nrSimplDone << " simplicial cones done" << endl;
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
  
  BigRat VM;

  if(do_leadCoeff){
    VM=I*factorial(deg(F)+rank-1);
    C.getIntData().setVirtualMultiplicity(mpq(VM));
    C.getIntData().setLeadCoef(mpq(I));
  }
  else
    C.getIntData().setIntegral(mpq(I));

   if(verbose_INT){
    verboseOutput() << "********************************************" << endl;
    verboseOutput() << result << " is " << endl << I << endl;
    verboseOutput() << "********************************************" << endl;
   }
   
   if(do_leadCoeff && verbose_INT){
    verboseOutput() << "Virtual multiplicity  is " << endl << VM << endl;
    verboseOutput() << "********************************************" << endl;
   }  
   
    verbose_INT=verbose_INTsave; 
}

CyclRatFunct evaluateFaceClasses(const vector<vector<CyclRatFunct> >& GFP,
                                    map<vector<long>,RingElem>& faceClasses){
// computes the generating rational functions
// for the denominator classes collected from proper faces and returns the sum

    SparsePolyRing R=owner(faceClasses.begin()->second);
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

    SparsePolyRing R=owner(denomClass.second[0]);
    
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

template<typename Integer>
void readInEx(Cone<Integer>& C, vector<pair<boost::dynamic_bitset<>, long> >& inExCollect, const size_t nrGen){

    size_t inExSize=C.getInclusionExclusionData().size(), keySize;   
    long mult;
    boost::dynamic_bitset<> indicator(nrGen);
    for(size_t i=0;i<inExSize;++i){
        keySize=C.getInclusionExclusionData()[i].first.size();
        indicator.reset();
        for(size_t j=0;j<keySize;++j){
            indicator.set(C.getInclusionExclusionData()[i].first[j]);
        }
        mult=C.getInclusionExclusionData()[i].second;
        inExCollect.push_back(pair<boost::dynamic_bitset<>, long>(indicator,mult));       
    }
}

template<typename Integer>
void readDecInEx(Cone<Integer>& C, const long& dim, list<STANLEYDATA_INT>& StanleyDec,
                vector<pair<boost::dynamic_bitset<>, long> >& inExCollect, const size_t nrGen){
// rads Stanley decomposition and InExSata from C
    
    if(C.isComputed(ConeProperty::InclusionExclusionData)){
        readInEx(C, inExCollect,nrGen);
    }

    STANLEYDATA_INT newSimpl;
    long i=0,j,det;
    newSimpl.key.resize(dim);
    
    long test;
    
    auto SD=C.getStanleyDec().begin();

    for(;SD!=C.getStanleyDec().end();++SD){
 
        test=0;
        for(i=0;i<dim;++i){
            newSimpl.key[i]=SD->key[i]+1;
            if(newSimpl.key[i]<=test){
                throw FatalException("Key of simplicial cone not ascending or out of range");
            }
            test=newSimpl.key[i];
        }
        
        det=SD->offsets.nr_of_rows();
        newSimpl.offsets.resize(det);
        for(i=0;i<det;++i)
            newSimpl.offsets[i].resize(dim);
        for(i=0;i<det;++i)
            for(j=0;j<dim;++j)
                convert(newSimpl.offsets[i][j],SD->offsets[i][j]);
        StanleyDec.push_back(newSimpl);
    }    
}

template<typename Integer>
void generalizedEhrhartSeries(Cone<Integer>& C){
  GlobalManager CoCoAFoundations;

  bool verbose_INTsave=verbose_INT;
  verbose_INT=C.get_verbose();
  
  if(verbose_INT){
    cout << "==========================================================" << endl;
    cout << "Generalized Ehrhart series " << endl;
    cout << "==========================================================" << endl << endl;
  }
  
  long i,j;
  
  vector<long> grading;
  convert(grading,C.getGrading());
  long gradingDenom;
  convert(gradingDenom,C.getGradingDenom());
  long rank=C.getRank();
  
  vector<vector<long> > gens;
  readGens(C,gens,grading,true);
  if(verbose_INT)
    cout << "Generators read" << endl;
  long dim=C.getEmbeddingDim();
  long maxDegGen=v_scalar_product(gens[gens.size()-1],grading)/gradingDenom; 
  
  list<STANLEYDATA_INT> StanleyDec;
  vector<pair<boost::dynamic_bitset<>, long> > inExCollect;
  readDecInEx(C,rank,StanleyDec,inExCollect,gens.size());
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
  
  SparsePolyRing R=NewPolyRing_DMPI(RingQQ(),dim+1,lex);
  SparsePolyRing RZZ=NewPolyRing_DMPI(RingZZ(),PPM(R)); // same indets and ordering as R
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

  bool homogeneous;
  RingElem F=processInputPolynomial(C.getPolynomial(),R,RZZ,primeFactors, primeFactorsNonhom,
                multiplicities,remainingFactor,homogeneous,false);
                
  vector<BigInt> Factorial(deg(F)+dim); // precomputed values
  for(i=0;i<deg(F)+dim;++i)
      Factorial[i]=factorial(i);
  
  ourFactorization FF(primeFactors,multiplicities,remainingFactor); // assembeles the data
  ourFactorization FFNonhom(primeFactorsNonhom,multiplicities,remainingFactor); // for output

  long nf=FF.myFactors.size();
  if(verbose_INT){
    cout <<"Factorization" << endl;  // we show the factorization so that the user can check
    for(i=0;i<nf;++i)
        cout << FFNonhom.myFactors[i] << "  mult " << FF.myMultiplicities[i] << endl;
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
        degree_b=v_scalar_product(degrees,S->offsets[i]);
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
  
  /* string outputName;
  if(pnm==pureName(project))
    outputName=project;
  else
    outputName=project+"."+pnm;
  
  if(output_dir!="")
      outputName=output_dir+pureName(outputName);
    
  writeGenEhrhartSeries(outputName, FFNonhom,HS,deg(F)+rank-1,commonDen); */
  
     
   verbose_INT=verbose_INTsave; 
}


#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template void integrate(Cone<long>& C, const bool do_leadCoeff);
#endif // NMZ_MIC_OFFLOAD
template void integrate(Cone<long long>& C, const bool do_leadCoeff);
template void integrate(Cone<mpz_class>& C, const bool do_leadCoeff);

#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template void generalizedEhrhartSeries<long>(Cone<long>& C);
#endif // NMZ_MIC_OFFLOAD
template void generalizedEhrhartSeries<long long>(Cone<long long>& C);
template void generalizedEhrhartSeries<mpz_class>(Cone<mpz_class>& C);



} // namespace libnormaliz
