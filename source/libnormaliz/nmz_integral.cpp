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

void writeIntegral(const string& project, const ourFactorization& FF,
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
        out << FF.myFactors[i] << "  mult " << FF.myMultiplicities[i] << endl;
    out << "Remaining factor " << FF.myRemainingFactor << endl << endl;
    
    if(do_leadCoeff){
        out << "Virtual leading coefficient: " << I << endl;
        out << endl << "Virtual multiplicity: " << I*factorial(virtDeg) << endl;
    }
    else
       out << "Integral: " << I << endl;
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
void integrate(Cone<Integer>& C, const bool do_leadCoeff, bool& homogeneous) {
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
  long dim=gens[0].size();

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

  C.setIntegral(mpq(I));
  
  BigRat VM;

  if(do_leadCoeff){
    VM=I*factorial(deg(F)+rank-1);
    C.setVirtualMultiplicity(mpq(VM));
    C.setLeadCoef(mpq(I));
  }

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

#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template void integrate(Cone<long>& C, const bool do_leadCoeff, bool& homogeneous);
#endif // NMZ_MIC_OFFLOAD
template void integrate(Cone<long long>& C, const bool do_leadCoeff, bool& homogeneous);
template void integrate(Cone<mpz_class>& C, const bool do_leadCoeff, bool& homogeneous);

} // namespace libnormaliz
