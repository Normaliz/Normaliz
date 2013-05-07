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

RingElem binomial(const RingElem& f, long k)
// computes binomial coefficient (f choose k)
{
    SparsePolyRing P=AsSparsePolyRing(owner(f));
    RingElem g(P);
    g=1;
    for(int i=0;i<k;i++)
        g*=(f-i)/(i+1);
    return(g);
}

RingElem ascFact(const RingElem& f, long k)
// computes (f+1)*...*(f+k)
{
    SparsePolyRing P=AsSparsePolyRing(owner(f));
    RingElem g(P);
    g=1;
    for(int i=0;i<k;i++)
        g*=(f+i+1);
    return(g);
}

RingElem descFact(const RingElem& f, long k)
// computes f*(f-1)*...*(f-k+1)
{
    SparsePolyRing P=AsSparsePolyRing(owner(f));
    RingElem g(P);
    g=1;
    for(int i=0;i<k;i++)
        g*=(f-i);
    return(g);
}

vector<RingElem> ourCoeffs(const RingElem& F, const long j){
// our version of expanding a poly nomial wrt to indeterminate j
// The return value is the vector of coefficients of x[j]^i
    vector<RingElem> c;
    SparsePolyRing P=AsSparsePolyRing(owner(F));
    RingElem x=indets(P)[j];
    if(F==0){
        c.push_back(zero(P));
        return(c);
    }

    vector<long> v(NumIndets(P));
    long k,cs;

    SparsePolyIter i=BeginIter(F);
    for (; !IsEnded(i); ++i){
        exponents(v,PP(i));
        k=v[j];
        cs=c.size();
        if(k>cs-1)
            c.resize(k+1,zero(P));
        v[j]=0;
        // c[k]+=monomial(P,coeff(i),v);
        PushBack(c[k],coeff(i),v);
    }
    return(c);
}

vector<long> MxV(const vector<vector<long> >& M, vector<long> V){
// matrix*vector
    vector<long> P(M.size());
    for(size_t i=0;i< M.size();++i){
        long s=0;
        for(size_t j=0;j<V.size();++j)
            s+=M[i][j]*V[j];
        P[i]=s;
    }
    return(P);
}

long scalProd(const vector<long>& a, const vector<long>& b){
    long s=0;
    for(size_t i=0;i<a.size();++i)
        s+=a[i]*b[i];
    return(s);
}

vector<RingElem> VxM(const vector<RingElem>& V, const vector<vector<long> >& M){
// vector*matrix
    SparsePolyRing R=AsSparsePolyRing(owner(V[0]));
    RingElem s(zero(R));
    vector<RingElem> P(M[0].size(),zero(R));
    for(size_t j=0;j<M[0].size();++j){
        s=0;
        for(size_t i=0;i<M.size();++i)
            s+=V[i]*M[i][j];
        P[j]=s;
    }
    return(P);
}


/*
RingElem affineLinearSubstitution(const RingElem& F,const vector<vector<long> >& A,
                     const vector<long>& b, const long& denom){
// NOT IN USE
    size_t i;
    SparsePolyRing R=AsSparsePolyRing(owner(F));
    size_t m=A.size();
    // long n=A[0].size();
    vector<RingElem> v(m,zero(R));
    RingElem g(zero(R));
    
    for(i=0;i<m;i++)
    {
        g=b[i];
        g=g/denom;
        v[i]=g+indets(R)[i+1];
    }
    vector<RingElem> w=VxM(v,A);
    vector<RingElem> w1(w.size()+1,zero(R));
    w1[0]=indets(R)[0];
    for(i=1;i<w1.size();++i)
        w1[i]=w[i-1];
    
    RingHom phi=PolyAlgebraHom(R,R,w1);
    RingElem G=phi(F);  
    return(G);   
}
*/

RingElem orderExpos(const RingElem& F, const vector<long>& degrees);
BigRat IntegralUnitSimpl(const RingElem& F, const vector<BigInt>& Factorial, const long& dim);



RingElem affineLinearSubstitutionFL(const factorization<RingElem>& FF,const vector<vector<long> >& A,
                     const vector<long>& b, const long& denom, const SparsePolyRing& R, const vector<long>& degrees, const BigInt& lcmDets){
// applies linar substitution y --> lcmDets*A(y+b/denom) to all factors in FF 
// and returns the product of the modified factorsafter ordering the exponent vectors

    size_t i;
    size_t m=A.size();
    vector<RingElem> v(m,zero(R));
    RingElem g(zero(R));
    
    for(i=0;i<m;i++)
    {
        g=b[i]*(lcmDets/denom);
        v[i]=g+lcmDets*indets(R)[i+1];
    }
    vector<RingElem> w=VxM(v,A);
    vector<RingElem> w1(w.size()+1,zero(R));
    w1[0]=one(R);
    w1[0]*=lcmDets;
    for(i=1;i<w1.size();++i) 
        w1[i]=w[i-1]; 
    
    RingHom phi=PolyAlgebraHom(R,R,w1);
    RingElem G(one(R));
    for(i=0;i<FF.myFactors.size();++i){
        if(FF.myExponents[i]==1)
            G*=phi(FF.myFactors[i]);
        else           
            G*=power(phi(FF.myFactors[i]),FF.myExponents[i]);
    }
    
    return(orderExpos(G,degrees));   
}

BigRat substituteAndIntegrate(const factorization<RingElem>& FF,const vector<vector<long> >& A,
                     const vector<long>& degrees, const SparsePolyRing& R, const vector<BigInt>& Factorial, 
                     const BigInt& lcmDegs){
// we need F to define the ring
// applies linar substitution y --> lcmDegs*(A/degrees)y to all factors in FF 
// where row A[i] is divided by degrees[i]
// and returns the product of the modified factors


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
    w1[0]=lcmDegs*one(R);
    for(i=1;i<w1.size();++i) // we have to shift w since the zeroth variable
        w1[i]=w[i-1];        // must not be touched
        
    
    RingHom phi=PolyAlgebraHom(R,R,w1);
    
    RingElem G(one(R));
    for(i=0;i<FF.myFactors.size();++i){
        G*=power(phi(FF.myFactors[i]),FF.myExponents[i]);
    }
    return(IntegralUnitSimpl(orderExpos(G,degrees),Factorial,rank));
}

vector<RingElem> homogComps(const RingElem& F){
// returns the vector of homogeneous components of F
// w.r.t. standard grading

    SparsePolyRing P=AsSparsePolyRing(owner(F));
    long dim=NumIndets(P);
    vector<long> v(dim);
    vector<RingElem> c(deg(F)+1,zero(P));
    long j,k;

    //TODO there is a leading_term() function coming in cocoalib
    //TODO maybe there will be even a "splice_leading_term"
    SparsePolyIter i=BeginIter(F);
    for (; !IsEnded(i); ++i){
        exponents(v,PP(i));
        k=0;
        for(j=0;j<dim;j++)
            k+=v[j];
        PushBack(c[k],coeff(i),v);
    }
    return(c);
}

RingElem homogenize(const RingElem& F){
// homogenizes F wrt the zeroth variable and returns the 
// homogenized polynomial

    SparsePolyRing P=AsSparsePolyRing(owner(F));
    int d=deg(F);    
    vector<RingElem> c(d+1,zero(P));
    c=homogComps(F);
    RingElem h(zero(P));
    for(int i=0;i<=d;++i)
        h+=c[i]*power(indets(P)[0],d-i);
    return(h);
}

RingElem makeZZCoeff(const RingElem& F, const SparsePolyRing RZZ){
// F is a polynomial over RingQQ with integral coefficients
// This function converts it into a polynomial over RingZZ

    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    RingElem G(zero(RZZ));
    BigRat QQCoeff;
    vector<long> v;
    for (; !IsEnded(mon); ++mon){
        IsRational(QQCoeff,coeff(mon));
        exponents(v,PP(mon));
        PushBack(G,LC(num(QQCoeff)*one(RZZ)),v);
    }
    return(G);
}

RingElem makeQQCoeff(const RingElem& F, const SparsePolyRing R){
// F is a polynomial over RingZZ
// This function converts it into a polynomial over RingQQ

    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    RingElem G(zero(R));
    BigInt ZZCoeff;
    vector<long> v;
    for (; !IsEnded(mon); ++mon){
        IsInteger(ZZCoeff,coeff(mon));
        exponents(v,PP(mon));
        PushBack(G,LC(ZZCoeff*one(R)),v);
    }
    return(G);
}
 

// the next routine reads the input polynomail in factiored form,
// factorizes the factors, makes them integral,
// passes to the highest homogeneous if necessary
// and homogenizes for substitution with integral coefficients
// the return value is the polynomial defined by the input
 
RingElem processInputPolynomial(const string& project, const SparsePolyRing& R, const SparsePolyRing& RZZ,
                vector<RingElem>& resPrimeFactors, vector<RingElem>& resPrimeFactorsNonhom, vector<long>& resMultiplicities,
                RingElem& remainingFactor, bool& homogeneous,const bool& do_leadCoeff){
   
  long i,j;
  vector<RingElem> factorsRead=readFactorList(project,R); // read factors of F
  vector<long> multiplicities;                            

  vector<RingElem> primeFactors; // for use in this routine
  vector<RingElem> primeFactorsNonhom; // return results will go into the "res" parameters for output 
  
  if(verbose_INT)
    cout << "Polynomial read" << endl;
    
  homogeneous=true;                                     // we factor the polynomials read
  for(i=0;i< (long) factorsRead.size();++i){ // and make them integral this way
                                             // they must further be homogenized
                                             // and converted to polynomials with ZZ 
                                             // coefficients (instead of inegral QQ)
                                             // The homogenization is necessary to allow
                                             // substitutions over ZZ
      RingElem G(factorsRead[i]);
      if(deg(G)==0){         
        remainingFactor*=G;  // constants go into remainingFactor
        continue;       
      }
      
    vector<RingElem> compsG= homogComps(G);
                             // we test for homogeneity. In case do_leadCoeff==true, polynomial
                             // is replaced by highest homogeneous component
    if(G!=compsG[compsG.size()-1]){            
       if(verbose_INT && homogeneous && do_leadCoeff) 
           cout << "Polynomial is inhomogeneous. Replacing it by highest hom. comp." << endl;
       homogeneous=false;
       if(do_leadCoeff){
           G=compsG[compsG.size()-1];
           factorsRead[i]=G;  // it may no longer be the factor read from input
       }    
    }
     
    factorization<RingElem> FF=factor(G);              // now the factorization and transfer to integer coefficients
    for(j=0;j< (long) FF.myFactors.size();++j){
        primeFactorsNonhom.push_back(FF.myFactors[j]); // these are the factors of the polynomial to be integrated
        primeFactors.push_back(makeZZCoeff(homogenize(FF.myFactors[j]),RZZ)); // these are the homogenized factors
        multiplicities.push_back(FF.myExponents[j]);                          // homogenized for substitution !
      }
      remainingFactor*=FF.myRemainingFactor;
  }
  
  // cout << "remainingFactor " << remainingFactor  << endl; 
  
 // it remains to collect multiple factors that can come from different input factors
  
  for(i=0;i< (long) primeFactors.size();++i)
    if(primeFactors[i]!=0)
        for(j=i+1;j< (long) primeFactors.size();++j)
            if(primeFactors[j]!=0 && primeFactors[i]==primeFactors[j]){
                primeFactors[j]=0;
                multiplicities[i]++;
            }
            
  for(i=0;i< (long) primeFactors.size();++i)  // now everything is transferred to the return parameters
    if(primeFactors[i]!=0){
        resPrimeFactorsNonhom.push_back(primeFactorsNonhom[i]);
        resPrimeFactors.push_back(primeFactors[i]);
        resMultiplicities.push_back(multiplicities[i]);
  }
  
  RingElem F(one(R));
  for(i=0;i< (long) factorsRead.size();++i)  // multiply the factors
        F*=factorsRead[i];
  
  // cout << "F= " << F << endl; 
    
  return(F);
}
