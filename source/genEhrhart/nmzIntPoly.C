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

bool compareLength(const RingElem& p, const RingElem& q){
    return(NumTerms(p)>NumTerms(q));
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

bool DDD=false;

vector<long> shiftVars(const vector<long>& v, const vector<long>& key){
// selects components of v and reorders them according to key
    vector<long> w(v.size(),0);
    for(size_t i=0;i<key.size();++i){
        w[i]=v[key[i]];
    }
    return(w);
}

void  makeLocalDegreesAndKey(const boost::dynamic_bitset<>& indicator,const vector<long>& degrees, vector<long>& localDeg,vector<long>& key){

    localDeg.clear();
    key.clear();
    key.push_back(0);
    for(size_t i=0;i<indicator.size();++i)
        if(indicator.test(i))
            key.push_back(i+1);
    for(size_t i=0;i<key.size()-1;++i)
        localDeg.push_back(degrees[key[i+1]-1]);  
}

void makeStartEnd(const vector<long>& localDeg, vector<long>& St, vector<long>& End){

    vector<long> denom=degrees2denom(localDeg); // first we must find the blocks of equal degree
    St.push_back(1);
    for(size_t i=0;i<denom.size();++i)
        if(denom[i]!=0){
            End.push_back(St[St.size()-1]+denom[i]-1);
            if(i<denom.size()-1)
                St.push_back(End[End.size()-1]+1);
        }
}

vector<long> orderExposInner(vector<long>& vin, const vector<long>& St, vector<long>& End){

    vector<long> v=vin;
    long p,s,pend,pst;
    bool ordered;
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
    return(v);
}

RingElem orderExpos(const RingElem& F, const vector<long>& degrees, const boost::dynamic_bitset<>& indicator, bool compactify){
 // orders the exponent vectors v of the terms of F
 // the exponents v[i] and v[j], i < j,  are swapped if
 // (1) degrees[i]==degrees[j] and (2) v[i] < v[j]
 // so that the exponents are descending in each degree block
 // the ordered exponent vectors are inserted into a map
 // and their coefficients are added
 // at the end the polynomial is rebuilt from the map
 // If compactify==true, the exponents will be shifted to the left in order to keep the correspondence
 // of variables to degrees
 // compactification not used at present (occurs only in restrictToFaces) 

    SparsePolyRing P=AsSparsePolyRing(owner(F));
    vector<long> v(NumIndets(P));
    vector<long> key,localDeg;
    key.reserve(v.size()+1);
    localDeg.reserve(degrees.size()+1);
    
    if(compactify){
        makeLocalDegreesAndKey(indicator,degrees,localDeg,key);  
    }
    else{
        localDeg=degrees;
    }
    
    vector<long> St,End;
    makeStartEnd(localDeg,St,End);

    // now the main job

    map<vector<long>,RingElem> orderedMons;  // will take the ordered exponent vectors
    map<vector<long>,RingElem>::iterator ord_mon;

    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    for (; !IsEnded(mon); ++mon){
      exponents(v,PP(mon)); // this function gives the exponent vector back as v
      if(compactify)
          v=shiftVars(v,key);
      v=orderExposInner(v,St,End);
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
 //JAA   cout << "Loop start " << orderedMons.size() <<  endl;
 //JAA   size_t counter=0;
    for(ord_mon=orderedMons.begin();ord_mon!=orderedMons.end();++ord_mon){
 //JAA       cout << counter++ << ord_mon->first << endl;
 //JAA       try {
        PushFront(r,ord_mon->second,ord_mon->first);
//JAA        }
 //JAA       catch(const std::exception& exc){cout << "Caught exception: " << exc.what() << endl;}
    }
//JAA    cout << "Loop end" << endl;
    return(r);
}


void restrictToFaces(const RingElem& G,RingElem& GOrder, vector<RingElem>& GRest,const vector<long> degrees, const vector<SIMPLINEXDATA_INT>& inExSimplData){
// Computesd the restrictions of G to the faces in inclusion-exclusion.
// All terms are simultaneously compactified and exponentwise ordered
// Polynomials returned in GRest
// Ordering is also applied to G itself, returned in GOrder
// Note: degrees are given for the full simplex. Therefore "local" degreees must be made 
// (depend only on face and not on offset, but generation here is cheap)

    SparsePolyRing P=AsSparsePolyRing(owner(G));

    vector<long> v(NumIndets(P));
    vector<long> w(NumIndets(P));
    vector<long> localDeg;
    localDeg.reserve(v.size());
    size_t dim=NumIndets(P)-1;
 
    // first we make the facewise data that are needed for the compactification and otrdering
    // of exponent vectors    
    vector<vector<long> > St(inExSimplData.size()),End(inExSimplData.size()),key(inExSimplData.size());
    vector<long> active;
    for(size_t i=0;i<inExSimplData.size();++i)
        if(!inExSimplData[i].done){
            active.push_back(i);
            makeLocalDegreesAndKey(inExSimplData[i].GenInFace ,degrees,localDeg,key[i]);
            makeStartEnd(localDeg,St[i],End[i]);
        }
    
    // now the same for the full simplex (localDeg=degrees)
    boost::dynamic_bitset<> fullSimpl(dim);
    fullSimpl.set();
    vector<long> StSimpl,EndSimpl;
    makeStartEnd(degrees,StSimpl,EndSimpl);
    
    vector<map<vector<long>,RingElem> > orderedMons(inExSimplData.size());  // will take the ordered exponent vectors
    map<vector<long>,RingElem> orderedMonsSimpl; 
    map<vector<long>,RingElem>::iterator ord_mon;


    boost::dynamic_bitset<> indicator(dim);

    // now we go over the terms of G
    SparsePolyIter term=BeginIter(G);
    //PPMonoid TT = PPM(AsSparsePolyRing(owner(G)));
    for (; !IsEnded(term); ++term){
        //PPMonoidElem mon(PP(term));
        exponents(v,PP(term));
        w=v;
        indicator.reset();
        for(size_t j=0;j<dim;++j)
            if(v[j+1]!=0)               // we must add 1 since the 0-th indeterminate is irrelevant here
                indicator.set(j);
        for(size_t i=0;i<active.size();++i){
            int j=active[i];
            if(indicator.is_subset_of(inExSimplData[j].GenInFace)){
                w=shiftVars(v,key[j]);
                w=orderExposInner(w,St[j],End[j]);
                // w=shiftVars(v,key[j]);
                ord_mon=orderedMons[j].find(w); // insert into map or add coefficient
                if(ord_mon!=orderedMons[j].end()){
                    ord_mon->second+=coeff(term);
                }
                else{
                    orderedMons[j].insert(pair<vector<long>,RingElem>(w,coeff(term)));
                }
             }
        } // for i
        
        v=orderExposInner(v,StSimpl,EndSimpl);
        ord_mon=orderedMonsSimpl.find(v); // insert into map or add coefficient
        if(ord_mon!=orderedMonsSimpl.end()){
            ord_mon->second+=coeff(term);
        }
        else{
            orderedMonsSimpl.insert(pair<vector<long>,RingElem>(v,coeff(term)));
        }
    } // loop over term
    
    // now we must make the resulting polynomials from the maps
    
    for(size_t i=0;i<active.size();++i){
        int j=active[i];
        for(ord_mon=orderedMons[j].begin();ord_mon!=orderedMons[j].end();++ord_mon){
            PushFront(GRest[j],ord_mon->second,ord_mon->first);
        }
        // cout << "GRest[j] " << j << " " << NumTerms(GRest[j]) << endl;
    }
            
    for(ord_mon=orderedMonsSimpl.begin();ord_mon!=orderedMonsSimpl.end();++ord_mon){
        PushFront(GOrder,ord_mon->second,ord_mon->first);
    }
    
}

long nrActiveFaces=0; 
long nrActiveFacesOld=0;

void all_contained_faces(const RingElem& G, RingElem& GOrder,const vector<long>& degrees, boost::dynamic_bitset<>& indicator, long Deg, 
                     vector<SIMPLINEXDATA_INT>& inExSimplData, vector<deque<pair<vector<long>,RingElem> > >& facePolys){
                     
    const SparsePolyRing R=AsSparsePolyRing(owner(G));
    vector<RingElem> GRest;
    // size_t dim=indicator.size();
    for(size_t i=0;i<inExSimplData.size();++i){
        GRest.push_back(zero(R));
        
        if(!indicator.is_subset_of(inExSimplData[i].GenInFace))  
            inExSimplData[i].done=true;       // done if face cannot contribute to result for this offset
        else
            inExSimplData[i].done=false;       // not done otherwise
    }
    restrictToFaces(G,GOrder,GRest,degrees,inExSimplData);
    int tn;
    if(omp_get_level()==0)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(1);
    for(size_t j=0;j<inExSimplData.size();++j){
        if(inExSimplData[j].done)
            continue;
        #pragma omp atomic
        nrActiveFaces++;
        // cout << "Push back " << NumTerms(GRest[j]);
        GRest[j]=power(indets(R)[0],Deg)*inExSimplData[j].mult*GRest[j];  // shift by degree of offset amd multiply by mult of face
        facePolys[tn].push_back(pair<vector<long>,RingElem>(inExSimplData[j].degrees,GRest[j]));
        // cout << " Now " << facePolys[tn].size() << endl;
    }
}
        
RingElem affineLinearSubstitutionFL(const ourFactorization& FF,const vector<vector<long> >& A,
                     const vector<long>& b, const long& denom, const SparsePolyRing& R, const vector<long>& degrees, const BigInt& lcmDets,
                     vector<SIMPLINEXDATA_INT>& inExSimplData,vector<deque<pair<vector<long>,RingElem> > >& facePolys){
// applies linar substitution y --> lcmDets*A(y+b/denom) to all factors in FF 
// and returns the product of the modified factorsafter ordering the exponent vectors

    size_t i;
    size_t m=A.size();
    size_t dim=m;                    // TO DO: eliminate this duplication
    vector<RingElem> v(m,zero(R));
    RingElem g(zero(R));
    
    for(i=0;i<m;i++)
    {
        g=b[i]*(lcmDets/denom);
        v[i]=g+lcmDets*indets(R)[i+1];
    }
    vector<RingElem> w=VxM(v,A);
    vector<RingElem> w1(w.size()+1,zero(R));
    w1[0]=RingElem(R,lcmDets);
    for(i=1;i<w1.size();++i) 
        w1[i]=w[i-1]; 
    
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
    
    if(inExSimplData.size()==0){    // not really necesary, but a slight shortcut
        boost::dynamic_bitset<> dummyInd;
        return(orderExpos(G,degrees,dummyInd,false));
    }
    
    // if(inExSimplData.size()!=0){
        long Deg=0;
        boost::dynamic_bitset<> indicator(dim); // indicates the non-zero components of b
        indicator.reset();
        for(size_t i=0;i<dim;++i)
            if(b[i]!=0){
                indicator.set(i);
                Deg+=degrees[i]*b[i];
            }
        Deg/=denom;
        RingElem Gorder(zero(R));            
        all_contained_faces(G,Gorder,degrees,indicator, Deg, inExSimplData,facePolys);
        return(Gorder);
    // }  
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
    for (; !IsEnded(mon); ++mon){
         PushBack(G,num(coeff(mon)),PP(mon));
    }
    return(G);
}


RingElem makeQQCoeff(const RingElem& F, const SparsePolyRing R){
// F is a polynomial over RingZZ
// This function converts it into a polynomial over RingQQ
    SparsePolyIter mon=BeginIter(F); // go over the given polynomial
    RingElem G(zero(R));
    for (; !IsEnded(mon); ++mon){
        PushBack(G,RingElem(RingQQ(),coeff(mon)),PP(mon));  
    }
    return(G);
}
 

// the next routine reads the input polynomail in factiored form,
// factorizes the factors, makes them integral,
// passes to the highest homogeneous component if necessary
// and homogenizes for substitution with integral coefficients if necessary
// 0-th variable used for homogenization
// the return value is the polynomial defined by the input
 
RingElem processInputPolynomial(const string& project, const SparsePolyRing& R, const SparsePolyRing& RZZ,
                vector<RingElem>& resPrimeFactors, vector<RingElem>& resPrimeFactorsNonhom, vector<long>& resMultiplicities,
                RingElem& remainingFactor, bool& homogeneous,const bool& do_leadCoeff){
// "res" stands for "result"
// resPrimeFactors are homogenized, the "nonhom" come from the original polynomial
   
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
        continue;            // this extra treatment would not be necessary      
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
           factorsRead[i]=G;  // though it may no longer be the factor read from input
       }    
    }
     
    factorization<RingElem> FF=factor(G);              // now the factorization and transfer to integer coefficients
    for(j=0;j< (long) FF.myFactors().size();++j){
        primeFactorsNonhom.push_back(FF.myFactors()[j]); // these are the factors of the polynomial to be integrated
        primeFactors.push_back(makeZZCoeff(homogenize(FF.myFactors()[j]),RZZ)); // the homogenized factors with ZZ coeff
        multiplicities.push_back(FF.myMultiplicities()[j]);                          // homogenized for substitution !
      }
      remainingFactor*=FF.myRemainingFactor();
  }

  
 // it remains to collect multiple factors that come from different input factors
  
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
  
  RingElem F(one(R));                        //th polynomial to be integrated
  for(i=0;i< (long) factorsRead.size();++i)  // with QQ coefficients
        F*=factorsRead[i]; 
    
  return(F);
}
