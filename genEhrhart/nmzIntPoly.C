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

RingElem affineLinearSubstitutionFL(const factorization<RingElem>& FF,const vector<vector<long> >& A,
                     const vector<long>& b, const long& denom, const RingElem& F){
// we need F to define the ring
// applies linar substitution y --> A(y+b/denom) to all factors in FF 
// and returns the product of the modified factors

    size_t i;
    SparsePolyRing R=AsSparsePolyRing(owner(F));
    size_t m=A.size();
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
    for(i=1;i<w1.size();++i) // we have to shift w since the zeroth variable
        w1[i]=w[i-1];        // must not be touched
    
    RingHom phi=PolyAlgebraHom(R,R,w1);
    RingElem G(phi(FF.myRemainingFactor));
    for(i=0;i<FF.myFactors.size();++i){
        G*=power(phi(FF.myFactors[i]),FF.myExponents[i]);
    }   
    return(G);   
}

RingElem homogeneousLinearSubstitutionFL(const factorization<RingElem>& FF,const vector<vector<long> >& A,
                     const vector<long>& degrees, const RingElem& F){
// we need F to define the ring
// applies linar substitution y --> (A/degress)y to all factors in FF 
// where row A[i] is divided by degrees[i]
// and returns the product of the modified factors

    size_t i;
    SparsePolyRing R=AsSparsePolyRing(owner(F));
    size_t m=A.size();
    vector<RingElem> v(m,zero(R));
    
    for(i=0;i<m;i++)
        v[i]=indets(R)[i+1]/degrees[i];
    vector<RingElem> w=VxM(v,A);
    vector<RingElem> w1(w.size()+1,zero(R));
    w1[0]=indets(R)[0];
    for(i=1;i<w1.size();++i) // we have to shift w since the zeroth variable
        w1[i]=w[i-1];        // must not be touched
        
    
    RingHom phi=PolyAlgebraHom(R,R,w1);
    RingElem G(phi(FF.myRemainingFactor));
    for(i=0;i<FF.myFactors.size();++i){
        G*=power(phi(FF.myFactors[i]),FF.myExponents[i]);
    }   
    return(G);   
}

vector<RingElem> homogComps(const RingElem& F){
// returns the vector of homogeneous components of F
// w.r.t. standard grading

    SparsePolyRing P=AsSparsePolyRing(owner(F));
    long dim=NumIndets(P);
    vector<long> v(dim);
    vector<RingElem> c(deg(F)+1,zero(P));
    long j,k;

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
