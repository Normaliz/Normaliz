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

class CyclRatFunct {
// class for rational functions whose denominator is a product
// of cyclotomic polynomials
// We work with denominators that are products of factors 1-t^i
// which is of course equivalent
// the numerator is a polynomial in its ring
// the denominator is an integer vector that at index i
// gives the multiplicity of 1-t^i in the denominator
// (the entry at index 0 is not used and must always be equal to 0)
public:

    RingElem num;
    vector<long> denom;

    void extendDenom(const vector<long>& target);
    void addCRF(const CyclRatFunct& r);
    void multCRF(const CyclRatFunct& r);
    void simplifyCRF();
    void set2(const RingElem& f, const vector<long>& d);
    void set2(const RingElem& f);
    void showCRF();
    void showCoprimeCRF();
    CyclRatFunct(const RingElem& c);
    CyclRatFunct(const RingElem& c,const vector<long>& d);

};
//class end *****************************************************************


void CyclRatFunct::extendDenom(const vector<long>& target)
// extends the denominator to target
// by multiplying the numrerator with the remaining factor
{
    RingElem t=indets(AsSparsePolyRing(owner(num)))[0];
    long i,ns=target.size(),nf=denom.size();
    for(i=1;i<ns;++i){
        if(i>nf-1)
            num*=power(1-power(t,i),(target[i]));
        else
            if(target[i]>denom[i])
                num*=power(1-power(t,i),(target[i]-denom[i]));
    }
    denom=target;
}

vector<long> lcmDenom(const vector<long>& df, const vector<long>& dg){
// computes the lcm of ztwo denominators as used in CyclRatFunct
// (1-t^i and 1-t^j, i != j, are considered as coprime)
    size_t nf=df.size(),ng=dg.size(),i;
    size_t n=max(nf,ng);
    vector<long> dh=df;
    dh.resize(n);
    for(i=1;i<n;++i)
        if(i<ng && dh[i]<dg[i])
            dh[i]=dg[i];
    return(dh);
}


vector<long> prodDenom(const vector<long>& df, const vector<long>& dg){
// as above, but computes the profduct
    size_t nf=df.size(),ng=dg.size(),i;
    size_t n=max(nf,ng);
    vector<long> dh=df;
    dh.resize(n);
    for(i=1;i<n;++i)
        if(i<ng)
            dh[i]+=dg[i];
    return(dh);
}

vector<long> degrees2denom(const vector<long>& d){
// converts a vector of degrees to a "denominator"
// listing at position i the multiplicity of i in d
    long m=0;
    size_t i;
    for(i=0;i<d.size();++i)
        m=max(m,d[i]);
    vector<long> e(m+1);
    for(i=0;i<d.size();++i)
        e[d[i]]++;
    return(e);
}

vector<long> denom2degrees(const vector<long>& d){
// the converse operation
    vector<long> denomDeg;
    for(size_t i=0;i<d.size();++i)
        for(long j=0;j<d[i];++j)
           denomDeg.push_back(i);
    return(denomDeg);
}

RingElem denom2poly(const SparsePolyRing& P, const vector<long>& d){
// converts a denominator into a real polynomial
// the variable for the denominator is x[0]
    RingElem t=indets(P)[0];
    RingElem f(one(P));
    for(size_t i=1;i<d.size();++i)
        f*=power(1-power(t,i),d[i]);
    return(f);
}

vector<long> makeDenom(long k,long n)
// makes the denominator (1-t^k)^n
{
    vector<long> d(k+1);
    d[k]=n;
    return(d);
}

void CyclRatFunct::addCRF(const CyclRatFunct& r){
// adds r to *this, r is preserved in its given form
    CyclRatFunct s(zero(AsSparsePolyRing(owner(num))));
    const vector<long> lcmden(lcmDenom(denom,r.denom));
    s=r;
    s.extendDenom(lcmden);
    extendDenom(lcmden);
    num+=s.num;
}

void CyclRatFunct::multCRF(const CyclRatFunct& r){
// nmultiplies *this by r
    num*=r.num;
    denom=prodDenom(denom,r.denom);
}

void CyclRatFunct::showCRF(){
    if(!verbose_INT)
        return;

    cout << num << endl;
    for(size_t i=1;i<denom.size();++i)
        cout << denom[i] << " ";
    cout << endl;
}

void CyclRatFunct::showCoprimeCRF(){
// shows *this also with coprime numerator and denominator
// makes only sense if only x[0] appears in the numerator (not checked)

    if(!verbose_INT)
        return;

    cout << "--------------------------------------------" << endl << endl;
    cout << "Given form" << endl << endl;
    showCRF();
    cout << endl;
    SparsePolyRing R=AsSparsePolyRing(owner(num));
    SparsePolyRing P=NewPolyRing_DMPI(RingQQ(),symbols("t"));
    vector<RingElem> Im(NumIndets(R),zero(P));
    Im[0]=indets(P)[0];
    RingHom phi=PolyAlgebraHom(R,P,Im);
    RingElem f(phi(num));
    RingElem g(denom2poly(P,denom));
    RingElem h=gcd(f,g);
    f/=h;
    g/=h;
    cout << "Coprime numerator (for denom with remaining factor 1)" << endl <<endl;
    factorization<RingElem> gf=factor(g);
    cout << f/gf.myRemainingFactor << endl << endl << "Factorization of denominator" << endl << endl;
    size_t nf=gf.myFactors.size();
    for(size_t i=0;i<nf;++i)
        cout << gf.myFactors[i] << "  mult " << gf.myMultiplicities[i] << endl;
    cout << "--------------------------------------------" << endl;

}

void CyclRatFunct::simplifyCRF(){
// cancels factors 1-t^i from the denominator that appear there explicitly
// (and not just as factors of 1-t^j for some j)

    SparsePolyRing R=AsSparsePolyRing(owner(num));
    long nd=denom.size();
    for(long i=1;i<nd;i++)
    {
        while(denom[i]>0)
        {
            if(!IsDivisible(num,1-power(indets(R)[0],i)))
                break;
            num/=1-power(indets(R)[0],i);
            denom[i]--;
        }
    }
}

void CyclRatFunct::set2(const RingElem& f, const vector<long>& d)
{
    num=f;
    denom=d;
}

void CyclRatFunct::set2(const RingElem& f)
{
    num=f;
    denom.resize(1,0);
}

CyclRatFunct::CyclRatFunct(const RingElem& c):num(c)
// constructor starting from a RingElem
// initialization necessary because RingElem has no default
// constructor
{
    denom.resize(1,0);
}

CyclRatFunct::CyclRatFunct(const RingElem& c,const vector<long>& d):num(c),denom(d){
}
