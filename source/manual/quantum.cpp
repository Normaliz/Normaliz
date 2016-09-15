#include <stdlib.h>
#include <vector>
#include <map>
#include <set>

#include <list>

#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <omp.h>
using namespace std;

#include "libnormaliz.h"
#include "cone.h"
#include "vector_operations.h"
#include "cone_property.h"
#include "integer.h"
//#include "libnormaliz/libnormaliz.cpp"
using namespace libnormaliz;

typedef long Integer;

// *********************** Write to database ********************

void SaveMat(const vector<vector<Integer> > Mat){
        size_t dim=Mat[0].size();
        string dim_str = static_cast<ostringstream*>( &(ostringstream() << dim-1) )->str();
        string database="Max_Polytopes_";
        database+=dim_str;
        string typ="log";
        Matrix<Integer>(Mat).print_append(database,typ);
}


// *********************** Mailing ********************

void run_pgm(const char *pgm, const char *par1, const char *par2, const char *par3)
// Baut einen Befehl auf, der dann extern ausgefuehrt wird
// Setzt sich zusammen aus den 4 Parametern, wobei hinter pgm
// ein Blank eingefuegt wird
{
    char *FullCommand;

    fflush(stdout);
    FullCommand= (char*) calloc (strlen(pgm)+1+strlen(par1)+
          strlen(par2)+strlen(par3)+3,sizeof(char));
    strcpy(FullCommand,pgm);
    strcat(FullCommand," ");
    strcat(FullCommand,par1);
    strcat(FullCommand,par2);
    strcat(FullCommand,par3);
    // CallZ++;
    // printf("Rufe %s %d\n",FullCommand,CallZ);
    // printf("Calling %s\n",FullCommand);
    // fflush(stdout);
    int i;
    if(system(NULL))
        i=system(FullCommand);
    else
    {
        printf("Nix geht mehr\n");
        exit(1);
    }
    i+=i;
    fflush(stdout);
}

template<typename Integer>
void Matrix<Integer>::print(const string& name,const string& suffix) const{
    string file_name = name+"."+suffix;
    const char* file = file_name.c_str();
    ofstream out(file);
    print(out);
    out.close();
}


void FileMail(const char *Addr, const char *FileName)
{
    string prog="mail -s \"QuantumXXXX\" ";
    string read_from=" < "; 
    run_pgm(prog.c_str(), Addr,read_from.c_str(),FileName);
}


template<typename Integer>
void MatMail(vector<vector<Integer> >& Mat){

    string name="Mat";
    string typ="mail";
    string addr="wbruns@uos.de";
    string file_name="Mat.mail";
    Matrix<Integer>(Mat).print(name,typ);
    FileMail(addr.c_str(),file_name.c_str());
}


// *************** Global variables that save parameters ***********************

int ChooseStart=0;
int ChoosePotential=0;
bool BreakReg=false;  // helps to break regular growth
bool FromSmall=false, HeightOneFirst=false, RandomExt=false, Print=false;
vector<vector<Integer> > Origin;
bool Verify=false;
bool AllExtensions=false;
bool CheckIllumination=false;
Integer MinHeight;

template<typename Integer>
class Polytope{

public:

    vector<vector<Integer> > Gen;
    vector<vector<Integer> > Supp;
    vector<vector<Integer> > LattPoints;
    vector<vector<Integer> > HB;
    vector<vector<Integer> > ExtRays;
    vector<Integer> FacetAreas;
    vector<Integer> Widths;
    vector<bool> Illu;
};

// *************** Helpers ***********************


template<typename Integer>
void scramble(vector<vector<Integer> >& Mat){
    for(size_t i=0;i<2*Mat.size();++i)
        swap(Mat[rand()%Mat.size()],Mat[rand()%Mat.size()]);
}

template<typename Integer>
vector<Integer> v_add_1(vector<Integer>& a,vector<Integer>& b){
// adds two vectors and returns the sum
   assert(a.size() == b.size());
   /* if (test_arithmetic_overflow) {  // does arithmetic tests
       return(v_add_overflow_check(a,b));
   } */
    size_t i,s=a.size();
    vector<Integer> d(s);
    for (i = 0; i <s; i++) {
        d[i]=a[i]+b[i];
    }
    return d;
}


template<typename Integer>
// checks whether point is in the cone defined by support hyperplanes
bool InCone(vector<Integer>& point, vector<vector<Integer> >& cone){

    for(size_t i=0;i<cone.size();++i)
        if(v_scalar_product(cone[i],point)<0)
            return false;
    return true;

}

bool prime(int n){
    for(int i=2;;++i){
        if(i*i>n)
            return true;
        if( n%i == 0)
            return false;
    }
}

template<typename Integer>
vector<int> get3primes(){

    int p,q,r;

    do{
        p=rand()%18+2;
    }while(!prime(p));

    do{
        q=rand()%18+2;
    }while(!(prime(q) && p!=q));

    do{
        r=rand()%18+2;
    }while(!(prime(r) && p!=r && q!=r));

    vector<int>  Primes(3);
    Primes[0]=p;
    Primes[1]=q;
    Primes[2]=r;
    // cout << Primes;
    sort(Primes.begin(),Primes.end());
    return Primes;

}

template<typename Integer>
bool is_reducible(vector<Integer>& c, vector<vector<Integer> >& Supp, vector<vector<Integer> >& Values){

    vector<Integer> v(Supp.size());
    for(size_t k=0;k<Supp.size();++k)
        v[k]=v_scalar_product(Supp[k],c);

    for(size_t m=0;m<Values.size();++m){
        bool reduced_by_this=true;
        for(size_t k=0; k<Values.size();++k)
            if(v[k]<Values[m][k]){
                reduced_by_this=false;
                break;
            }
        if(reduced_by_this){
            return true;
        }
    }
    return false;
}

/*template<typename Integer>
Integer Iabs(Integer n){
    if(n>=0)
        return n;
    return -n;
}*/

template<typename Integer>
int direction(vector<Integer>& v1, vector<Integer>& v2){

    bool First=true;
    int dir=0;
    Integer DiffMax;
    int dim=v1.size();
    Integer fact=v1[dim-1];
    for(size_t i=0;i<dim-1;++i)
        if(First || Iabs(fact*v2[i]-v1[i])>DiffMax){
            First=false;
            DiffMax=Iabs(fact*v2[i]-v1[i]);
            dir=2*i;
            if(fact*v2[i]-v1[i]<0)
                dir++;
        }
    return dir;

}



// ************** Access to libnormaliz ********************

template<typename Integer>
bool  Normalize(vector<vector<Integer> >& Gen, vector<vector<Integer> >& HB){


    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::HilbertBasis);
    Cone<Integer> StartCone(Gen, Type::integral_closure);
    StartCone.compute(StartPolyWanted);

    HB=StartCone.getHilbertBasis();
    sort(HB.begin(),HB.end());

    size_t k=0;
    int dim=Gen[0].size();

    for(;k<HB.size();++k)
        if(HB[k][dim-1]!=1)
                break;
    return(k==HB.size() && StartCone.getRank()==dim);
}

template<typename Integer>
bool  IsNormal(vector<vector<Integer> >& Gen){
    
    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::HilbertBasis);
    // StartPolyWanted.set(ConeProperty::DualMode);
    Cone<Integer> StartCone(Gen, Type::integral_closure);
    StartCone.compute(StartPolyWanted);
    return StartCone.isIntegrallyClosed();
}

template<typename Integer>
bool  IsNormalPolytope(vector<vector<Integer> >& Gen){
    
    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::HilbertBasis);
    // StartPolyWanted.set(ConeProperty::DualMode);
    Cone<Integer> StartCone(Gen, Type::integral_closure);
    StartCone.compute(StartPolyWanted);
    return StartCone.isDeg1HilbertBasis();
}

template<typename Integer>
bool  NormalizeAndExtRays(vector<vector<Integer> >& Gen, vector<vector<Integer> >& HB,vector<vector<Integer> >& ExtremeRays){

    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::HilbertBasis);
    Cone<Integer> StartCone(Gen, Type::integral_closure);
    StartCone.compute(StartPolyWanted);

    HB=StartCone.getHilbertBasis();
    ExtremeRays=StartCone.getExtremeRays();

    sort(HB.begin(),HB.end());
    sort(ExtremeRays.begin(),ExtremeRays.end());

    size_t k=0;
    int dim=Gen[0].size();

    for(;k<HB.size();++k)
        if(HB[k][dim-1]!=1)
                break;
    return(k==HB.size() && StartCone.getRank()==dim);
}

template<typename Integer>
Integer Multiplicity(vector<vector<Integer> >& Gen){

    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::Multiplicity);
    Cone<Integer> StartCone(Gen, Type::integral_closure);
    StartCone.compute(StartPolyWanted);
    mpz_class M=StartCone.getMultiplicity().get_num();
    return explicit_cast_to_long(M);

}

template<typename Integer>
vector<vector<Integer> > ExtremeRays(vector<vector<Integer> >& Gen){

    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::SupportHyperplanes);
    Cone<Integer> StartCone = libnormaliz::Cone<Integer>(Gen,Type::integral_closure);
    StartCone.compute(StartPolyWanted);
    return StartCone.getExtremeRays();
}

template<typename Integer>
vector<vector<Integer> > SuppHyps(vector<vector<Integer> >& Gen){

    ConeProperties StartPolyWanted;
    StartPolyWanted.set(ConeProperty::SupportHyperplanes);
    Cone<Integer> StartCone = libnormaliz::Cone<Integer>(Gen,Type::integral_closure);
    StartCone.compute(StartPolyWanted);
    return StartCone.getSupportHyperplanes();
}

template<typename Integer>
vector<vector<Integer> > Deg1FromSupps(vector<vector<Integer> >& SuppHyps){

    int dim=SuppHyps[0].size();

    vector< vector<Integer > > StartPolySupp, Grading;
    Grading.resize(1);  //make grading = last coordinate
    Grading[0].resize(dim);
    for(int p=0;p<dim;++p)
        Grading[0][p]=0;
    Grading[0][dim-1]=1;

    map <InputType, vector< vector<Integer> > > CandPolytope;
    CandPolytope[Type::grading]=Grading;
    CandPolytope[Type::inequalities]=SuppHyps;

    ConeProperties CandPolyWantedDeg1;
    CandPolyWantedDeg1.set(ConeProperty::Deg1Elements);
    CandPolyWantedDeg1.set(ConeProperty::ApproximateRatPolytope);

    Cone<Integer> Cand(CandPolytope);
    Cand.compute( CandPolyWantedDeg1);
    vector<vector<Integer> > Deg1Elements=Cand.getDeg1Elements();
    sort(Deg1Elements.begin(),Deg1Elements.end());

    return Deg1Elements;
}

template<typename Integer>
vector<vector<Integer> > Deg1FromSupps_mpz(vector<vector<Integer> >& CandPolySupp){

    int dim=CandPolySupp[0].size();

    ConeProperties CandPolyWantedDeg1;
    CandPolyWantedDeg1.set(ConeProperty::Deg1Elements);
    CandPolyWantedDeg1.set(ConeProperty::ApproximateRatPolytope);


    vector<vector<mpz_class> > mpzCandPolySupp;
    mpzCandPolySupp.resize(CandPolySupp.size());
    for(size_t j=0;j<CandPolySupp.size();++j){
        mpzCandPolySupp[j].resize(dim);
        for(size_t k=0;k<dim;++k)
            mpzCandPolySupp[j][k]=to_mpz(CandPolySupp[j][k]);
    }
    vector<vector<mpz_class> > mpzGrading;
    mpzGrading.resize(1);
    mpzGrading[0].resize(dim);
    long dummy=0;
    for(size_t i=0;i<mpzGrading.size();++i)
        mpzGrading[0][i]=to_mpz(dummy);
    dummy=1;
    mpzGrading[0][dim-1]=to_mpz(dummy);

    map <InputType, vector< vector<mpz_class> > > mpzCandPolytope;
    mpzCandPolytope[Type::grading]=mpzGrading;
    mpzCandPolytope[Type::inequalities]=mpzCandPolySupp;
    Cone<mpz_class> mpzCandCone = libnormaliz::Cone<mpz_class>(mpzCandPolytope);

    // verbose=true;
    mpzCandCone.compute(CandPolyWantedDeg1);
    // verbose=false;

    vector<vector<mpz_class> > mpzCandJumps=mpzCandCone.getDeg1Elements();
    vector<vector<Integer> > CandJumps;
    CandJumps.resize(mpzCandJumps.size());
    for(size_t j=0;j<CandJumps.size();++j){
        CandJumps[j].resize(dim);
        for(size_t k=0;k<dim;++k)
            CandJumps[j][k]=explicit_cast_to_long<mpz_class>(mpzCandJumps[j][k]);
    }
    sort(CandJumps.begin(),CandJumps.end());
    return CandJumps;

}

template<typename Integer>
vector<vector<Integer> > Deg1FromGens(vector<vector<Integer> >& Gen){

    ConeProperties CandPolyWantedDeg1;
    CandPolyWantedDeg1.set(ConeProperty::Deg1Elements);

    Cone<Integer> Cand(Gen, Type::integral_closure);
    Cand.compute( CandPolyWantedDeg1);
    vector<vector<Integer> > Deg1Elements=Cand.getDeg1Elements();
    // sort(Deg1Elements.begin(),Deg1Elements.end());

    return Deg1Elements;
}

template<typename Integer>
bool  IsVeryAmple(vector<vector<Integer> >& Gen){

    size_t dim=Gen[0].size();

    ConeProperties CandPolyWantedDeg1Elements;
    CandPolyWantedDeg1Elements.set(ConeProperty::Deg1Elements);
    
    // Matrix<Integer>(Gen).pretty_print(cout);    
    
    Cone<Integer> Cand(Gen, Type::integral_closure);
    Cand.compute( CandPolyWantedDeg1Elements);
    vector<vector<Integer> > ExtRays=Cand.getExtremeRays();
    vector<vector<Integer> > LattP=Cand.getDeg1Elements();
    
    // cout << "Gitterpunkte berechnet" << endl;
    
    vector<vector<Integer> > CornerCone;
    CornerCone.reserve(LattP.size()-1);
    
    size_t i;
    
    for(i=0;i<ExtRays.size();++i){
    
        CornerCone.clear();
    
        for(size_t j=0;j<LattP.size();++j){
            if(LattP[j]==ExtRays[i])
                continue;
            vector<Integer> diff(dim);
            for(size_t k=0;k<dim;++k)
                diff[k]=LattP[j][k]-ExtRays[i][k];
            CornerCone.push_back(diff);    
        }
        
        // Matrix<Integer>(CornerCone).pretty_print(cout);
        
        if(!IsNormal(CornerCone))
            break;               
    }
    return i==ExtRays.size();

}
// *************** Heights and potentials ***********************

template<typename Integer>
vector<Integer> FacetAreas(Polytope<Integer>& Poly){

    vector<vector<Integer> > GensInFacet;
    vector<Integer> Areas;
    if(Poly.Supp.size()==0)
        Poly.Supp=SuppHyps(Poly.Gen);
    GensInFacet.reserve(Poly.Gen.size());
    for(size_t i=0;i<Poly.Supp.size();++i){
        GensInFacet.clear();
        for(size_t j=0;j<Poly.Gen.size();++j)
            if(v_scalar_product(Poly.Gen[j],Poly.Supp[i])==0)
                GensInFacet.push_back(Poly.Gen[j]);
        Areas.push_back(Multiplicity(GensInFacet));
    }
    return Areas;
}

template<typename Integer>
Integer SurfaceArea(Polytope<Integer>& Poly){

    Integer Area=0;
    if(Poly.Supp.size()==0)
        Poly.Supp=SuppHyps(Poly.Gen);
    vector<Integer> Areas=FacetAreas(Poly);
    for(size_t i=0;i<Areas.size();++i)
        Area+=Areas[i];
    // cout << "Area " << Area << endl;
    return Area;
}

template<typename Integer>
float EuclDist(vector<Integer>& p1, vector<Integer>& p2){
// computes the Euclidean distance of two vectors

    float dist=0;
    for(size_t i=0;i<p1.size();++i)
        dist+=pow(p1[i]-p2[i],2);
    return sqrt(dist);

}

template<typename Integer>
float DistPotential(vector<Integer>& p1,Polytope<Integer>& Poly){
// computes the distance based potential
    float pot=0;
    for(size_t i=0;i<Poly.Gen.size();++i){
        pot+=1.0/pow(EuclDist(p1,Poly.Gen[i]),p1.size()-3);
    }
    return pot;
}

template<typename Integer>
Integer HeightOverPolytope(vector<Integer>& p1, Polytope<Integer>& Poly){
// returns the maximal height over the polytope
    Integer height=0;
    for(size_t i=0;i<Poly.Supp.size();++i){
        Integer h=v_scalar_product(Poly.Supp[i],p1);
        if(h<0 && -h >height)
            height=-h;
    }
    return height;
}

template<typename Integer>
float HeightPotential(vector<Integer>& p1, Polytope<Integer>& Poly){
// computes potential based on all heights: 11/sum of negative heights
    Integer pot=0;
    for(size_t i=0;i<Poly.Supp.size();++i){
        Integer h=v_scalar_product(Poly.Supp[i],p1);
        if(h<0)
            pot-=h;
    }
    // cout << "pot " << pot << endl;
    return 1.0/pot;
}

template<typename Integer>
float MaxHeightPotential(vector<Integer>& p1, Polytope<Integer>& Poly){
// computes potential based on maximal heights

    Integer h=HeightOverPolytope(p1,Poly);
    if(h==0){
        cout << "Height 0 error " << endl;
        exit(1);
    }
    return 1.0/h;
}

template<typename Integer>
float MinHeightPotential(vector<Integer>& p1,Polytope<Integer>& Poly){
// computes potential based on minimal height
    Integer height=0;
    for(size_t i=0;i<Poly.Supp.size();++i){
        Integer h=v_scalar_product(Poly.Supp[i],p1);
        if(h<0 && (-h <height || height==0))
            height=-h;
    }
    // cout << "pot " << pot << endl;
    return 1.0/height;
}

template<typename Integer>
float VolumePotential(vector<Integer>& cand,Polytope<Integer>& Poly){

    // cout << "cand " << cand;
    if(Poly.FacetAreas.size()==0)
        Poly.FacetAreas=FacetAreas(Poly);    
    Integer vol=0;
    for(size_t i=0;i<Poly.Supp.size();++i){
        Integer height=v_scalar_product(cand,Poly.Supp[i]);
        // cout << -height << " ";
        if(height<0)
            vol+=(-height)*Poly.FacetAreas[i];
    }
    // cout << "vol " << vol << endl;
    return 1.0/vol;

}

template<typename Integer>
float AreaPotential(Polytope<Integer>& Poly){

    Integer Area=SurfaceArea(Poly);
    return 1.0/(float) Area;
}

template<typename Integer>
float AveragePotential(Polytope<Integer>& Poly){

    Integer Area=SurfaceArea(Poly);    
    return (float) Poly.Supp.size() /(float) Area;
}

template<typename Integer>
float NrSuppPotential(Polytope<Integer>& Poly){

    if(Poly.Supp.size()==0)
        Poly.Supp=SuppHyps(Poly.Gen);
    return (float) Poly.Supp.size();
}


template<typename Integer>
float   ComputePotential(Polytope<Integer> CandPoly, Polytope<Integer>& Poly){

    vector<Integer> cand=CandPoly.Gen.back();
    float pot=1.0;

    switch(ChoosePotential){
        case 0:
            pot=MaxHeightPotential(cand,Poly);
            break;
        case 1:
            pot=MinHeightPotential(cand,Poly);
            break;
        case 2:
            pot=HeightPotential(cand,Poly);
            break;
        case 3:
            pot=DistPotential(cand,Poly);
            break;
        case 4:
            pot=VolumePotential(cand,Poly);
            break;
        case 5:
            pot=AreaPotential(CandPoly);
            break;
        case 6:
            pot=AveragePotential(CandPoly);
            break;
        case 7:
            pot=NrSuppPotential(CandPoly); 
            break; 
    }
    if(FromSmall)
        return 1.0/pot;
    else
        return pot; 
}

// *************** Makers ***********************

template<typename Integer>
vector<vector<Integer> >  ReadMat(string project){
// reads one matrix from .in file

    string name_in=project;
    const char* file_in=name_in.c_str();
    ifstream in;
    in.open(file_in,ifstream::in);
    if (in.is_open()==false){
        cerr << "Cannot find input file" << endl;
        exit(1);
    }

    int nrows,ncols;
    in >> nrows;
    in >> ncols;

    if(nrows==0 || ncols==0){
        cerr << "Matrix empty" << endl;
        exit(1);
    }


    int i,j,entry;
    vector<vector<Integer> > result(nrows);

    for(i=0;i<nrows;++i)
        for(j=0;j<ncols;++j){
            in >> entry;
            result[i].push_back(entry);
        }

    vector<vector<Integer> > HB;
    if(!Normalize(result,HB)){
        cout << HB;
        cout << "============" << endl;
        cout << "Input polytope not normal" << endl;
        exit(1);
    }
    return HB;
}

template<typename Integer>
vector<vector<Integer> > UnitSimplex(int dim){ //cone dimension
    vector<vector<Integer> > StartPolyGen;       // make unit simplex homogenized
    StartPolyGen.resize(dim);
    for(int k=0;k<dim;++k){
        StartPolyGen[k].resize(dim);
        for(size_t p=0;p<dim;++p)
            StartPolyGen[k][p]=0;
        StartPolyGen[k][k]=1;
        StartPolyGen[k][dim-1]=1;
    }
    return StartPolyGen;
}

template<typename Integer>
vector<vector<Integer> > Ball(size_t& m){ //

    vector<vector<Integer> > Gen;
    while(m<50)
    {
    Gen.clear();
    int n=m;
    int s=n;
    for(int i=-s;i<=s;++i)
        for(int j=-s;j<=s;++j)
            for(int k=-s;k<=s;++k)
              for(int q=-s; q<=s;++q)
                if(i*i+j*j+k*k+q*q<=s){
                    vector<Integer> point(5);
                    point[0]=i;
                    point[1]=j;
                    point[2]=k;
                    point[3]=q;
                    point[4]=1;
                    Gen.push_back(point);
                }
    // cout << Gen;
    // cout << "-------------" << endl;
    vector<vector<Integer> > HB;
    if(!Normalize(Gen,HB)){
        cout << "radius " << m << " with " << Gen.size() << " points not normal" << endl;
        ++m;
    }
    else
    cout << "radius^2 " << s << " with " << Gen.size() << " points NORMAL" << endl;
        break;

    } // while
    return Gen;
} 

/* template<typename Integer>
vector<vector<Integer> > Ellipsoid(){ //
	vector<vector<Integer> > Gen;	
	bool Normal=false;
	while(!Normal){
		
		Integer a=1+rand()%7, b=1+rand()%7, c=1+rand()%7, d=1+rand()%7;
	
		Integer A=b*b*c*c*d*d;
		Integer B=a*a*c*c*d*d;
		Integer C=a*a*b*b*d*d;
		Integer D=a*a*b*b*c*c;
		Integer E=a*a*b*b*c*c*d*d;
	
		Gen.clear();
		for(int i=-a;i<=a;++i)
			for(int j=-b;j<=b;++j)
				for(int k=-c;k<=c;++k)
					for(int q=-d; q<=d;++q)
						if(A*i*i+B*j*j+C*k*k+D*q*q<=E){
							vector<Integer> point(5);
							point[0]=i;
							point[1]=j;
							point[2]=k;
							point[3]=q;
							point[4]=1;
							Gen.push_back(point);
						}
						// cout << Gen;
						// cout << "-------------" << endl;
						vector<vector<Integer> > HB;
						if(!Normalize(Gen,HB)){
							cout << "semi axes " << a << " " << b << " " << c << " " << d << " not normal" << endl;
						}
						else
							cout << "semi axes " << a << " " << b << " " << c << " " << d << " NORMAL" << endl;
						break;
						
	} // while
	return Gen;
} */

/* template<typename Integer>
vector<vector<Integer> > Ellipsoid(){ //
	vector<vector<Integer> > Gen;	
	bool Normal=false;
	while(!Normal){
		
		Integer a=1+rand()%3, b=1+rand()%4, c=1+rand()%5, d=1+rand()%7, e=1+rand()%7;
	
		Integer A=b*b*c*c*d*d*e*e;
		Integer B=a*a*c*c*d*d*e*e;
		Integer C=a*a*b*b*d*d*e*e;
		Integer D=a*a*b*b*c*c*e*e;
        Integer E=a*a*b*b*c*c*d*d;
		Integer F=a*a*b*b*c*c*d*d*e*e;
	
		Gen.clear();
		for(int i=-a-1;i<=a+1;++i)
			for(int j=-b-1;j<=b+1;++j)
				for(int k=-c-1;k<=c+1;++k)
					for(int q=-d-1; q<=d+1;++q)
                      for(int r=-e-1;r<=e+1;++r)
						if(A*i*i+B*j*j+C*k*k+D*q*q+E*r*r<=F){
							vector<Integer> point(6);
							point[0]=i;
							point[1]=j;
							point[2]=k;
							point[3]=q;
                            point[4]=r;
							point[5]=1;
							Gen.push_back(point);
						}
						// cout << Gen;
						// cout << "-------------" << endl;
						vector<vector<Integer> > HB;
						if(!Normalize(Gen,HB)){
							cout << "semi axes " << a << " " << b << " " << c << " " << d << " " << e << " not normal" << endl;
						}
						else
							cout << "semi axes " << a << " " << b << " " << c << " " << d << " " << e << " NORMAL" << endl;
						break;
						
	} // while
	return Gen;
} */

// 1426674743


template<typename Integer>
vector<vector<Integer> > Ellipsoid(){ //
	vector<vector<Integer> > Gen;	
	bool Normal=false;
	
	float den=7.0;
	while(!Normal){
		
		Integer a=1+rand()%3, b=1+rand()%4, c=1+rand()%5, d=1+rand()%7, e=1+rand()%7;
		Integer ca=rand()%7, cb=rand()%7, cc=rand()%7, cd=rand()%7, ce=rand()%7;
		
		Integer A=b*b*c*c*d*d*e*e;
		Integer B=a*a*c*c*d*d*e*e;
		Integer C=a*a*b*b*d*d*e*e;
		Integer D=a*a*b*b*c*c*e*e;
		Integer E=a*a*b*b*c*c*d*d;
		Integer F=a*a*b*b*c*c*d*d*e*e;
		
		Gen.clear();
		for(int i=-a-1;i<=a+1;++i)
			for(int j=-b-1;j<=b+1;++j)
				for(int k=-c-1;k<=c+1;++k)
					for(int q=-d-1; q<=d+1;++q)
						for(int r=-e-1;r<=e+1;++r){
							
							float ir=(float) i- (float) ca/den;
							float jr=(float) j- (float) cb/den;
							float kr=(float) k- (float) cc/den;
							float qr=(float) q- (float) cd/den;
							float rr=(float) r- (float) ce/den;
							
							// cout << ir << " " << jr <<  " " << kr << " "  << qr << " " << rr << endl;
							
							if(A*ir*ir+B*jr*jr+C*kr*kr+D*qr*qr+E*rr*rr<=F){
								vector<Integer> point(6);
								point[0]=i;
								point[1]=j;
								point[2]=k;
								point[3]=q;
								point[4]=r;
								point[5]=1;
								Gen.push_back(point);
							}
							
						}
							// cout << Gen;
							// cout << "-------------" << endl;
							// vector<vector<Integer> > HB;
							if(!IsNormal(Gen)){
								cout << "semi axes " << a << " " << b << " " << c << " " << d << " " << e << " not normal" << endl;
								cout << "center " << ca << " " << cb << " " << cc << " " << cd << " " << ce << endl;
								cout << Gen;
								cout << "-----------" << endl;
								// cout << HB;
								exit(0);
							}
							else{
								cout << "semi axes " << a << " " << b << " " << c << " " << d << " " << e << " NORMAL" << endl;
								cout << "center " << ca << " " << cb << " " << cc << " " << cd << " " << ce << endl;
							}
							// break;
							
	} // while
	return Gen;
} 



template<typename Integer>
vector<vector<Integer> > CrossPoly(int n){ // makes 3-dim cross polytope with axes n, n+1, n^2+n+1

    vector<vector<Integer> > CrossGen(6);
    for(size_t i=0;i<6;++i){
        CrossGen[i].resize(4,0);
        CrossGen[i][3]=1;
    }
    CrossGen[0][0]=n; // 2;//n;
    CrossGen[1][0]=-CrossGen[0][0];
    CrossGen[2][1]=n+1;// 2; // 3; // n+1;
    CrossGen[3][1]=-CrossGen[2][1];
    CrossGen[4][2]=n*n+n+1; // 2*n+1; // 6*n+1; // n*n+n+1;
    CrossGen[5][2]=-CrossGen[4][2];

    cout << CrossGen;

    cout << "----------------" << endl;

    vector<vector<Integer> > StartPolyGen;
    Normalize(CrossGen,StartPolyGen); // is automatically normal
    return StartPolyGen;
}



template<typename Integer>
vector<vector<Integer> > RectSimp(){ // makes 3-dim rectangular simplex with prime axes

    vector<vector<Integer> > CrossGen(4);
    for(size_t i=0;i<4;++i){
        CrossGen[i].resize(4,0);
        CrossGen[i][3]=1;
    }

    bool Normal;

    vector<vector<Integer> > StartPolyGen;

    do{
    vector<int> Primes=get3primes<int>();

    CrossGen[1][0]=Primes[0];
    CrossGen[2][1]=Primes[1];
    CrossGen[3][2]=Primes[2];

    cout << CrossGen;

    Normal=Normalize(CrossGen, StartPolyGen);
    cout << "----------------" << endl;
    if(!Normal){
        cout << "Not normal"<<endl;
        cout << "----------------" << endl;
    }

    }while(!Normal);

    // cout << StartPolyGen;

    // cout << "==============" << endl;
    return StartPolyGen;
}

template<typename Integer>
vector<vector<Integer> > RandPoly(int dim,int AddGen,int bound){ //cone dimension
// makes a random normal polytope

    vector<vector<Integer> > StartPolyGen;
    size_t nrgen=dim+rand()%(AddGen+1);      // dim <= nr_gen <= dim+AddGen
    StartPolyGen.resize(nrgen);

    bool Normal=false;
    vector<vector<Integer> > HB;

    do{
        for(size_t k=0;k<nrgen;++k){
            StartPolyGen[k].resize(dim);
            for(size_t p=0;p<dim-1;++p)
                StartPolyGen[k][p]=rand() % (bound+1);
            StartPolyGen[k][dim-1]=1;
        }

    Normal=Normalize(StartPolyGen,HB);
    }while(!Normal);

    return HB;
}

template<typename Integer>
vector<vector<Integer> > MakePara(int Dim, Integer DetBound){

    vector<vector<Integer> > PreErz(Dim-1);
    for(size_t i=0;i<PreErz.size();++i)
        PreErz[i].resize(Dim,0);

    Integer p;
    p=rand()%5+1;
    // for(p=rand()%5+1;!((p%2==1) || (p==2));p=rand()%5+1);

    Integer Det=p;
    PreErz[0][0]=p;

    for(size_t i=1;i<PreErz.size();++i){
        for(size_t j=0;j<i;++j)
            PreErz[i][j]=rand()%14-7;

        p=DetBound/Det;
        if(p==0) p=1; if(p>10) p=10;
        PreErz[i][i]=rand()%p+1;

        v_make_prime(PreErz[i]);
        Det*=PreErz[i][i];
    }

    vector<vector<Integer> > ConeGen(1);
    ConeGen[0].resize(Dim,0);
    vector<Integer> v(Dim);
    ConeGen[0][Dim-1]=1;
    for(size_t j=0;j<PreErz.size();++j){
        size_t kk=ConeGen.size();
        for(size_t k=0;k<kk;++k)
            ConeGen.push_back(v_add_1(PreErz[j],ConeGen[k]));
    }
    return ConeGen;

}



template<typename Integer>
bool VertexShrink(vector<vector<Integer> >& StartGen, vector<vector<Integer> >& Result){

    int dim=StartGen[0].size();

    vector<vector<Integer> > ExtremeRays, TestExtRays, HilbertBasis;

    if(!NormalizeAndExtRays(StartGen,HilbertBasis,ExtremeRays)){
        cout << "Warning: polytope for vertex shrink not normal or not of full rank" << endl;
        return false;
    }


    if(HilbertBasis.size()==dim){  // unimodular simplex, considered failure
        Result=HilbertBasis;  // but return the outcome
        return false;
    }

    // cout << "First test done" << endl;

    vector<vector<Integer> > TestHB;
    while(1){
        size_t nrExt=ExtremeRays.size();
        vector<Integer> v;
        size_t i=0;
        for(;i<nrExt;++i){
            v=ExtremeRays[0];
            for(size_t j=0;j<nrExt-1;++j)
                swap(ExtremeRays[j],ExtremeRays[j+1]);
            ExtremeRays.pop_back(); // if destructive will be reinserted below

            if(NormalizeAndExtRays(ExtremeRays,TestHB,TestExtRays))
                break; // non-destructive

            ExtremeRays.push_back(v); // restore destructive corner
        }
        if(i==nrExt){ //tight cone (cannot be unimodular simplex at this point)
            Result=HilbertBasis;
            return true;
        }

        // found a nondestructive corner

        if(TestHB.size()==dim){ // unimodular simplex
            Result=TestHB;
            return false;
        }

         HilbertBasis=TestHB; // proceed with smaller Hilbert basis
         ExtremeRays=TestExtRays; // and new extreme rays
    }
}


template<typename Integer>
bool shrink(vector<vector<Integer> >& StartGen, vector<vector<Integer> >& Result){

    int dim=StartGen[0].size();

    vector<vector<Integer> > ExtremeRays, TestExtRays, HilbertBasis;

    if(!NormalizeAndExtRays(StartGen,HilbertBasis,ExtremeRays)){
        cout << "Warning: polytope for shrink not normal or not of full rank" << endl;
        return false;
    }


    if(HilbertBasis.size()==dim){  // unimodular simplex, considered failure
        Result=HilbertBasis;  // but return the outcome
        return false;
    }

    // cout << "First test done" << endl;

    vector<vector<Integer> > TestHB;
    while(1){
        size_t nrExt=ExtremeRays.size();
        size_t i=0;
        for(;i<nrExt;++i){
            vector<vector<Integer> > ReducedHB;
            for(size_t k=0;k<HilbertBasis.size();++k)
                if(HilbertBasis[k]!=ExtremeRays[i])
                    ReducedHB.push_back(HilbertBasis[k]);
            if(NormalizeAndExtRays(ReducedHB,TestHB,TestExtRays))
                break; // non-destructive
        }
        if(i==nrExt){ //tight cone (cannot be unimodular simplex at this point)
            Result=HilbertBasis;
            return true;
        }

        // found a nondestructive corner

        if(TestHB.size()==dim){ // unimodular simplex
            Result=TestHB;
            return false;
        }

         HilbertBasis=TestHB; // proceed with smaller Hilbert basis
         ExtremeRays=TestExtRays; // and new extreme rays
    }
}


template<typename Integer>
vector<vector<Integer> > ShrunkPara(int dim, Integer DetBound){

    while(1){
        vector<vector<Integer> > ParaExtRays=MakePara(dim,DetBound);
        vector<vector<Integer> > StartGen;
        bool success=VertexShrink(ParaExtRays,StartGen);
        if(success)
            return(StartGen);
    }
}

// ***************************** Stop extension along line ***********************************

template<typename Integer>
bool CheckPointsOnLine(vector<vector<Integer> >& Gen){

    size_t nr_gen=Gen.size();
    if(nr_gen<11)
        return false;
    size_t dim=Gen[0].size();
    vector<Integer> diff(dim);
    for(size_t j=0;j<dim;++j)
        diff[j]=Gen[nr_gen-1][j]-Gen[nr_gen-2][j];
    for(long i=nr_gen-2;i>nr_gen-11;--i)
        for(size_t j=0;j<dim;++j)
            if(Gen[i][j]-Gen[i-1][j]!=diff[j])
                return false;
    return true;
}

// *************** Pick extension ***********************

template<typename Integer>
vector<Integer> PickJump(vector<vector<Integer> > Jumps, Polytope<Integer>& StartPoly, bool RandomExt){

    assert(Jumps.size()!=0);
    
    vector<Integer> MinPotentialJump;
    float Potential, MinPotential=0;
    bool first=true;
    
    if(RandomExt){
        size_t pick=rand()%Jumps.size();
        MinPotentialJump=Jumps[pick];
    }
    else{
        
        size_t MinPotentialIndex;
        
        Polytope<Integer> TestPoly;
        size_t ts=StartPoly.Gen.size()+1;
        TestPoly.Gen=StartPoly.Gen;
        TestPoly.Gen.resize(ts);  
                
        for(size_t i=0;i<Jumps.size();++i){
            TestPoly.Gen[ts-1]=Jumps[i];  
            Potential=ComputePotential(TestPoly,StartPoly);
            if(first || Potential < MinPotential){
                first=false;
                MinPotential=Potential;
                MinPotentialIndex=i;                
            }
        }
        
        MinPotentialJump=Jumps[MinPotentialIndex];
    }
        
    if(HeightOverPolytope(MinPotentialJump,StartPoly)>1 || Print){
	Integer HOP=HeightOverPolytope(MinPotentialJump,StartPoly);
        cout << "ext ht " << HOP << " pot " << MinPotential;

        cout << " by " << MinPotentialJump;
        if(Print) cout << "Now " << StartPoly.Gen.size()+1 << " lattice oints" << endl;
        cout << "hts ";
        for(size_t k=0;k<StartPoly.Supp.size();++k){
            Integer test= v_scalar_product(MinPotentialJump,StartPoly.Supp[k]);
            cout << -test  << " ";
        }
        cout << endl;
    if(FromSmall && HOP > MinHeight)
	   cout << "Nonminimal height extension" << endl;
	
        cout << "----" << endl;       
    }
    
    return MinPotentialJump;
} 


// *************** Ht1 extensions ***********************

template<typename Integer>
void  DoHeightOneExtensions(Polytope<Integer>& StartPoly, int TooLarge){

    int dim=StartPoly.Gen[0].size();

    while((int) StartPoly.Gen.size()<=TooLarge){  // we take heught 1 extensions by random choice as ling as possible

        StartPoly.Supp=SuppHyps(StartPoly.Gen);
        vector<vector<Integer> > CandPolySupp=StartPoly.Supp;

        for(size_t i=0;i< CandPolySupp.size();++i) // move hyperplanes outward
            CandPolySupp[i][dim-1]++;

        vector<vector<Integer> > Ht1Points=Deg1FromSupps(CandPolySupp);

        vector<vector<Integer> > ReallyHt1;  // discard points in the start polytope
        for(size_t i=0;i<Ht1Points.size();++i){
            if(!InCone<Integer>(Ht1Points[i], StartPoly.Supp))
                ReallyHt1.push_back(Ht1Points[i]);
        }

        if(ReallyHt1.empty()){
            // cout << "Empty" << endl;
            break;
        }

        vector<Integer> Jump=PickJump(ReallyHt1,StartPoly,true);
        StartPoly.Gen.push_back(Jump); // extend by picked vector
    }  // while
}


// *************** Higher extensions ***********************


template<typename Integer>
void MakeJumpCandidates(list<vector<Integer> >& CandJumps1, Polytope<Integer>& StartPoly, bool Verification){

    int dim=StartPoly.Gen[0].size();

    StartPoly.Supp=SuppHyps(StartPoly.Gen);

    vector<Integer> Widths(StartPoly.Supp.size(),0);

    // compute widths

    for(size_t i=0;i<StartPoly.Supp.size();++i){
      for(size_t g=0;g<StartPoly.Gen.size();++g){
          Integer test=v_scalar_product(StartPoly.Supp[i],StartPoly.Gen[g]);
          if(test>Widths[i])
              Widths[i]=test;
      }
    }

    if(Print) cout << "Widths " << Widths;
    

    // compute offsets of hyperplanes
    
    Integer factor=1;
    if(Verification)
        factor=dim-3; // dim is cone dimension !!!!!

    for(size_t i=0;i<Widths.size();++i){
        Widths[i]*=factor;
        Widths[i]++;
    }

    if(dim>=7 &&!Verification){
        for(size_t i=0;i<StartPoly.Supp.size();++i)
            Widths[i]=5+Widths[i]/2;
    }
    
    if(AllExtensions)
        for(size_t i=0;i<StartPoly.Supp.size();++i)
            Widths[i]=Widths[i]/3+1;
        

    if(Print) cout << "Distances " << Widths;
    
    // adjust hyperplanes

    vector<vector<Integer> > CandPolySupp=StartPoly.Supp; 
    for(size_t i=0;i< CandPolySupp.size();++i){
        CandPolySupp[i][dim-1]+=Widths[i];
    }

    // cout << CandPolySupp;

    vector<vector<Integer> > CandJumps;

    if(true){ // (dim<=5){
        CandJumps=Deg1FromSupps(CandPolySupp);
    }
    else{
        CandJumps=Deg1FromSupps_mpz(CandPolySupp);
    }

    scramble(CandJumps);  // is reproducable --- will be sorted later if necessary

    // cout << CandJumps;

    if(CandJumps.size()>100000 || Print)
        cout << "CSize " << CandJumps.size() << endl;

    MinHeight=0;

    for(size_t i=0;i<CandJumps.size();++i){
        if(!InCone<Integer>(CandJumps[i], StartPoly.Supp)){
            CandJumps1.push_back(CandJumps[i]);
            if(MinHeight==0 || HeightOverPolytope(CandJumps[i], StartPoly) < MinHeight)
                MinHeight=HeightOverPolytope(CandJumps[i], StartPoly);
        }
    }

    if(MinHeight>1 || Print)
        cout << "MinHeight " << MinHeight << endl;
}
  

template<typename Integer>
void FindExtension(list<vector<Integer> >& CandJumps1, Polytope<Integer>& StartPoly, bool Verification, bool& Cotight){

    Cotight=false;

    typename list<vector<Integer> >::iterator q;
    Polytope<Integer> TestPoly;

    // for(q=CandJumps1.begin();q!=CandJumps1.end();++q)
    //    cout << *q;
        
    StartPoly.ExtRays=ExtremeRays(StartPoly.Gen);
    TestPoly.Gen=StartPoly.Gen;
    vector<Integer> center=TestPoly.Gen[0];
    for(size_t i=1;i<TestPoly.Gen.size();++i)
        center=v_add_1(center,TestPoly.Gen[i]);
        
    size_t ts=TestPoly.Gen.size();
    TestPoly.Gen.resize(ts+1);
    int dim=StartPoly.Gen[0].size();

    
    bool ExtensionFound=false, NoJumps;
    vector<Integer> MinPotentialJump;
    

    if(CheckIllumination){
	   StartPoly.Illu.resize(StartPoly.Supp.size(),false);
    }
    
    vector<vector<Integer> > VAExt;    
    vector<vector<Integer> > Jumps;
    vector<vector<Integer> > NormalExt;
       
    for(q=CandJumps1.begin();q!=CandJumps1.end();++q){
    
        TestPoly.Gen[ts]=*q;
        
        bool OneMore=false, Normal=false, VA=false;

        TestPoly.LattPoints=Deg1FromGens(TestPoly.Gen);
        if(TestPoly.LattPoints.size()==StartPoly.Gen.size()+1)
            OneMore=true;
            
        if(!AllExtensions && !OneMore)
            continue;

        if(IsNormalPolytope(TestPoly.Gen)){
            Normal=true;
            ExtensionFound=true;
        }
        
        if(FromSmall && !Normal && HeightOverPolytope(*q,StartPoly)==MinHeight){
            cout << "Min height, but not jump "<< *q;
            vector<vector<Integer> > TestHB;
            Normalize(TestPoly.Gen,TestHB);
            bool Deg3=false;
            for(size_t i=0;i<TestHB.size();++i)
                if(TestHB[i][dim-1]>2)
                    Deg3=true;
            if(!Deg3)
                cout << "No witness degree > 2" << endl;            
        }
        
        if(!Normal){
            if(IsVeryAmple(TestPoly.Gen)){
                // VAExtensionFound=true;
                VA=true;            
            }
        }
        
        /* if(FromSmall && !Normal)
            cout << "One more LP, but not jump "<< *q; */
            
        if(OneMore && Normal){ 
            Jumps.push_back(*q);
	        if(CheckIllumination){
			   for(size_t i=0;i<StartPoly.Supp.size();++i)
			    if(v_scalar_product(StartPoly.Supp[i],*q)<0)
				  StartPoly.Illu[i]=true;
		   }	
	   }

        if(OneMore && VA)
            VAExt.push_back(*q);
            
       if(AllExtensions && Normal)
            NormalExt.push_back(*q);
                          
    } // loop over CandJumps
    
    NoJumps=Jumps.empty();
    
    if(ExtensionFound){
        if(AllExtensions)
            MinPotentialJump=PickJump(NormalExt,StartPoly,RandomExt);
        else        
            MinPotentialJump=PickJump(Jumps,StartPoly,RandomExt);
    }
    
    if(CheckIllumination){
		size_t DarkFacets=0;
		for(size_t i=0;i<StartPoly.Illu.size();++i)
		    if(!StartPoly.Illu[i]){
			cout << "Facet " << StartPoly.Supp[i] << " not illuminated" << endl;
			DarkFacets++;
		    }
		  if(DarkFacets >= dim-1){
			for(size_t i=0;i<StartPoly.ExtRays.size();++i){
			    size_t j;
			    for(j=0;j<StartPoly.Supp.size();++j)
				if(v_scalar_product(StartPoly.Supp[j],StartPoly.ExtRays[i])==0 && StartPoly.Illu[j])
				    break;
			    if(j==StartPoly.Supp.size())
				cout << "Vertex " << StartPoly.ExtRays[i] << " not illuminated" << endl;
		       }	
		}
		if(DarkFacets>0){
		    cout << StartPoly.ExtRays;
		    cout << "*****" << endl;
		}
    }	
    
    if(NoJumps){
    
        Cotight=true;
        
        if(Verification){
            cout << "Verification successful" << endl << flush;
            
            cout << "Checking minimality" << endl;
            vector<vector<Integer> > Result;
            
            shrink(StartPoly.Gen,Result);
            if(Result.size()==StartPoly.Gen.size())
                cout << "Isolated state found" << endl;
            else
                cout << "C<an be shrunk" << endl;
            return;
            
            if(!VAExt.empty()){
                cout << "Very ample extension by " << endl;
                cout << VAExt;
            }
                    
        }
        
        cout << "Cotight" << endl;
        cout << StartPoly.Gen << endl;
        
        cout <<"Starting from " << endl;
        cout << Origin << endl << flush;
        
        cout << endl << "Now verifying" << endl;
        return;
    }
    
    if(Verification){
        cout << "Verification failed" << endl << endl << flush;
        return;    
    }
    

    StartPoly.Gen.push_back(MinPotentialJump);
    if(AllExtensions){
        Normalize(StartPoly.Gen,StartPoly.Gen);
        cout << "Now " << StartPoly.Gen.size() << " LP" << endl;
    }
}


// *************** Main  ***********************

int main(int argc, char* argv[])
{
    
    time_t ticks;
    
    if(argc<4){
        cout << "Not eneough parameters" << endl;
        exit(1);
    }
    
    Polytope<Integer> StartPoly;
    bool ReadTicks=false;

    string options(argv[1]);
    for(size_t i=0;i<options.size();++i){
        switch (options[i]) {
            case 'u':
                break;
            case 'r':
                ChooseStart=1;
                break;
            case 'i':
                ChooseStart=3;
                StartPoly.Gen=ReadMat<Integer>(string(argv[2]));
                break;
            case 's': 
                ChooseStart=2;  
                break;
            case 'c': 
                ChooseStart=4;  
                break;
            case 'S': 
                ChooseStart=5;  
                break;
            case 'b': 
                ChooseStart=6;  
                break;
			case 'e': 
				ChooseStart=7;  
				break;
            case 'B': 
                BreakReg=true;  
                break;
            case '-': 
                FromSmall=true;
                break;
            case '1':
                HeightOneFirst=true;
                break;
            case 'R':
                RandomExt=true;
                break;
            case 'P':
                Print=true;
                break;
            case 'm':
                ChoosePotential=1;
                break;
            case 'E':
                ChoosePotential=3;
                break;
            case 'H':
                ChoosePotential=2;
                break;
            case 'v':
                ChoosePotential=4;
                break; 
            case 'A':
                ChoosePotential=5;
                break;            
            case 'a':
                ChoosePotential=6;
                    break;
            case 'f':
                ChoosePotential=7;
                break;
            case 't':
                ReadTicks=true;
                break;
            case 'V':
                Verify=true;
                break;
            case 'X':
                AllExtensions=true;
                break;  
	    case 'I':
                CheckIllumination=true;
                break;                                                                                             
            default:
                cout<<"Fatal error: Unknown option -"<<options[i]<<endl;
                exit(1);
        }
    }

    
    if(ReadTicks){
        sscanf(argv[argc-1],"%ld",&ticks);    
        cout << "input ticks " << ticks << endl;
        srand(ticks);   
    }
    else{   
        srand(time(&ticks));
        cout << "time ticks " << ticks << endl;
    }
    
    if(FromSmall)
        HeightOneFirst=true;
    
    cout << "Polytope generation method " << ChooseStart << " From small " << FromSmall 
        << " HeightOneFirst " << HeightOneFirst << " Random " << RandomExt << " BreakReg " << BreakReg 
                  << " Chosen Potential " << ChoosePotential << " Print " << Print << endl;

    int dim1, dim;
    if(ChooseStart==3){
     dim1=StartPoly.Gen[0].size()-1;
    }    
    else{    
        sscanf(argv[2],"%d",&dim1); // polytope dimension
    }
    cout << "Dimension " << dim1 << endl;
    dim=dim1+1; // cone dimension
    
    // if(dim1==3)
    //     AllExtensions=true;
    
    // omp_set_num_threads(dim1);
    omp_set_num_threads(1);
    
    int TooLarge; // stop extension if number of lattice points exceeds this number
    sscanf(argv[3],"%d",&TooLarge); 
    cout << "TooLarge " << TooLarge << endl;
    
    int AddGen;  // random polytope has <= dim+AddGen generators    
    int bound; // bounds the random entries
    if(ChooseStart==1){
        if(argc<6){
            cout << "Not eneough parameters" << endl;
            exit(1);            
        }
        sscanf(argv[4],"%d",&bound);
        cout << "Bound " << bound << endl;
        sscanf(argv[5],"%d",&AddGen);
        cout << "Add Gen " << AddGen << endl;
    }
    
    int intDetBound;    

    if(ChooseStart==2){
        if(argc<5){
            cout << "Not eneough parameters" << endl;
            exit(1);            
        }
        sscanf(argv[4],"%d",&intDetBound);
        
        cout << "DetBound " << intDetBound << endl;
    }
    
    if(BreakReg)
        AllExtensions=false;

    if(CheckIllumination){
	AllExtensions=false;
	BreakReg=false;
    }
        
    
    Integer DetBound=intDetBound;
    
    // vector< vector<Integer > > StartPolySupp;
    
// ************************************ Now the main loop in which start polytopes are generated

    size_t PolyCount=0;
    
    while(1){
        
        PolyCount++; // cout << "Poly " << PolyCount << endl;
        
        if(ChooseStart==3 && PolyCount>=2)  // in case of input polytope there is only one
            break;
        
        if(PolyCount%500==0)
            cout << "****** Polytope " << PolyCount << endl;
        
        switch(ChooseStart){
            case 0:
                StartPoly.Gen=UnitSimplex<Integer>(dim);              
                break;
            case 1:
                StartPoly.Gen=RandPoly<Integer>(dim,AddGen,bound);
                TooLarge=StartPoly.Gen.size();
                break;
            case 2:
                StartPoly.Gen=ShrunkPara<Integer>(dim,DetBound);
                break;
            case 3:  // read from file
                break;
            case 4:
                StartPoly.Gen=CrossPoly<Integer>(PolyCount);
                TooLarge=StartPoly.Gen.size();
                break;
            case 5:
                StartPoly.Gen=RectSimp<Integer>();
                TooLarge=StartPoly.Gen.size();
                break;
            case 6:
                StartPoly.Gen=Ball<Integer>(PolyCount);
                TooLarge=StartPoly.Gen.size();
                break;
			case 7:
				StartPoly.Gen=Ellipsoid<Integer>();
				TooLarge=StartPoly.Gen.size();
				break;
            default:
                cout<<"Fatal error: Option nnot yet implemented" <<endl;
                exit(1);
        }
                        
        // cout << StartPoly.Gen;
        // cout << "-------------" << endl;
        
        if(ChooseStart==0){        
            int r=rand()%10+1;  // we first extend by a number of ht1 extension to get randomness
            for(int j=0;j<r;++j){
                DoHeightOneExtensions(StartPoly,dim1+r);
            }
        } 
     
        
        if(!HeightOneFirst)
                cout << "+++ Exrending Poly " << PolyCount << " #LP " << StartPoly.Gen.size() << endl;
                
                        
        if(!(HeightOneFirst || RandomExt)){
            cout << ExtremeRays(StartPoly.Gen);        
        }
        
        Origin=StartPoly.Gen;
        
        
        // ************************************ Now the inner loop in which the start polytopes is extended
        
        while((int) StartPoly.Gen.size()<=TooLarge){ // stop when too large 
                
        
            if(HeightOneFirst){  // try height 1 extensions first    
                DoHeightOneExtensions(StartPoly,TooLarge);
            } 
             
            
            if((int) StartPoly.Gen.size()>TooLarge){ // maximal extension reached by ht 1 extensions
                // cout << "Too large" << endl;
                break;
            }
            
            if(HeightOneFirst)
                cout << "No ht 1 ext Poly " << PolyCount << " #LP " << StartPoly.Gen.size() << endl;            
                
            if(Print) cout << "Poly " << PolyCount << " #LP " << StartPoly.Gen.size() << endl;
            
            list<vector<Integer> >   CandJumps1;
            
            // vector<vector<Integer> > StartPolySupp;
            
            bool Cotight, WeaklyCotight;            
            MakeJumpCandidates(CandJumps1, StartPoly,false);            
            FindExtension(CandJumps1,StartPoly,false,WeaklyCotight);
            
            if(WeaklyCotight){            
                list<vector<Integer> > CandJumpsVer; 
                bool SavePrint=Print;
                Print=true;
                bool SaveBreakReg=BreakReg;
                BreakReg=false;
                bool SaveAllExtensions=AllExtensions;
                AllExtensions=false;           
                MakeJumpCandidates(CandJumpsVer, StartPoly,true);            
                FindExtension(CandJumpsVer,StartPoly,true,Cotight);
                Print=SavePrint;
                AllExtensions=SaveAllExtensions;
                BreakReg=SaveBreakReg;
                MatMail(StartPoly.Gen);
                if(Cotight)
                    SaveMat(StartPoly.Gen);
            }
            
            if(WeaklyCotight)
                break;
            
            bool single_direction=CheckPointsOnLine(StartPoly.Gen);
            
            if(single_direction){
                cout << "STOPPED: ext along line" << endl;
                break;
            }
                
            
        } // while over extension
        
        // if(PolyCount >=1) break;
        
    } // while over polytopes        
    
    exit(0);
}
