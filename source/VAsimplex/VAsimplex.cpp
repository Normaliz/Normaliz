#include <stdlib.h>
#include <vector>
#include <fstream>
#include <omp.h>
using namespace std;

#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/cone.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/cone_property.h"
#include "libnormaliz/integer.h"
using namespace libnormaliz;

typedef long long Integer;

Cone<Integer> nonnormal_rand_simplex(size_t dim, long bound){
    
    vector<vector<Integer> > vertices(dim+1,vector<Integer> (dim));    
    while(true){  // an eternal loop ...
        for(size_t i=0;i<=dim;++i){
            for(size_t j=0;j<dim;++j)
                vertices[i][j]=rand()%(bound+1);            
        }
        Cone<Integer> Simplex(Type::polytope,vertices);
        // we must check the rank and normality
        if(Simplex.getRank()==dim+1 && !Simplex.isDeg1HilbertBasis())
            return Simplex;        
    }
    vector<vector<Integer> > dummy_gen(1,vector<Integer>(1,1)); // to make the compiler happy
    return Cone<Integer>(Type::cone,dummy_gen); 
}

template<typename Integer>
bool  IsVeryAmple(const vector<vector<Integer> >& Gen, size_t& nr_normal_corners){

    size_t dim=Gen[0].size();

    ConeProperties CandPolyWantedDeg1Elements;
    CandPolyWantedDeg1Elements.set(ConeProperty::Deg1Elements);
    
    // Matrix<Integer>(Gen).pretty_print(cout);    
    vector<vector<Integer> > GenCopy=Gen;    
    Cone<Integer> Cand(Type::integral_closure, GenCopy);
    Cand.compute( CandPolyWantedDeg1Elements);
    vector<vector<Integer> > ExtRays=Cand.getExtremeRays();
    vector<vector<Integer> > LattP=Cand.getDeg1Elements();
    
    // cout << "Gitterpunkte berechnet" << endl;
    
    vector<vector<Integer> > CornerCone;
    CornerCone.reserve(LattP.size()-1);
    
    size_t i;
    // cout << "========================" << endl;
    
    nr_normal_corners=0;
    
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
        // cout <<"---------" << endl;
        

        Cone<Integer> Corner(Type::cone,CornerCone);
        if(Corner.isIntegrallyClosed())
            nr_normal_corners++;               
    }
    return nr_normal_corners==ExtRays.size();

}


int main(int argc, char* argv[]){

    time_t ticks;
    srand(time(&ticks));
    cout << "Seed " <<ticks << endl;  // we may want to reproduce the run

    size_t polytope_dim;
    sscanf(argv[1],"%uld",&polytope_dim);
    size_t cone_dim=polytope_dim+1;
    long bound;
    sscanf(argv[2],"%d",&bound);
    vector<Integer> grading(cone_dim,0); // at some points we need the explicit grading
    grading[polytope_dim]=1;
    
    size_t nr_simplex=0; // for the progress report
    size_t nr_normal_corners;
    vector<size_t> corner_stat(polytope_dim+2,0);
    
    while(true){
        omp_set_num_threads(1);
        Cone<Integer> Candidate=nonnormal_rand_simplex(polytope_dim,bound);
        nr_simplex++;
        if(nr_simplex%10000 ==0)
                cout << "simplex " << nr_simplex << " " << corner_stat;
        if(!IsVeryAmple(Candidate.getExtremeRays(),nr_normal_corners)){
            corner_stat[nr_normal_corners]++;
            continue;
        }
        cout << "very ample simplex found" << endl;
        cout << "Vertices" << endl;
        Candidate.getExtremeRaysMatrix().pretty_print(cout);
        cout << "Number of lattice points = " << Candidate.getNrDeg1Elements();
        cout << " Multiplicity = " << Candidate.getMultiplicity() << endl; 
        break;        
    } // end while
}  //end main
