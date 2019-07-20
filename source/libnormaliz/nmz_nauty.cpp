/*
 * Normaliz
 * Copyright (C) 2007-2019  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

//---------------------------------------------------------------------------

#include <boost/dynamic_bitset.hpp>
#include<map>

#include "libnormaliz/integer.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/nmz_nauty.h"
#include "libnormaliz/vector_operations.h"

#ifdef NMZ_NAUTY

// #define MAXN 5000    /* Define this before including nauty.h */
// we use dynamic allocation

#include <nauty/nauty.h>

namespace libnormaliz {
using namespace std;

vector<vector<long> > CollectedAutoms;

void getmyautoms(int count, int *perm, int *orbits,
               int numorbits, int stabvertex, int n){
    int i;
    vector<long> this_perm(n);
    for (i = 0; i < n; ++i) this_perm[i] = perm[i];
    CollectedAutoms.push_back(this_perm);
}

template<typename Integer>
void makeMM(BinaryMatrix& MM, const Matrix<Integer>& Generators,
                const Matrix<Integer>& LinForms, bool zero_one){
    
    key_t i,j,k;
    size_t mm=Generators.nr_of_rows();
    size_t nn=LinForms.nr_of_rows();
    Matrix<long> MVal(mm,nn);

    long new_val=0;
    Integer val;
    map<Integer,long> Values;
    for(i=0;i<mm; ++i){
        for(j=0;j<nn;++j){
            val=v_scalar_product(Generators[i],LinForms[j]);
            if(zero_one && val!=0)
                val=1;
            auto v=Values.find(val);
            if(v!=Values.end()){
                MVal[i][j]=v->second;
            }
            else{
                Values[val]=new_val;
                MVal[i][j]=new_val;
                new_val++;
            }
        }
    }
    
    MM.set_offset((long) 0);
    for(i=0;i<mm; ++i){
        for(j=0;j<nn;++j)
            MM.insert(MVal[i][j],i,j);
    }
}

template<typename Integer>
void makeMMFromGensOnly_inner(BinaryMatrix& MM, const Matrix<Integer>& Generators, const Matrix<Integer>& SpecialLinForms){
    
    size_t mm=Generators.nr_of_rows();
    size_t dim=Generators.nr_of_columns();    
    
    Matrix<Integer> ScalarProd(dim,dim);
    
    for(size_t i=0;i<mm;++i){
        for(size_t j=0;j<dim;++j){
            for(size_t k=0;k<dim;++k){
                ScalarProd[j][k]+=Generators[i][j]*Generators[i][k]; 
            }
        }        
    }
    
    Integer dummy;    
    Matrix<Integer> SPInv=ScalarProd.invert(dummy);    
    Matrix<Integer> LinForms=Generators.multiplication(SPInv);    
    LinForms.append(SpecialLinForms);    
    
    makeMM(MM,Generators,LinForms,false);    
}

template<typename Integer>
void makeMMFromGensOnly(BinaryMatrix& MM, const Matrix<Integer>& Generators,
                        const Matrix<Integer>& SpecialLinForms){
    
    Matrix<mpz_class> Generators_mpz;      // we go through mpz_class since taking inverse matrices
    convert(Generators_mpz,Generators);    // is extremely critical 
    Matrix<mpz_class> SpecialLinForms_mpz;
    convert(SpecialLinForms_mpz,SpecialLinForms);
    makeMMFromGensOnly_inner(MM,Generators_mpz,SpecialLinForms_mpz);
}

template<>
void makeMMFromGensOnly(BinaryMatrix& MM, const Matrix<renf_elem_class>& Generators, 
                        const Matrix<renf_elem_class>& SpecialLinForms){
    
    makeMMFromGensOnly_inner(MM,Generators,SpecialLinForms);    
}


template<typename Integer>
nauty_result compute_automs_by_nauty_Gens_LF(const Matrix<Integer>& Generators, size_t nr_special_gens, 
                            const Matrix<Integer>& LinForms, const size_t nr_special_linforms, bool zero_one){
    
    CollectedAutoms.clear();
    
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    
    options.userautomproc = getmyautoms;
    options.getcanon = TRUE;
    
    int n,m;
    
    options.writeautoms = FALSE;
    options.defaultptn = FALSE;
    
    size_t mm=Generators.nr_of_rows();
    size_t nn=LinForms.nr_of_rows();
    
    BinaryMatrix MM(mm,nn);
    makeMM(MM,Generators,LinForms,zero_one);
        
    size_t ll=MM.nr_layers();
        
    size_t layer_size=mm+nn;
    n=ll*layer_size;
    m = SETWORDSNEEDED(n);
    
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,n,m,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    
    EMPTYGRAPH(g,m,n);
    
    key_t i,j,k;
    
    for(i=0;i<layer_size;++i){   // make vertical edges over all layers
        for(k=1;k<ll;++k)
            ADDONEEDGE(g,(k-1)*layer_size+i,k*layer_size+i,m);
    }
    
    for(i=0;i<mm;++i){   // make horizontal edges layer by layer
        for(j=0;j<nn;++j){
            for(k=0;k<ll;++k){
                if(MM.test(i,j,k))  // k is the number of layers below the current one
                    ADDONEEDGE(g,k*layer_size+i,k*layer_size+mm+j,m);
            }
        }
    }           
    
    for(int ii=0;ii<n;++ii){ // prepare partitions
        lab[ii]=ii;
        ptn[ii]=1;
    }
    
    for(k=0;k<ll;++k){ // make partitions layer by layer
        ptn[k*layer_size+ mm-1]=0; // row vertices in one partition
        // for(size_t s=0; s< nr_special_gens;++s) // speciall generators in extra partitions (makes them fixed points)
        //    ptn[k*layer_size+s]=0; unclear 
        ptn[(k+1)*layer_size-1]=0; // column indices in the next
        for(size_t s=0; s< nr_special_linforms;++s) // special linear forms in extra partitions
            ptn[(k+1)*layer_size-2-s]=0;            
    } 

    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    
    // vector<vector<long> > AutomsAndOrbits(2*CollectedAutoms.size());
    // AutomsAndOrbits.reserve(2*CollectedAutoms.size()+3);
    
    nauty_result result;

    for(k=0;k<CollectedAutoms.size();++k){
        vector<key_t> GenPerm(mm);
        for(i=0;i<mm;++i)
            GenPerm[i]=CollectedAutoms[k][i];
        result.GenPerms.push_back(GenPerm);
        vector<key_t> LFPerm(nn-nr_special_linforms);  // we remove the special linear forms here
        for(i=mm;i<mm+nn-nr_special_linforms;++i)
            LFPerm[i-mm]=CollectedAutoms[k][i]-mm;
        result.LinFormPerms.push_back(LFPerm);
    }    
    
    vector<key_t> GenOrbits(mm);
    for(i=0;i<mm;++i)
        GenOrbits[i]=orbits[i];
    result.GenOrbits=GenOrbits;
    
    vector<key_t> LFOrbits(nn-nr_special_linforms); // we remove the special linear forms here
    for(i=0;i<nn-nr_special_linforms;++i)
        LFOrbits[i]=orbits[i+mm]-mm;
    result.LinFormOrbits=LFOrbits;
 
    result.order=mpz_class(stats.grpsize1);
    
    vector<key_t> row_order(mm), col_order(nn);
    for(key_t i=0;i<mm;++i)
        row_order[i]=lab[i];
    for(key_t i=0;i<nn;++i)
        col_order[i]=lab[mm+i]-mm;
    
    result.CanLabellingGens=row_order;
    
    result.CanType=MM.reordered(row_order,col_order);
    
    nauty_freedyn();
    
    return result;
        
}

//====================================================================

template<typename Integer>
nauty_result compute_automs_by_nauty_FromGensOnly(const Matrix<Integer>& Generators,  size_t nr_special_gens,
            const Matrix<Integer>& SpecialLinForms){
    
    size_t mm=Generators.nr_of_rows();
    
    size_t nr_special_linforms=SpecialLinForms.nr_of_rows();
    
    // LinForms.append(SpecialLinForms);
    
    CollectedAutoms.clear();
    
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    
    options.userautomproc = getmyautoms;
    options.getcanon = TRUE;
    
    int n,m;
    
    options.writeautoms = FALSE;
    options.defaultptn = FALSE;
    
    
    BinaryMatrix MM(mm,mm+nr_special_linforms);
    makeMMFromGensOnly(MM,Generators,SpecialLinForms);
        
    size_t ll=MM.nr_layers();
        
    size_t layer_size=mm+nr_special_linforms;
    n=ll*layer_size;        // total number of vertices
    m = SETWORDSNEEDED(n);
    
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC2(graph,cg,cg_sz,n,m,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    
    EMPTYGRAPH(g,m,n);
    
    key_t i,j,k;
    
    for(i=0;i<layer_size;++i){   // make vertical edges over all layers
        for(k=1;k<ll;++k)
            ADDONEEDGE(g,(k-1)*layer_size+i,k*layer_size+i,m);
    }
    
    for(i=0;i<mm;++i){   // make horizontal edges layer by layer
        for(j=0;j<=i;++j){  // take lower triangularr matrix inclcudung diagonal
            for(k=0;k<ll;++k){
                if(MM.test(i,j,k))  // k is the number of layers below the current one
                    ADDONEEDGE(g,k*layer_size+i,k*layer_size+j,m);
            }
        }
    }  
    
    // we add the edges that connect generators and special linear forms
    for(i=mm;i<mm+nr_special_linforms;++i){
        for(j=0;j<mm;++j){
            if(MM.test(j,i,k))  // here we use that the special linear forms appear in columns: i <--> j
                ADDONEEDGE(g,k*layer_size+i,k*layer_size+j,m);
        }
        
    }
    
    for(int ii=0;ii<n;++ii){ // prepare partitions
        lab[ii]=ii;  // label of vertex
        ptn[ii]=1;   // indicatorvector for partitions: 0 indicates end of partition
    }
    
    for(k=0;k<ll;++k){ // make partitions layer by layer
        ptn[k*layer_size+ mm-1]=0; // row vertices in one partition
        // for(size_t s=0; s< nr_special_gens;++s) // speciall generators in extra partitions (makes them fixed points)
        //     ptn[k*layer_size+s]=0;   correct ??????
        // ptn[(k+1)*layer_size-1]=0; // column indices in the next NO COLIMN INDICES IN THIS VARIANT
        for(size_t s=0; s< nr_special_linforms;++s) // special linear forms in extra partitions
            ptn[(k+1)*layer_size-2-s]=0;            
    } 

    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,cg);
    
    nauty_result result;

    for(k=0;k<CollectedAutoms.size();++k){
        vector<key_t> GenPerm(mm);
        for(i=0;i<mm;++i)  // remove special lion forms
            GenPerm[i]=CollectedAutoms[k][i];
        result.GenPerms.push_back(GenPerm);
    }    
    
    vector<key_t> GenOrbits(mm);
    for(i=0;i<mm;++i)
        GenOrbits[i]=orbits[i]; // remove special lion forms
    result.GenOrbits=GenOrbits;
 
    result.order=mpz_class(stats.grpsize1);
    
    vector<key_t> row_order(mm);
    for(key_t i=0;i<mm;++i)
        row_order[i]=lab[i];
    
    result.CanLabellingGens=row_order;
    
    nauty_freedyn();
    
    // CanType=MM.reordered(row_order,col_order);
    
    cout << "ORDER " << result.order << endl;
    
    return result;
    
}


#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template nauty_result compute_automs_by_nauty_Gens_LF(const Matrix<long>& Generators, size_t nr_special_gens, 
                        const Matrix<long>& LinForms,const size_t nr_special_linforms, bool zero_one);

template nauty_result compute_automs_by_nauty_FromGensOnly(const Matrix<long>& Generators,  size_t nr_special_gens,
            const Matrix<long>& SpecialLinForms);
#endif // NMZ_MIC_OFFLOAD
template nauty_result compute_automs_by_nauty_Gens_LF(const Matrix<long long>& Generators, size_t nr_special_gens, 
                        const Matrix<long long>& LinForms,const size_t nr_special_linforms, bool zero_one);

template nauty_result compute_automs_by_nauty_FromGensOnly(const Matrix<long long>& Generators,  size_t nr_special_gens,
            const Matrix<long long>& SpecialLinForms);

template nauty_result compute_automs_by_nauty_Gens_LF(const Matrix<mpz_class>& Generators, size_t nr_special_gens, 
                        const Matrix<mpz_class>& LinForms,const size_t nr_special_linforms, bool zero_one);

template nauty_result compute_automs_by_nauty_FromGensOnly(const Matrix<mpz_class>& Generators,  size_t nr_special_gens,
            const Matrix<mpz_class>& SpecialLinForms);
#ifdef ENFNORMALIZ
template nauty_result compute_automs_by_nauty_Gens_LF(const Matrix<renf_elem_class>& Generators, 
                        size_t nr_special_gens, const Matrix<renf_elem_class>& LinForms,
                        const size_t nr_special_linforms,  bool zero_one);

template nauty_result compute_automs_by_nauty_FromGensOnly(const Matrix<renf_elem_class>& Generators,  size_t nr_special_gens,
            const Matrix<renf_elem_class>& SpecialLinForms);
#endif

} // namespace

#endif // NMZ_NAUTY

