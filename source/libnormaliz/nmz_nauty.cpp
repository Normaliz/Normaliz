/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
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

#include "libnormaliz/integer.h"
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
vector<vector<long> > compute_automs_by_nauty(const vector<vector<Integer> >& Generators, const vector<vector<Integer> >& Support_Hyperplanes){
  
    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    
    options.userautomproc = getmyautoms;
    
    int n,m;
    
    /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
     *   Here we change those options that we want to be different from the
     *   defaults.  writeautoms=TRUE causes automorphisms to be written.     */
    
    options.writeautoms = FALSE;
    options.defaultptn = FALSE;
    
    size_t mm=Generators.size();
    size_t nn=Support_Hyperplanes.size();
    
    vector<vector<Integer> >  MM(mm,vector<Integer> (nn));
    size_t i,j,k;
    
    for(i=0;i<mm; ++i){
        for(j=0;j<nn;++j)
            MM[i][j]=v_scalar_product(Generators[i],Support_Hyperplanes[j]);      
    }
        
    bool first=true;
    Integer mini=0;
    Integer maxi=0;
    
    for(i=0;i<mm;++i){
        for(j=0;j<nn;++j){
            if(first || MM[i][j]< mini)
                mini= MM[i][j];
            if(first || MM[i][j] > maxi)
                maxi= MM[i][j];
            first=FALSE;
        }
    }
    
    // printf("max %d\n",maxi);
    // printf("min %d\n",mini);
    cout << "max " << maxi << " min " << mini << endl;
    
    for(i=0;i<mm;++i)
        for(j=0;j<nn;++j)
            MM[i][j]-=mini;
    maxi-=mini;
    
    Integer test=1;
    long ll=0;
    for(i=0;;++i){
        if(test>maxi)
            break;
        ll++;
        test*=2;        
    }
    
    cout << "log " << ll << endl;
    vector<long> bin_exp(ll);
    
    size_t layer_size=mm+nn;
    n=ll*layer_size;
    m = SETWORDSNEEDED(n);
    
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    
    DYNALLOC2(graph,g,g_sz,m,n,"malloc");
    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    
    EMPTYGRAPH(g,m,n);
    
    for(i=0;i<layer_size;++i){   // make vertical edges for all vertices
        for(k=1;k<ll;++k)
            ADDONEEDGE(g,(k-1)*layer_size+i,k*layer_size+i,m);
    }
    
    for(i=0;i<mm;++i){   // make horizontal edges layer by layer
        for(j=0;j<nn;++j){
            test=MM[i][j];
            for(k=0;k<ll;++k){ // binary expansion of matrix entry
                bin_exp[k]=convertTo<long>(test%2);
                test/=2;
            }
            for(k=0;k<ll;++k){
                if(bin_exp[k]==1)  // k is the number of layers below the current one
                    ADDONEEDGE(g,k*layer_size+i,k*layer_size+mm+j,m);
            }
        }
    }           
    
    for(i=0;i<n;++i){ // prepare partitions
        lab[i]=i;
        ptn[i]=1;
    }
    
    for(k=0;k<ll;++k){ // make partitions layer by layer
        ptn[k*layer_size+ mm-1]=0; // tow vertices in one partition
        ptn[(k+1)*layer_size-1]=0; // column indices in the next
    } 
    
    printf("Generators for Aut(C[%d]):\n",n);

    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
    
    printf("Automorphism group size = ");
    writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
    printf("\n");
    
    printf("\n===================\n");
    
    vector<vector<long> > AutomsAndOrbits(2*CollectedAutoms.size());
    AutomsAndOrbits.reserve(2*CollectedAutoms.size()+2);

    for(k=0;k<CollectedAutoms.size();++k){
        vector<long> GenPerm(mm);
        for(i=0;i<mm;++i)
            GenPerm[i]=CollectedAutoms[k][i];
        AutomsAndOrbits[k]=GenPerm;
        vector<long> LFPerm(nn);
        for(i=mm;i<mm+nn;++i)
            LFPerm[i-mm]=CollectedAutoms[k][i]-mm;
        AutomsAndOrbits[k+CollectedAutoms.size()]=LFPerm;        
    }
    
    vector<long> GenOrbits(mm);
    for(i=0;i<mm;++i)
        GenOrbits[i]=orbits[i];
    AutomsAndOrbits.push_back(GenOrbits);
    
    vector<long> LFOrbits(nn);
    for(i=0;i<nn;++i)
        LFOrbits[i]=orbits[i+mm]-mm;
    AutomsAndOrbits.push_back(LFOrbits);
    
    
    for(k=0;k<AutomsAndOrbits.size();++k){
        cout << AutomsAndOrbits[k] << endl;
    }
	
	nauty_freedyn();
    
	return AutomsAndOrbits;
        
}

#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template vector<vector<long> > compute_automs_by_nauty(const vector<vector<long> >& Generators, const vector<vector<long> >& Support_Hyperplanes);
#endif // NMZ_MIC_OFFLOAD
template vector<vector<long> > compute_automs_by_nauty(const vector<vector<long long> >& Generators, const vector<vector<long long> >& Support_Hyperplanes);
template vector<vector<long> > compute_automs_by_nauty(const vector<vector<mpz_class> >& Generators, const vector<vector<mpz_class> >& Support_Hyperplanes);

} // namespace

#endif // NMZ_NAUTY
