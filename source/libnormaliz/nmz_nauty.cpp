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

#include "libnormaliz/nmz_nauty.h"
#include "libnormaliz/vector_operations.h"

#ifdef NMZ_NUTY

#define MAXN 5000    /* Define this before including nauty.h */

#include <nauty/nauty.h>

namespace libnormaliz {
using namespace std;

vector<vector<long> > CollectedAutoms;

void getmyautoms(int count, int *perm, int *orbits,
               int numorbits, int stabvertex, int n)
{
    int i;
    vector<long> this_perm(n);
    for (i = 0; i < n; ++i) this_perm[i] = perm[i];
    CollectedAutoms.push_back(this_perm);
}

template<typename Integer>
vector<vector<long> > compute_automs(const vector<vector<Integer> >& Generators, const vector<vector<Integer> >& Support_Hyperplanes){
  
    graph g[MAXN*MAXM];
    int lab[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    
    options.userautomproc = getmyautoms;
    
    int n,m,v;
    
    /* Default options are set by the DEFAULTOPTIONS_GRAPH macro above.
     *   Here we change those options that we want to be different from the
     *   defaults.  writeautoms=TRUE causes automorphisms to be written.     */
    
    options.writeautoms = TRUE;
    options.defaultptn = FALSE;
    
    size_t mm=Generators.size();
    size_t nn=Support_Hyperplanes.size();
    
    vector<vector<Integer> >  MM(mm,vector<Integer> (nn));
    int dummy,i,j,k;
    
    for(i=0;i<mm; ++i){
	for(j=0;j<nn;++j)
	  MM[i][j]=v_scalar_product(Generators[i],Support_Hyperplanes[j]);
      
    }
        
        boolean first=TRUE;
        int mini=0;
        int maxi=0;
        
        for(i=0;i<mm;++i){
            for(j=0;j<nn;++j){
                if(first || MM[i][j]< mini)
                    mini= MM[i][j];
                if(first || MM[i][j] > maxi)
                    maxi= MM[i][j];
                first=FALSE;
            }
        }
        
        printf("max %d\n",maxi);
        printf("min %d\n",mini);
        
        for(i=0;i<mm;++i)
            for(j=0;j<nn;++j)
                MM[i][j]-=mini;
            maxi-=mini;
        
        int test=1;
        int ll=0;
        for(i=0;;++i){
            if(test>maxi)
                break;
            ll++;
            test*=2;        
        }
        
        printf("log %d\n",ll);
        int bin_exp[ll];
        
        int layer_size=mm+nn;
        n=ll*layer_size;
        
        if (n > MAXN)
        {
            printf("n must be in the range 1..%d\n",MAXN);
            exit(1);
        }
        
        /* The nauty parameter m is a value such that an array of
         *       m setwords is sufficient to hold n bits.  The type setword
         *       is defined in nauty.h.  The number of bits in a setword is
         *       WORDSIZE, which is 16, 32 or 64.  Here we calculate
         *       m = ceiling(n/WORDSIZE).                                  */
        
        m = SETWORDSNEEDED(n);
        
        /* The following optional call verifies that we are linking
         *       to compatible versions of the nauty routines.            */
        
        nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
        
        /* Now we create the cycle.  First we zero the graph, than for
         *       each v, we add the edge (v,v+1), where values are mod n. */
        
        EMPTYGRAPH(g,m,n);
        // for (v = 0; v < n; ++v)  ADDONEEDGE(g,v,(v+1)%n,m);
        
        for(i=0;i<layer_size;++i){   // make vertical edges for all vertices
            for(k=1;k<ll;++k)
                ADDONEEDGE(g,(k-1)*layer_size+i,k*layer_size+i,m);
        }
        
        for(i=0;i<mm;++i){   // make horizontal edges layer by layer
            for(j=0;j<nn;++j){
                test=MM[i][j];
                for(k=0;k<ll;++k){ // binary expansion of matrix entry
                    bin_exp[k]=test%2;
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
        
        /* Since we are not requiring a canonical labelling, the last
         *       parameter to densenauty() is not required and can be NULL. */
        
        densenauty(g,lab,ptn,orbits,&options,&stats,m,n,NULL);
        
        /* The size of the group is returned in stats.grpsize1 and
         *       stats.grpsize2. */
        
        printf("Automorphism group size = ");
        writegroupsize(stdout,stats.grpsize1,stats.grpsize2);
        printf("\n");
        
        printf("\n===================\n");
        
        auto it=CollectedAutoms.begin();
        //for(;it!=CollectedAutoms.end();++it)
         //   cout << *it;
        // }
	
	return CollectedAutoms;
        
}

#ifndef NMZ_MIC_OFFLOAD  //offload with long is not supported
template vector<vector<long> > compute_automs(const vector<vector<long> >& Generators, const vector<vector<long> >& Support_Hyperplanes);
#endif // NMZ_MIC_OFFLOAD
template vector<vector<long> > compute_automs(const vector<vector<long long> >& Generators, const vector<vector<long long> >& Support_Hyperplanes);
// template vector<vector<long> > compute_automs(const vector<vector<mpz_class> >& Generators, const vector<vector<mpz_class> >& Support_Hyperplanes);

} // namespace

#endif // NMZ_NAUTY
