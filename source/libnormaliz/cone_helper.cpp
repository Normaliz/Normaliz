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

#ifdef NMZ_MIC_OFFLOAD
#pragma offload_attribute (push, target(mic))
#endif

#include "libnormaliz/cone_helper.h"
#include "libnormaliz/vector_operations.h"
#include "libnormaliz/my_omp.h"

namespace libnormaliz {
using std::vector;

//---------------------------------------------------------------------------

// determines the maximal subsets in a vector of subsets given by their indicator vectors
// result returned in is_max_subset -- must be initialized outside
// only set to false in this routine
// if a set occurs more than once, only the last instance is recognized as maximal
void maximal_subsets(const vector<vector<bool> >& ind, vector<bool>& is_max_subset) {

    if(ind.size()==0)
        return;

    size_t nr_sets=ind.size();
    size_t card=ind[0].size();
    vector<key_t> elem(card);

    for (size_t i = 0; i <nr_sets; i++) {
        if(!is_max_subset[i])  // already known to be non-maximal
            continue;

        size_t k=0; // counts the number of elements in set with index i
        for (size_t j = 0; j <card; j++) {
            if (ind[i][j]) {
                elem[k]=j;
                k++;
            }
        }

        for (size_t j = 0; j <nr_sets; j++) {
            if (i==j || !is_max_subset[j] ) // don't compare with itself or something known not to be maximal
                continue;
            size_t t;
            for (t = 0; t<k; t++) {
                if (!ind[j][elem[t]])
                    break; // not a superset
            }
            if (t==k) { // found a superset
                is_max_subset[i]=false;
                break; // the loop over j
            }
        }
    }
}

//---------------------------------------------------------------------------
// computes c1*v1-c2*v2
template<typename Integer>
vector<Integer> FM_comb(Integer c1, const vector<Integer>& v1,Integer c2, const vector<Integer>& v2, bool& is_zero){
 
    size_t dim=v1.size();
    vector<Integer> new_supp(dim);
    is_zero=false;
    size_t k=0;
    for(;k<dim;++k){
        new_supp[k]=c1*v1[k]-c2*v2[k];
        if(!check_range(new_supp[k]))
            break;    
    }
    Integer g=0;
    if(k==dim)
        g=v_make_prime(new_supp);
    else{ // redo in GMP if necessary
        #pragma omp atomic
        GMP_hyp++;
        vector<mpz_class> mpz_neg(dim), mpz_pos(dim), mpz_sum(dim);
        convert(mpz_neg, v1);
        convert(mpz_pos, v2);
        for (k = 0; k <dim; k++)
            mpz_sum[k]=convertTo<mpz_class>(c1)*mpz_neg[k]-
                    convertTo<mpz_class>(c2)*mpz_pos[k];
        mpz_class GG=v_make_prime(mpz_sum);
        convert(new_supp, mpz_sum);
        convert(g,GG);
    }
    if(g==0)
        is_zero=true;
    return new_supp;
}

template<typename IntegerPL>
vector<size_t> order_supps(const Matrix<IntegerPL>& Supps){
    
    assert(Supps.nr_of_rows()>0);
    size_t dim=Supps.nr_of_columns();

    vector<pair<IntegerPL,size_t> > NewPos,NewNeg; // to record the order of the support haperplanes
    for(size_t i=0;i<Supps.nr_of_rows();++i){
        if(Supps[i][dim-1] >= 0)
            NewPos.push_back(make_pair(-Supps[i][dim-1],i));
        else
            NewNeg.push_back(make_pair(Supps[i][dim-1],i));
    }
    sort(NewPos.begin(),NewPos.end());
    sort(NewNeg.begin(),NewNeg.end());
    
    size_t min_length=NewNeg.size();
    if(NewPos.size()<min_length)
        min_length=NewPos.size();
    
    vector<size_t> Order;
    
    for(size_t i=0;i<min_length;++i){
        Order.push_back(NewPos[i].second);
        Order.push_back(NewNeg[i].second);
    }
    for(size_t i=min_length;i<NewPos.size();++i)
        Order.push_back(NewPos[i].second);
    for(size_t i=min_length;i<NewNeg.size();++i)
        Order.push_back(NewNeg[i].second);
    
    assert(Order.size()==Supps.nr_of_rows());
    
    /* for(size_t i=0;i<Order.size();++i)
        cout << Supps[Order[i]][dim-1] << " ";
    cout << endl;*/
    
    return Order;
}

//---------------------------------------------------------------------------
template<typename IntegerPL>
void compute_projections(vector<Matrix<IntegerPL> >& AllSupps, vector<vector<size_t> >& AllOrders, size_t dim,
                         vector< boost::dynamic_bitset<> >& Ind, size_t rank, bool verbose){
    
    INTERRUPT_COMPUTATION_BY_EXCEPTION
    
    if(dim==1)
        return;
    
    const Matrix<IntegerPL> & Supps=AllSupps[dim];
    size_t dim1=dim-1;
    
    if(verbose)
        verboseOutput() << "embdim " << dim  << " inequalities " << Supps.nr_of_rows() << endl; 
    
    // cout << "SSS" << Ind.size() << " " << Ind;

    // We now augment the given cone by the last basis vector and its negative
    // Afterwards we project modulo the subspace spanned by them
    
    vector<key_t> Neg, Pos; // for the Fourier-Motzkin elimination of inequalities    
    Matrix<IntegerPL> SuppsProj(0,dim); // for the support hyperplanes of the projection
    Matrix<IntegerPL> EqusProj(0,dim); // for the equations (both later minimized)
    
    // First we make incidence vectors with the given generators    
    vector< boost::dynamic_bitset<> > NewInd; // for the incidence vectors of the new hyperplanes
    
     boost::dynamic_bitset<> TRUE(Ind[0].size());
     TRUE.set();
     
     vector<bool> IsEquation(Supps.nr_of_rows());
     
    bool rank_goes_up=false; // if we add the last unit vector
    size_t PosEquAt=0; // we memorize the positions of pos/neg equations if rank goes up
    size_t NegEquAt=0;

    for(size_t i=0;i<Supps.nr_of_rows();++i){
        if(Ind[i]==TRUE)
            IsEquation[i]=true;
        
        if(Supps[i][dim1]==0){  // already independent of last coordinate
            if(IsEquation[i])
                EqusProj.append(Supps[i]); // is equation
            else{
                SuppsProj.append(Supps[i]); // neutral support hyperplane
                NewInd.push_back(Ind[i]);
            }
            continue;
        }
        if(IsEquation[i])
            rank_goes_up=true;
        if(Supps[i][dim1]>0){
            if(IsEquation[i])
                PosEquAt=i;
            Pos.push_back(i);
            continue;
        }
        Neg.push_back(i);
        if(IsEquation[i])
            NegEquAt=i;
    }
    
    // cout << "Nach Pos/Neg " << EqusProj.nr_of_rows() << " " << Pos.size() << " " << Neg.size() << endl;
    
    // now the elimination, matching Pos and Neg
    
    // cout << "rank_goes_up " << rank_goes_up << endl;
    
        bool skip_remaining;
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif
    
    if(rank_goes_up){
        for(size_t i=0;i<Pos.size();++i){ // match pos and neg equations
            size_t p=Pos[i];
            if(!IsEquation[p])
                continue;
            IntegerPL PosVal=Supps[p][dim1];
            for(size_t j=0;j<Neg.size();++j){
                size_t n=Neg[j];
                if(!IsEquation[n])
                    continue;
                IntegerPL NegVal=Supps[n][dim1];
                bool is_zero;
                // cout << Supps[p];
                // cout << Supps[n];
                vector<IntegerPL> new_equ=FM_comb(PosVal,Supps[n],NegVal,Supps[p],is_zero);
                // cout << "zero " << is_zero << endl;
                // cout << "=====================" << endl;
                if(is_zero)
                    continue;
                EqusProj.append(new_equ);
            }
        }
        
        for(size_t i=0;i<Pos.size();++i){ // match pos inequalities with a negative equation
            size_t p=Pos[i];
            if(IsEquation[p])
                continue;
            IntegerPL PosVal=Supps[p][dim1];
            IntegerPL NegVal=Supps[NegEquAt][dim1];
            vector<IntegerPL> new_supp(dim);
            bool is_zero;
            new_supp=FM_comb(PosVal,Supps[NegEquAt],NegVal,Supps[p],is_zero);
            /* cout << Supps[NegEquAt];
            cout << Supps[p];
            cout << new_supp;
            cout << "zero " << is_zero << endl;
            cout << "+++++++++++++++++++++" << endl; */
            if(is_zero) // cannot happen, but included for analogy
                continue;
            SuppsProj.append(new_supp);
            NewInd.push_back(Ind[p]);            
        }
        
        for(size_t j=0;j<Neg.size();++j){ // match neg inequalities with a posizive equation
            size_t n=Neg[j];
            if(IsEquation[n])
                continue;
            IntegerPL PosVal=Supps[PosEquAt][dim1];
            IntegerPL NegVal=Supps[n][dim1];
            vector<IntegerPL> new_supp(dim);
            bool is_zero;
            new_supp=FM_comb(PosVal,Supps[n],NegVal,Supps[PosEquAt],is_zero);
            /* cout << Supps[PosEquAt];
            cout << Supps[n];
            cout << new_supp;
            cout << "zero " << is_zero << endl;
            cout << "=====================" << endl;*/
            
            if(is_zero) // cannot happen, but included for analogy
                continue;
            SuppsProj.append(new_supp);
            NewInd.push_back(Ind[n]);            
        }
    }
    
    // cout << "Nach RGU " << EqusProj.nr_of_rows() << " " << SuppsProj.nr_of_rows() << endl;
        
    if(!rank_goes_up){ // must match pos and neg hyperplanes
        
        skip_remaining=false;
        
        size_t min_nr_vertices=rank-2;
        /*if(rank>=3){
            min_nr_vertices=1;
            for(long i=0;i<(long) rank -3;++i)
                min_nr_vertices*=2;
                
        }*/
        
        #pragma omp parallel for schedule(dynamic)
        for(size_t i=0;i<Pos.size();++i){
            
            if (skip_remaining) continue;
        
#ifndef NCATCH
        try {
#endif
 
            INTERRUPT_COMPUTATION_BY_EXCEPTION
            
            size_t p=Pos[i];
            IntegerPL PosVal=Supps[p][dim1];
            vector<key_t> PosKey;
            for(size_t k=0;k<Ind[i].size();++k)
                if(Ind[p][k])
                    PosKey.push_back(k);
            
            for(size_t j=0;j<Neg.size();++j){
                size_t n=Neg[j];
                // // to give a facet of the extended cone
                // match incidence vectors
                boost::dynamic_bitset<> incidence(TRUE.size());
                size_t nr_match=0;
                vector<key_t> CommonKey;
                for(size_t k=0;k<PosKey.size();++k)
                    if(Ind[n][PosKey[k]]){
                        incidence[PosKey[k]]=true;
                        CommonKey.push_back(PosKey[k]);
                        nr_match++;
                    }
                if(rank>=2 && nr_match<min_nr_vertices) // cannot make subfacet of augmented cone
                    continue; 

                bool IsSubfacet=true;
                for(size_t k=0;k<Supps.nr_of_rows();++k){
                    if(k==p || k==n || IsEquation[k])
                        continue;
                    bool contained=true;
                    for(size_t j=0;j<CommonKey.size();++j){
                        if(!Ind[k][CommonKey[j]]){
                            contained=false;
                            break;                           
                        }                        
                    }
                    if(contained){
                        IsSubfacet=false;
                        break;
                    }
                }
                if(!IsSubfacet)
                    continue;
                //}
                
                IntegerPL NegVal=Supps[n][dim1];
                vector<IntegerPL> new_supp(dim);
                bool is_zero;
                new_supp=FM_comb(PosVal,Supps[n],NegVal,Supps[p],is_zero);                
                if(is_zero) // linear combination is 0
                    continue;
                
                if(nr_match==TRUE.size()){ // gives an equation
                    #pragma omp critical(NEWEQ)
                    EqusProj.append(new_supp);
                    continue;
                }
                #pragma omp critical(NEWSUPP)
                {
                SuppsProj.append(new_supp);
                NewInd.push_back(incidence);
                }
            }
            
#ifndef NCATCH
            } catch(const std::exception& ) {
                tmp_exception = std::current_exception();
                skip_remaining = true;
                #pragma omp flush(skip_remaining)
            }
#endif
        }
        
#ifndef NCATCH
        if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
        
    } // !rank_goes_up
    
    // cout << "Nach FM " << EqusProj.nr_of_rows() << " " << SuppsProj.nr_of_rows() << endl;
    
    Ind.clear(); // no longer needed
    
    EqusProj.resize_columns(dim1); // cut off the trailing 0     
    SuppsProj.resize_columns(dim1); // project hyperplanes
 
    // Equations have not yet been appended to support hypwerplanes
    EqusProj.row_echelon(); // reduce equations
    // cout << "Nach eche " << EqusProj.nr_of_rows() << endl;
    /* for(size_t i=0;i<EqusProj.nr_of_rows(); ++i)
        cout << EqusProj[i]; */
    SuppsProj.append(EqusProj); // append them as pairs of inequalities
    EqusProj.scalar_multiplication(-1);
    SuppsProj.append(EqusProj);
    
    // Now we must make the new indicator matrix    
    // We must add indictor vectors for the equations
    for(size_t i=0;i<2*EqusProj.nr_of_rows();++i)
        NewInd.push_back(TRUE);    

    AllOrders[dim1]=order_supps(SuppsProj);
    swap(AllSupps[dim1],SuppsProj);
    
    size_t new_rank=dim1-EqusProj.nr_of_rows();
    
    compute_projections(AllSupps,AllOrders,dim-1,NewInd, new_rank,verbose);
}

//---------------------------------------------------------------------------
template<typename IntegerPL,typename IntegerRet>
bool fiber_interval(IntegerRet& MinInterval, IntegerRet& MaxInterval, const vector<IntegerRet>& base_point, 
                    const Matrix<IntegerPL>& Supps, const vector<size_t> Order){

        bool FirstMin=true, FirstMax=true;
        vector<IntegerPL> LiftedGen;
        convert(LiftedGen,base_point);
        // cout << LiftedGen;
        size_t check_supps=Supps.nr_of_rows();
        if(check_supps>1000)
            check_supps=1000;
        for(size_t j=0;j<check_supps;++j){
            IntegerPL Den=Supps[Order[j]].back();
            if(Den==0)
                continue;
            IntegerPL Num= -v_scalar_product_vectors_unequal_lungth(LiftedGen,Supps[Order[j]]);
            // cout << "Num " << Num << endl;
            IntegerRet Quot;
            bool frac=int_quotient(Quot,Num,Den);
            IntegerRet Bound=0;
            //frac=(Num % Den !=0);
            if(Den>0){ // we must produce a lower bound of the interval
                if(Num>=0){  // true quot >= 0
                    Bound=Quot;
                    if(frac)
                        Bound++;
                }
                else // true quot < 0
                    Bound=-Quot;
                if(FirstMin || Bound > MinInterval){
                    MinInterval=Bound;
                    FirstMin=false;
                }
            }
            if(Den<0){ // we must produce an upper bound of the interval
                if(Num >= 0){ // true quot <= 0
                    Bound=-Quot;
                    if(frac)
                        Bound--;                    
                }
                else // true quot > 0
                    Bound=Quot;
                if(FirstMax || Bound < MaxInterval){
                    MaxInterval=Bound;
                    FirstMax=false;
                }
            }
            if(!FirstMax && !FirstMin && MaxInterval<MinInterval)
                return false; // interval empty
    }
    return true; // interval nonempty
}

///---------------------------------------------------------------------------
template<typename IntegerPL,typename IntegerRet>
void lift_points_to_this_dim(Matrix<IntegerRet>& Deg1, const Matrix<IntegerRet>& Deg1Proj,
                             const vector<Matrix<IntegerPL> >& AllSupps, const vector<vector<size_t> >& AllOrders){

    size_t dim=Deg1Proj.nr_of_columns()+1;
    size_t dim1=dim-1;
    vector<Matrix<IntegerRet> > Deg1Thread(omp_get_max_threads());
    for(size_t i=0;i<Deg1Thread.size();++i)
        Deg1Thread[i].resize(0,dim);
    
    bool skip_remaining;
#ifndef NCATCH
    std::exception_ptr tmp_exception;
#endif
    
    skip_remaining=false;
    
    #pragma omp parallel
    {
    int tn;
    if(omp_get_level()==0)
        tn=0;
    else    
        tn = omp_get_ancestor_thread_num(1);
 
    #pragma omp for schedule(dynamic)
    for(size_t i=0;i<Deg1Proj.nr_of_rows();++i){
        
        if (skip_remaining) continue;
        
#ifndef NCATCH
        try {
#endif

        IntegerRet MinInterval=0, MaxInterval=0; // the fiber over Deg1Proj[i] is an interval -- 0 to make gcc happy        
        fiber_interval(MinInterval, MaxInterval, Deg1Proj[i], AllSupps[dim],AllOrders[dim]);
        // cout << "Min " << MinInterval << " Max " << MaxInterval << endl;
        for(IntegerRet k=MinInterval;k<=MaxInterval;++k){
            
            INTERRUPT_COMPUTATION_BY_EXCEPTION
            
            vector<IntegerRet> NewPoint(dim);
            for(size_t j=0;j<dim1;++j)
                NewPoint[j]=Deg1Proj[i][j];
            NewPoint[dim1]=k;
            Deg1Thread[tn].append(NewPoint);
        }
        
#ifndef NCATCH
        } catch(const std::exception& ) {
            tmp_exception = std::current_exception();
            skip_remaining = true;
            #pragma omp flush(skip_remaining)
        }
#endif

    } // lifting
    } // pararllel
    
#ifndef NCATCH
    if (!(tmp_exception == 0)) std::rethrow_exception(tmp_exception);
#endif
    
    for(size_t i=0;i<Deg1Thread.size();++i)
        Deg1.append(Deg1Thread[i]);
    

    /* Deg1.pretty_print(cout);
    cout << "*******************" << endl; */
}

///---------------------------------------------------------------------------
template<typename IntegerPL,typename IntegerRet>
void lift_points_by_generation(Matrix<IntegerRet>& Deg1, const vector<Matrix<IntegerPL> >& AllSupps,
                                const vector<vector<size_t> >& AllOrders, const IntegerRet GD, bool verbose){
    if(verbose)
        verboseOutput() << "Lifting" << endl;
    
    size_t dim=AllSupps.size()-1;
    assert(dim>=2);
    
    Matrix<IntegerRet> Deg1Proj(0,1);
    vector<IntegerRet> One(1,GD);
        Deg1Proj.append(One);
    
    for(size_t i=2; i<=dim;++i){
        Deg1=Matrix<IntegerRet>(0,i);
        lift_points_to_this_dim(Deg1,Deg1Proj,AllSupps,AllOrders);
        if(verbose)    
            verboseOutput() <<  "embdim " << i << " Deg1Elements " << Deg1.nr_of_rows() << endl;
        if(i<dim){
            swap(Deg1,Deg1Proj);
        }               
    }
    if(verbose)
        verboseOutput() << "Lifting done" << endl;
}

///---------------------------------------------------------------------------
template<typename IntegerPL,typename IntegerRet>
void lift_point_recursively(vector<IntegerRet>& finaL_latt_point, const vector<IntegerRet>& latt_point_proj, 
                            const vector<Matrix<IntegerPL> >& AllSupps, const vector<vector<size_t> >& AllOrders, 
                            const vector<IntegerRet>& excluded_point){
 
    size_t dim1=latt_point_proj.size();
    size_t dim=dim1+1;
    size_t final_dim=AllSupps.size()-1;
    
    IntegerRet MinInterval=0, MaxInterval=0; // the fiber over Deg1Proj[i] is an interval -- 0 to make gcc happy        
    fiber_interval(MinInterval, MaxInterval, latt_point_proj, AllSupps[dim],AllOrders[dim]);
    for(IntegerRet k=MinInterval;k<=MaxInterval;++k){
        
        INTERRUPT_COMPUTATION_BY_EXCEPTION
        
        vector<IntegerRet> NewPoint(dim);
        for(size_t j=0;j<dim1;++j)
            NewPoint[j]=latt_point_proj[j];
        NewPoint[dim1]=k;
        if(dim==final_dim && NewPoint!=excluded_point){
            finaL_latt_point=NewPoint;
            break;
        }
        if(dim<final_dim){
            lift_point_recursively(finaL_latt_point, NewPoint,AllSupps,AllOrders,excluded_point);
            if(finaL_latt_point.size()>0)
                break;
        }
    }
}

///---------------------------------------------------------------------------
template<typename IntegerPL,typename IntegerRet>
void find_single_point(Matrix<IntegerRet>& Deg1, const vector<Matrix<IntegerPL> >& AllSupps,
                        const vector<vector<size_t> >& AllOrders, 
                        const IntegerRet& GD, const vector<IntegerRet>& excluded_point){
    
    size_t dim=AllSupps.size()-1;
    assert(dim>=2);
    
    vector<IntegerRet> start(1,GD);
    vector<IntegerRet> final_latt_point;
    lift_point_recursively(final_latt_point,start,AllSupps,AllOrders,excluded_point);
    if(final_latt_point.size()>0)
        Deg1.append(final_latt_point);
}

//---------------------------------------------------------------------------
template<typename IntegerPL,typename IntegerRet>
void project_and_lift_inner(Matrix<IntegerRet>& Deg1, const Matrix<IntegerPL>& Supps,
                            vector< boost::dynamic_bitset<> >& Ind, const IntegerRet& GD, size_t rank, 
                            bool verbose, bool all_points, const vector<IntegerRet>& excluded_point){

// Project-and-lift for lattice points in a polytope. 
// The first coordinate is homogenizing. Its value for polytope points ism set by GD so that
// a grading denominator 1=1 can be accomodated.
// We need only the support hyperplanes Supps and the facet-vertex incidence matrix Ind.
// Its rows correspond to facets.

    size_t dim=Supps.nr_of_columns(); // our local embedding dimension    
    vector<Matrix<IntegerPL> > AllSupps(dim+1);
    vector<vector<size_t> > AllOrders(dim+1);

    AllSupps[dim]=Supps;
    AllOrders[dim]=order_supps(Supps);
    compute_projections(AllSupps, AllOrders,dim, Ind, rank, verbose);
    if(all_points)
        lift_points_by_generation(Deg1,AllSupps,AllOrders,GD,verbose);
    else
        find_single_point(Deg1,AllSupps,AllOrders,GD,excluded_point);
}

//---------------------------------------------------------------------------

#ifdef NMZ_MIC_OFFLOAD
#pragma offload_attribute (pop)
#endif

} //end namespace libnormaliz
