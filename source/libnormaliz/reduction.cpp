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

#include "reduction.h"

namespace libnormaliz {
using namespace std;

//---------------------------------------------------------------------------

template<typename Integer>
Candidate<Integer>::Candidate(const vector<Integer>& v, const vector<Integer>& val, long sd){
    cand(v);
    values(val);
    sort_deg(sd);
    reducible(true);
    original_generator(false);
    // in_HB(false);   
}

//---------------------------------------------------------------------------

template<typename Integer>
Candidate<Integer>::Candidate(const vector<Integer>& v, const Full_Cone<Integer>& C){
    cand=v;
    values.resize(C.nrSupport_Hyperplanes);
    typename list<vector<Integer> >::const_iterator h;
    size_t i=0;
    for(h=C.Support_Hyperplanes.begin();h!=C.Support_Hyperplanes.end();++h){
        values[i]=v_scalar_product(v,*h);
        ++i;
    }
    sort_deg=explicit_cast_to_long<Integer>(v_scalar_product(v,C.Sorting));
    original_generator=false;
    // in_HB=false;
}

//---------------------------------------------------------------------------

template<typename Integer>
Candidate<Integer>::Candidate(const vector<Integer>& v, size_t max_size){
    cand=v;
    values.resize(max_size,0);
    sort_deg=0;
    reducible=true;
    original_generator=false;
    // in_HB=false;   
}

//---------------------------------------------------------------------------

template<typename Integer>
Candidate<Integer>::Candidate(size_t cand_size, size_t val_size){
    // cand=v;
    values.resize(val_size,0);
    cand.resize(cand_size,0);
    sort_deg=0;
    reducible=true;
    original_generator=false;
    // in_HB=false; 
}



//---------------------------------------------------------------------------

template<typename Integer>
CandidateList<Integer>::CandidateList(){

}


//---------------------------------------------------------------------------

template<typename Integer>
CandidateList<Integer>::CandidateList(bool dual_algorithm){

    dual=dual_algorithm;  
}

//---------------------------------------------------------------------------

/*template<typename Integer>
void CandidateList<Integer>::insert(const vector<Integer>& v, Full_Cone<Integer>& C){
    insert(v,C.Hyperplanes,C.Sorting);
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::insert(const vector<Integer>& v, const list<vector<Integer> >& SuppHyps, 
            const size_t& nrSuppHyps, const vector<Integer>& Sorting){

    typename list<vector<Integer> >::const_iterator h;
    Integer sd;
    vector<Integer> val(nrSuppHyps);
    size_t i=0;
    for(h=SuppHyps.begin();h!=SuppHyps.end();++h){
        val[i]=v_scalar_product(v,*h);
        ++i;
    }
    sd=explicit_cast_to_long<Integer>(v_scalar_product(*v,Sorting));
    Candidates.push_back(Candidate(v,val,sd));
} */

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::is_reducible(const vector<Integer>& values, const long sort_deg) const {
 
    long sd;
    /* if(dual)
        sd=sort_deg;
    else */
        sd=sort_deg/2;
    size_t kk=0;
    typename list<Candidate<Integer> >::const_iterator r;
    for(r=Candidates.begin();r!=Candidates.end();++r){
        if(sd < r->sort_deg){
            return(false);
        }
        size_t i=0;
        if(values[kk]<r->values[kk])
                continue;
        for(;i<values.size();++i)
            if(values[i]<r->values[i]){
                kk=i;
                break;
            }
        if(i==values.size()){
            return(true);
        }
   }   
   return(false);    
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::is_reducible_last_hyp(const vector<Integer>& values, const long sort_deg) const {
 
    long sd;
    /* if(dual)
        sd=sort_deg;
    else */
        sd=sort_deg/2;
    size_t kk=0;
    typename list<Candidate<Integer> >::const_iterator r;
    for(r=Candidates.begin();r!=Candidates.end();++r){
        if(sd < r->sort_deg){
            return(false);
        }
        size_t i=0;
        
        if(values[last_hyp]<r->values[last_hyp])
                continue;
        
        if(values[kk]<r->values[kk])
                continue;
        for(;i<values.size();++i)
            if(values[i]<r->values[i]){
                kk=i;
                break;
            }
        if(i==values.size()){
            return(true);
        }
   }   
   return(false);    
}


//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::is_reducible_last_hyp(Candidate<Integer>& c) const {

    /*if(dual && c.in_HB)
        c.reducible=false;
    else */
        c.reducible=is_reducible_last_hyp(c.values, c.sort_deg);
    return(c.reducible);
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::is_reducible(Candidate<Integer>& c) const {

    /*if(dual && c.in_HB)
        c.reducible=false;
    else */
        c.reducible=is_reducible(c.values, c.sort_deg);
    return(c.reducible);
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::is_reducible(vector<Integer> v,Candidate<Integer>& cand, const Full_Cone<Integer>& C) const {
    cand=Candidate<Integer>(v,C);
    return((*this).is_reducible(cand));
}

//---------------------------------------------------------------------------

// Fourth version with parallelization and tables
template<typename Integer>
void CandidateList<Integer>::reduce_by(CandidateList<Integer>& Reducers){

        typename list<Candidate<Integer> >::iterator c;
        size_t cpos,csize=Candidates.size();
        
        CandidateTable<Integer> ReducerTable(Reducers);
        
        #pragma omp parallel private(c,cpos) firstprivate(ReducerTable)
        {
        
        c=Candidates.begin();
        cpos=0;
        
        #pragma omp for schedule(dynamic)
        for (size_t k=0; k<csize; ++k) {
            for(;k > cpos; ++cpos, ++c) ;
            for(;k < cpos; --cpos, --c) ;
        
            ReducerTable.is_reducible(*c);
        }
        
        }// end parallel
        
        // erase reducibles
        for(c=Candidates.begin();c!=Candidates.end();){
            if((*c).reducible)
                c=Candidates.erase(c);
            else // continue
                ++c;
        }      
}

//---------------------------------------------------------------------------

/*template<typename Integer>
void CandidateList<Integer>::auto_reduce(){
cout << "Size " << Candidates.size() << endl;
    reduce_by(*this);
}*/

template<typename Integer>
void CandidateList<Integer>::auto_reduce(){

    if(empty())
        return;
        
    sort_by_deg();
    auto_reduce_sorted();
}

template<typename Integer>
void CandidateList<Integer>::auto_reduce_sorted(){
// uses generations defined by degrees

    if(empty())
        return;

    CandidateList<Integer> Irreducibles(dual), CurrentReducers(dual);
    long irred_degree;
    size_t cs=Candidates.size();
    if(verbose && cs > 1000){
            verboseOutput() << "auto-reduce " << cs << " candidates, degrees <= "; 
    }
    
    typename list<Candidate<Integer> >::iterator c;
    while(!Candidates.empty()){
        irred_degree=Candidates.begin()->sort_deg*2-1;
        if(verbose && cs > 1000){
            verboseOutput() << irred_degree << " " << flush;
        }
        for(c=Candidates.begin();c!=Candidates.end() && c->sort_deg <=irred_degree;++c); // find location for splicing
        CurrentReducers.Candidates.splice(CurrentReducers.Candidates.begin(),Candidates,Candidates.begin(),c);
        reduce_by(CurrentReducers);
        Irreducibles.Candidates.splice(Irreducibles.Candidates.end(),CurrentReducers.Candidates);
    }
    if(verbose && cs > 1000){
            verboseOutput() << endl;
    }
    Candidates.splice(Candidates.begin(),Irreducibles.Candidates);
}

/*
//---------------------------------------------------------------------------
template<typename Integer>
void CandidateList<Integer>::unique_auto_reduce(bool only_unique){

    unique_vectors();
    if(only_unique) // in this case we have only to make unique
            return;
    auto_reduce();        
}
*/

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::reduce_by_and_insert(Candidate<Integer>& cand, const CandidateList<Integer>& Reducers){
    bool irred=!Reducers.is_reducible(cand);
    if(irred)
        Candidates.push_back(cand);
    return irred;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::reduce_by_and_insert(const vector<Integer>& v, Full_Cone<Integer>& C, CandidateList<Integer>& Reducers){
    Candidate<Integer> cand(v,C);
    return reduce_by_and_insert(cand,Reducers);
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::unique_vectors(){

    assert(dual);

    if(empty())
        return;
        
    // sort_by_val();

    typename list<Candidate<Integer> >::iterator h,h_start,prev;
    h_start=Candidates.begin();

    h_start++;    
    for(h=h_start;h!=Candidates.end();){
        prev=h;
        prev--;
        if(h->values==prev->values)  // since cone may not be pointed in the dual , vectors
            h=Candidates.erase(h);   // must be made unique modulo the unit group
        else                         // values gives standard embedding
            ++h;
    }
}


//---------------------------------------------------------------------------
/*
template<typename Integer>
void CandidateList<Integer>::select_HB(size_t guaranteed_HB_deg){

    typename list<Candidate<Integer> >::iterator h;
    for(h=Candidates.begin(); h!=Candidates.end();++h)
        if(h->old_tot_deg<=guaranteed_HB_deg)
            h->in_HB=true;

}
*/

//---------------------------------------------------------------------------

template<typename Integer>
bool deg_compare(const Candidate<Integer>& a, const Candidate<Integer>& b){
    return(a.sort_deg < b.sort_deg);
}

//---------------------------------------------------------------------------

template<typename Integer>
bool val_compare(const Candidate<Integer>& a, const Candidate<Integer>& b){
    if(a.sort_deg<b.sort_deg)
        return(true);
    if(a.sort_deg==b.sort_deg)
        return(a.values < b.values);
    return false;
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::sort_by_deg(){

    Candidates.sort(deg_compare<Integer>);

}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::sort_by_val(){

    Candidates.sort(val_compare<Integer>);

}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::clear(){
    Candidates.clear();
}

//---------------------------------------------------------------------------

template<typename Integer>
size_t CandidateList<Integer>::size(){
    return Candidates.size();
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateList<Integer>::empty(){
    return Candidates.empty();
}


//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::merge(CandidateList<Integer>& NewCand){
    Candidates.merge(NewCand.Candidates,deg_compare<Integer>);
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::merge_by_val(CandidateList<Integer>& NewCand){
    Candidates.merge(NewCand.Candidates,val_compare<Integer>);
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::push_back(const Candidate<Integer>& cand){
    // cout << cand;
    Candidates.push_back(cand);
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::extract(list<vector<Integer> >& V_List){
    typename list<Candidate<Integer> >::iterator c;
    for(c=Candidates.begin();c!=Candidates.end();++c)
    V_List.push_back(c->cand);
                
}

//---------------------------------------------------------------------------

template<typename Integer>
void CandidateList<Integer>::splice(CandidateList<Integer>& NewCand){
    Candidates.splice(Candidates.begin(),NewCand.Candidates);
}

//---------------------------------------------------------------------------

template<typename Integer>
CandidateTable<Integer>::CandidateTable(CandidateList<Integer>& CandList){
    typename list<Candidate<Integer> >::iterator c;
    for(c=CandList.Candidates.begin();c!=CandList.Candidates.end();++c)
        CandidatePointers.push_back(&(*c));
        dual=CandList.dual;
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateTable<Integer>::is_reducible(Candidate<Integer>& c){
    c.reducible=is_reducible(c.values, c.sort_deg);
    return(c.reducible);
}

//---------------------------------------------------------------------------

template<typename Integer>
bool CandidateTable<Integer>::is_reducible(const vector<Integer>& values, const long sort_deg) {

    long sd;
    /* if(dual)
        sd=sort_deg;
    else */
        sd=sort_deg/2;
    size_t kk=0;
    typename list<Candidate<Integer>* >::iterator r;
    for(r=CandidatePointers.begin();r!=CandidatePointers.end();++r){
        if(sd < (*r)->sort_deg){
            return(false);
        }
        size_t i=0;
        if(values[kk]<(*r)->values[kk])
                continue;
        for(;i<values.size();++i)
            if(values[i]<(*r)->values[i]){
                kk=i;
                break;
            }
        if(i==values.size()){
            CandidatePointers.splice(CandidatePointers.begin(),CandidatePointers,r);
            return(true);
        }
   }   
   return(false);    
}

 
} // namespace
