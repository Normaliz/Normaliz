/*
 * Normaliz
 * Copyright (C) 2007-2022  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

#include <fstream>
#include <sstream>

#include "libnormaliz/cone.h"
#include "libnormaliz/induction.h"

namespace libnormaliz{
using std::vector;
using std::string;
using std::list;
using std::ifstream;
using std::ofstream;

template<typename Integer>
Induction<Integer>::Induction(){

}

template<typename Integer>
Induction<Integer>::Induction(const vector<Integer>& fus_type, const vector<key_t>& fus_duality , const vector<Integer>& FusRing){

    FusBasic = FusionBasic();
    FusBasic.fusion_rank = fus_type.size();
    fusion_rank = FusBasic.fusion_rank;
    /* for(auto& f:fus_type){
        fusion_type_string += to_string(f) + " ";
    }
    FusBasic.fusion_type_string = fusion_type_string; */
    duality = fus_duality;
    FusBasic.duality = duality;
    FusBasic.fusion_type = fusion_coincidence_pattern(fus_type);
    FusComp = FusionComp<Integer>(FusBasic);
    FusComp.make_CoordMap();
    Tables = FusComp.make_all_data_tables(FusRing);
    cout << FusRing;
    for(auto& T: Tables)
        T.debug_print();

    FPdim = 0;
    for(auto& f: fus_type)
        FPdim += f*f;

    for(Integer t = 1; t <= FPdim; ++t){

        if(FPdim % t == 0)
            divisors.push_back(t);
    }

    cout << FPdim << " " << divisors;

    EVMat.resize(fusion_rank, fusion_rank);
    for(size_t s =0; s < fusion_rank; ++s){
        for(size_t l = 0; l < fusion_rank; ++l){
            Integer S = 0;
            for(size_t t = 0; t < fusion_rank; ++t){
                for(size_t k = 0; k < fusion_rank; ++k)
                    S += N(t, duality[t], k)*N(k,l,s);
            }
            EVMat[s][l] = S;
        }
    }
    EVMat.debug_print();

    cout << "-----------------" << endl;

    for(size_t i = 0; i< divisors.size(); ++i){
        cout << divisors[i] << " " << EVMat.mult_of_eigenvalue(divisors[i]) << endl;
    }
}

template<>
Induction<renf_elem_class>::Induction(const vector<renf_elem_class>& fus_type, const vector<key_t>& fus_duality , const vector<renf_elem_class>& FusRing){

    assert(false);
}

template<typename Integer>
Integer Induction<Integer>::N(const key_t i, const key_t j, const key_t k){
    return Tables[i][j][k];
}

template class Induction<mpz_class>;
template class Induction<long long>;
template class Induction<long>;
#ifdef ENFNORMALIZ
template class Induction<renf_elem_class>;
#endif



} // namespace
