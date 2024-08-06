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

#ifndef LIBNORMALIZ_INDUCTION_H_
#define LIBNORMALIZ_INDUCTION_H_

#include <vector>
#include <list>

#include "libnormaliz/general.h"
#include "libnormaliz/matrix.h"
#include "libnormaliz/dynamic_bitset.h"
#include "libnormaliz/nmz_polynomial.h"
#include "libnormaliz/fusion.h"

namespace libnormaliz {
using std::vector;

template <typename Integer>
class Induction {

public:

    Matrix<Integer> F;

    size_t fusion_rank;
    vector<Integer> fusion_type; // to be made from fusion_type_string
    string fusion_type_string;
    vector<key_t> duality;

    vector<Integer> FusRing;

    Integer FPdim;

    FusionBasic FusBasic;
    FusionComp<Integer> FusComp;
    vector<Matrix<Integer> >  Tables;

    vector<Integer> divisors;
    Matrix<Integer> EVMat;

    Integer N(const key_t i, const key_t j, const key_t k);

    Induction();
    Induction(const vector<Integer>& fus_type, const vector<key_t>& fus_duality , const vector<Integer>& FusRing);



}; // class Induction end

} // namespace

#endif /* LIBNORMALIZ_INDUCTION_H */
