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

    bool verbose;

    bool mult_of_ev_ok;

    Matrix<Integer> F;

    size_t fusion_rank;
    vector<Integer> fusion_type; // to be made from fusion_type_string
    string fusion_type_string;
    vector<key_t> duality;

    vector<Integer> ImageRing;

    Integer FPdim;
    Integer FPSquare;

    FusionBasic FusBasic;
    FusionComp<Integer> FusComp;
    vector<Matrix<Integer> >  Tables;

    vector<Integer> divisors;
    vector<size_t> MultEV;
    Matrix<Integer> EVMat;

    vector<Integer> low_m;  // the m_i for i < r

    vector<Matrix<Integer> > InductionMatrices;
    vector<Matrix<Integer> > LowParts;
    // vector<Matrix<Integer> > LowPartsBounds;

    size_t iupper_bound;


    Integer N(const key_t i, const key_t j, const key_t k);

    map< Integer, Matrix<Integer > > LowRepresentations;  // F_ij for i <= r (counting from 1), F_i1 = 1
    Matrix<Integer > HighRepresentations; // F_ij for i > 1 (counting from 1), F_i1 = 0
    Matrix<Integer> Bounds;

    Induction();
    Induction(const vector<Integer>& fus_type, const vector<key_t>& fus_duality , const vector<Integer>& FusRing, bool verb);

    void start_low_parts();
    void build_low_matrices(Matrix<Integer> matrix_so_far, Matrix<Integer> bounds_so_far);
    void from_low_to_full(const Matrix<Integer>& ThisLowPart);
    void augment_induction_matrices();

    Matrix<Integer> make_allowed_transpositions(Matrix<Integer> FusionMap);

    void compute();

    //void extend_matrix(Matrix<Integer> matrix_so_far, key_t rep_index, Matrix<Integer> bounds_so_far, Integer FPdim_so_far);

}; // class Induction end

} // namespace

#endif /* LIBNORMALIZ_INDUCTION_H */
