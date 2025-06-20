/*
 * Normaliz
 * Copyright (C) 2007-2025  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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

#ifndef LIBNORMALIZ_FAVE_LATTICE_H_
#define LIBNORMALIZ_FAVE_LATTICE_H_

#include <vector>
#include <set>
#include <list>
#include <map>

#include <libnormaliz/general.h>
#include <libnormaliz/matrix.h>
#include "libnormaliz/dynamic_bitset.h"

namespace libnormaliz {
using std::map;
using std::pair;
using std::vector;

template <typename Integer>
class FaceLattice {
    bool verbose;
    bool inhomogeneous;

    size_t nr_supphyps;
    size_t nr_extr_rec_cone;
    size_t nr_vert;
    size_t nr_gens;

    size_t dim;  // we assume pointed!

    Matrix<Integer> SuppHyps;  // local storage for support hypeplanes
    Matrix<Integer> VerticesOfPolyhedron;
    Matrix<Integer> ExtremeRaysRecCone;

    map<dynamic_bitset, int> FaceLat;
    vector<dynamic_bitset> SuppHypInd;
    vector<size_t> f_vector;

    vector<vector<key_t> > SuppHypPermutations;
    vector<vector<key_t> > ExtRayPermutations;

   public:
    FaceLattice(Matrix<Integer>& SupportHyperplanes,
                const Matrix<Integer>& Vert,
                const Matrix<Integer>& ExtRaysRC,
                const bool cone_inhomogeneous,
                bool swap_allowed = true);
    FaceLattice();
    void compute(const long face_codim_bound, const bool verbose, bool change_integer_type, const bool only_f_vector);
    void compute_orbits(const long face_codim_bound, const bool verbose, bool change_integer_type, const bool only_f_vector);
    vector<size_t> getFVector();
    void get(map<dynamic_bitset, int>& FaceLatticeOutput);
    void get(vector<dynamic_bitset>& SuppHypIndOutput);
    void set_supphyp_permutations(const vector<vector<key_t> >& gens, const bool verb);
    void set_extray_permutations(const vector<vector<key_t> >& gens, const bool verb);
    dynamic_bitset normal_form(const dynamic_bitset& arg, const vector<vector<key_t> >& Perms);
};

}  // namespace libnormaliz

#endif /* LIBNORMALIZ_FAVE_LATTICE_H__ */
