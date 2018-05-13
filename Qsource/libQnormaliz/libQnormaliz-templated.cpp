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

#include "libQnormaliz/Qlibnormaliz.cpp"
#include "libQnormaliz/Qinteger.cpp"
#include "libQnormaliz/Qvector_operations.cpp"
#include "libQnormaliz/Qmatrix.cpp"
// #include "libQnormaliz/Qsimplex.cpp"
#include "libQnormaliz/Qlist_operations.cpp"
#include "libQnormaliz/Qsublattice_representation.cpp"
//#include "libQnormaliz/Qreduction.cpp"
#include "libQnormaliz/Qproject_and_lift.cpp"
#include "libQnormaliz/Qfull_cone.cpp"
// #include "libQnormaliz/Qcone_dual_mode.cpp"
#include "libQnormaliz/Qcone.cpp"

namespace libQnormaliz {

template class Cone<mpq_class>;
template class Matrix<mpq_class>;
template class Sublattice_Representation<mpq_class>;
template class Full_Cone<mpq_class>;

#ifdef ENFNORMALIZ
template class Cone<renf_elem_class>;
template class Matrix<renf_elem_class>;
template class Sublattice_Representation<renf_elem_class>;
template class Full_Cone<renf_elem_class>;
#endif

}

