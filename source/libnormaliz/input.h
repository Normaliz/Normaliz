/*
 * Normaliz
 * Copyright (C) 2007-2021  W. Bruns, B. Ichim, Ch. Soeger, U. v. d. Ohe
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

#include <iostream>
#include <cctype>  // std::isdigit
#include <limits>  // numeric_limits

#include "libnormaliz/options.h"
#include "libnormaliz/input_type.h"
#include "libnormaliz/list_and_map_operations.h"
#include "libnormaliz/cone_property.h"

#ifndef NORMALIZ_INPUT_H
#define NORMALIZ_INPUT_H


namespace libnormaliz {
    
template <typename Number>
map<Type::InputType, vector<vector<Number> > > readNormalizInput(istream& in,
                                                                 OptionsHandler& options,
                                                                 map<NumParam::Param, long>& num_param_input,
                                                                 string& polynomial,
                                                                 renf_class_shared& number_field);
} // namespace

#endif
