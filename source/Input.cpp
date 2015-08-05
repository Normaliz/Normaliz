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

#include "Options.h"
#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/map_operations.h"

template <typename Integer>
map <Type::InputType, vector< vector<Integer> > > readNormalizInput (istream& in) {

    string type_string;
    long i,j;
    long nr_rows,nr_columns;;
    InputType input_type = Type::integral_closure;
    Integer number;
    map<Type::InputType, vector< vector<Integer> > > input_map;
    typename map<Type::InputType, vector< vector<Integer> > >::iterator it;

    while (in.good()) {
        in >> nr_rows;
        if(in.fail())
            break;
        in >> nr_columns;
        if((nr_rows <0) || (nr_columns < 0)){
            cerr << "Error while reading a "<<nr_rows<<"x"<<nr_columns<<" matrix form the input!" << endl;
            throw BadInputException();
        }
        vector< vector<Integer> > M(nr_rows,vector<Integer>(nr_columns));
        for(i=0; i<nr_rows; i++){
            for(j=0; j<nr_columns; j++) {
                in >> number;
                M[i][j] = number;
            }
        }

        in>>type_string;

        if ( in.fail() ) {
            cerr << "Error while reading a "<<nr_rows<<"x"<<nr_columns<<" matrix form the input!" << endl;
            throw BadInputException();
        }

        input_type = to_type(type_string);

        //check if this type already exists
        if(exists_element(input_map, input_type)){
            cerr << "Error: Multiple inputs of type \"" << type_string <<"\" are not allowed!" << endl;
            throw BadInputException();
        }
        input_map[input_type] = M;
    }

    return input_map;
}
