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

#include <iostream>
#include <cctype>       // std::isdigit
#include <limits>       // numeric_limits

#include "Options.h"
#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/map_operations.h"
#include "libnormaliz/cone_property.h"

// eats up a comment, stream must start with "/*", eats everything until "*/"
void skip_comment(istream& in) {
    int i = in.get();
    int j = in.get();
    if (i != '/' || j != '*') {
        cerr << "Error: Bad comment start!" << endl;
        throw BadInputException();
    }
    while (in.good()) {
        in.ignore(numeric_limits<streamsize>::max(), '*'); //ignore everything until next '*'
        i = in.get();
        if (in.good() && i == '/') return; // successfully skipped comment
    }
    cerr << "Error: Incomplete comment!" << endl;
    throw BadInputException();
}

template<typename Integer>
void save_matrix(map<Type::InputType, vector<vector<Integer> > >& input_map,
        InputType input_type, const string& type_string, const vector<vector<Integer> >& M) {
    //check if this type already exists
    if (exists_element(input_map, input_type)) {
        cerr << "Error: Multiple inputs of type \"" << type_string
                << "\" are not allowed!" << endl;
        throw BadInputException();
    }
    input_map[input_type] = M;
}

template <typename Integer>
map <Type::InputType, vector< vector<Integer> > > readNormalizInput (istream& in, OptionsHandler& options) {

    string type_string;
    long i,j;
    long nr_rows,nr_columns;
    InputType input_type;
    Integer number;
    ConeProperty::Enum cp;

    map<Type::InputType, vector< vector<Integer> > > input_map;
    typename map<Type::InputType, vector< vector<Integer> > >::iterator it;

    in >> std::ws;  // eat up any leading white spaces
    int c = in.peek();
    if ( c == EOF ) {
        cerr << "Error: Empty input file!" << endl;
        throw BadInputException();
    }
    bool new_input_syntax = !std::isdigit(c);

    if (new_input_syntax) {
        long dim;
        while (in.peek() == '/') {
            skip_comment(in);
            in >> std::ws;
        }
        in >> type_string;
        if (!in.good() || type_string != "amb_space") {
            cerr << "Error: First entry must be \"amb_space\"!" << endl;
            throw BadInputException();
        }
        in >> dim;
        if (!in.good() || dim <= 0) {
            cerr << "Error: Bad amb_space value!" << endl;
            throw BadInputException();
        }
        while (in.good()) {
            in >> std::ws;  // eat up any leading white spaces
            c = in.peek();
            if (c == EOF) break;
            if (c == '/') {
                skip_comment(in);
            } else {
                in >> type_string;
                if (in.fail()) {
                    cerr << "Error: Could not read type string!" << endl;
                    throw BadInputException();
                }
                if (std::isdigit(c)) {
                    cerr << "Error: Unexpected number "<< type_string << " when expecting a type !" << endl;
                    throw BadInputException();
                }
                if (isConeProperty(cp, type_string)) {
                    options.activateInputFileConeProperty(cp);
                    continue;
                }
                if (type_string == "BigInt") {
                    options.activateInputFileBigInt();
                    continue;
                }
                if (type_string == "total_degree") {
                    input_type = Type::grading;
                    save_matrix(input_map, input_type, type_string, vector< vector<Integer> >(1,vector<Integer>(dim+type_nr_columns_correction(input_type),1)));
                    continue;
                }
                if (type_string == "nonnegative") {
                    input_type = Type::signs;
                    save_matrix(input_map, input_type, type_string, vector< vector<Integer> >(1,vector<Integer>(dim+type_nr_columns_correction(input_type),1)));
                    continue;
                }


                input_type = to_type(type_string);

                if (type_is_vector(input_type)) {
                    nr_rows = 1;
                    in >> std::ws;  // eat up any leading white spaces
                    c = in.peek();
                    if (!std::isdigit(c) && c != '-') {
                        string vec_kind;
                        in >> vec_kind;
                        if (vec_kind == "unit_vector") {
                            long pos = 0;
                            in >> pos;
                            if (in.fail()) {
                                cerr << "Error while reading " << type_string << " as a unit_vector form the input!" << endl;
                                throw BadInputException();
                            }

                            vector< vector<Integer> > e_i = vector< vector<Integer> >(1,vector<Integer>(dim+type_nr_columns_correction(input_type),0));
                            if (pos < 1 || pos > static_cast<long>(e_i[0].size())) {
                                cerr << "Error while reading " << type_string << " as a unit_vector "<< pos <<" form the input!" << endl;
                                throw BadInputException();
                            }
                            pos--; // in input file counting starts from 1
                            e_i[0].at(pos) = 1;
                            save_matrix(input_map, input_type, type_string, e_i);
                            continue;
                        }
                    }
                } else {
                    in >> nr_rows;
                }
                nr_columns = dim + type_nr_columns_correction(input_type);
                if(in.fail() || nr_rows < 0) {
                    cerr << "Error while reading " << type_string << " (a "<<nr_rows<<"x"<<nr_columns<<" matrix) form the input!" << endl;
                    throw BadInputException();
                }
                vector< vector<Integer> > M(nr_rows,vector<Integer>(nr_columns));
                for(i=0; i<nr_rows; i++){
                    for(j=0; j<nr_columns; j++) {
                        in >> M[i][j];
                    }
                }
                save_matrix(input_map, input_type, type_string, M);
            }
            if (in.fail()) {
                cerr << "Error while reading " << type_string << " (a "<<nr_rows<<"x"<<nr_columns<<" matrix) form the input!" << endl;
                throw BadInputException();
            }
        }
    } else {
        // old input syntax
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

            in >> type_string;

            if ( in.fail() ) {
                cerr << "Error while reading a "<<nr_rows<<"x"<<nr_columns<<" matrix form the input!" << endl;
                throw BadInputException();
            }

            input_type = to_type(type_string);

            //check if this type already exists
            save_matrix(input_map, input_type, type_string, M);
        }
    }
    return input_map;
}
