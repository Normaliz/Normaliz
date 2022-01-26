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
#include <fstream>

#include "libnormaliz/options.h"
#include "libnormaliz/input_type.h"
#include "libnormaliz/list_and_map_operations.h"
#include "libnormaliz/cone_property.h"
#include "libnormaliz/matrix.h"

#ifndef NORMALIZ_INPUT_H
#define NORMALIZ_INPUT_H


namespace libnormaliz {

template <typename Number>
map<Type::InputType, vector<vector<Number> > > readNormalizInput(istream& in,
                                                                 OptionsHandler& options,
                                                                 map<NumParam::Param, long>& num_param_input,
                                                                 string& polynomial,
                                                                 renf_class_shared& number_field);

// here defined for use in interfaces
void read_number_field_strings(istream& in, string& mp_string, string& indet, string& emb_string);

//---------------------------------------------------------------------------
//                     Number input
//---------------------------------------------------------------------------

inline mpq_class mpq_read(istream& in) {
    const string numeric = "+-0123456789/.e";
    in >> std::ws;
    string s;
    char c;
    bool is_float = false;
    while (in.good()) {
        c = in.peek();
        size_t pos = numeric.find(c);
        if (pos == string::npos)
            break;
        if (pos > 12)
            is_float = true;
        in >> c;
        s += c;
    }

    if (s == "") {
        string t;
        t += c;
        throw BadInputException("Empty number string preceding character " + t +
                                ". Most likely mismatch of amb_space and matrix format or forgotten keyword.");
    }

    // cout << "t " << s << " f " << is_float << endl;

    if (s[0] == '+')
        s = s.substr(1);  // must suppress + sign for mpq_class

    try {
        if (!is_float) {
            return mpq_class(s);
        }
        else
            return dec_fraction_to_mpq(s);
    } catch (const std::exception& e) {
        cerr << e.what() << endl;
        throw BadInputException("Illegal number string " + s + " in input, Exiting.");
    }
}

// To be used in input.cpp
inline void string2coeff(mpq_class& coeff, istream& in, const string& s) {  // in here superfluous parameter

    stringstream sin(s);
    coeff = mpq_read(sin);
    // coeff=mpq_class(s);
}

// To be used from other sources
inline void string2coeff(mpq_class& coeff, const string& s) {

    // cout << "SSSSSS " << s << endl;

    const string numeric = "+-0123456789/.e "; // must allow blank
    for(auto& c: s){
        size_t pos = numeric.find(c);
        if(pos == string::npos)
            throw BadInputException("Illegal character in numerical string");
    }


    stringstream sin(s);
    coeff = mpq_read(sin);
    // coeff=mpq_class(s);
}

inline void read_number(istream& in, mpq_class& number) {
    number = mpq_read(in);
}

inline void read_number(istream& in, long& number) {
    in >> number;
}

inline void read_number(istream& in, long long& number) {
    in >> number;
}

inline void read_number(istream& in, nmz_float& number) {
    in >> number;
}

inline void read_number(istream& in, mpz_class& number) {
    in >> number;
}

#ifdef ENFNORMALIZ

inline void string2coeff(renf_elem_class& coeff, istream& in, const string& s) {  // we need in to access the renf

    try {
        coeff = renf_elem_class(*renf_class::get_pword(in), s);
    } catch (const std::exception& e) {
        cerr << e.what() << endl;
        throw BadInputException("Illegal number string " + s + " in input, Exiting.");
    }
}

inline void read_number(istream& in, renf_elem_class& number) {
    // in >> number;

    char c;

    in >> ws;
    c = in.peek();
    if (c != '(' && c != '\'' && c != '\"') {  // rational number
        mpq_class rat = mpq_read(in);
        number = renf_elem_class(rat);
        return;
    }

    // now we have a proper field element

    in >> c;  // read (

    string num_string;
    bool skip = false;
    while (in.good()) {
        c = in.peek();
        if (c == ')' || c == '\'' || c == '\"') {
            in >> c;
            break;
        }
        if (c == '~' || c == '=' || c == '[')  // skip the approximation
            skip = true;
        in.get(c);
        if (in.fail())
            throw BadInputException("Error in reading number: field element not terminated");
        if (!skip)
            num_string += c;
    }
    string2coeff(number, in, num_string);
}
#endif

template <typename Integer>
Matrix<Integer> readMatrix(const string project) {
    // reads one matrix from file with name project
    // format: nr of rows, nr of colimns, entries
    // all separated by white space

    string name_in = project;
    const char* file_in = name_in.c_str();
    ifstream in;
    in.open(file_in, ifstream::in);
    if (in.is_open() == false)
        throw BadInputException("readMatrix cannot find file");
    int nrows, ncols;
    in >> nrows;
    in >> ncols;

    if (nrows == 0 || ncols == 0)
        throw BadInputException("readMatrix finds matrix empty");

    int i, j;
    Matrix<Integer> result(nrows, ncols);

    for (i = 0; i < nrows; ++i)
        for (j = 0; j < ncols; ++j) {
            read_number(in, result[i][j]);
            if (in.fail())
                throw BadInputException("readMatrix finds matrix corrupted");
        }
    return result;
}

} // namespace

#endif
