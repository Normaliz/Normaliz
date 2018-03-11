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

#include "Qoptions.h"
#include "libQnormaliz/libQnormaliz.h"
#include "libQnormaliz/Qmap_operations.h"
#include "libQnormaliz/Qcone_property.h"

// eats up a comment, stream must start with "/*", eats everything until "*/"
void skip_comment(istream& in) {
    int i = in.get();
    int j = in.get();
    if (i != '/' || j != '*') {
        throw BadInputException("Bad comment start!");
    }
    while (in.good()) {
        in.ignore(numeric_limits<streamsize>::max(), '*'); //ignore everything until next '*'
        i = in.get();
        if (in.good() && i == '/') return; // successfully skipped comment
    }
    throw BadInputException("Incomplete comment!");
}

template<typename Number>
void save_matrix(map<Type::InputType, vector<vector<Number> > >& input_map,
        InputType input_type, const vector<vector<Number> >& M) {
    //check if this type already exists
    if (exists_element(input_map, input_type)) {
        /*throw BadInputException("Multiple inputs of type \"" + type_string
                + "\" are not allowed!");*/
	input_map[input_type].insert(input_map[input_type].end(),M.begin(),M.end());
        return;
    }
    input_map[input_type] = M;
}

template<typename Number>
void save_empty_matrix(map<Type::InputType, vector<vector<Number> > >& input_map,
        InputType input_type){
    
    vector<vector<Number> > M;
    save_matrix(input_map, input_type, M);   
}

template <typename Number>
vector<vector<Number> > transpose_mat(const vector<vector<Number> >& mat){

    if(mat.size()==0 || mat[0].size()==0)
        return vector<vector<Number> >(0);
    size_t m=mat[0].size();
    size_t n=mat.size();
    vector<vector<Number> > transpose(m,vector<Number> (n,0));
    for(size_t i=0;i<m;++i)
        for(size_t j=0;j<n;++j)
            transpose[i][j]=mat[j][i];
    return transpose;
}

template <typename Number>
void append_row(const vector<Number> row, map <Type::InputType, vector< vector<Number> > >& input_map,
                    Type::InputType input_type) {
    
    vector<vector<Number> > one_row(1,row);
    save_matrix(input_map,input_type,one_row); 
}

template <typename Number>
void process_constraint(const string& rel, const vector<Number>& left, Number right, const Number modulus, 
                        map <Type::InputType, vector< vector<Number> > >& input_map, bool forced_hom) {
    
    vector<Number> row=left;
    bool inhomogeneous=false;
    if(right!=0 || rel=="<" || rel==">")
        inhomogeneous=true;
    string modified_rel=rel;
    bool strict_inequality=false;
    if(rel=="<"){
        strict_inequality=true;
        right-=1;
        modified_rel="<=";
        
    }
    if(rel==">"){
        strict_inequality=true;
        right+=1;
        modified_rel=">=";
    }
    if(strict_inequality && forced_hom){
            throw BadInputException("Strict inequality not allowed in hom_constraints!");
    }
    if(inhomogeneous || forced_hom)
        row.push_back(-right); // rhs --> lhs
    if(modified_rel=="<="){ // convert <= to >=
        for(size_t j=0; j<row.size();++j)
            row[j]=-row[j];
        modified_rel=">=";
    }
    if(rel=="~")
        row.push_back(modulus);

    if(inhomogeneous && !forced_hom){
        if(modified_rel=="="){
            append_row(row,input_map,Type::inhom_equations);
            return;
        }
        if(modified_rel==">="){
            append_row(row,input_map,Type::inhom_inequalities);
            return;
        }
        if(modified_rel=="~"){
            append_row(row,input_map,Type::inhom_congruences);
            return;
        }
    }
    else {
        if(modified_rel=="="){
            append_row(row,input_map,Type::equations);
            return;
        }
        if(modified_rel==">="){
            append_row(row,input_map,Type::inequalities);
            return;
        }
        if(modified_rel=="~"){
            append_row(row,input_map,Type::congruences);
            return;
        }                
    }
    throw BadInputException("Illegal constrint type "+rel+" !");
}

template <typename Number>
bool read_modulus(istream& in, Number& modulus) {

    in >> std::ws;  // gobble any leading white space
    char dummy;
    in >> dummy;
    if(dummy != '(')
      return false;
    in >> modulus;
    if(in.fail() || modulus==0)
        return false;
    in >> std::ws;  // gobble any white space before closing
    in >> dummy;
    if(dummy != ')')
        return false;
    return true;
}

template <typename Number>
void read_constraints(istream& in, long dim, map <Type::InputType, vector< vector<Number> > >& input_map, bool forced_hom) {

    long nr_constraints;
    in >> nr_constraints;
    
    if(in.fail() || nr_constraints < 0) {
        throw BadInputException("Cannot read "
        + toString(nr_constraints) + " constraints!");
    }
    long hom_correction=0;
    if(forced_hom)
        hom_correction=1;
    for(long i=0;i< nr_constraints; ++i) {
        vector<Number> left(dim-hom_correction);
        for(long j=0;j<dim-hom_correction;++j){
            in >> left[j];
        }
        string rel, modulus_str;
        Number right, modulus=0;
        in >> rel;
        in >> right;
        if(rel=="~") {
            if(!read_modulus(in,modulus))
                throw BadInputException("Error while reading modulus of congruence!");
        }
        if (in.fail()) {
            throw BadInputException("Error while reading constraint!");
        }
        process_constraint(rel,left,right,modulus,input_map,forced_hom);        
    }
}

template <typename Number>
bool read_sparse_vector(istream& in, vector<Number>& input_vec, long length){
    
    input_vec=vector<Number> (length,0);
    char dummy;
    
    while(true){
        in >> std::ws;
        int c = in.peek();
        if(c==';'){
            in >> dummy; // swallow ;
            return true;
        }
        long pos;
        in >> pos;
        if(in.fail())
            return false;
        pos--;
        if(pos<0 or pos>=length)
            return false;
        in >> std::ws;
        c=in.peek();
        if(c!=':')
            return false;
        in >> dummy; // skip :
        Number value;
        in >> value;
        if(in.fail())
            return false;
        input_vec[pos]=value;        
    }
}

template <typename Number>
bool read_formatted_vector(istream& in, vector<Number>& input_vec) {

    input_vec.clear();
    in >> std::ws;
    char dummy;
    in >> dummy; // read first proper character
    if(dummy!='[')
        return false;
    bool one_more_entry_required=false;
    while(true){
        in >> std::ws;
        if(!one_more_entry_required && in.peek()==']'){
            in >> dummy;
            return true;
        }
        Number number;
        in >> number;
        if(in.fail())
            return false;
        input_vec.push_back(number);
        in >> std::ws;
        one_more_entry_required=false;
        if(in.peek()==',' || in.peek()==';'){  // skip potential separator
            in >> dummy;
            one_more_entry_required=true;
        }
    }
}

template <typename Number>
bool read_formatted_matrix(istream& in, vector<vector<Number> >& input_mat, bool transpose) {
    input_mat.clear();
    in >> std::ws;
    char dummy;
    in >> dummy; // read first proper character
    if(dummy!='[')
        return false;
    bool one_more_entry_required=false;
    while(true){
        in >> std::ws;
        if(!one_more_entry_required && in.peek()==']'){ // closing ] found
            in >> dummy;
            if(transpose)
                input_mat=transpose_mat(input_mat);
            return true;
        }
        vector<Number> input_vec;
        if(!read_formatted_vector(in,input_vec)){
            throw BadInputException("Error in reading input vector!");
        }
        if(input_mat.size()>0 && input_vec.size()!=input_mat[0].size()){
            throw BadInputException("Rows of input matrix have unequal lengths!");
        }            
        input_mat.push_back(input_vec);
        in >> std::ws;
        one_more_entry_required=false;
        if(in.peek()==',' || in.peek()==';'){ // skip potential separator
            in >> dummy;
            one_more_entry_required=true;
        }
    }
}

template <typename Number>
void read_number_field(istream &in)
{
    throw NumberFieldInputException();
}

#ifdef ENFNORMALIZ

template<>
void read_number_field<renf_elem_class>(istream &in)
{
    renf_class nf;
    in >> nf;
    if (in.fail()) {
        throw BadInputException("Could not read number field!");
    }
    in >> set_renf(nf);
}
#endif


template <typename Number>
map <Type::InputType, vector< vector<Number> > > readNormalizInput (istream& in, OptionsHandler& options) {

    string type_string;
    long i,j;
    long nr_rows,nr_columns,nr_rows_or_columns;
    InputType input_type;
    Number number;
    ConeProperty::Enum cp;

    map<Type::InputType, vector< vector<Number> > > input_map;
    typename map<Type::InputType, vector< vector<Number> > >::iterator it;

    in >> std::ws;  // eat up any leading white spaces
    int c = in.peek();
    if ( c == EOF ) {
        throw BadInputException("Empty input file!");
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
            throw BadInputException("First entry must be \"amb_space\"!");
        }
        bool dim_known=false;
        in >> std::ws;
        c=in.peek();
        if(c=='a'){
            string dummy;
            in >> dummy;
            if(dummy!="auto"){
                throw BadInputException("Bad amb_space value!");
            }
        }
        else{            
            in >> dim;
            if (!in.good() || dim <= 0) {
                throw BadInputException("Bad amb_space value!");
            }
            dim_known=true;
        }
        while (in.good()) {    //main loop
            
            bool transpose=false;
            in >> std::ws;  // eat up any leading white spaces
            c = in.peek();
            if (c == EOF) break;
            if (c == '/') {
                skip_comment(in);
            } else {
                in >> type_string;
                if (in.fail()) {
                    throw BadInputException("Could not read type string!");
                }
                if (std::isdigit(c)) {
                    throw BadInputException("Unexpected number " + type_string
                            + " when expecting a type!");
                }
                if (isConeProperty(cp, type_string)) {
                    options.activateInputFileConeProperty(cp);
                    continue;
                }
                /* if (type_string == "BigInt") {
                    options.activateInputFileBigInt();
                    continue;
                } */
                if (type_string == "LongLong") {
                    options.activateInputFileLongLong();
                    continue;
                }
                if (type_string == "number_field") {
                    read_number_field<Number>(in);
                    continue;
                }
                if (type_string == "total_degree") {
                    if(!dim_known){
                        throw BadInputException("Ambient space must be known for "+type_string+"!");
                    }
                    input_type = Type::grading;
                    save_matrix(input_map, input_type, vector< vector<Number> >(1,vector<Number>(dim+type_nr_columns_correction(input_type),1)));
                    continue;
                }
                if (type_string == "nonnegative") {
                    if(!dim_known){
                        throw BadInputException("Ambient space must be known for "+type_string+"!");
                    }
                    input_type = Type::signs;
                    save_matrix(input_map, input_type, vector< vector<Number> >(1,vector<Number>(dim+type_nr_columns_correction(input_type),1)));
                    continue;
                }
                if(type_string == "constraints") {
                    if(!dim_known){
                        throw BadInputException("Ambient space must be known for "+type_string+"!");
                    }
                    read_constraints(in,dim,input_map,false);
                    continue;
                }
                if(type_string == "hom_constraints") {
                    if(!dim_known){
                        throw BadInputException("Ambient space must be known for "+type_string+"!");
                    }
                    read_constraints(in,dim,input_map,true);
                    continue;
                }


                input_type = to_type(type_string);
                if(dim_known)
                    nr_columns = dim + type_nr_columns_correction(input_type);

                if (type_is_vector(input_type)) {
                    nr_rows_or_columns = nr_rows = 1;
                    in >> std::ws;  // eat up any leading white spaces
                    c = in.peek();
                    if (c=='u') { // must be unit vector
                        string vec_kind;
                        in >> vec_kind;
                        if (vec_kind != "unit_vector") {
                            throw BadInputException("Error while reading "
                            + type_string 
                            + ": unit_vector expected!");                            
                        }

                        long pos = 0;
                        in >> pos;
                        if (in.fail()) {
                            throw BadInputException("Error while reading "
                                    + type_string 
                                    + " as a unit_vector!");
                        }
                        
                        if(!dim_known){
                            throw BadInputException("Ambient space must be known for unit vector "+type_string+"!");
                        }

                        vector< vector<Number> > e_i = vector< vector<Number> >(1,vector<Number>(nr_columns,0));
                        if (pos < 1 || pos > static_cast<long>(e_i[0].size())) {
                            throw BadInputException("Error while reading "
                                    + type_string + " as a unit_vector "
                                    + toString(pos) + "!");
                        }
                        pos--; // in input file counting starts from 1
                        e_i[0].at(pos) = 1;
                        save_matrix(input_map, input_type, e_i);
                        continue;
                    } // end unit vector
                    
                    if(c=='s'){   // must be "sparse"
                        string vec_kind;
                        in >> vec_kind;
                        if (vec_kind != "sparse") {
                            throw BadInputException("Error while reading "
                            + type_string 
                            + ": sparse vector expected!");                            
                        }
                        
                        if(!dim_known){
                            throw BadInputException("Ambient space must be known for sparse vector "+type_string+"!");
                        }

                        vector<Number> sparse_vec;
                        nr_columns = dim + type_nr_columns_correction(input_type);
                        bool success = read_sparse_vector(in,sparse_vec,nr_columns);
                        if(!success){
                            throw BadInputException("Error while reading "
                            + type_string 
                            + " as a sparse vector!");
                        }
                        save_matrix(input_map, input_type, vector<vector<Number> > (1,sparse_vec)); 
                        continue;                        
                    }
                    
                    if (c == '[') { // must be formatted vector
                        vector<Number> formatted_vec;
                        bool success = read_formatted_vector(in,formatted_vec);
                        if(!dim_known){
                            dim=formatted_vec.size()- type_nr_columns_correction(input_type);
                            dim_known=true;
                            nr_columns = dim + type_nr_columns_correction(input_type);
                        }
                        if(!success || (long) formatted_vec.size()!=nr_columns){
                            throw BadInputException("Error while reading "
                            + type_string 
                            + " as a formatted vector!");
                        }
                        save_matrix(input_map, input_type, vector<vector<Number> > (1,formatted_vec)); 
                        continue;
                    } // end formatted vector

                } else {  // end vector, it is a matrix. Plain vector read as a one row matrix later on
                    in >> std::ws;
                    c = in.peek();
                    
                    if(c!='[' && !std::isdigit(c)){ // must be transpose
                        string transpose_str;
                        in >> transpose_str;
                        if(transpose_str!="transpose"){
                                throw BadInputException("Illegal keyword "+transpose_str+" following matrix type!");
                        }
                        transpose=true;
                        in >> std::ws;
                        c = in.peek();                                               
                    }
                    if(c=='['){ // it is a formatted matrix
                        vector<vector<Number> > formatted_mat;
                        bool success=read_formatted_matrix(in,formatted_mat, transpose);
                        if(!success){
                            throw BadInputException("Error while reading formatted matrix "
                            + type_string + "!");    
                        }
                        if(formatted_mat.size() ==0){ // empty matrix
                            input_type = to_type(type_string);
                            save_empty_matrix(input_map, input_type);
                            continue;
                        }
                        if(!dim_known){
                            dim=formatted_mat[0].size()- type_nr_columns_correction(input_type);
                            dim_known=true;
                            nr_columns = dim + type_nr_columns_correction(input_type);
                        }
                        
                        if((long) formatted_mat[0].size()!=nr_columns){
                            throw BadInputException("Error while reading formatted matrix "
                                + type_string + "!");    
                        }
                        
                        save_matrix(input_map, input_type, formatted_mat);
                        continue;
                    }  // only plain matrix left
                    
                    in >> nr_rows_or_columns; // is number of columns if transposed
                    nr_rows=nr_rows_or_columns; // most of the time
                }
                
                if(!dim_known){
                    throw BadInputException("Ambient space must be known for plain matrix or vector "+type_string+"!");
                }
                
                if(transpose)
                    swap(nr_rows,nr_columns);

                if(in.fail() || nr_rows_or_columns < 0) {
                    throw BadInputException("Error while reading "
                            + type_string + " (a " + toString(nr_rows)
                            + "x" + toString(nr_columns)
                            + " matrix) !");
                }
                if(nr_rows==0){
                    input_type = to_type(type_string);
                    save_empty_matrix(input_map, input_type);
                    continue;
                }
                
                vector< vector<Number> > M(nr_rows);
                in >> std::ws;
                c=in.peek();
                if(c=='s'){ // must be sparse
                    string sparse_test;
                    in >> sparse_test;
                    if (sparse_test!= "sparse") {
                        throw BadInputException("Error while reading "
                        + type_string 
                        + ": sparse matrix expected!");                            
                    }
                    for(long i=0;i<nr_rows;++i){
                        bool success=read_sparse_vector(in,M[i],nr_columns);
                        if(!success){
                            throw BadInputException("Error while reading "
                            + type_string 
                            + ": corrupted sparse matrix");                        
                        }
                        
                    }
                } else{ // dense matrix
                    for(i=0; i<nr_rows; i++){
                        M[i].resize(nr_columns);
                        for(j=0; j<nr_columns; j++) {
                            in >> M[i][j];
                        }
                    }
                }
                if(transpose)
                    M=transpose_mat(M);
                save_matrix(input_map, input_type, M);
            }
            if (in.fail()) {
                throw BadInputException("Error while reading " + type_string
                        + " (a " + toString(nr_rows) + "x"
                        + toString(nr_columns) + " matrix) !");
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
                throw BadInputException("Error while reading a "
                        + toString(nr_rows) + "x" + toString(nr_columns)
                        + " matrix !");
            }
            vector< vector<Number> > M(nr_rows,vector<Number>(nr_columns));
            for(i=0; i<nr_rows; i++){
                for(j=0; j<nr_columns; j++) {
                    in >> number;
                    M[i][j] = number;
                }
            }

            in >> type_string;

            if ( in.fail() ) {
                throw BadInputException("Error while reading a " 
                        + toString(nr_rows) + "x" + toString(nr_columns)
                        + " matrix!");
            }

            input_type = to_type(type_string);

            //check if this type already exists
            save_matrix(input_map, input_type, M);
        }
    }
    return input_map;
}
