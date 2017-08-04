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

#include "options.h"
#include "libnormaliz/libnormaliz.h"
#include "libnormaliz/map_operations.h"
#include "libnormaliz/cone_property.h"

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

void save_matrix(map<Type::InputType, vector<vector<mpq_class> > >& input_map,
        InputType input_type, const vector<vector<mpq_class> >& M) {
    //check if this type already exists
    if (exists_element(input_map, input_type)) {
        /*throw BadInputException("Multiple inputs of type \"" + type_string
                + "\" are not allowed!");*/
	input_map[input_type].insert(input_map[input_type].end(),M.begin(),M.end());
        return;
    }
    input_map[input_type] = M;
}

void save_empty_matrix(map<Type::InputType, vector<vector<mpq_class> > >& input_map,
        InputType input_type){
    
    vector<vector<mpq_class> > M;
    save_matrix(input_map, input_type, M);   
}


vector<vector<mpq_class> > transpose_mat(const vector<vector<mpq_class> >& mat){

    if(mat.size()==0 || mat[0].size()==0)
        return vector<vector<mpq_class> >(0);
    size_t m=mat[0].size();
    size_t n=mat.size();
    vector<vector<mpq_class> > transpose(m,vector<mpq_class> (n,0));
    for(size_t i=0;i<m;++i)
        for(size_t j=0;j<n;++j)
            transpose[i][j]=mat[j][i];
    return transpose;
}


void append_row(const vector<mpq_class> row, map <Type::InputType, vector< vector<mpq_class> > >& input_map,
                    Type::InputType input_type) {
    
    vector<vector<mpq_class> > one_row(1,row);
    save_matrix(input_map,input_type,one_row); 
}


void process_constraint(const string& rel, const vector<mpq_class>& left, mpq_class right, const mpq_class modulus, 
                        map <Type::InputType, vector< vector<mpq_class> > >& input_map, bool forced_hom) {
    
    vector<mpq_class> row=left;
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

mpq_class mpq_read(istream& in){
    const string numeric="+-0123456789/.e";
    in >> std::ws;
    string s;
    bool is_float=false;
    while(true){
        char c = in.peek();
        size_t pos=numeric.find(c);
        if(pos==string::npos)
            break;
        if(pos>12)
            is_float=true;
        in >> c;
            s+=c;
    }
    
    if(s=="")
        throw BadInputException("Error in input file. Most lekely mismatch of amb_space and matrix format.");
    
    // cout << "t " << s << " f " << is_float << endl;
    
    if(!is_float)
        return mpq_class(s);
    
    return dec_fraction_to_mpq(s);
}


bool read_modulus(istream& in, mpq_class& modulus) {

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


void read_symbolic_constraint(istream& in, string& rel, vector<mpq_class>& left, mpq_class& right, mpq_class& modulus, bool forced_hom) {
    
    bool congruence=false;
    bool modulus_read=false;
    mpq_class side=1,sign;
    right=0;
    long hom_correction=0;
    if(forced_hom)
        hom_correction=1;
    
    in >> std::ws;
    char c = in.peek();
    
    while(true){
        if(c=='('){   
            if(modulus_read || !congruence || !read_modulus(in,modulus))
                    throw BadInputException("Error while reading modulus of congruence!");
            modulus_read=true;
            in >> std::ws;
            c = in.peek();
        }
        if(modulus_read && c!=';')
            throw BadInputException("Error while reading modulus of congruence!");
        if(c==';'){
            if(rel=="")
                throw BadInputException("Missing relation in constraint");
            in >> c;
            if(congruence && !modulus_read)
                throw BadInputException("Modulus missing in congrruence");
            // cout << "LLLLL " << left << " " << rel << " RRR " << right << endl;
            return;
        }
        
        bool rel_read=false;
        
        if(c=='~' || c=='<' || c=='>' || c=='='){
            rel_read=true;
            if(rel!="")
                throw BadInputException("Error while reading relation in constraint!");                
            if(c=='~')
                congruence=true;
            in >> c;
            rel+=c;
        }
        c = in.peek();
        if(rel!="" && (rel=="<" || rel==">") && c=='='){
            in >> c;
            rel+=c;
        }
        in >> std::ws;
        c = in.peek();     
        if(rel_read){
            side=-1;
            continue;
        }
        sign=1;
        if(c=='-'){
            sign=-1;
            in >> c;            
        }
        if(c=='+'){
            in >> c;            
        }
        mpq_class entry=1;
        in >> std::ws;
        c = in.peek();
        if(c!='x'){
            if(c=='+' || c=='-')
                throw BadInputException("Double sign in constraint");
            entry=mpq_read(in);
            if(in.fail())
                throw BadInputException("Error while reading coefficient in constraint");
            in >> std::ws;
            c = in.peek();
        }
        if(c!='x'){
            right-=side*sign*entry;
            continue;            
        }
        in >> c;
        in >> std::ws;
        c = in.peek();
        if(c!='[')
            throw BadInputException("Error while reading index in constraint");
        in >> c;
        long index;
        in >> index;
        if(in.fail() || index <1 || index+hom_correction> (long) left.size())
            throw BadInputException("Error while reading index in constraint");
        index-=1;
        left[index]+=side*sign*entry;
        in >> std::ws;
        c = in.peek();
        if(c!=']')
            throw BadInputException("Error while reading index in constraint");
        in >> c;
        in >> std::ws;
        c = in.peek();
        continue;
    }
    
}



void read_constraints(istream& in, long dim, map <Type::InputType, vector< vector<mpq_class> > >& input_map, bool forced_hom) {

    long nr_constraints;
    in >> nr_constraints;
    
    if(in.fail() || nr_constraints < 0) {
        throw BadInputException("Cannot read "
        + toString(nr_constraints) + " constraints!");
    }
    if(nr_constraints==0)
        return;
    
    bool symbolic=false;
    
    in >> std::ws;
    int c = in.peek();
    if(c=='s'){
        string dummy;
        in >> dummy;
        if(dummy!="symbolic")
            throw BadInputException("Illegal keyword " + dummy
                                + " in input!");
        symbolic=true;
    }  
    
    long hom_correction=0;
    if(forced_hom)
        hom_correction=1;
    for(long i=0;i< nr_constraints; ++i) {
        
        vector<mpq_class> left(dim-hom_correction);
        string rel, modulus_str;
        mpq_class right, modulus=0;
        
        if(symbolic){
            read_symbolic_constraint(in,rel,left,right,modulus,forced_hom);            
        }
        else{ // ordinary constraint read here
            for(long j=0;j<dim-hom_correction;++j){
                left[j]=mpq_read(in);
            }
            in >> rel;
            right=mpq_read(in);
            if(rel=="~") {
                if(!read_modulus(in,modulus))
                    throw BadInputException("Error while reading modulus of congruence!");
            }
            if (in.fail()) {
                throw BadInputException("Error while reading constraint!");
            }
        }
        process_constraint(rel,left,right,modulus,input_map,forced_hom);        
    }
}

void read_polynomial(istream& in, string& polynomial) {

    char c;
    while(true){
        in >> c;
        if(in.fail())
                throw BadInputException("Error while reading polynomial!");
        if(c==';'){
            if(polynomial.size()==0)
                throw BadInputException("Error while reading polynomial!");
            return;
        }
        polynomial+=c;
    }
}


bool read_sparse_vector(istream& in, vector<mpq_class>& input_vec, long length){
    
    input_vec=vector<mpq_class> (length,0);
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
        if(pos<0 || pos>=length)
            return false;
        in >> std::ws;
        c=in.peek();
        if(c!=':')
            return false;
        in >> dummy; // skip :
        mpq_class value;
        // in >> value;
        value=mpq_read(in);
        if(in.fail())
            return false;
        input_vec[pos]=value;        
    }
}


bool read_formatted_vector(istream& in, vector<mpq_class>& input_vec) {

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
        mpq_class number;
        number=mpq_read(in);
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


bool read_formatted_matrix(istream& in, vector<vector<mpq_class> >& input_mat, bool transpose) {
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
        vector<mpq_class> input_vec;
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
    


map <Type::InputType, vector< vector<mpq_class> > > readNormalizInput (istream& in, OptionsHandler& options, string& polynomial, long& nr_coeff_quasipol) {

    string type_string;
    long i,j;
    long nr_rows,nr_columns,nr_rows_or_columns;
    InputType input_type;
    mpq_class number;
    ConeProperty::Enum cp;
    bool we_have_a_polynomial=false;
    bool we_have_nr_coeff=false;

    map<Type::InputType, vector< vector<mpq_class> > > input_map;
    typename map<Type::InputType, vector< vector<mpq_class> > >::iterator it;

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
                if (type_string == "NoExtRaysOutput") {
                    options.activateNoExtRaysOutput();
                    continue;
                }
                if (type_string == "total_degree") {
                    if(!dim_known){
                        throw BadInputException("Ambient space must be known for "+type_string+"!");
                    }
                    input_type = Type::grading;
                    save_matrix(input_map, input_type, vector< vector<mpq_class> >(1,vector<mpq_class>(dim+type_nr_columns_correction(input_type),1)));
                    continue;
                }
                if (type_string == "nonnegative") {
                    if(!dim_known){
                        throw BadInputException("Ambient space must be known for "+type_string+"!");
                    }
                    input_type = Type::signs;
                    save_matrix(input_map, input_type, vector< vector<mpq_class> >(1,vector<mpq_class>(dim+type_nr_columns_correction(input_type),1)));
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
                
                if(type_string == "polynomial") {
                    if(we_have_a_polynomial)
                        throw BadInputException("Only one polynomial allowed");
                    read_polynomial(in,polynomial);
                    we_have_a_polynomial=true;
                    continue;
                }
                if(type_string=="nr_coeff_quasipol"){
                    if(we_have_nr_coeff)
                        throw BadInputException("Only one nr_coeff_quasipol allowed");
                    in >> nr_coeff_quasipol;
                    we_have_nr_coeff=true;
                    if(in.fail())
                        throw BadInputException("Error while reading nr_coeff_quasipol");
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

                        vector< vector<mpq_class> > e_i = vector< vector<mpq_class> >(1,vector<mpq_class>(nr_columns,0));
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

                        vector<mpq_class> sparse_vec;
                        nr_columns = dim + type_nr_columns_correction(input_type);
                        bool success = read_sparse_vector(in,sparse_vec,nr_columns);
                        if(!success){
                            throw BadInputException("Error while reading "
                            + type_string 
                            + " as a sparse vector!");
                        }
                        save_matrix(input_map, input_type, vector<vector<mpq_class> > (1,sparse_vec)); 
                        continue;                        
                    }
                    
                    if (c == '[') { // must be formatted vector
                        vector<mpq_class> formatted_vec;
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
                        save_matrix(input_map, input_type, vector<vector<mpq_class> > (1,formatted_vec)); 
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
                        vector<vector<mpq_class> > formatted_mat;
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
                
                vector< vector<mpq_class> > M(nr_rows);
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
                            M[i][j]=mpq_read(in);
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
            vector< vector<mpq_class> > M(nr_rows,vector<mpq_class>(nr_columns));
            for(i=0; i<nr_rows; i++){
                for(j=0; j<nr_columns; j++) {
                    number=mpq_read(in);
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