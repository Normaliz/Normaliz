#include "output.h"

InputType to_type(string& type_string) {
    if (type_string=="0"||type_string=="integral_closure") {
        return Type::integral_closure;
    }
    if (type_string=="1"||type_string=="normalization") {
        return Type::normalization;
    }
    if (type_string=="2"||type_string=="polytope") {
        return Type::polytope;
    }
    if (type_string=="3"||type_string=="rees_algebra") {
        return Type::rees_algebra;
    }
    if (type_string=="4"||type_string=="hyperplanes") {
        return Type::hyperplanes;
    }
    if (type_string=="5"||type_string=="equations") {
        return Type::equations;
    }
    if (type_string=="6"||type_string=="congruences") {
        return Type::congruences;
    }
    if (type_string=="10"||type_string=="lattice_ideal") {
        return Type::lattice_ideal;
    }
    if (type_string=="grading") {
        return Type::grading;
    }
    
    cerr<<"Warning: Unknown type \""<<type_string<<"\"! Will be replaced with type integral_closure."<<endl;
    return Type::integral_closure;
}


template <typename Integer>
map <Type::InputType, vector< vector<Integer> > > readNormalizInput (istream& in, Output<Integer>& O) {

    string type_string;
    size_t i,j;
    size_t nr_rows,nr_columns;;
    InputType input_type = Type::integral_closure;
    Integer number;
    map<Type::InputType, vector< vector<Integer> > > input_map;
    typename map<Type::InputType, vector< vector<Integer> > >::iterator it;

    while (in.good()) {
        in >> nr_rows;
        if(in.fail())
            break;
        in >> nr_columns;
        Matrix<Integer> M(nr_rows,nr_columns);
        for(i=1; i<=nr_rows; i++){
            for(j=1; j<=nr_columns; j++) {
                in >> number;
                M.write(i,j,number);
            }
        }

        in>>type_string;

        if ( in.fail() ) {
            throw BadInputException();
        }

        input_type = to_type(type_string);
        if (input_type == Type::polytope)
            O.set_type(OT_POLYTOPE);
        if (input_type == Type::rees_algebra)
            O.set_type(OT_REES);

        //check if this type already exists and merge data then
        it = input_map.find(input_type);
        if (it == input_map.end()) {
            input_map.insert(make_pair(input_type, M.get_elements()));
        } else { //in this case we merge the data
            vector< vector<Integer> > v = M.get_elements();
            it->second.insert(it->second.end(), v.begin(), v.end());
            //TODO shoud not be necessary to create copy v
            //M.get_elements().begin(), M.get_elements.end());
        }
    }

    return input_map;
}
