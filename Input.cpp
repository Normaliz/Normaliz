#include "libnormaliz/libnormaliz.h"

template <typename Integer>
map <Type::InputType, vector< vector<Integer> > > readNormalizInput (istream& in) {

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
        vector< vector<Integer> > M(nr_rows,vector<Integer>(nr_columns));
        for(i=0; i<nr_rows; i++){
            for(j=0; j<nr_columns; j++) {
                in >> number;
                M[i][j] = number;
            }
        }

        in>>type_string;

        if ( in.fail() ) {
            throw BadInputException();
        }

        input_type = to_type(type_string);

        //check if this type already exists and merge data then
        it = input_map.find(input_type);
        if (it == input_map.end()) {
            input_map.insert(make_pair(input_type, M));
        } else { //in this case we merge the data
            it->second.insert(it->second.end(), M.begin(), M.end());
        }
    }

    return input_map;
}
