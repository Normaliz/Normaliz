#include "libnormaliz.cpp"

namespace libnormaliz {

template class Cone<long long int>;
template class Cone<mpz_class>;

template class Matrix<long long int>;
template class Matrix<mpz_class>;

template class Sublattice_Representation<long long int>;
template class Sublattice_Representation<mpz_class>;

template class Lineare_Transformation<long long int>;
template class Lineare_Transformation<mpz_class>;

template int decimal_length<long long int>(long long int);
template int decimal_length<mpz_class>(mpz_class);

}
