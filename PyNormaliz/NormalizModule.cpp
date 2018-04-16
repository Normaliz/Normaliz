/***************************************************************************
 * 
 * Include
 * 
 ***************************************************************************/

#include <Python.h>
using namespace std;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <string>
using std::string;

#include <libnormaliz/cone.h>
#include <libnormaliz/map_operations.h>

using libnormaliz::Cone;
//using libnormaliz::ConeProperty;
using libnormaliz::ConeProperties;
using libnormaliz::Sublattice_Representation;
using libnormaliz::Type::InputType;

#include <vector>
using std::map;
using std::vector;
using std::pair;

#include<csignal>

typedef int py_size_t;

/***************************************************************************
 * 
 * Macros for exception handling
 * 
 ***************************************************************************/

#define FUNC_BEGIN try {

#define FUNC_END \
    } catch (libnormaliz::InterruptException& e ) {\
        PyOS_setsig(SIGINT,current_interpreter_sigint_handler);\
        libnormaliz::nmz_interrupted = false; \
        PyErr_SetInterrupt(); \
        PyErr_CheckSignals(); \
        return NULL; \
    } catch (libnormaliz::NormalizException& e) { \
        PyOS_setsig(SIGINT,current_interpreter_sigint_handler);\
        PyErr_SetString( NormalizError, e.what() ); \
        return NULL; \
    } catch( exception& e ) { \
        PyOS_setsig(SIGINT,current_interpreter_sigint_handler);\
        PyErr_SetString( PyNormaliz_cppError, e.what() ); \
        return NULL; \
    }

/***************************************************************************
 * 
 * Signal handling
 * 
 ***************************************************************************/

void signal_handler( int signal ){
    libnormaliz::nmz_interrupted = true;
}

/***************************************************************************
 * 
 * Static objects
 * 
 ***************************************************************************/


static PyObject * NormalizError;
static PyObject * PyNormaliz_cppError;
static const char* cone_name = "Cone";
static const char* cone_name_long = "Cone<long long>";
static string cone_name_str( cone_name );
static string cone_name_str_long( cone_name_long );

static PyOS_sighandler_t current_interpreter_sigint_handler;

static PyObject * RationalHandler = NULL;
static PyObject * VectorHandler = NULL;
static PyObject * MatrixHandler = NULL;

/***************************************************************************
 * 
 * Call func on one argument
 * 
 ***************************************************************************/

PyObject* CallPythonFuncOnOneArg( PyObject* function, PyObject* single_arg ){
    PyObject* single_arg_tuple = PyTuple_Pack(1,single_arg);
    PyObject* return_obj = PyObject_CallObject(function,single_arg_tuple);
    Py_DecRef(single_arg);
    Py_DecRef(single_arg_tuple);
    return return_obj;
}

/***************************************************************************
 * 
 * Compiler version control
 * 
 ***************************************************************************/

#if PY_MAJOR_VERSION >= 3
#define string_check PyUnicode_Check
#else
#define string_check PyString_Check
#endif

// Hacky 64-bit check. Works for windows and gcc, probably not clang.
// Check windows
#if _WIN32 || _WIN64
#if _WIN64
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif

// Check GCC
#if __GNUC__
#if __x86_64__ || __ppc64__ || __aarch64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif

#ifndef NMZ_RELEASE
    static_assert(false,
       "Your Normaliz version (unknown) is to old! Update to 3.5.4 or newer.");
#endif
#if NMZ_RELEASE < 30504
    static_assert(false, "Your Normaliz version is to old! Update to 3.5.4 or newer.");
#endif

/***************************************************************************
 * 
 * Python-C data conversion functions
 * 
 ***************************************************************************/

string PyUnicodeToString( PyObject* in ){
#if PY_MAJOR_VERSION >= 3
  string out = "";
  int length = PyUnicode_GET_SIZE( in );
  for( int i = 0; i < length; i++ ){
      out += PyUnicode_READ_CHAR( in, i );
  }
  return out;
#else
  char* out = PyString_AsString( in );
  return string(out);
#endif
}

PyObject * StringToPyUnicode( string in ){
#if PY_MAJOR_VERSION >= 3
  return PyUnicode_FromString( in.c_str() );
#else
  return PyString_FromString( in.c_str() );
#endif
}

// Boolean conversion

inline PyObject* BoolToPyBool( bool in ){
  return in ? Py_True : Py_False;
}

// Converting MPZ's to PyLong and back via strings. Worst possible solution ever.

bool PyNumberToNmz( PyObject * in, mpq_class& out ){
  if( PyFloat_Check( in ) ){
      mpq_class temp(PyFloat_AsDouble(in));
      out.swap(temp);
      return true;
  }
  PyObject * in_as_string = PyObject_Str( in );
  const char* in_as_c_string = PyUnicodeToString( in_as_string ).c_str();
  out.set_str( in_as_c_string, 10 );
  return true;
}

PyObject* NmzToPyNumber( mpz_class in ){
  string mpz_as_string = in.get_str();
  char* mpz_as_c_string = const_cast<char*>(mpz_as_string.c_str());
  char * pend;
  PyObject* ret_val = PyLong_FromString( mpz_as_c_string, &pend, 10 );
  return ret_val;
}

PyObject* NmzToPyList( mpq_class in ){
    PyObject* out_list = PyList_New( 2 );
    PyList_SetItem( out_list, 0, NmzToPyNumber( in.get_num() ) );
    PyList_SetItem( out_list, 1, NmzToPyNumber( in.get_den() ) );
    if(RationalHandler!=NULL)
        out_list = CallPythonFuncOnOneArg(RationalHandler,out_list);
    return out_list;
}

bool PyNumberToNmz( PyObject* in, long long & out ){
  
  int overflow;
  out = PyLong_AsLongLongAndOverflow( in, &overflow );
  if( overflow == -1 )
    return false;
  return true;
  
}

PyObject* NmzToPyNumber( long long in ){
  
  return PyLong_FromLongLong( in );
  
}

PyObject* NmzToPyNumber( libnormaliz::key_t in ){
  
  return PyLong_FromLong( in );
  
}

#ifdef ENVIRONMENT64
PyObject* NmzToPyNumber( size_t in ){
  
  return PyLong_FromLong( in );
  
}
#endif

PyObject* NmzToPyNumber( long in ){
  
  return PyLong_FromLong( in );
  
}

PyObject* NmzToPyNumber( double in ){
  
  return PyFloat_FromDouble( in );
  
}

template<typename Integer>
bool PyNumberToNmz(Integer& x, Integer &out){
  
  return Integer::unimplemented_function;
  
}

template<typename Integer>
PyObject* NmzToPyNumber(Integer &in){
  
  return Integer::unimplemented_function;
  
}

template<typename Integer>
static bool PyListToNmz( vector<Integer>& out, PyObject* in ){
  if (!PyList_Check(in))
        return false;
    const int n = PyList_Size(in);
    out.resize(n);
    for (int i = 0; i < n; ++i) {
        PyObject* tmp = PyList_GetItem(in, i);
        if (!PyNumberToNmz(tmp, out[i]))
            return false;
    }
    return true;
}

template<typename Integer>
static bool PyIntMatrixToNmz( vector<vector<Integer> >& out, PyObject* in ){
  if (!PyList_Check( in ) )
        return false;
    const int nr = PyList_Size( in );
    out.resize(nr);
    for (int i = 0; i < nr; ++i) {
        bool okay = PyListToNmz(out[i], PyList_GetItem(in, i));
        if (!okay)
            return false;
    }
    return true;
}

template<typename Integer>
static bool PyInputToNmz( vector<vector<Integer> >& out, PyObject* in ){
    bool check_input;
    check_input = PyIntMatrixToNmz( out, in );
    if(check_input)
        return true;
    out.resize(1);
    check_input = PyListToNmz( out[0], in );
    if(check_input){
        return true;
    }
    return false;
}

template<typename Integer>
PyObject* NmzVectorToPyList(const vector<Integer>& in)
{
    PyObject* vector;
    const size_t n = in.size();
    vector = PyList_New(n);
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(vector, i, NmzToPyNumber(in[i]));
    }
    if(VectorHandler!=NULL)
        vector = CallPythonFuncOnOneArg(VectorHandler,vector);
    return vector;
}

PyObject* NmzBoolVectorToPyList(const vector<bool>& in)
{
    PyObject* vector;
    const size_t n = in.size();
    vector = PyList_New(n);
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(vector, i, BoolToPyBool(in[i]));
    }
    if(VectorHandler!=NULL)
        vector = CallPythonFuncOnOneArg(VectorHandler,vector);
    return vector;
}

PyObject* NmzBoolMatrixToPyList(const vector< vector<bool> >& in)
{
    PyObject* matrix;
    const size_t n = in.size();
    matrix = PyList_New( n );
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(matrix, i, NmzBoolVectorToPyList(in[i]));
    }
    if(MatrixHandler!=NULL)
        matrix = CallPythonFuncOnOneArg(MatrixHandler,matrix);
    return matrix;
}

template<typename Integer>
PyObject* NmzMatrixToPyList(const vector< vector<Integer> >& in)
{
    PyObject* matrix;
    const size_t n = in.size();
    matrix = PyList_New( n );
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(matrix, i, NmzVectorToPyList(in[i]));
    }
    if(MatrixHandler!=NULL)
        matrix = CallPythonFuncOnOneArg(MatrixHandler,matrix);
    return matrix;
}

PyObject* NmzHilbertSeriesToPyList(const libnormaliz::HilbertSeries& HS, bool is_HSOP)
{   
    PyObject* return_list = PyList_New( 3 );
    if(is_HSOP){
        PyList_SetItem(return_list, 0, NmzVectorToPyList(HS.getHSOPNum()));
        PyList_SetItem(return_list, 1, NmzVectorToPyList(libnormaliz::to_vector(HS.getHSOPDenom())));
        PyList_SetItem(return_list, 2, NmzToPyNumber(HS.getShift()));
    }else{
        PyList_SetItem(return_list, 0, NmzVectorToPyList(HS.getNum()));
        PyList_SetItem(return_list, 1, NmzVectorToPyList(libnormaliz::to_vector(HS.getDenom())));
        PyList_SetItem(return_list, 2, NmzToPyNumber(HS.getShift()));
    }
    return return_list;
}

template<typename Integer>
PyObject* NmzWeightedEhrhartSeriesToPyList(const std::pair<libnormaliz::HilbertSeries,Integer>& HS)
{   
    PyObject* return_list = PyList_New( 4 );
    PyList_SetItem(return_list, 0, NmzVectorToPyList(HS.first.getNum()));
    PyList_SetItem(return_list, 1, NmzVectorToPyList(libnormaliz::to_vector(HS.first.getDenom())));
    PyList_SetItem(return_list, 2, NmzToPyNumber(HS.first.getShift()));
    PyList_SetItem(return_list, 3, NmzToPyNumber(HS.second) );
    return return_list;
}

template<typename Integer>
PyObject* NmzHilbertQuasiPolynomialToPyList(const libnormaliz::HilbertSeries& HS)
{
    vector< vector<Integer> > HQ = HS.getHilbertQuasiPolynomial();
    const size_t n = HS.getPeriod();
    PyObject* return_list = PyList_New(n+1);
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(return_list, i, NmzVectorToPyList(HQ[i]));
    }
    PyList_SetItem(return_list, n, NmzToPyNumber(HS.getHilbertQuasiPolynomialDenom()));
    return return_list;
}

template<typename Integer>
PyObject* NmzWeightedEhrhartQuasiPolynomialToPyList(const libnormaliz::IntegrationData& int_data)
{
    vector< vector<Integer> > ehrhart_qp = int_data.getWeightedEhrhartQuasiPolynomial();
    const size_t n = ehrhart_qp.size();
    PyObject* return_list = PyList_New(n+1);
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(return_list, i, NmzVectorToPyList(ehrhart_qp[i]));
    }
    PyList_SetItem(return_list, n, NmzToPyNumber(int_data.getWeightedEhrhartQuasiPolynomialDenom()));
    return return_list;
}

template<typename Integer>
PyObject* NmzTriangleListToPyList(const vector< pair<vector<libnormaliz::key_t>, Integer> >& in)
{
    const size_t n = in.size();
    PyObject* M = PyList_New( n );
    for (size_t i = 0; i < n; ++i) {
        // convert the pair
        PyObject* pair = PyList_New(2);
        PyList_SetItem(pair, 0, NmzVectorToPyList<libnormaliz::key_t>(in[i].first));
        PyList_SetItem(pair, 1, NmzToPyNumber(in[i].second));
        PyList_SetItem(M, i, pair);
    }
    return M;
}

template<typename Integer>
PyObject* NmzStanleyDataToPyList(const libnormaliz::STANLEYDATA<Integer>& StanleyData)
{
    PyObject* pair = PyList_New(2);
    PyList_SetItem(pair, 0, NmzVectorToPyList<libnormaliz::key_t>(StanleyData.key));
    PyList_SetItem(pair, 1, NmzMatrixToPyList(StanleyData.offsets.get_elements()));
    return pair;
}

template<typename Integer>
PyObject* NmzStanleyDecToPyList(const list<libnormaliz::STANLEYDATA<Integer> >& StanleyDec)
{
    const size_t n = StanleyDec.size();
    PyObject* M = PyList_New( n );
    typename list<libnormaliz::STANLEYDATA<Integer> >::const_iterator S = StanleyDec.begin();
    for (size_t i = 0; i < n; ++i) {        
        PyList_SetItem(M, i,NmzStanleyDataToPyList(*S) );
        ++S;
    }
    return M;
}

template<typename Integer>
static PyObject* _NmzBasisChangeIntern(Cone<Integer>* C)
{
    Sublattice_Representation<Integer> bc = C->getSublattice();

    PyObject* res = PyList_New( 3 );
    PyList_SetItem(res, 0, NmzMatrixToPyList(bc.getEmbedding()));
    PyList_SetItem(res, 1, NmzMatrixToPyList(bc.getProjection()));
    PyList_SetItem(res, 2, NmzToPyNumber(bc.getAnnihilator()));
    // Dim, Rank, Equations and Congruences are already covered by special functions
    // ditto ExternalIndex
    return res;
}

/***************************************************************************
 * 
 * PyCapsule handler functions
 * 
 ***************************************************************************/

void delete_cone_mpz( PyObject* cone ){
  Cone<mpz_class> * cone_ptr = reinterpret_cast<Cone<mpz_class>* >( PyCapsule_GetPointer( cone, cone_name ) );
  delete cone_ptr;
}

void delete_cone_long( PyObject* cone ){
  Cone<long long> * cone_ptr = reinterpret_cast<Cone<long long>* >( PyCapsule_GetPointer( cone, cone_name_long ) );
  delete cone_ptr;
}

Cone<long long>* get_cone_long( PyObject* cone ){
  return reinterpret_cast<Cone<long long>*>( PyCapsule_GetPointer( cone, cone_name_long ) );
}

Cone<mpz_class>* get_cone_mpz( PyObject* cone ){
  return reinterpret_cast<Cone<mpz_class>*>( PyCapsule_GetPointer( cone, cone_name ) );
}


PyObject* pack_cone( Cone<mpz_class>* C ){
  return PyCapsule_New( reinterpret_cast<void*>( C ), cone_name, &delete_cone_mpz );
}

PyObject* pack_cone( Cone<long long>* C ){
  return PyCapsule_New( reinterpret_cast<void*>( C ), cone_name_long, &delete_cone_long );
}

bool is_cone( PyObject* cone ){
  if( PyCapsule_CheckExact( cone ) ){
    // compare as string
    return cone_name_str == string(PyCapsule_GetName( cone )) || cone_name_str_long == string(PyCapsule_GetName( cone ));
  }
  return false;
}

/***************************************************************************
 * 
 * Cone property list
 * 
 ***************************************************************************/

/*
@Name NmzListConeProperties
@Arguments none
@Description
Returns two lists of strings.
The first list are all cone properties that define compute
goals in Normaliz (see Normaliz manual for details)
The second list are all cone properties that define internal
control flow control in Normaliz, and which should not be used
to get results of computations.
All entries of the first list can be passed to NmzResult
to get the result of a normaliz computation.
All entries of the second list can be passed to NmzCompute
to set different options for Normaliz computations.
*/
static PyObject* NmzListConeProperties(PyObject* args)
{
    FUNC_BEGIN
    
    PyObject* return_list = PyList_New( 2 );
    
    ConeProperties props;
    for(int i=0; i < libnormaliz::ConeProperty::EnumSize;i++){
        props.set( static_cast<libnormaliz::ConeProperty::Enum>(i) );
    }
    
    ConeProperties goals = props.goals();
    ConeProperties options = props.options();
    
    int number_goals = goals.count();
    int number_options = options.count();
    
    PyObject* goal_list = PyList_New( number_goals );
    PyObject* option_list = PyList_New( number_options );

    PyList_SetItem( return_list, 0, goal_list );
    PyList_SetItem( return_list, 1, option_list );
    
    int list_position = 0;
    for(int i=0; i < libnormaliz::ConeProperty::EnumSize;i++){
      if(goals.test(static_cast<libnormaliz::ConeProperty::Enum>(i))){
          string name = libnormaliz::toString(static_cast<libnormaliz::ConeProperty::Enum>(i));
          PyList_SetItem( goal_list, list_position, StringToPyUnicode( name ) );
          list_position++;
      }
    }
    
    list_position = 0;
    for(int i=0; i < libnormaliz::ConeProperty::EnumSize;i++){
      if(options.test(static_cast<libnormaliz::ConeProperty::Enum>(i))){
          string name = libnormaliz::toString(static_cast<libnormaliz::ConeProperty::Enum>(i));
          PyList_SetItem( option_list, list_position, StringToPyUnicode( name ) );
          list_position++;
      }
    }
    
    return return_list;
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * NmzCone
 * 
 ***************************************************************************/

template<typename Integer>
static PyObject* _NmzConeIntern(PyObject * args, PyObject* kwargs)
{
    map <InputType, vector< vector<mpq_class> > > input;
    
    PyObject* input_list;
    
    bool grading_polynomial = false;
    string polynomial;
    
    if( PyTuple_Size(args)==1 ){
        input_list = PyTuple_GetItem( args, 0 );
        if( ! PyList_Check( input_list ) ){
            PyErr_SetString( PyNormaliz_cppError, "Single argument must be a list" );
            return NULL;
        }
        input_list = PyList_AsTuple( input_list );
    }else{
        input_list = args;
    }
    
    const int n = PyTuple_Size(input_list);
    if (n&1) {
        PyErr_SetString( PyNormaliz_cppError, "Number of arguments must be even" );
        return NULL;
    }
    for (int i = 0; i < n; i += 2) {
        PyObject* type = PyTuple_GetItem(input_list, i);
        if (!string_check(type)) {
            PyErr_SetString( PyNormaliz_cppError, "Odd entries must be strings" );
            return NULL;
        }
        
        string type_str = PyUnicodeToString( type );
        
        if( type_str.compare( "polynomial" ) == 0 ){
            PyObject* M = PyTuple_GetItem(input_list, i+1);
            polynomial =  PyUnicodeToString( M );
            grading_polynomial = true;
            continue;
        }
        
        PyObject* M = PyTuple_GetItem(input_list, i+1);
        if(M==Py_None)
            continue;
        vector<vector<mpq_class> > Mat;
        bool okay = PyInputToNmz(Mat, M);
        if (!okay) {
            PyErr_SetString( PyNormaliz_cppError, "Even entries must be matrices" );
            return NULL;
        }

        input[libnormaliz::to_type(type_str)] = Mat;
    }

    if(kwargs!=NULL){
        PyObject* keys = PyDict_Keys(kwargs);
        PyObject* values = PyDict_Values(kwargs);
        const int length = PyList_Size(keys);
        for(int i = 0; i<length; i++ ){
            string type_string = PyUnicodeToString( PyList_GetItem( keys, i ) );
            PyObject* current_value = PyList_GetItem( values, i );
            if(current_value==Py_None)
                continue;
            if( type_string.compare( "polynomial" ) == 0 ){
                polynomial = PyUnicodeToString( current_value );
                grading_polynomial = true;
                continue;
            }
            vector<vector<mpq_class> > Mat;
            bool okay = PyInputToNmz(Mat, current_value);
            if (!okay) {
                PyErr_SetString( PyNormaliz_cppError, "Even entries must be matrices" );
                return NULL;
            }
            input[libnormaliz::to_type(type_string)] = Mat;
        }
    }

    Cone<Integer>* C = new Cone<Integer>(input);
    
    if( grading_polynomial ){
        C->setPolynomial( polynomial );
    }
    
    PyObject* return_container = pack_cone( C );
    
    return return_container;
}

/*
@Name NmzCone
@Arguments <keywords>
@Description
Constructs a normaliz cone object. The keywords must be 
Normaliz input types, and the values for the keys matrices
(consisting of either Longs, Floats, or strings for rationals),
lists for single vector input types, or bools for boolean input type.
Special cases are a string describing a polynomial for the polynomial
input type, and the CreateAsLongLong keyword to restrict normaliz computations
to machine integers instead of arbitrary precision numbers.
*/
PyObject* _NmzCone(PyObject* self, PyObject* args, PyObject* kwargs)
{
    FUNC_BEGIN
    
    static const char* string_for_keyword_argument = "CreateAsLongLong";
    PyObject* create_as_long_long;
    
#if PY_MAJOR_VERSION >= 3
    PyObject* key = PyUnicode_FromString( string_for_keyword_argument );
#else
    PyObject* key = PyString_FromString( const_cast<char*>(string_for_keyword_argument) );
#endif
    
    if( kwargs != NULL && PyDict_Contains( kwargs, key ) == 1 ){
        create_as_long_long = PyDict_GetItem( kwargs, key );
        PyDict_DelItem( kwargs, key );
    }else{
        create_as_long_long = Py_False;
    }
    
    if( create_as_long_long!=Py_True ){
        return _NmzConeIntern<mpz_class>(args,kwargs);
    }else{
        return _NmzConeIntern<long long>(args,kwargs);
    }

    FUNC_END
}

/*
@Name NmzConeCopy
@Arguments Cone
@Description
Returns a copy of the cone.
*/
PyObject* _NmzConeCopy( PyObject* self, PyObject* args )
{
    FUNC_BEGIN
    PyObject* cone = PyTuple_GetItem( args, 0 );
    if( !is_cone(cone) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }

    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        Cone<mpz_class>* new_cone = new Cone<mpz_class>(*cone_ptr);
        return pack_cone(new_cone);
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        Cone<long long>* new_cone = new Cone<long long>(*cone_ptr);
        return pack_cone(new_cone);
    }
    FUNC_END
}

/***************************************************************************
 * 
 * NmzHilbertSeries
 * 
 ***************************************************************************/

template<typename Integer>
PyObject* NmzHilbertSeries(Cone<Integer>* C, PyObject* args)
{
    FUNC_BEGIN
    
    const int arg_len = PyTuple_Size(args);
    
    if(arg_len==1){
        bool is_HSOP = C->isComputed(libnormaliz::ConeProperty::HSOP);
        return NmzHilbertSeriesToPyList(C->getHilbertSeries(),is_HSOP);
    }
    
    PyObject* is_HSOP = PyTuple_GetItem( args, 1 );
    
    if( is_HSOP == Py_True ){
        if (!C->isComputed(libnormaliz::ConeProperty::HSOP)) C->compute(libnormaliz::ConeProperty::HSOP);
        return NmzHilbertSeriesToPyList(C->getHilbertSeries(),true);
    }else{
        return NmzHilbertSeriesToPyList(C->getHilbertSeries(),false);
    }
    
    FUNC_END
}


PyObject* NmzHilbertSeries_Outer(PyObject* self, PyObject* args){
  
  FUNC_BEGIN
  
  PyObject* cone = PyTuple_GetItem( args, 0 );
  
  if( !is_cone(cone) ){
      PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
      return NULL;
  }
  
  current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
  
  if( cone_name_str == string(PyCapsule_GetName(cone)) ){
      Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
      PyObject* return_value = NmzHilbertSeries(cone_ptr, args);
      PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
      return return_value;
  }else{
      Cone<long long>* cone_ptr = get_cone_long(cone);
      PyObject* return_value = NmzHilbertSeries(cone_ptr,args);
      PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
      return return_value;
  }
  
  FUNC_END
  
}

/***************************************************************************
 * 
 * NmzCompute
 * 
 ***************************************************************************/


template<typename Integer>
PyObject* _NmzCompute(Cone<Integer>* C, PyObject* args)
{
    FUNC_BEGIN
    
    const int arg_len = PyTuple_Size(args);
    
    PyObject* to_compute;
    
    if(arg_len==2){
        PyObject* first_arg = PyTuple_GetItem(args,1);
        if(PyList_CheckExact( first_arg )){
            to_compute = first_arg;
        }else{
            to_compute = PyList_New( 1 );
            int result = PyList_SetItem( to_compute, 0, first_arg );
            if(result!=0){
                PyErr_SetString( PyNormaliz_cppError, "List could not be created" );
                return NULL;
            }
        }
    }else{
        to_compute = PyList_New( arg_len - 1 );
        for( int i = 1;i<arg_len;i++){
            PyList_SetItem( to_compute, i, PyTuple_GetItem( args, i ) );
        }
    }

    ConeProperties propsToCompute;
    const int n = PyList_Size(to_compute);
    
    for (int i = 0; i < n; ++i) {
        PyObject* prop = PyList_GetItem(to_compute, i);
        if (!string_check(prop)) {
            PyErr_SetString( PyNormaliz_cppError, "All elements must be strings" );
            return NULL;
        }
        string prop_str(PyUnicodeToString(prop));
        propsToCompute.set( libnormaliz::toConeProperty(prop_str) );
    }
    
    ConeProperties notComputed = C->compute(propsToCompute);
    
    // Cone.compute returns the not computed properties
    // we return a bool, true when everything requested was computed
    return notComputed.none() ? Py_True : Py_False;
    FUNC_END
}


PyObject* _NmzCompute_Outer(PyObject* self, PyObject* args){
  
  FUNC_BEGIN
  
  current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
  
  PyObject* cone = PyTuple_GetItem( args, 0 );
  
  PyObject* result;
  
  if( !is_cone(cone) ){
      PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
      return NULL;
  }
  
  if( cone_name_str == string(PyCapsule_GetName(cone)) ){
      Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
      result = _NmzCompute(cone_ptr, args);
  }else{
      Cone<long long>* cone_ptr = get_cone_long(cone);
      result = _NmzCompute(cone_ptr,args);
  }
  
  PyOS_setsig(SIGINT,current_interpreter_sigint_handler);
  
  return result;
  
  FUNC_END
  
}

/***************************************************************************
 * 
 * NmzIsComputed
 * 
 ***************************************************************************/

/*
@Name NmzIsComputed
@Arguments <cone>, <property_string>
@Desctiption
Returns if the cone property <property_string> is computed in the cone <cone>.
*/
template<typename Integer>
PyObject* NmzIsComputed(Cone<Integer>* C, PyObject* prop)
{
    FUNC_BEGIN
    
    libnormaliz::ConeProperty::Enum p = libnormaliz::toConeProperty(PyUnicodeToString( prop ) );
  
    return C->isComputed(p) ? Py_True : Py_False;

    FUNC_END
}

PyObject* NmzIsComputed_Outer(PyObject* self, PyObject* args)
{
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    PyObject* to_compute = PyTuple_GetItem( args, 1 );
    
    if( !is_cone(cone) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        return NmzIsComputed(cone_ptr, to_compute);
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        return NmzIsComputed(cone_ptr,to_compute);
    }
    
    FUNC_END
}

/***************************************************************************
 * 
 * NmzResult
 * 
 ***************************************************************************/

/*
@Name NmzResult
@Arguments <cone>,<cone property string>,<keys>
@Description
Returns the cone property belonging to the string <cone property string> of
cone <cone>. Please see the Normaliz manual for details on which cone properties are available.
Here are some special outputs that might differ from Normaliz:
* HilbertSeries and WeightedEhrhartSeries
  The returned object is a list with three entries: The first one describes the
  numerator of the hilbert series, the second one the denominator, and the last one
  is the shift. If you pass the HSOP option, output will be done in HSOP format.
* Grading
  Returns a list with two entries. First is the grading, second one is the grading denominator.
* Sublattice
  Returns a list with three entries. First is the embedding of the sublattice, second is the projection
  third is the annihilator.
* IntegerHull and ProjectCone return new cones.
* StanleyDec
  Returns a list containing the Stanley decomposition. All entries are 2-tuples. First entry in the tuple is the
  key, second the decomposition data.
*/
template<typename Integer>
PyObject* _NmzResultImpl(Cone<Integer>* C, PyObject* prop_obj)
{
    
    string prop = PyUnicodeToString( prop_obj );

    libnormaliz::ConeProperty::Enum p = libnormaliz::toConeProperty(prop);
    
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    ConeProperties notComputed = C->compute(ConeProperties(p));
    PyOS_setsig(SIGINT, current_interpreter_sigint_handler );
    
    if (notComputed.any()) {
        return Py_None;
    }

    switch (p) {
    case libnormaliz::ConeProperty::Generators:
        return NmzMatrixToPyList(C->getGenerators());

    case libnormaliz::ConeProperty::ExtremeRays:
        return NmzMatrixToPyList(C->getExtremeRays());

    case libnormaliz::ConeProperty::VerticesOfPolyhedron:
        return NmzMatrixToPyList(C->getVerticesOfPolyhedron());

    case libnormaliz::ConeProperty::SupportHyperplanes:
        return NmzMatrixToPyList(C->getSupportHyperplanes());

    case libnormaliz::ConeProperty::SuppHypsFloat:
        return NmzMatrixToPyList(C->getSuppHypsFloat());

    case libnormaliz::ConeProperty::TriangulationSize:
        return NmzToPyNumber(C->getTriangulationSize());

    case libnormaliz::ConeProperty::TriangulationDetSum:
        return NmzToPyNumber(C->getTriangulationDetSum());

    case libnormaliz::ConeProperty::Triangulation:
        return NmzTriangleListToPyList<Integer>(C->getTriangulation());

    case libnormaliz::ConeProperty::Multiplicity:
        return NmzToPyList(C->getMultiplicity());
    
    case libnormaliz::ConeProperty::Integral:
        return NmzToPyList(C->getIntegral());
    
    case libnormaliz::ConeProperty::VirtualMultiplicity:
        return NmzToPyList(C->getVirtualMultiplicity());

    case libnormaliz::ConeProperty::RecessionRank:
        return NmzToPyNumber(C->getRecessionRank());

    case libnormaliz::ConeProperty::AffineDim:
        return NmzToPyNumber(C->getAffineDim());

    case libnormaliz::ConeProperty::ModuleRank:
        return NmzToPyNumber(C->getModuleRank());

    case libnormaliz::ConeProperty::HilbertBasis:
        return NmzMatrixToPyList(C->getHilbertBasis());

    case libnormaliz::ConeProperty::MaximalSubspace:
        return NmzMatrixToPyList(C->getMaximalSubspace());

    case libnormaliz::ConeProperty::ModuleGenerators:
        return NmzMatrixToPyList(C->getModuleGenerators());

    case libnormaliz::ConeProperty::Deg1Elements:
        return NmzMatrixToPyList(C->getDeg1Elements());

    case libnormaliz::ConeProperty::HilbertSeries:
    case libnormaliz::ConeProperty::EhrhartSeries:
        {
        bool is_HSOP = C->isComputed(libnormaliz::ConeProperty::HSOP);
        return NmzHilbertSeriesToPyList(C->getHilbertSeries(),is_HSOP);
        }

    case libnormaliz::ConeProperty::WeightedEhrhartSeries:
        return NmzWeightedEhrhartSeriesToPyList(C->getWeightedEhrhartSeries());
        

    case libnormaliz::ConeProperty::Grading:
        {
        vector<Integer> grad = C->getGrading();
        Integer denom = C->getGradingDenom();
        PyObject * return_list = PyList_New(2);
        PyList_SetItem( return_list, 0, NmzVectorToPyList(grad) );
        PyList_SetItem( return_list, 1, NmzToPyNumber( denom ) );
        return return_list;
        }

    case libnormaliz::ConeProperty::IsPointed:
        return BoolToPyBool(C->isPointed());

    case libnormaliz::ConeProperty::IsDeg1ExtremeRays:
        return BoolToPyBool(C->isDeg1ExtremeRays());

    case libnormaliz::ConeProperty::IsDeg1HilbertBasis:
        return BoolToPyBool(C->isDeg1HilbertBasis());

    case libnormaliz::ConeProperty::IsIntegrallyClosed:
        return BoolToPyBool(C->isIntegrallyClosed());

    case libnormaliz::ConeProperty::OriginalMonoidGenerators:
        return NmzMatrixToPyList(C->getOriginalMonoidGenerators());

    case libnormaliz::ConeProperty::IsReesPrimary:
        return BoolToPyBool(C->isReesPrimary());

    case libnormaliz::ConeProperty::ReesPrimaryMultiplicity:
        return NmzToPyNumber(C->getReesPrimaryMultiplicity());

    case libnormaliz::ConeProperty::StanleyDec:
        return NmzStanleyDecToPyList(C->getStanleyDec());

    case libnormaliz::ConeProperty::ExcludedFaces:
        return NmzMatrixToPyList(C->getExcludedFaces());

    case libnormaliz::ConeProperty::Dehomogenization:
        return NmzVectorToPyList(C->getDehomogenization());

    case libnormaliz::ConeProperty::InclusionExclusionData:
        return NmzTriangleListToPyList<long>(C->getInclusionExclusionData());

    case libnormaliz::ConeProperty::ClassGroup:
        return NmzVectorToPyList(C->getClassGroup());
    
    case libnormaliz::ConeProperty::IsInhomogeneous:
        return BoolToPyBool(C->isInhomogeneous());
    
    /* Sublattice properties */
    
    case libnormaliz::ConeProperty::Equations:
        return NmzMatrixToPyList(C->getSublattice().getEquations());
    
    case libnormaliz::ConeProperty::Congruences:
        return NmzMatrixToPyList(C->getSublattice().getCongruences());
    
    case libnormaliz::ConeProperty::EmbeddingDim:
        return NmzToPyNumber(C->getEmbeddingDim());
    
    case libnormaliz::ConeProperty::Rank:
        return NmzToPyNumber(C->getRank());
    
    case libnormaliz::ConeProperty::Sublattice:
        return _NmzBasisChangeIntern(C);
    
    case libnormaliz::ConeProperty::ExternalIndex:
        return NmzToPyNumber(C->getSublattice().getExternalIndex());
    
    case libnormaliz::ConeProperty::InternalIndex:
        return NmzToPyNumber(C->getIndex());
    
    case libnormaliz::ConeProperty::WitnessNotIntegrallyClosed:
        return NmzVectorToPyList(C->getWitnessNotIntegrallyClosed());
    
    
    /* New stuff */
    
    case libnormaliz::ConeProperty::GradingDenom:
        return NmzToPyNumber(C->getGradingDenom());
    
    case libnormaliz::ConeProperty::UnitGroupIndex:
        return NmzToPyNumber(C->getUnitGroupIndex());
    
    case libnormaliz::ConeProperty::ModuleGeneratorsOverOriginalMonoid:
        return NmzMatrixToPyList(C->getModuleGeneratorsOverOriginalMonoid());
    
    case libnormaliz::ConeProperty::IntegerHull:
    {
        Cone<Integer>* hull = new Cone<Integer>( C->getIntegerHullCone() );
        return pack_cone( hull ); 
    }
    
    case libnormaliz::ConeProperty::ProjectCone:
    {
        Cone<Integer>* projection = new Cone<Integer>(C->getProjectCone());
        return pack_cone( projection );
    }
    
    case libnormaliz::ConeProperty::HilbertQuasiPolynomial:
        return NmzHilbertQuasiPolynomialToPyList<mpz_class>(C->getHilbertSeries()); //FIXME: Why is this return value not parametrized, but mpz_class only?
    
    case libnormaliz::ConeProperty::WeightedEhrhartQuasiPolynomial:
        return NmzWeightedEhrhartQuasiPolynomialToPyList<mpz_class>(C->getIntData());
    
    case libnormaliz::ConeProperty::IsTriangulationNested:
        return BoolToPyBool(C->isTriangulationNested());
        
    case libnormaliz::ConeProperty::IsTriangulationPartial:
        return BoolToPyBool(C->isTriangulationPartial());
        
    case libnormaliz::ConeProperty::ConeDecomposition:
        return NmzBoolMatrixToPyList(C->getOpenFacets());
    
    case libnormaliz::ConeProperty::IsGorenstein:
        return BoolToPyBool(C->isGorenstein());
        
    case libnormaliz::ConeProperty::GeneratorOfInterior:
        return NmzVectorToPyList(C->getGeneratorOfInterior());
        
    case libnormaliz::ConeProperty::VerticesFloat:
        return NmzMatrixToPyList(C->getVerticesFloat());
        
    case libnormaliz::ConeProperty::Volume:
        return NmzToPyList(C->getVolume());

    case libnormaliz::ConeProperty::EuclideanVolume:
        return NmzToPyNumber(C->getEuclideanVolume());

//  the following properties are compute options and do not return anything
    case libnormaliz::ConeProperty::DualMode:
    case libnormaliz::ConeProperty::DefaultMode:
    case libnormaliz::ConeProperty::Approximate:
    case libnormaliz::ConeProperty::BottomDecomposition:
    case libnormaliz::ConeProperty::KeepOrder:
    case libnormaliz::ConeProperty::NoBottomDec:
    case libnormaliz::ConeProperty::PrimalMode:
    case libnormaliz::ConeProperty::Symmetrize:
    case libnormaliz::ConeProperty::NoSymmetrization:
    case libnormaliz::ConeProperty::BigInt:
    case libnormaliz::ConeProperty::NoNestedTri:
    case libnormaliz::ConeProperty::HSOP:
    case libnormaliz::ConeProperty::Projection:
    case libnormaliz::ConeProperty::NoProjection:
    case libnormaliz::ConeProperty::ProjectionFloat:
    case libnormaliz::ConeProperty::SCIP:
    case libnormaliz::ConeProperty::NoPeriodBound:
    case libnormaliz::ConeProperty::NoLLL:
    case libnormaliz::ConeProperty::NoRelax:
        PyErr_SetString( PyNormaliz_cppError, "ConeProperty is input-only" );
        return NULL;
#if NMZ_RELEASE >= 30200
    case libnormaliz::ConeProperty::NoSubdivision:
        PyErr_SetString( PyNormaliz_cppError, "ConeProperty is input-only" );
        return NULL;
#endif
    default:
        PyErr_SetString( PyNormaliz_cppError, "Unknown cone property" );
        return NULL;
        break;
    }

    return Py_None;
}

PyObject* _NmzResult( PyObject* self, PyObject* args, PyObject* kwargs ){
    
    FUNC_BEGIN

    PyObject* cone = PyTuple_GetItem( args, 0 );
    PyObject* prop = PyTuple_GetItem( args, 1 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    if( !string_check( prop ) ){
        PyErr_SetString( PyNormaliz_cppError, "Second argument must be a unicode string" );
        return NULL;
    }
    
    if(kwargs){
        RationalHandler = PyDict_GetItemString(kwargs,"RationalHandler");
        VectorHandler = PyDict_GetItemString(kwargs,"VectorHandler");
        MatrixHandler = PyDict_GetItemString(kwargs,"MatrixHandler");
    }

    PyObject* result;

    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        result = _NmzResultImpl(cone_ptr, prop);
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        result = _NmzResultImpl(cone_ptr, prop);
    }

    RationalHandler = NULL;
    VectorHandler = NULL;
    MatrixHandler = NULL;

    return result;


    
    FUNC_END
}

/***************************************************************************
 * 
 * Python verbosity
 * 
 ***************************************************************************/

PyObject* NmzSetVerboseDefault( PyObject* self, PyObject* args)
{
    FUNC_BEGIN
    PyObject * value = PyTuple_GetItem( args, 0 );
    if (value != Py_True && value != Py_False){
        PyErr_SetString( PyNormaliz_cppError, "Argument must be True or False" );
        return NULL;
    }
    return BoolToPyBool(libnormaliz::setVerboseDefault(value == Py_True));
    FUNC_END
}

template<typename Integer>
PyObject* NmzSetVerbose(Cone<Integer>* C, PyObject* value)
{
    FUNC_BEGIN
    bool old_value;
    old_value = C->setVerbose(value == Py_True);
    return BoolToPyBool(old_value);
    FUNC_END
}

PyObject* NmzSetVerbose_Outer(PyObject* self, PyObject* args)
{
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    PyObject* value = PyTuple_GetItem( args, 1 );
    if (value != Py_True && value != Py_False){
        PyErr_SetString( PyNormaliz_cppError, "Second argument must be True or False" );
        return NULL;
    }
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        return NmzSetVerbose(cone_ptr, value);
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        return NmzSetVerbose(cone_ptr, value);
    }
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * Get Polynomial
 * 
 ***************************************************************************/

PyObject* NmzGetPolynomial(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        PyObject* return_value = StringToPyUnicode( (cone_ptr->getIntData()).getPolynomial() );
        PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
        return return_value;
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        PyObject* return_value = StringToPyUnicode( (cone_ptr->getIntData()).getPolynomial() );
        PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
        return return_value;
    }
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * NrCoeffQuasiPol
 * 
 ***************************************************************************/

PyObject* NmzSetNrCoeffQuasiPol( PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    PyObject* bound_py = PyTuple_GetItem( args, 1 );
    
    int overflow;
    long bound = PyLong_AsLongLongAndOverflow( bound_py, &overflow );
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        cone_ptr->setNrCoeffQuasiPol(bound);
        return Py_True;
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        cone_ptr->setNrCoeffQuasiPol(bound);
        return Py_True;
    }
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * Get Symmetrized cone
 * 
 ***************************************************************************/

PyObject* NmzSymmetrizedCone(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        Cone<mpz_class>* symm_cone = &(cone_ptr->getSymmetrizedCone());
        PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
        if( symm_cone==0 ){
            return Py_None;
        }
        symm_cone = new Cone<mpz_class>( *symm_cone );
        return pack_cone( symm_cone );
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        Cone<long long>* symm_cone = &(cone_ptr->getSymmetrizedCone());
        PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
        if( symm_cone==0 ){
            return Py_None;
        }
        symm_cone = new Cone<long long>( *symm_cone );
        return pack_cone( symm_cone );
    }
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * Get euclidian volume
 * 
 ***************************************************************************/

PyObject* NmzGetEuclideanVolume(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        PyObject* return_value = NmzToPyNumber(cone_ptr->getEuclideanVolume());
        PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
        return return_value;
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        PyObject* return_value = NmzToPyNumber(cone_ptr->getEuclideanVolume());
        PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
        return return_value;
    }
    PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * Get expanded hilbert series
 * 
 ***************************************************************************/

PyObject* NmzGetHilbertSeriesExpansion(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    PyObject* py_degree = PyTuple_GetItem( args, 1 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    if( !PyLong_Check( py_degree ) ){
        PyErr_SetString( PyNormaliz_cppError, "Second argument must be a long" );
        return NULL;
    }
    
    long degree = PyLong_AsLong( py_degree );
    libnormaliz::HilbertSeries HS;
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        HS = cone_ptr->getHilbertSeries();
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        HS = cone_ptr->getHilbertSeries();
    }
    
    HS.set_expansion_degree(degree);
    PyObject* result = NmzVectorToPyList( HS.getExpansion() );
    PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
    
    return result;
    
    FUNC_END
    
}


PyObject* NmzGetWeightedEhrhartSeriesExpansion(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    PyObject* py_degree = PyTuple_GetItem( args, 1 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    if( !PyLong_Check( py_degree ) ){
        PyErr_SetString( PyNormaliz_cppError, "Second argument must be a long" );
        return NULL;
    }
    
    long degree = PyLong_AsLong( py_degree );
    pair<libnormaliz::HilbertSeries,mpz_class> ES;
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    
    if( cone_name_str == string(PyCapsule_GetName(cone)) ){
        Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
        ES = cone_ptr->getWeightedEhrhartSeries();
    }else{
        Cone<long long>* cone_ptr = get_cone_long(cone);
        ES = cone_ptr->getWeightedEhrhartSeries();
    }
    
    ES.first.set_expansion_degree(degree);
    PyObject* result = NmzVectorToPyList( ES.first.getExpansion() );
    PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
    
    return result;
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * Set number of threads
 * 
 ***************************************************************************/

PyObject* NmzSetNumberOfNormalizThreads(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* num_treads = PyTuple_GetItem( args, 0 );
    
    if( !PyLong_Check( num_treads ) ){
        PyErr_SetString( PyNormaliz_cppError, "First argument must be an integer" );
        return NULL;
    }
    
    long num_threads_long = PyLong_AsLong( num_treads );
    
    num_threads_long = libnormaliz::set_thread_limit( num_threads_long );
    
    return PyLong_FromLong( num_threads_long );
    
    FUNC_END
    
}

/***************************************************************************
 * 
 * Python init stuff
 * 
 ***************************************************************************/

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject * error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef PyNormaliz_cppMethods[] = {
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {"NmzCone",  (PyCFunction)_NmzCone, METH_VARARGS|METH_KEYWORDS,
     "Create a cone"},
    {"NmzConeCopy",  (PyCFunction)_NmzConeCopy,METH_VARARGS,
     "Copy an existing cone" },
    {"NmzCompute", (PyCFunction)_NmzCompute_Outer, METH_VARARGS,
     "Compute some stuff"},
    {"NmzIsComputed", (PyCFunction)NmzIsComputed_Outer, METH_VARARGS,
     "Check if property is computed "},
    {"NmzResult", (PyCFunction)_NmzResult, METH_VARARGS|METH_KEYWORDS,
      "Return cone property" },
    { "NmzSetVerboseDefault", (PyCFunction)NmzSetVerboseDefault, METH_VARARGS,
      "Set verbosity" },
    { "NmzSetVerbose", (PyCFunction)NmzSetVerbose_Outer, METH_VARARGS,
      "Set verbosity of cone" },
    { "NmzListConeProperties", (PyCFunction)NmzListConeProperties,METH_NOARGS,
      "List all available properties" },
    { "NmzHilbertSeries", (PyCFunction)NmzHilbertSeries_Outer, METH_VARARGS,
      "Returns Hilbert series, either HSOP or not" },
    { "NmzGetPolynomial", (PyCFunction)NmzGetPolynomial, METH_VARARGS,
      "Returns grading polynomial" },
    { "NmzSymmetrizedCone", (PyCFunction)NmzSymmetrizedCone, METH_VARARGS,
      "Returns symmetrized cone" },
    { "NmzSetNumberOfNormalizThreads", (PyCFunction)NmzSetNumberOfNormalizThreads, METH_VARARGS,
      "Sets the Normaliz thread limit" },
    { "NmzSetNrCoeffQuasiPol", (PyCFunction)NmzSetNrCoeffQuasiPol, METH_VARARGS,
      "Sets the period bound for the quasi-polynomial" },
    { "NmzGetEuclideanVolume", (PyCFunction)NmzGetEuclideanVolume, METH_VARARGS,
      "Returns euclidean volume of cone as float" },
    { "NmzGetHilbertSeriesExpansion", (PyCFunction)NmzGetHilbertSeriesExpansion, METH_VARARGS,
      "Returns expansion of the hilbert series" },
    { "NmzGetWeightedEhrhartSeriesExpansion", (PyCFunction)NmzGetWeightedEhrhartSeriesExpansion, METH_VARARGS,
      "Returns expansion of the weighted Ehrhart series" },
    {NULL, }        /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3

static int PyNormaliz_cpp_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int PyNormaliz_cpp_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "PyNormaliz_cpp",
        NULL,
        sizeof(struct module_state),
        PyNormaliz_cppMethods,
        NULL,
        PyNormaliz_cpp_traverse,
        PyNormaliz_cpp_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC PyInit_PyNormaliz_cpp(void)

#else
#define INITERROR return

extern "C" void initPyNormaliz_cpp(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("PyNormaliz_cpp", PyNormaliz_cppMethods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException(const_cast<char*>("PyNormaliz_cpp.INITError"), NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    
    NormalizError = PyErr_NewException(const_cast<char*>("Normaliz.error"), NULL, NULL );
    Py_INCREF( NormalizError );
    PyNormaliz_cppError = PyErr_NewException(const_cast<char*>("Normaliz.interface_error"), NULL, NULL );
    Py_INCREF( PyNormaliz_cppError );
    
    PyModule_AddObject( module, "error", NormalizError );
    PyModule_AddObject( module, "error", PyNormaliz_cppError );
    
    current_interpreter_sigint_handler = PyOS_getsig( SIGINT );

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

