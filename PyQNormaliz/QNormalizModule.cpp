/***************************************************************************
 * 
 * Include
 * 
 ***************************************************************************/

#include <Python.h>
using namespace std;

#include <e-antic/renfxx.h>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <string>
using std::string;

#include <sstream>
using std::istringstream;

#include <libQnormaliz/Qcone.h>
#include <libQnormaliz/Qmap_operations.h>



using libQnormaliz::Cone;
//using libQnormaliz::QConeProperty;
using libQnormaliz::ConeProperties;
using libQnormaliz::Sublattice_Representation;
using libQnormaliz::QType::InputType;

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
    } catch (libQnormaliz::InterruptException& e ) {\
        PyOS_setsig(SIGINT,current_interpreter_sigint_handler);\
        libQnormaliz::nmz_interrupted = false; \
        PyErr_SetInterrupt(); \
        PyErr_CheckSignals(); \
        return NULL; \
    } catch (libQnormaliz::NormalizException& e) { \
        PyOS_setsig(SIGINT,current_interpreter_sigint_handler);\
        PyErr_SetString( QNormalizError, e.what() ); \
        return NULL; \
    } catch( exception& e ) { \
        PyOS_setsig(SIGINT,current_interpreter_sigint_handler);\
        PyErr_SetString( PyQNormaliz_cppError, e.what() ); \
        return NULL; \
    }

/***************************************************************************
 * 
 * Signal handling
 * 
 ***************************************************************************/

void signal_handler( int signal ){
    libQnormaliz::nmz_interrupted = true;
}

/***************************************************************************
 * 
 * Static objects
 * 
 ***************************************************************************/


static PyObject * QNormalizError;
static PyObject * PyQNormaliz_cppError;
static const char* cone_name = "Cone";
static string cone_name_str( cone_name );

static PyOS_sighandler_t current_interpreter_sigint_handler;

static PyObject * RationalHandler = NULL;
static PyObject * NumberfieldElementHandler = NULL;
static PyObject * VectorHandler = NULL;
static PyObject * MatrixHandler = NULL;


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

// #ifndef NMZ_RELEASE
//     static_assert(false,
//        "Your Normaliz version (unknown) is to old! Update to 3.5.2 or newer.");
// #endif
// #if NMZ_RELEASE < 30502
//     static_assert(false, "Your Normaliz version is to old! Update to 3.5.2 or newer.");
// #endif

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

// Do conversion for number field poly elements
bool PyPolyStringToNmz( PyObject* in, vector<mpq_class>& out ){
    string c_string = PyUnicodeToString(in);
    vector<mpq_class> temp = poly_components(c_string);
    out.swap(temp);
    return true;
}

// Converting MPZ's to PyLong and back via strings. Worst possible solution ever.
bool PyNumberToNmz( PyObject *, mpz_class& );

bool PyNumberToNmz( PyObject * in, mpq_class& out ){
    if( PyFloat_Check( in ) ){
        mpq_class temp(PyFloat_AsDouble(in));
        out.swap(temp);
        return true;
    }
    if( PyList_Check( in ) ){
        PyObject* py_num = PyList_GetItem(in,0);
        PyObject* py_denom = PyList_GetItem(in,1);
        mpz_class num;
        if(!PyNumberToNmz(py_num,num)){
            return false;
        }
        mpz_class denom;
        if(!PyNumberToNmz(py_denom,denom)){
            return false;
        }
        mpq_class temp(num,denom);
        out.swap(temp);
        return true;
    }
    PyObject * in_as_string = PyObject_Str( in );
    const char* in_as_c_string = PyUnicodeToString( in_as_string ).c_str();
    out.set_str( in_as_c_string, 10 );
    return true;
}

bool PyNumberToNmz( PyObject * in, mpz_class& out ){
  if( !PyLong_Check( in ) ){
      cerr << "here1" << endl;
      return false;
  }
  int overflow;
  long input_long = PyLong_AsLongAndOverflow( in, &overflow );
  if( overflow != -1 ){
      mpz_class temp(input_long);
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

PyObject* NmzToPyNumber( mpq_class in ){
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

PyObject* NmzToPyNumber( libQnormaliz::key_t in ){
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
    if (!PyList_Check(in)){
            return false;
    }
    const int n = PyList_Size(in);
    out.resize(n);
    for (int i = 0; i < n; ++i) {
        PyObject* tmp = PyList_GetItem(in, i);
        if (!PyNumberToNmz(tmp, out[i])){
            return false;
        }
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

template<typename NumberField, typename NumberFieldElem>
bool prepare_nf_input( vector< vector<NumberFieldElem> >& out, PyObject* in, NumberField* nf ){
    if (!PyList_Check( in ) )
        return false;
    const int nr = PyList_Size( in );
    out.resize(nr);
    for (int i = 0; i < nr; ++i) {
        PyObject * current_row = PyList_GetItem(in,i);
        int current_length = PyList_Size(current_row);
        for(int j=0; j < current_length;j++){
            PyObject * current_element = PyList_GetItem(current_row,j);
            vector<mpq_class> current_vector;
            bool current_res;
            if(PyList_Check(current_element)){
                current_res = PyListToNmz(current_vector,current_element);
            }
            else if(string_check(current_element)) {
                current_res = PyPolyStringToNmz(current_element,current_vector);
            } else {
                mpq_class temp;
                current_res = PyNumberToNmz(current_element,temp);
                current_vector.push_back(temp);
            }
            if(!current_res)
                return false;
            fmpq_poly_t current_poly;
            fmpq_poly_init(current_poly);
            vector2fmpq_poly(current_poly,current_vector);
            NumberFieldElem current_elem(nf->get_renf());
            current_elem = current_poly;
            out[i].push_back(current_elem);
        }
    }
    return true;
}

template<typename Integer>
PyObject* NmzVectorToPyList(const vector<Integer>& in, bool do_callback=true)
{
    PyObject* vector;
    const size_t n = in.size();
    vector = PyList_New(n);
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(vector, i, NmzToPyNumber(in[i]));
    }
    if(do_callback&&VectorHandler!=NULL)
        vector = CallPythonFuncOnOneArg(VectorHandler,vector);
    return vector;
}

PyObject* NmzToPyNumber( renf_elem_class in ){
    fmpq_poly_t current;
    fmpq_poly_init(current);
    in.get_fmpq_poly(current);
    vector<mpq_class> output;
    fmpq_poly2vector(output,current);
    PyObject * out_list = NmzVectorToPyList(output,false);
    if(NumberfieldElementHandler!=NULL)
        out_list = CallPythonFuncOnOneArg(NumberfieldElementHandler,out_list);
    return out_list;
}

PyObject* NmzBoolVectorToPyList(const vector<bool>& in)
{
    PyObject* vector;
    const size_t n = in.size();
    vector = PyList_New(n);
    for (size_t i = 0; i < n; ++i) {
        PyList_SetItem(vector, i, BoolToPyBool(in[i]));
    }
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

template<typename Integer>
PyObject* NmzTriangleListToPyList(const vector< pair<vector<libQnormaliz::key_t>, Integer> >& in)
{
    const size_t n = in.size();
    PyObject* M = PyList_New( n );
    for (size_t i = 0; i < n; ++i) {
        // convert the pair
        PyObject* pair = PyList_New(2);
        PyList_SetItem(pair, 0, NmzVectorToPyList<libQnormaliz::key_t>(in[i].first));
        PyList_SetItem(pair, 1, NmzToPyNumber(in[i].second));
        PyList_SetItem(M, i, pair);
    }
    return M;
}

// template<typename Integer>
// static PyObject* _NmzBasisChangeIntern(Cone<Integer>* C)
// {
//     Sublattice_Representation<Integer> bc = C->getSublattice();

//     PyObject* res = PyList_New( 3 );
//     PyList_SetItem(res, 0, NmzMatrixToPyList(bc.getEmbedding()));
//     PyList_SetItem(res, 1, NmzMatrixToPyList(bc.getProjection()));
//     PyList_SetItem(res, 2, NmzToPyNumber(bc.getAnnihilator()));
//     // Dim, Rank, Equations and Congruences are already covered by special functions
//     // ditto ExternalIndex
//     return res;
// }

/***************************************************************************
 * 
 * PyCapsule handler functions
 * 
 ***************************************************************************/

struct NumberFieldCone{
    renf_class* nf;
    Cone<renf_elem_class>* cone;
};

void delete_cone_renf( PyObject* cone ){
  NumberFieldCone * cone_ptr = reinterpret_cast<NumberFieldCone*>( PyCapsule_GetPointer( cone, cone_name ) );
  delete cone_ptr->cone;
  delete cone_ptr->nf;
}

Cone<renf_elem_class>* get_cone_renf( PyObject* cone ){
  NumberFieldCone* cone_ptr = reinterpret_cast<NumberFieldCone*>( PyCapsule_GetPointer( cone, cone_name ) );
  return cone_ptr->cone;
}

PyObject* pack_cone( Cone<renf_elem_class>* C, renf_class* nf ){
  NumberFieldCone* cone_ptr = new NumberFieldCone();
  cone_ptr->nf = nf;
  cone_ptr->cone = C;
  return PyCapsule_New( reinterpret_cast<void*>( cone_ptr ), cone_name, &delete_cone_renf );
}

bool is_cone( PyObject* cone ){
  if( PyCapsule_CheckExact( cone ) ){
    // compare as string
    return cone_name_str == string(PyCapsule_GetName( cone ));
  }
  return false;
}

/***************************************************************************
 * 
 * Cone property list
 * 
 ***************************************************************************/

static PyObject* NmzListConeProperties(PyObject* args)
{
    FUNC_BEGIN
    
    PyObject* return_list = PyList_New( 2 );
    
    ConeProperties props;
    for(int i=0; i < libQnormaliz::QConeProperty::EnumSize;i++){
        props.set( static_cast<libQnormaliz::QConeProperty::Enum>(i) );
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
    for(int i=0; i < libQnormaliz::QConeProperty::EnumSize;i++){
      if(goals.test(static_cast<libQnormaliz::QConeProperty::Enum>(i))){
          string name = libQnormaliz::toString(static_cast<libQnormaliz::QConeProperty::Enum>(i));
          PyList_SetItem( goal_list, list_position, StringToPyUnicode( name ) );
          list_position++;
      }
    }
    
    list_position = 0;
    for(int i=0; i < libQnormaliz::QConeProperty::EnumSize;i++){
      if(options.test(static_cast<libQnormaliz::QConeProperty::Enum>(i))){
          string name = libQnormaliz::toString(static_cast<libQnormaliz::QConeProperty::Enum>(i));
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

template<typename NumberField, typename NumberFieldElem>
static PyObject* _NmzQCone_internal(PyObject * args, PyObject* kwargs)
{

    NumberField * renf = new NumberField();

    PyObject* number_field_data = PyDict_GetItemString(kwargs,"number_field");
    istringstream number_field_data_stream(PyUnicodeToString(number_field_data));
    number_field_data_stream >> *renf;
    // Error handling
    // number_field_data_stream >> set_renf(renf);
    // Error handling

    map <InputType, vector< vector<NumberFieldElem> > > input;
    
    PyDict_DelItemString(kwargs,"number_field");
    if(kwargs!=NULL){
        PyObject* keys = PyDict_Keys(kwargs);
        PyObject* values = PyDict_Values(kwargs);
        const int length = PyList_Size(keys);
        for(int i = 0; i<length; i++ ){
            string type_string = PyUnicodeToString( PyList_GetItem( keys, i ) );
            PyObject* current_value = PyList_GetItem( values, i );
            if(current_value==Py_None)
                continue;
            vector<vector<NumberFieldElem> > Mat;
            bool okay = prepare_nf_input(Mat, current_value,renf);
            if (!okay) {
                PyErr_SetString( PyQNormaliz_cppError, "Could not read matrix" );
                return NULL;
            }
            input[libQnormaliz::to_type(type_string)] = Mat;
        }
    }

    Cone<NumberFieldElem>* C = new Cone<NumberFieldElem>(input);
    
    PyObject* return_container = pack_cone( C, renf );
    
    return return_container;
}

PyObject* _NmzQCone( PyObject* self, PyObject* args, PyObject* kwargs){
    return _NmzQCone_internal<renf_class,renf_elem_class>(args,kwargs);
}

// /***************************************************************************
//  * 
//  * NmzCompute
//  * 
//  ***************************************************************************/

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
                PyErr_SetString( PyQNormaliz_cppError, "List could not be created" );
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
            PyErr_SetString( PyQNormaliz_cppError, "All elements must be strings" );
            return NULL;
        }
        string prop_str(PyUnicodeToString(prop));
        propsToCompute.set( libQnormaliz::toConeProperty(prop_str) );
    }

    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    ConeProperties notComputed = C->compute(propsToCompute);
    PyOS_setsig(SIGINT,current_interpreter_sigint_handler);

    // Cone.compute returns the not computed properties
    // we return a bool, true when everything requested was computed
    return notComputed.none() ? Py_True : Py_False;
    FUNC_END
}


PyObject* _NmzCompute_Outer(PyObject* self, PyObject* args){
  
  FUNC_BEGIN
  
  PyObject* cone = PyTuple_GetItem( args, 0 );
  
  PyObject* result;
  
  if( !is_cone(cone) ){
      PyErr_SetString( PyQNormaliz_cppError, "First argument must be a cone" );
      return NULL;
  }
  
  Cone<renf_elem_class>* cone_ptr = get_cone_renf(cone);
  result = _NmzCompute(cone_ptr, args);
  
  return result;
  
  FUNC_END
  
}

/***************************************************************************
 * 
 * NmzIsComputed
 * 
 ***************************************************************************/

template<typename Integer>
PyObject* NmzIsComputed(Cone<Integer>* C, PyObject* prop)
{
    FUNC_BEGIN
    
    libQnormaliz::QConeProperty::Enum p = libQnormaliz::toConeProperty(PyUnicodeToString( prop ) );
  
    return C->isComputed(p) ? Py_True : Py_False;

    FUNC_END
}

PyObject* NmzIsComputed_Outer(PyObject* self, PyObject* args)
{
    FUNC_BEGIN
    
    PyObject* cone = PyTuple_GetItem( args, 0 );
    PyObject* to_compute = PyTuple_GetItem( args, 1 );
    
    if( !is_cone(cone) ){
        PyErr_SetString( PyQNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    Cone<renf_elem_class>* cone_ptr = get_cone_renf(cone);
    return NmzIsComputed(cone_ptr, to_compute);
    
    FUNC_END
}

/***************************************************************************
 * 
 * NmzResult
 * 
 ***************************************************************************/

template<typename Integer>
PyObject* _NmzResultImpl(Cone<Integer>* C, PyObject* prop_obj)
{
    
    FUNC_BEGIN

    string prop = PyUnicodeToString( prop_obj );

    libQnormaliz::QConeProperty::Enum p = libQnormaliz::toConeProperty(prop);
    current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    ConeProperties notComputed = C->compute(ConeProperties(p));
    PyOS_setsig(SIGINT,current_interpreter_sigint_handler);
    
    if (notComputed.any()) {
        return Py_None;
    }

    switch (p) {
    case libQnormaliz::QConeProperty::Generators:
        return NmzMatrixToPyList(C->getGenerators());

    case libQnormaliz::QConeProperty::ExtremeRays:
        return NmzMatrixToPyList(C->getExtremeRays());

    case libQnormaliz::QConeProperty::VerticesOfPolyhedron:
        return NmzMatrixToPyList(C->getVerticesOfPolyhedron());

    case libQnormaliz::QConeProperty::SupportHyperplanes:
        return NmzMatrixToPyList(C->getSupportHyperplanes());

//     case libQnormaliz::QConeProperty::SuppHypsFloat:
//         return NmzMatrixToPyList(C->getSuppHypsFloat());

//     case libQnormaliz::QConeProperty::TriangulationSize:
//         return NmzToPyNumber(C->getTriangulationSize());

//     case libQnormaliz::QConeProperty::TriangulationDetSum:
//         return NmzToPyNumber(C->getTriangulationDetSum());

    case libQnormaliz::QConeProperty::Triangulation:
        return NmzTriangleListToPyList<Integer>(C->getTriangulation());

//     case libQnormaliz::QConeProperty::Multiplicity:
//         return NmzToPyList(C->getMultiplicity());
    
//     case libQnormaliz::QConeProperty::Integral:
//         return NmzToPyList(C->getIntegral());
    
//     case libQnormaliz::QConeProperty::VirtualMultiplicity:
//         return NmzToPyList(C->getVirtualMultiplicity());

//     case libQnormaliz::QConeProperty::RecessionRank:
//         return NmzToPyNumber(C->getRecessionRank());

    case libQnormaliz::QConeProperty::AffineDim:
        return NmzToPyNumber(C->getAffineDim());

//     case libQnormaliz::QConeProperty::ModuleRank:
//         return NmzToPyNumber(C->getModuleRank());

//     case libQnormaliz::QConeProperty::HilbertBasis:
//         return NmzMatrixToPyList(C->getHilbertBasis());

    case libQnormaliz::QConeProperty::MaximalSubspace:
        return NmzMatrixToPyList(C->getMaximalSubspace());

//     case libQnormaliz::QConeProperty::ModuleGenerators:
//         return NmzMatrixToPyList(C->getModuleGenerators());

//     case libQnormaliz::QConeProperty::Deg1Elements:
//         return NmzMatrixToPyList(C->getDeg1Elements());

//     case libQnormaliz::QConeProperty::HilbertSeries:
//         {
//         bool is_HSOP = C->isComputed(libQnormaliz::QConeProperty::HSOP);
//         return NmzHilbertSeriesToPyList(C->getHilbertSeries(),is_HSOP);
//         }

//     case libQnormaliz::QConeProperty::WeightedEhrhartSeries:
//         return NmzWeightedEhrhartSeriesToPyList(C->getWeightedEhrhartSeries());
        

//     case libQnormaliz::QConeProperty::Grading:
//         {
//         vector<Integer> grad = C->getGrading();
//         Integer denom = C->getGradingDenom();
//         PyObject * return_list = PyList_New(2);
//         PyList_SetItem( return_list, 0, NmzVectorToPyList(grad) );
//         PyList_SetItem( return_list, 1, NmzToPyNumber( denom ) );
//         return return_list;
//         }

    case libQnormaliz::QConeProperty::IsPointed:
        return BoolToPyBool(C->isPointed());

//     case libQnormaliz::QConeProperty::IsDeg1ExtremeRays:
//         return BoolToPyBool(C->isDeg1ExtremeRays());

//     case libQnormaliz::QConeProperty::IsDeg1HilbertBasis:
//         return BoolToPyBool(C->isDeg1HilbertBasis());

//     case libQnormaliz::QConeProperty::IsIntegrallyClosed:
//         return BoolToPyBool(C->isIntegrallyClosed());

//     case libQnormaliz::QConeProperty::OriginalMonoidGenerators:
//         return NmzMatrixToPyList(C->getOriginalMonoidGenerators());

//     case libQnormaliz::QConeProperty::IsReesPrimary:
//         return BoolToPyBool(C->isReesPrimary());

//     case libQnormaliz::QConeProperty::ReesPrimaryMultiplicity:
//         return NmzToPyNumber(C->getReesPrimaryMultiplicity());

//     case libQnormaliz::QConeProperty::StanleyDec:
//         return NmzStanleyDecToPyList(C->getStanleyDec());

//     case libQnormaliz::QConeProperty::ExcludedFaces:
//         return NmzMatrixToPyList(C->getExcludedFaces());

    case libQnormaliz::QConeProperty::Dehomogenization:
        return NmzVectorToPyList(C->getDehomogenization());

//     case libQnormaliz::QConeProperty::InclusionExclusionData:
//         return NmzTriangleListToPyList<long>(C->getInclusionExclusionData());

//     case libQnormaliz::QConeProperty::ClassGroup:
//         return NmzVectorToPyList(C->getClassGroup());
    
    case libQnormaliz::QConeProperty::IsInhomogeneous:
        return BoolToPyBool(C->isInhomogeneous());
    
//     /* Sublattice properties */
    
    case libQnormaliz::QConeProperty::Equations:
        return NmzMatrixToPyList(C->getSublattice().getEquations());
    
//     case libQnormaliz::QConeProperty::Congruences:
//         return NmzMatrixToPyList(C->getSublattice().getCongruences());
    
    case libQnormaliz::QConeProperty::EmbeddingDim:
        return NmzToPyNumber(C->getEmbeddingDim());
    
    case libQnormaliz::QConeProperty::Rank:
        return NmzToPyNumber(C->getRank());
    
//     case libQnormaliz::QConeProperty::Sublattice:
//         return _NmzBasisChangeIntern(C);
    
//     case libQnormaliz::QConeProperty::ExternalIndex:
//         return NmzToPyNumber(C->getSublattice().getExternalIndex());
    
//     case libQnormaliz::QConeProperty::InternalIndex:
//         return NmzToPyNumber(C->getIndex());
    
//     case libQnormaliz::QConeProperty::WitnessNotIntegrallyClosed:
//         return NmzVectorToPyList(C->getWitnessNotIntegrallyClosed());
    
    
//     /* New stuff */
    
//     case libQnormaliz::QConeProperty::GradingDenom:
//         return NmzToPyNumber(C->getGradingDenom());
    
//     case libQnormaliz::QConeProperty::UnitGroupIndex:
//         return NmzToPyNumber(C->getUnitGroupIndex());
    
//     case libQnormaliz::QConeProperty::ModuleGeneratorsOverOriginalMonoid:
//         return NmzMatrixToPyList(C->getModuleGeneratorsOverOriginalMonoid());
    
//     case libQnormaliz::QConeProperty::IntegerHull:
//     {
//         Cone<Integer>* hull = new Cone<Integer>( C->getIntegerHullCone() );
//         return pack_cone( hull ); 
//     }
    
//     case libQnormaliz::QConeProperty::ProjectCone:
//     {
//         Cone<Integer>* projection = new Cone<Integer>(C->getProjectCone());
//         return pack_cone( projection );
//     }
    
//     case libQnormaliz::QConeProperty::HilbertQuasiPolynomial:
//         return NmzHilbertQuasiPolynomialToPyList<mpz_class>(C->getHilbertSeries()); //FIXME: Why is this return value not parametrized, but mpz_class only?
    
//     case libQnormaliz::QConeProperty::WeightedEhrhartQuasiPolynomial:
//         return NmzWeightedEhrhartQuasiPolynomialToPyList<mpz_class>(C->getIntData());
    
//     case libQnormaliz::QConeProperty::IsTriangulationNested:
//         return BoolToPyBool(C->isTriangulationNested());
        
//     case libQnormaliz::QConeProperty::IsTriangulationPartial:
//         return BoolToPyBool(C->isTriangulationPartial());
        
    case libQnormaliz::QConeProperty::ConeDecomposition:
        return NmzBoolMatrixToPyList(C->getOpenFacets());
    
//     case libQnormaliz::QConeProperty::IsGorenstein:
//         return BoolToPyBool(C->isGorenstein());
        
//     case libQnormaliz::QConeProperty::GeneratorOfInterior:
//         return NmzVectorToPyList(C->getGeneratorOfInterior());
        
//     case libQnormaliz::QConeProperty::VerticesFloat:
//         return NmzMatrixToPyList(C->getVerticesFloat());
        
//     case libQnormaliz::QConeProperty::Volume:
//         return NmzToPyList(C->getVolume());

//     case libQnormaliz::QConeProperty::EuclideanVolume:
//         return NmzToPyNumber(C->getEuclideanVolume());

//  the following properties are compute options and do not return anything
    case libQnormaliz::QConeProperty::DualMode:
    case libQnormaliz::QConeProperty::DefaultMode:
    case libQnormaliz::QConeProperty::Approximate:
    case libQnormaliz::QConeProperty::BottomDecomposition:
    case libQnormaliz::QConeProperty::KeepOrder:
    case libQnormaliz::QConeProperty::NoBottomDec:
    case libQnormaliz::QConeProperty::PrimalMode:
    case libQnormaliz::QConeProperty::Symmetrize:
    case libQnormaliz::QConeProperty::NoSymmetrization:
    case libQnormaliz::QConeProperty::BigInt:
    case libQnormaliz::QConeProperty::HSOP:
        PyErr_SetString( PyQNormaliz_cppError, "QConeProperty is input-only" );
        return NULL;
// #if NMZ_RELEASE >= 30200
//     case libQnormaliz::QConeProperty::NoSubdivision:
//         PyErr_SetString( PyQNormaliz_cppError, "QConeProperty is input-only" );
//         return NULL;
// #endif
    default:
        PyErr_SetString( PyQNormaliz_cppError, "Unknown cone property" );
        return NULL;
        break;
    }

    return Py_None;

    FUNC_END
}

PyObject* get_dict_elem_or_null(PyObject* dict, const char* key){
    PyObject* return_obj = PyDict_GetItemString( dict, key );
    return return_obj;
}

PyObject* _NmzResult( PyObject* self, PyObject* args, PyObject* kwargs ){
    
    FUNC_BEGIN
    PyObject* cone = PyTuple_GetItem( args, 0 );
    PyObject* prop = PyTuple_GetItem( args, 1 );
    
    if( !is_cone( cone ) ){
        PyErr_SetString( PyQNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    if( !string_check( prop ) ){
        PyErr_SetString( PyQNormaliz_cppError, "Second argument must be a unicode string" );
        return NULL;
    }
    
    Cone<renf_elem_class>* cone_ptr = get_cone_renf(cone);


    //Set handlers
    if(kwargs){
        RationalHandler = PyDict_GetItemString(kwargs,"RationalHandler");
        NumberfieldElementHandler = get_dict_elem_or_null(kwargs,"NumberfieldElementHandler");
        VectorHandler = get_dict_elem_or_null(kwargs,"VectorHandler");
        MatrixHandler = get_dict_elem_or_null(kwargs,"MatrixHandler");
    }
    
    PyObject* result = _NmzResultImpl(cone_ptr, prop);

    RationalHandler = NULL;
    NumberfieldElementHandler = NULL;
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
        PyErr_SetString( PyQNormaliz_cppError, "Argument must be True or False" );
        return NULL;
    }
    return BoolToPyBool(libQnormaliz::setVerboseDefault(value == Py_True));
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
        PyErr_SetString( PyQNormaliz_cppError, "First argument must be a cone" );
        return NULL;
    }
    
    PyObject* value = PyTuple_GetItem( args, 1 );
    if (value != Py_True && value != Py_False){
        PyErr_SetString( PyQNormaliz_cppError, "Second argument must be True or False" );
        return NULL;
    }
    
    Cone<renf_elem_class>* cone_ptr = get_cone_renf(cone);
    return NmzSetVerbose(cone_ptr, value);
    
    FUNC_END
    
}

// /***************************************************************************
//  * 
//  * Get euclidian volume
//  * 
//  ***************************************************************************/

// PyObject* NmzGetEuclideanVolume(PyObject* self, PyObject* args ){
    
//     FUNC_BEGIN
    
//     PyObject* cone = PyTuple_GetItem( args, 0 );
    
//     if( !is_cone( cone ) ){
//         PyErr_SetString( PyQNormaliz_cppError, "First argument must be a cone" );
//         return NULL;
//     }
    
//     current_interpreter_sigint_handler = PyOS_setsig(SIGINT,signal_handler);
    
//     if( cone_name_str == string(PyCapsule_GetName(cone)) ){
//         Cone<mpz_class>* cone_ptr = get_cone_mpz(cone);
//         PyObject* return_value = NmzToPyNumber(cone_ptr->getEuclideanVolume());
//         PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
//         return return_value;
//     }else{
//         Cone<long long>* cone_ptr = get_cone_long(cone);
//         PyObject* return_value = NmzToPyNumber(cone_ptr->getEuclideanVolume());
//         PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
//         return return_value;
//     }
//     PyOS_setsig( SIGINT, current_interpreter_sigint_handler );
    
//     FUNC_END
    
// }

/***************************************************************************
 * 
 * Set number of threads
 * 
 ***************************************************************************/

PyObject* NmzSetNumberOfNormalizThreads(PyObject* self, PyObject* args ){
    
    FUNC_BEGIN
    
    PyObject* num_treads = PyTuple_GetItem( args, 0 );
    
    if( !PyLong_Check( num_treads ) ){
        PyErr_SetString( PyQNormaliz_cppError, "First argument must be an integer" );
        return NULL;
    }
    
    long num_threads_long = PyLong_AsLong( num_treads );
    
    num_threads_long = libQnormaliz::set_thread_limit( num_threads_long );
    
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

static PyMethodDef PyQNormaliz_cppMethods[] = {
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {"NmzCone",  (PyCFunction)_NmzQCone, METH_VARARGS|METH_KEYWORDS,
     "Create a cone"},
    // {"NmzConeCopy",  (PyCFunction)_NmzConeCopy,METH_VARARGS,
    //  "Copy an existing cone" },
    {"NmzCompute", (PyCFunction)_NmzCompute_Outer, METH_VARARGS,
     "Compute some stuff"},
    {"NmzIsComputed", (PyCFunction)NmzIsComputed_Outer, METH_VARARGS,
     "Check if property is computed "},
    // {"NmzSetGrading", (PyCFunction)NmzSetGrading, METH_VARARGS,
    //  "Reset the grading of a cone"},
    {"NmzResult", (PyCFunction)_NmzResult, METH_VARARGS|METH_KEYWORDS,
      "Return cone property" },
    { "NmzSetVerboseDefault", (PyCFunction)NmzSetVerboseDefault, METH_VARARGS,
      "Set verbosity" },
    { "NmzSetVerbose", (PyCFunction)NmzSetVerbose_Outer, METH_VARARGS,
      "Set verbosity of cone" },
    { "NmzListConeProperties", (PyCFunction)NmzListConeProperties,METH_NOARGS,
      "List all available properties" },
    // { "NmzHilbertSeries", (PyCFunction)NmzHilbertSeries_Outer, METH_VARARGS,
    //   "Returns Hilbert series, either HSOP or not" },
    // { "NmzGetPolynomial", (PyCFunction)NmzGetPolynomial, METH_VARARGS,
    //   "Returns grading polynomial" },
    { "NmzSetNumberOfNormalizThreads", (PyCFunction)NmzSetNumberOfNormalizThreads, METH_VARARGS,
      "Sets the Normaliz thread limit" },
    // { "NmzSetNrCoeffQuasiPol", (PyCFunction)NmzSetNrCoeffQuasiPol, METH_VARARGS,
    //   "Sets the period bound for the quasi-polynomial" },
    // { "NmzGetEuclideanVolume", (PyCFunction)NmzGetEuclideanVolume, METH_VARARGS,
    //   "Returns euclidean volume of cone as float" },
    // { "NmzGetHilbertSeriesExpansion", (PyCFunction)NmzGetHilbertSeriesExpansion, METH_VARARGS,
    //   "Returns expansion of the hilbert series" },
    // { "NmzGetWeightedEhrhartSeriesExpansion", (PyCFunction)NmzGetWeightedEhrhartSeriesExpansion, METH_VARARGS,
    //   "Returns expansion of the weighted Ehrhart series" },
    {NULL, }        /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3

static int PyQNormaliz_cpp_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int PyQNormaliz_cpp_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "PyQNormaliz_cpp",
        NULL,
        sizeof(struct module_state),
        PyQNormaliz_cppMethods,
        NULL,
        PyQNormaliz_cpp_traverse,
        PyQNormaliz_cpp_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC PyInit_PyQNormaliz_cpp(void)

#else
#define INITERROR return

extern "C" void initPyQNormaliz_cpp(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("PyQNormaliz_cpp", PyQNormaliz_cppMethods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException(const_cast<char*>("PyQNormaliz_cpp.INITError"), NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }
    
    QNormalizError = PyErr_NewException(const_cast<char*>("QNormaliz.error"), NULL, NULL );
    Py_INCREF( QNormalizError );
    PyQNormaliz_cppError = PyErr_NewException(const_cast<char*>("QNormaliz.interface_error"), NULL, NULL );
    Py_INCREF( PyQNormaliz_cppError );
    
    PyModule_AddObject( module, "error", QNormalizError );
    PyModule_AddObject( module, "error", PyQNormaliz_cppError );
    
    current_interpreter_sigint_handler = PyOS_getsig( SIGINT );

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

