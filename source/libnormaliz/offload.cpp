#ifdef NMZ_MIC_OFFLOAD

#pragma offload_attribute (push, target(mic))
#include "offload.h"
#include "matrix.h"
#include "full_cone.h"
#include <iostream>

namespace libnormaliz {

using namespace std;

template<typename Integer>
void fill_vector(vector<Integer>& v, long size, Integer* data)
{
  for (long i=0; i<size; i++)
    v[i] = data[i];
}

template<typename Integer>
void fill_plain_vector(Integer* data, long size, const vector<Integer>& v)
{
  for (long i=0; i<size; i++)
      data[i] = v[i];

}

template<typename Integer>
void fill_matrix(Matrix<Integer>& M, long rows, long cols, Integer* data)
{
  for (long i=0; i<rows; i++)
    for (long j=0; j<cols; j++)
    {
      M[i][j] = data[i*cols+j];
    }
}

template<typename Integer>
void fill_plain(Integer* data, long rows, long cols, const Matrix<Integer>& M)
{
  for (long i=0; i<rows; i++)
    for (long j=0; j<cols; j++)
      data[i*cols+j] = M[i][j];

}
#pragma offload_attribute (pop)

//---------------------begin-class-implementation----------------------------

template<typename Integer>
OffloadHandler<Integer>::OffloadHandler(const Full_Cone<Integer>& fc, int mic_number)
  : mic_nr(mic_number)
{
  create_full_cone(fc);

  transfer_bools(fc);
  transfer_support_hyperplanes(fc);
  transfer_grading(fc);            // including truncation and shift
  //TODO transfer_triangulation_info(fc); // extreme rays, deg1_triangulation, Order_Vector

  //prepare_pyramid_evaluation();    //
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::create_full_cone(const Full_Cone<Integer>& fc)
{
  const Matrix<Integer>& M = fc.Generators;
  long nr = M.nr_of_rows();
  long nc = M.nr_of_columns();
  long size = nr*nc;
  Integer *data = new Integer[size];
  fill_plain(data, nr, nc, M);

  cout << "Offload data to mic, fc_ptr value on cpu " << fc_ptr << endl;
  // offload to mic, copy data and free it afterwards, but keep a pointer to the created C++ matrix
  #pragma offload target(mic:mic_nr) in(nr,nc) in(data: length(size) ONCE)
  {
    Matrix<Integer> gens(nr, nc);
    fill_matrix(gens, nr, nc, data);
    fc_ptr = new Full_Cone<Integer>(gens);
    cout << "fc_ptr value on mic " << fc_ptr << endl;
  }
  cout << "After offload fc_ptr value on cpu " << fc_ptr << endl;
  delete[] data;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_bools(const Full_Cone<Integer>& fc)
{
  cout << "transfer_bools" << endl;
  #pragma offload target(mic:mic_nr)
  {
    bool foo = fc_ptr->inhomogeneous;  // prevents segfault
    fc_ptr->inhomogeneous = fc.inhomogeneous;
    fc_ptr->do_Hilbert_basis = fc.do_Hilbert_basis;
    fc_ptr->do_h_vector = fc.do_h_vector;
    fc_ptr->keep_triangulation = fc.keep_triangulation;
    fc_ptr->do_multiplicity = fc.do_multiplicity;
    fc_ptr->do_determinants = fc.do_determinants;
    fc_ptr->do_triangulation = fc.do_triangulation;
    fc_ptr->do_deg1_elements = fc.do_deg1_elements;
    fc_ptr->do_Stanley_dec = fc.do_Stanley_dec;
    fc_ptr->do_approximation = fc.do_approximation;
    fc_ptr->do_default_mode = fc.do_default_mode;
  }
  cout << "transfer_bools done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_support_hyperplanes(const Full_Cone<Integer>& fc)
{
  cout << "transfer_support_hyperplanes" << endl;
  const Matrix<Integer>& M = fc.Support_Hyperplanes;
  long nr = M.nr_of_rows();
  long nc = M.nr_of_columns();
  long size = nr*nc;
  assert(size > 0); // make sure there are support hyperplanes computed
  Integer *data = new Integer[size];
  fill_plain(data, nr, nc, M);

  cout << "Offload data to mic, fc_ptr value on cpu " << fc_ptr << endl;
  // offload to mic, copy data and free it afterwards, but keep a pointer to the created C++ matrix
  #pragma offload target(mic:mic_nr) in(nr,nc) in(data: length(size) ONCE)
  {
    cout << "fc_ptr value on mic " << fc_ptr << endl;
    fc_ptr->Support_Hyperplanes = Matrix<Integer>(nr, nc);
    fill_matrix(fc_ptr->Support_Hyperplanes, nr, nc, data);
    fc_ptr->is_Computed.set(ConeProperty::SupportHyperplanes);
    fc_ptr->do_all_hyperplanes = false;
  }
  cout << "After offload fc_ptr value on cpu " << fc_ptr << endl;
  delete[] data;

  cout << "transfer_support_hyperplanes done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_grading(const Full_Cone<Integer>& fc)
{
  cout << "transfer_grading" << endl;
  long dim = fc.dim;
  if (fc.inhomogeneous) 
  {
    Integer *data = new Integer[dim];
    fill_plain_vector(data, dim, fc.Truncation);

    #pragma offload target(mic:mic_nr) in(dim) in(data: length(dim) ONCE)
    {
      fc_ptr->Truncation = vector<Integer>(dim);
      fill_vector(fc_ptr->Truncation, dim, data);
    }
    delete[] data;
  }

  if (fc.isComputed(ConeProperty::Grading))
  {
    Integer *data = new Integer[dim];
    fill_plain_vector(data, dim, fc.Grading);

    #pragma offload target(mic:mic_nr) in(dim) in(data: length(dim) ONCE)
    {
      fc_ptr->Grading = vector<Integer>(dim);
      fill_vector(fc_ptr->Grading, dim, data);
      fc_ptr->is_Computed.set(ConeProperty::Grading);
      fc_ptr->set_degrees();
    }
    delete[] data;
  }

  if (fc.isComputed(ConeProperty::Shift))
  {
    #pragma offload target(mic:mic_nr)
    {
      fc_ptr->shift = fc.shift;
      fc_ptr->is_Computed.set(ConeProperty::Shift);
    }
  }

  cout << "transfer_grading done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::print_on_mic() const
{
  cout << "Offloaded print" << endl;
  #pragma offload target(mic:mic_nr)
  {
    fc_ptr->Generators.pretty_print(cout);
  }
}

template<typename Integer>
void OffloadHandler<Integer>::compute_on_mic(long a, long b)
{
  cout << "Offloaded computation" << endl;
  #pragma offload target(mic:mic_nr)
  {
    cout << "Rank computed on mic " << fc_ptr->Generators.rank() << endl;
  }
}

template<typename Integer>
long OffloadHandler<Integer>::collect_data()
{
//  #pragma offload target(mic:mic_nr) out()
	return 0;
}

//---------------------------------------------------------------------------

template<typename Integer>
Matrix<Integer> OffloadHandler<Integer>::transfer_from_mic()
{
  cout << "Transfer data back to cpu" << endl;
  long a,b;
  #pragma offload target(mic:mic_nr) out(a) out(b)
  {
    a = fc_ptr->Generators.nr_of_rows();
    b = fc_ptr->Generators.nr_of_columns();
  }

  cout << "transfer matrix of size a = " << a << " by  b = " << b << endl;
  long size = a*b;
  Integer *ret = new Integer[size];
  #pragma offload target(mic:mic_nr) out(ret : length(size))
  {
    fill_plain(ret, a, b, fc_ptr->Generators);
  }
  Matrix<Integer> M (a, b);
  fill_matrix(M, a, b, ret);
  delete[] ret;
  return M;
}

template<typename Integer>
OffloadHandler<Integer>::~OffloadHandler()
{
  cout << "OffloadHandler destructor" << endl;
  #pragma offload target(mic:mic_nr)
  {
    delete fc_ptr;
  }
}


/***************** Instantiation for template parameter long long *****************/

template class OffloadHandler<long long int>;
//template class OffloadHandler<mpz_class>;


/***************** Offload test *****************/

void offload_test()
{
  typedef long long Integer;

  // initial offload for better timing comparisons of the following offloads
  cout << "initial offload for better timing comparisons of the following offloads"
       << endl;
  #pragma offload target(mic:0)
  { }
  cout << "done." << endl;

  int a = 4, b = 3;
  long SIZE = a*b;
  Integer data[SIZE];
  for (long i=0; i<SIZE; i++) data[i] = i+1;

  Matrix<Integer> m1(a,b);
  Matrix<Integer> m2(a+1,b+1);
  m1.random(10);
  m2.random(10);

  // offload the full cones
  Full_Cone<Integer> fc1(m1);
  fc1.get_supphyps_from_copy(true);          // from_scratch = true
  OffloadHandler<Integer> fc1_off(fc1);
  cout << "first offload completed" << endl;
  fc1_off.print_on_mic();

  Full_Cone<Integer> fc2(m2);
  fc2.get_supphyps_from_copy(true);          // from_scratch = true
  OffloadHandler<Integer> fc2_off(fc2);
  fc2_off.print_on_mic();

  // work with m1
  fc1_off.print_on_mic();
  fc1_off.compute_on_mic(1,2);
  fc1_off.compute_on_mic(0,2);

  // work with m2
  fc2_off.print_on_mic();
  fc2_off.compute_on_mic(1,2);
  fc2_off.compute_on_mic(2,1);

  // get results back
  Matrix<Integer> ret = fc1_off.transfer_from_mic();
  ret.read();

  // destructor of Offload handler will clear data on mic
}

} // end namespace libnormaliz


#else //NMZ_MIC_OFFLOAD

// no offloading available
namespace libnormaliz {
void offload_test() {}
} // end namespace libnormaliz

#endif //NMZ_MIC_OFFLOAD
