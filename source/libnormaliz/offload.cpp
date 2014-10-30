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
OffloadHandler<Integer>::OffloadHandler(Full_Cone<Integer>& fc, int mic_number)
  : mic_nr(mic_number),
    local_fc_ref(fc)
{
  create_full_cone();

  transfer_bools();
  transfer_support_hyperplanes();
  transfer_grading();            // including truncation and shift
  transfer_triangulation_info(); // extreme rays, deg1_triangulation, Order_Vector

  //prepare_pyramid_evaluation();    //
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::create_full_cone()
{
  const Matrix<Integer>& M = local_fc_ref.Generators;
  long nr = M.nr_of_rows();
  long nc = M.nr_of_columns();
  long size = nr*nc;
  Integer *data = new Integer[size];
  fill_plain(data, nr, nc, M);

  cout << "Offload data to mic, offload_fc_ptr value on cpu " << offload_fc_ptr << endl;
  // offload to mic, copy data and free it afterwards, but keep a pointer to the created C++ matrix
  #pragma offload target(mic:mic_nr) in(nr,nc) in(data: length(size) ONCE)
  {
    Matrix<Integer> gens(nr, nc);
    fill_matrix(gens, nr, nc, data);
    offload_fc_ptr = new Full_Cone<Integer>(gens);
    cout << "offload_fc_ptr value on mic " << offload_fc_ptr << endl;
  }
  cout << "After offload offload_fc_ptr value on cpu " << offload_fc_ptr << endl;
  delete[] data;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_bools()
{
  cout << "transfer_bools" << endl;
  Full_Cone<Integer>& foo_loc = local_fc_ref;  // prevents segfault
  //TODO segfaults should be resolved in intel compiler version 2015
  #pragma offload target(mic:mic_nr)
  {
    bool foo = offload_fc_ptr->inhomogeneous;  // prevents segfault
    offload_fc_ptr->inhomogeneous      = foo_loc.inhomogeneous;
    offload_fc_ptr->do_Hilbert_basis   = foo_loc.do_Hilbert_basis;
    offload_fc_ptr->do_h_vector        = foo_loc.do_h_vector;
    offload_fc_ptr->keep_triangulation = foo_loc.keep_triangulation;
    offload_fc_ptr->do_multiplicity    = foo_loc.do_multiplicity;
    offload_fc_ptr->do_determinants    = foo_loc.do_determinants;
    offload_fc_ptr->do_triangulation   = foo_loc.do_triangulation;
    offload_fc_ptr->do_deg1_elements   = foo_loc.do_deg1_elements;
    offload_fc_ptr->do_Stanley_dec     = foo_loc.do_Stanley_dec;
    offload_fc_ptr->do_approximation   = foo_loc.do_approximation;
    offload_fc_ptr->do_default_mode    = foo_loc.do_default_mode;
    // deg1_generated could be set more precise
    offload_fc_ptr->deg1_triangulation = foo_loc.deg1_generated;
  }
  cout << "transfer_bools done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_support_hyperplanes()
{
  cout << "transfer_support_hyperplanes" << endl;
  const Matrix<Integer>& M = local_fc_ref.Support_Hyperplanes;
  long nr = M.nr_of_rows();
  long nc = M.nr_of_columns();
  long size = nr*nc;
  assert(size > 0); // make sure there are support hyperplanes computed
  Integer *data = new Integer[size];
  fill_plain(data, nr, nc, M);

  cout << "Offload data to mic, offload_fc_ptr value on cpu " << offload_fc_ptr << endl;
  // offload to mic, copy data and free it afterwards, but keep a pointer to the created C++ matrix
  #pragma offload target(mic:mic_nr) in(nr,nc) in(data: length(size) ONCE)
  {
    cout << "offload_fc_ptr value on mic " << offload_fc_ptr << endl;
    offload_fc_ptr->Support_Hyperplanes = Matrix<Integer>(nr, nc);
    fill_matrix(offload_fc_ptr->Support_Hyperplanes, nr, nc, data);
    offload_fc_ptr->is_Computed.set(ConeProperty::SupportHyperplanes);
    offload_fc_ptr->do_all_hyperplanes = false;
  }
  cout << "After offload offload_fc_ptr value on cpu " << offload_fc_ptr << endl;
  delete[] data;

  cout << "transfer_support_hyperplanes done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_grading()
{
  cout << "transfer_grading" << endl;
  long dim = local_fc_ref.dim;
  if (local_fc_ref.inhomogeneous)
  {
    Integer *data = new Integer[dim];
    fill_plain_vector(data, dim, local_fc_ref.Truncation);

    #pragma offload target(mic:mic_nr) in(dim) in(data: length(dim) ONCE)
    {
      offload_fc_ptr->Truncation = vector<Integer>(dim);
      fill_vector(offload_fc_ptr->Truncation, dim, data);
    }
    delete[] data;
  }

  if (local_fc_ref.isComputed(ConeProperty::Grading))
  {
    Integer *data = new Integer[dim];
    fill_plain_vector(data, dim, local_fc_ref.Grading);

    #pragma offload target(mic:mic_nr) in(dim) in(data: length(dim) ONCE)
    {
      offload_fc_ptr->Grading = vector<Integer>(dim);
      fill_vector(offload_fc_ptr->Grading, dim, data);
      offload_fc_ptr->is_Computed.set(ConeProperty::Grading);
      offload_fc_ptr->set_degrees();
    }
    delete[] data;
  }

  if (local_fc_ref.isComputed(ConeProperty::Shift))
  {
    #pragma offload target(mic:mic_nr)
    {
      offload_fc_ptr->shift = local_fc_ref.shift;
      offload_fc_ptr->is_Computed.set(ConeProperty::Shift);
    }
  }

  cout << "transfer_grading done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_triangulation_info()
{
  cout << "transfer_triangulation_info" << endl;
  long dim = local_fc_ref.dim;
  long nr_gen = local_fc_ref.nr_gen;

  if (local_fc_ref.isComputed(ConeProperty::ExtremeRays))
  {
    bool *data = new bool[nr_gen];
    fill_plain_vector(data, nr_gen, local_fc_ref.Extreme_Rays);

    #pragma offload target(mic:mic_nr) in(nr_gen) in(data: length(nr_gen) ONCE)
    {

      offload_fc_ptr->Extreme_Rays = vector<bool>(nr_gen);
      fill_vector(offload_fc_ptr->Extreme_Rays, nr_gen, data);
    }
    delete[] data;
  }

  // always transfer the order vector  //TODO ensure it is computed!
  {
    Integer *data = new Integer[dim];
    fill_plain_vector(data, dim, local_fc_ref.Order_Vector);

    #pragma offload target(mic:mic_nr) in(dim) in(data: length(dim) ONCE)
    {
      offload_fc_ptr->Order_Vector = vector<Integer>(dim);
      fill_vector(offload_fc_ptr->Order_Vector, dim, data);
    }
    delete[] data;
  }

  if (local_fc_ref.isComputed(ConeProperty::Shift))
  {
    #pragma offload target(mic:mic_nr)
    {
      offload_fc_ptr->shift = local_fc_ref.shift;
      offload_fc_ptr->is_Computed.set(ConeProperty::Shift);
    }
  }

  cout << "transfer_triangulation_info done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::print_on_mic() const
{
  cout << "Offloaded print" << endl;
  #pragma offload target(mic:mic_nr)
  {
    offload_fc_ptr->Generators.pretty_print(cout);
  }
}

template<typename Integer>
void OffloadHandler<Integer>::compute_on_mic(long a, long b)
{
  cout << "Offloaded computation" << endl;
  #pragma offload target(mic:mic_nr)
  {
    cout << "Rank computed on mic " << offload_fc_ptr->Generators.rank() << endl;
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
    a = offload_fc_ptr->Generators.nr_of_rows();
    b = offload_fc_ptr->Generators.nr_of_columns();
  }

  cout << "transfer matrix of size a = " << a << " by  b = " << b << endl;
  long size = a*b;
  Integer *ret = new Integer[size];
  #pragma offload target(mic:mic_nr) out(ret : length(size))
  {
    fill_plain(ret, a, b, offload_fc_ptr->Generators);
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
    delete offload_fc_ptr;
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
  fc1.Order_Vector = vector<Integer>(b);
  OffloadHandler<Integer> fc1_off(fc1);
  cout << "first offload completed" << endl;
  fc1_off.print_on_mic();

  Full_Cone<Integer> fc2(m2);
  fc2.get_supphyps_from_copy(true);          // from_scratch = true
  fc2.Order_Vector = vector<Integer>(b+1);
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
