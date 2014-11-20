#ifdef NMZ_MIC_OFFLOAD

#pragma offload_attribute (push, target(mic))
#include "offload_handler.h"
#include "offload.h"  // offload system header
#include "matrix.h"
#include "full_cone.h"
#include "list_operations.h"
#include "vector_operations.h"
#include "my_omp.h"
#include "HilbertSeries.h"
#include <iostream>

namespace libnormaliz {

using namespace std;

// transfering vector
template<typename Integer>
void fill_vector(vector<Integer>& v, long size, Integer* data)
{
  for (long i=0; i<size; i++)
    v[i] = data[i];
}

template<typename Integer>
void fill_plain(Integer* data, long size, const vector<Integer>& v)
{
  for (long i=0; i<size; i++)
      data[i] = v[i];
}

// transfering Matrix
template<typename Integer>
void fill_matrix(Matrix<Integer>& M, long rows, long cols, Integer* data)
{
  for (long i=0; i<rows; i++)
    for (long j=0; j<cols; j++)
      M[i][j] = data[i*cols+j];
}

template<typename Integer>
void fill_plain(Integer* data, long rows, long cols, const Matrix<Integer>& M)
{
  for (long i=0; i<rows; i++)
    for (long j=0; j<cols; j++)
      data[i*cols+j] = M[i][j];
}

// transfering list<vector>
// the vectors may have different lengths
// fill_list_vector also creates the list entries and appendes them to the list!
template<typename Integer>
void fill_list_vector(list< vector<Integer> >& l, long plain_size, Integer* data)
{
  Integer* data_end = data + plain_size; // position after last entry
  while (data < data_end)
  {
    l.push_back(vector<Integer>(*data));
    fill_vector(l.back(), *data, data+1);
    data += *data + 1;
  }
}

template<typename Integer>
void fill_plain(Integer* data, long size, const list< vector<Integer> >& l)
{
  long v_size;
  typename list< vector<Integer> >::const_iterator it;
  for (it = l.begin(); it != l.end(); it++)
  {
    v_size = it->size();
    *data = v_size;
    fill_plain(++data, v_size, *it);
    data += v_size;
  }
}

template<typename Integer>
long plain_size(const list< vector<Integer> >& l)
{
  long size = 0;
  typename list< vector<Integer> >::const_iterator it;
  for (it = l.begin(); it != l.end(); it++)
    size += it->size() + 1;
  return size;
}
#pragma offload_attribute (pop)

//-------------------------- OffloadHandler ---------------------------------

template<typename Integer>
OffloadHandler<Integer>::OffloadHandler(Full_Cone<Integer>& fc, int mic_number)
  : mic_nr(mic_number),
    running(false),
    local_fc_ref(fc)
{
  create_full_cone();

  transfer_bools();
  transfer_support_hyperplanes();
  transfer_grading();            // including truncation and shift
  transfer_triangulation_info(); // extreme rays, deg1_triangulation, Order_Vector

  primal_algorithm_initialize();
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
  // offload to mic, copy data and free it afterwards, but keep a pointer to the created Full_Cone
  #pragma offload target(mic:mic_nr) in(nr,nc) in(data: length(size) ONCE)
  {
    //omp_set_num_threads(118);
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
    offload_fc_ptr->pointed = foo_loc.pointed;  // was locally computed in MicOffloader
    offload_fc_ptr->is_Computed.set(ConeProperty::IsPointed);
    verbose = true;
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
    offload_fc_ptr->nrSupport_Hyperplanes = nr;
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
    fill_plain(data, dim, local_fc_ref.Truncation);

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
    fill_plain(data, dim, local_fc_ref.Grading);

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
    fill_plain(data, nr_gen, local_fc_ref.Extreme_Rays);

    #pragma offload target(mic:mic_nr) in(nr_gen) in(data: length(nr_gen) ONCE)
    {

      offload_fc_ptr->Extreme_Rays = vector<bool>(nr_gen);
      fill_vector(offload_fc_ptr->Extreme_Rays, nr_gen, data);
      offload_fc_ptr->is_Computed.set(ConeProperty::ExtremeRays);
    }
    delete[] data;
  }

  // always transfer the order vector
  {
    Integer *data = new Integer[dim];
    fill_plain(data, dim, local_fc_ref.Order_Vector);

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

  if (!local_fc_ref.Comparisons.empty())
  {
    long size = local_fc_ref.Comparisons.size();

    size_t *data = new size_t[size];
    fill_plain(data, size, local_fc_ref.Comparisons);

    #pragma offload target(mic:mic_nr) in(size) in(data: length(size) ONCE)
    {
      offload_fc_ptr->Comparisons.resize(size);
      fill_vector(offload_fc_ptr->Comparisons, size, data);
      offload_fc_ptr->nrTotalComparisons = offload_fc_ptr->Comparisons[size-1];
    }
    delete[] data;
  }
  cout << "transfer_triangulation_info done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::primal_algorithm_initialize()
{
  #pragma offload target(mic:mic_nr) signal(&running)
  {
    //TODO handle also inhomogeneous, excluded_faces, approx, ...
    offload_fc_ptr->do_vars_check();
    offload_fc_ptr->primal_algorithm_initialize();

    cout << "create 4 mio empty simplices ..." << flush;
    SHORTSIMPLEX<Integer> simp;
    simp.key = vector<key_t>(offload_fc_ptr->dim);
    simp.height = 0;
    simp.vol = 0;
    offload_fc_ptr->FreeSimpl.insert(offload_fc_ptr->FreeSimpl.end(), 4000000, simp);
    cout << "done" << endl;
  }
  running = true;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::transfer_pyramids(const list< vector<key_t> >& pyramids)
{
  cout << "transfer_pyramids" << endl;
  long size = plain_size(pyramids);

  key_t *data = new key_t[size];
  fill_plain(data, size, pyramids);

  wait();
  #pragma offload target(mic:mic_nr) in(size) in(data: length(size) ONCE)
  {
    fill_list_vector(offload_fc_ptr->Pyramids[0], size, data);
    offload_fc_ptr->nrPyramids[0] = offload_fc_ptr->Pyramids[0].size();
  }
  delete[] data;

  cout << "transfer_pyramids done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::evaluate_pyramids()
{
  cout << "evaluate_pyramids" << endl;
  wait();
  #pragma offload target(mic:mic_nr) signal(&running)
  {
    offload_fc_ptr->evaluate_stored_pyramids(0);
  }
  running = true;
  cout << "evaluate_pyramids done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::complete_evaluation()
{
  cout << "complete_evaluation" << endl;
  wait();
  #pragma offload target(mic:mic_nr) signal(&running)
  {
    offload_fc_ptr->primal_algorithm_finalize();
  }
  running = true;
  cout << "complete_evaluation done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::print_on_mic() const
{
  cout << "Offloaded print" << endl;
  #pragma offload target(mic:mic_nr)
  {
    cout << offload_fc_ptr->Pyramids;
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

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::collect_data()
{
  wait();
  cout << "collect_data" << endl;
  collect_integers(); // TriangulationSize, DetSum, Multiplicity, ...
  collect_hilbert_series();
  collect_candidates(); // Hilbert basis, degree 1 elements
  cout << "collect_data done" << endl;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::collect_integers()
{
  {
    size_t col_totalNrSimplices, col_nrSimplicialPyr, col_totalNrPyr;
    Integer col_detSum;
    #pragma offload target(mic:mic_nr) out(col_totalNrSimplices) out(col_nrSimplicialPyr) out(col_totalNrPyr) out(col_detSum)
    {
      col_totalNrSimplices = offload_fc_ptr->totalNrSimplices;
      col_nrSimplicialPyr  = offload_fc_ptr->nrSimplicialPyr;
      col_totalNrPyr       = offload_fc_ptr->totalNrPyr;
      col_detSum           = offload_fc_ptr->detSum;
    }
    local_fc_ref.totalNrSimplices += col_totalNrSimplices;
    local_fc_ref.nrSimplicialPyr  += col_nrSimplicialPyr;
    local_fc_ref.totalNrPyr       += col_totalNrPyr;
    local_fc_ref.detSum           += col_detSum;
  }

  if (local_fc_ref.do_triangulation && local_fc_ref.do_evaluation
      && local_fc_ref.isComputed(ConeProperty::Grading))
  {
    cout << "collecting multiplicity ..." << endl;
    long size;
    std::string* str_ptr;
    #pragma offload target(mic:mic_nr) out(size) nocopy(str_ptr: length(0))
    {
      str_ptr = new std::string(offload_fc_ptr->multiplicity.get_str());
      size = str_ptr->length()+1;
    }

    char* c_str = new char[size];
    #pragma offload target(mic:mic_nr) in(size) out(c_str: length(size)) nocopy(str_ptr: length(0))
    {
      std::strcpy (c_str, str_ptr->c_str());
      delete str_ptr;
    }
    mpq_class coll_mult(c_str);
    delete c_str;
    local_fc_ref.multiplicity += coll_mult;
    cout << "collecting multiplicity done" << endl;
  }
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::collect_hilbert_series()
{

  if (local_fc_ref.do_h_vector)
  {
    cout << "collecting Hilbert series ..." << endl;
    long size;
    std::string* str_ptr;
    #pragma offload target(mic:mic_nr) out(size) nocopy(str_ptr: length(0))
    {
      str_ptr = new std::string(offload_fc_ptr->Hilbert_Series.to_string_rep());
      size = str_ptr->length()+1;
    }

    char* c_str = new char[size];
    #pragma offload target(mic:mic_nr) in(size) out(c_str: length(size)) nocopy(str_ptr: length(0))
    {
      std::strcpy (c_str, str_ptr->c_str());
      delete str_ptr;
    }
    HilbertSeries col_HS = HilbertSeries(string(c_str));
    delete c_str;
    local_fc_ref.Hilbert_Series += col_HS;
    cout << "collecting Hilbert series done" << endl;
  }
}

//---------------------------------------------------------------------------

#pragma offload_attribute (push, target(mic))
template<typename Integer>
bool is_ori_gen(const Candidate<Integer>& c)
{
  return c.original_generator;
}
#pragma offload_attribute (push, target(mic))

//---------------------------------------------------------------------------


template<typename Integer>
void OffloadHandler<Integer>::collect_candidates()
{
  if (local_fc_ref.do_Hilbert_basis)
  {
    cout << "collect Hilbert basis" << endl;
    long size;

    #pragma offload target(mic:mic_nr) out(size)
    {
      // remove all original generators, they are also inserted on the host
      offload_fc_ptr->OldCandidates.Candidates.remove_if(is_ori_gen<Integer>);
      // or in c++11 with lambda function
      // offload_fc_ptr->OldCandidates.Candidates.remove_if([](const Candidate<Integer>& c){ return c.original_generator; });
      offload_fc_ptr->OldCandidates.extract(offload_fc_ptr->Hilbert_Basis);
      offload_fc_ptr->OldCandidates.Candidates.clear();

      // using the same methods as for pyramids
      // handling list of vectors of possible different lenghts
      size = offload_fc_ptr->Hilbert_Basis.size() * (offload_fc_ptr->dim + 1);
    }
    if (size > 0) {
      Integer *data = new Integer[size];

      #pragma offload target(mic:mic_nr) in(size) out(data: length(size) ONCE)
      {
        fill_plain(data, size, offload_fc_ptr->Hilbert_Basis);
      }
//      fill_list_vector(local_fc_ref.Hilbert_Basis, size, data);
      list< vector<Integer> > coll_HB;
      fill_list_vector(coll_HB, size, data);
      delete[] data;
      CandidateList<Integer> cand_l;
      while (!coll_HB.empty())
      {
        cand_l.push_back(Candidate<Integer>(coll_HB.front(),local_fc_ref));
        coll_HB.pop_front();
      }
      cout << "CandidateList complete" << endl;
//      #pragma omp critical(CANDIDATES)
      local_fc_ref.NewCandidates.splice(cand_l);

      local_fc_ref.NewCandidates.reduce_by(local_fc_ref.OldCandidates);
      local_fc_ref.update_reducers();

    } // if (size > 0)
    cout << "collect Hilbert basis done" << endl;
  }

  if (local_fc_ref.do_deg1_elements)
  {
    cout << "collect degree 1 elements" << endl;
    long size;

    #pragma offload target(mic:mic_nr) out(size)
    {
      // using the same methods as for pyramids
      // handling list of vectors of possible different lenghts
      size = offload_fc_ptr->Deg1_Elements.size() * (offload_fc_ptr->dim + 1);
    }
    if (size > 0) {
      Integer *data = new Integer[size];

      #pragma offload target(mic:mic_nr) in(size) out(data: length(size) ONCE)
      {
        fill_plain(data, size, offload_fc_ptr->Deg1_Elements);
      }
      list< vector<Integer> > coll_Deg1;
      fill_list_vector(coll_Deg1, size, data);
      delete[] data;
      local_fc_ref.Deg1_Elements.splice(local_fc_ref.Deg1_Elements.end(),coll_Deg1);
    }
    cout << "collect degree 1 elements done" << endl;
  }
}

//---------------------------------------------------------------------------

template<typename Integer>
bool OffloadHandler<Integer>::is_running()
{
#ifndef __MIC__
  if (running)
  {
    running = ! _Offload_signaled(mic_nr, &running);
  }
#endif
  return running;
}

//---------------------------------------------------------------------------

template<typename Integer>
void OffloadHandler<Integer>::wait()
{
  if (is_running())
  {
    #pragma offload_wait target(mic:mic_nr) wait(&running)
    running = false;
  }
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

//-------------------------- MicOffloader -----------------------------------

template<typename Integer>
MicOffloader<Integer>::MicOffloader()
: is_init(false),
  nr_mic(0),
  handler_ptr(NULL)
{
}

//---------------------------------------------------------------------------

template<typename Integer>
MicOffloader<Integer>::~MicOffloader()
{
  if (is_init) delete handler_ptr;
}

//---------------------------------------------------------------------------

template<typename Integer>
void MicOffloader<Integer>::init(Full_Cone<Integer>& fc)
{
  if (!is_init)
  {
cout << "_Offload_get_device_number()" << _Offload_get_device_number() << endl;
cout << "_Offload_number_of_devices()" << _Offload_number_of_devices() << endl;
    //TODO check preconditions
    assert(fc.Order_Vector.size() == fc.dim);
    fc.get_supphyps_from_copy(false);          // (bool from_scratch)
    fc.check_pointed();

    // create handler
    handler_ptr = new OffloadHandler<Integer>(fc);
    is_init = true;
  }
}

//---------------------------------------------------------------------------

template<typename Integer>
void MicOffloader<Integer>::offload_pyramids(Full_Cone<Integer>& fc)
{
    if (!is_init) init(fc);

    //offload some pyramids //TODO move only a part
    list< vector<key_t> > pyrs;
    size_t nr_transfer = fc.nrPyramids[0]/2;
    //pyrs.splice(pyrs.end(), fc.Pyramids[0]);
    typename list< vector<key_t> >::iterator transfer_end(fc.Pyramids[0].begin());
    for (size_t i = 0; i < nr_transfer; ++i, ++transfer_end) ;
    pyrs.splice(pyrs.end(), fc.Pyramids[0], fc.Pyramids[0].begin(), transfer_end);
    fc.nrPyramids[0] -= nr_transfer;
    handler_ptr->transfer_pyramids(pyrs);
    cout << "offload: transfered " << pyrs.size() << " pyramids to mic." << endl;
    pyrs.clear();

    //compute on mics
    handler_ptr->evaluate_pyramids();
}


//---------------------------------------------------------------------------

template<typename Integer>
void MicOffloader<Integer>::complete_evaluation()
{
  if (is_init)
  {
    handler_ptr->complete_evaluation();
  }
}

//---------------------------------------------------------------------------

template<typename Integer>
void MicOffloader<Integer>::finalize()
{
  if (is_init)
  {
    handler_ptr->complete_evaluation();
    handler_ptr->collect_data();
    delete handler_ptr;
    handler_ptr = NULL;
    is_init = false;
  }
}


/***************** Instantiation for template parameter long long *****************/

template class MicOffloader<long long int>;
template class OffloadHandler<long long int>;

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
  m1.random(10);

  // offload the full cone
  Full_Cone<Integer> fc1(m1);
  fc1.get_supphyps_from_copy(true);          // from_scratch = true
  fc1.Order_Vector = vector<Integer>(b);
  OffloadHandler<Integer> fc1_off(fc1);
  cout << "first offload completed" << endl;
  fc1_off.print_on_mic();

  // work with fc1
  fc1_off.print_on_mic();
  fc1_off.compute_on_mic(1,2);
  fc1_off.compute_on_mic(0,2);

  // get results back
  Matrix<Integer> ret = fc1_off.transfer_from_mic();
  ret.read();

  // destructor of Offload handler will clear data on mic
}

} // end namespace libnormaliz




#endif //NMZ_MIC_OFFLOAD
