#ifndef OFFLOAD_H
#define OFFLOAD_H

#ifdef NMZ_MIC_OFFLOAD
#define ONCE alloc_if(1) free_if(1)
#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)

#include "libnormaliz.h"
#include <list>
#include <vector>

namespace libnormaliz {

// forward-declarations
template<typename Integer> class Full_Cone;
template<typename Integer> class Matrix;


template<typename Integer>
class OffloadHandler
{
public:
  // transfers Full_Cone to mic:mic_nr and keeps handle
  OffloadHandler(Full_Cone<Integer>&, int mic_number = 0);

  // destructor deletes Full_Cone on mic
  ~OffloadHandler();

  void transfer_pyramids(const std::list< std::vector<key_t> >& pyramids);
  void evaluate_pyramids();
  void finalize_evaluation();

  // prints the offloaded generators
  void print_on_mic() const;

  // does a test computation on the mic and prints the result
  void compute_on_mic(long a, long b);

  // transfer the generator matrix back
  Matrix<Integer> transfer_from_mic();

private:
  const int mic_nr;
  Full_Cone<Integer>& local_fc_ref;
  Full_Cone<Integer>* offload_fc_ptr;

  // init routines
  void create_full_cone();
  void transfer_bools();
  void transfer_support_hyperplanes();
  void transfer_grading();
  void transfer_triangulation_info();
  void primal_algorithm_initialize();

  // collect data routines
  void collect_data();
  void collect_integers(); // TriangulationSize, DetSum, Multiplicity
  void collect_hilbert_series();
  void collect_candidates(); // Hilbert basis, degree 1 elements

};

template<typename Integer>
class MicOffloader
{
public:
  MicOffloader();
  ~MicOffloader();

  void offload_pyramids(Full_Cone<Integer>& fc);

private:
  bool is_init;
  int nr_mic;
  OffloadHandler<Integer>* handler_ptr;

  void init(Full_Cone<Integer>& fc);
};

void offload_test();

} // end namespace libnormaliz

#endif //NMZ_MIC_OFFLOAD

#endif //OFFLOAD_H
