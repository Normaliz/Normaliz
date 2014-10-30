#ifndef OFFLOAD_H
#define OFFLOAD_H

#ifdef NMZ_MIC_OFFLOAD
#define ONCE alloc_if(1) free_if(1)
#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)
#endif //NMZ_MIC_OFFLOAD

namespace libnormaliz {

#ifdef NMZ_MIC_OFFLOAD

template<typename Integer> class Full_Cone; //forward-declaration
template<typename Integer> class Matrix; //forward-declaration

template<typename Integer>
class OffloadHandler
{
public:
  // transfers Full_Cone to mic:mic_nr and keeps handle
  OffloadHandler(Full_Cone<Integer>&, int mic_number = 0);

  // destructor deletes Full_Cone on mic
  ~OffloadHandler();

  // prints the offloaded generators
  void print_on_mic() const;

  // does a test computation on the mic and prints the result
  void compute_on_mic(long a, long b);

  // gets a result back from the mic
  long collect_data();

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

};

#endif //NMZ_MIC_OFFLOAD

void offload_test();

} // end namespace libnormaliz

#endif //OFFLOAD_H
