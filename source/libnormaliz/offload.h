#ifndef OFFLOAD_H
#define OFFLOAD_H

#ifdef NMZ_MIC_OFFLOAD
#define ONCE alloc_if(1) free_if(1)
#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)

#include "matrix.h"
typedef libnormaliz::Matrix<long> Matrix;

class OffloadHandler
{
public:
  OffloadHandler(const Matrix&, int mic_number = 0); // transfers Matrix to mic:mic_nr and keeps handle
  ~OffloadHandler();                              // destructor deletes Matrix on mic
  void print_on_mic() const;
  void compute_on_mic(long a, long b);
  long collect_data();

  Matrix transfer_from_mic();

private:
  const int mic_nr;
  Matrix* m_ptr;
};

#endif //NMZ_MIC_OFFLOAD

void offload_test();

#endif //OFFLOAD_H
