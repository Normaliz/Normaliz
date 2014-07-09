#ifdef NMZ_MIC_OFFLOAD


#pragma offload_attribute (push, target(mic))
#include "offload.h"
#include <iostream>
using namespace std;

void fill_matrix(Matrix &M, long rows, long cols, long* data)
{
  for (long i=0; i<rows; i++)
    for (long j=0; j<cols; j++)
    {
      M[i][j] = data[i*cols+j];
    }
}

void fill_plain(long* data, long rows, long cols, const Matrix &M)
{
  for (long i=0; i<rows; i++)
    for (long j=0; j<cols; j++)
      data[i*cols+j] = M[i][j];

}
#pragma offload_attribute (pop)



OffloadHandler::OffloadHandler(const Matrix& M, int mic_number): mic_nr(mic_number)
{
  long nr = M.nr_of_rows();
  long nc = M.nr_of_columns();
  long size = nr*nc;
  long *data = new long[size];
  fill_plain(data, nr, nc, M);

  cout << "Offload data to mic, m_ptr value on cpu " << m_ptr << endl;
  // offload to mic, copy data and free it afterwards, but keep a pointer to the created C++ matrix
  #pragma offload target(mic:mic_nr) in(nr,nc) in(data: length(size) ONCE)
  {
    m_ptr = new Matrix(nr, nc);
    fill_matrix(*m_ptr, nr, nc, data);
    cout << "m_ptr value on mic " << m_ptr << endl;
  }
  cout << "After offload m_ptr value on cpu " << m_ptr << endl;
  delete[] data;
}

void OffloadHandler::print_on_mic() const
{
  cout << "Offloaded print" << endl;
  #pragma offload target(mic:mic_nr)
  {
    m_ptr->pretty_print(cout);
  }
}


void OffloadHandler::compute_on_mic(long a, long b)
{
  cout << "Offloaded computation" << endl;
  #pragma offload target(mic:mic_nr)
  {
    cout << "Rank computed on mic " << m_ptr->rank() << endl;
  }
}


long OffloadHandler::collect_data()
{
//  #pragma offload target(mic:mic_nr) out()
	return 0;
}

Matrix OffloadHandler::transfer_from_mic()
{
  cout << "Transfer data back to cpu" << endl;
  long a,b;
  #pragma offload target(mic:mic_nr) out(a) out(b)
  {
    a = m_ptr->nr_of_rows();
    b = m_ptr->nr_of_columns();
  }

  cout << "a = " << a << "   b = " << b << endl;
  long size = a*b;
  long *ret = new long[size];
  #pragma offload target(mic:mic_nr) out(ret : length(size))
  {
    fill_plain(ret, a, b, *m_ptr);
  }
  Matrix M (a, b);
  fill_matrix(M, a, b, ret);
  delete[] ret;
  return M;
}

OffloadHandler::~OffloadHandler()
{
  cout << "OffloadHandler destructor" << endl;
  #pragma offload target(mic:mic_nr)
  {
    delete m_ptr;
  }
}



/***************** Offload test *****************/

void offload_test()
{
  // initial offload for better timing comparisons of the following offloads
  cout << "initial offload for better timing comparisons of the following offloads"
       << endl;
  #pragma offload target(mic:0)
  { }
  cout << "done." << endl;

  int a = 3, b = 4;
  long SIZE = a*b;
  long data[SIZE];
  for (long i=0; i<SIZE; i++) data[i] = i+1;

  Matrix m1(a,b);
  Matrix m2(b,a);
  m1.random();
  m2.random();

  // offload the matrices
  OffloadHandler m1_off(m1);
  m1_off.print_on_mic();

  OffloadHandler m2_off(m2);
  m2_off.print_on_mic();

  // work with m1
  m1_off.print_on_mic();
  m1_off.compute_on_mic(1,2);
  m1_off.compute_on_mic(0,2);

  // work with m2
  m2_off.print_on_mic();
  m2_off.compute_on_mic(1,2);
  m2_off.compute_on_mic(2,1);

  // get results back
  Matrix ret = m1_off.transfer_from_mic();
  ret.read();

  // destructor of Offload handler will clear data on mic
}
#else
// no offloading available
void offload_test() {}

#endif //NMZ_MIC_OFFLOAD
