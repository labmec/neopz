/**
 * @file TPZParallelUtils.h
 *
 * @brief Defines useful utility functions related to multithreading for use in the NeoPZ library
 *
 * @ingroup util
 * (Note: this needs exactly one @defgroup somewhere)
 * @date Mar 2021
 * @author Francisco Orlandini
 * Contact: francisco.orlandini@gmail.com
 *
 */
#ifndef PZPARALLELUTILS_H
#define PZPARALLELUTILS_H
#include <vector>
#include <thread>
#include <atomic>

namespace pzutils{
//!Sets number of mkl threads for current threads and return previous value
int SetNumThreadsLocalMKL(const int nt); //does nothing if not using mkl
/**based on: https://ideone.com/Z7zldb and netgen source code*/

/** This function aims to enable parallel for loops from first <= i < last
as in 
#pragma omp parallelfor.
It is the user responsability to use the necessary mutexes and locks.
*/
template<typename TFunc>
void ParallelFor( int ibeg, int iend, const TFunc & f )
{
  const auto nthreads = std::thread::hardware_concurrency();
  std::vector<std::thread> threadvec;
  const int sz = (iend-ibeg)/nthreads;
  for (int i=0; i<nthreads; i++)
  {
    const int myibeg = ibeg + i*sz;
    const int myiend = ibeg + (i+1)*sz;
    threadvec.push_back(
      std::thread( [myibeg,myiend,&f] ()
      {
        //we prevent mkl launching multiple threads inside this thread
        auto mklthreads = SetNumThreadsLocalMKL(1);
        for(int it = myibeg; it < myiend; it++)
          {f(it);}
        SetNumThreadsLocalMKL(mklthreads);
      }));
  }

  for (int i=0; i<nthreads; i++)
    {threadvec[i].join();}
  //possible remainder
  const int lastbeg = ibeg + nthreads*sz;
  for(int it = lastbeg; it < iend; it++){
    f(it);
  }
}


/*based on
https://stackoverflow.com/questions/48794815/replacing-pragma-omp-atomic-with-c-atomics
and netgen source code
*/

template<typename T>
inline std::atomic<T> & AtomicCast (T & d)
{
  return reinterpret_cast<std::atomic<T>&> (d);
}
/** This function aims to enable atomic addition operation
as in 
#pragma omp atomic
a += b
*/
template<typename T>
inline T AtomicAdd( T & sum, T val )
{
  std::atomic<T> & asum = AtomicCast(sum);
  T current = asum.load();
  T desired{0};
  do{
    desired = current+val;
  }
  while (!asum.compare_exchange_weak(current, desired));
  return current;
}
}//namespace
#endif
