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
  for (int i=0; i<nthreads; i++)
  {
    const int myibeg = ibeg + (iend-ibeg)*i/nthreads;
    const int myiend = ibeg + (iend-ibeg)*(i+1)/nthreads;
    threadvec.push_back(
      std::thread( [myibeg,myiend,&f] ()
      {
        for(int it = myibeg; it < myiend; it++)
          {f(it);}
      }));
  }

  for (int i=0; i<nthreads; i++)
    {threadvec[i].join();}
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
  while (!asum.compare_exchange_weak(current, current + val))
    ;
  return current;
}
}//namespace
#endif
