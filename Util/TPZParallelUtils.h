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

#endif
