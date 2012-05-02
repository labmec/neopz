/** 
 * @file 
 * @brief Contains the implementation of the methods to TPZFunction class. 
 */

#include "pzfunction.h"

template<class TVar>
TPZFunction<TVar>::TPZFunction()
{
}

template<class TVar>
TPZFunction<TVar>::~TPZFunction()
{
}

template class TPZFunction<float>;
template class TPZFunction<double>;
template class TPZFunction<long double>;

template class TPZFunction<std::complex<float> >;
template class TPZFunction<std::complex<double> >;
template class TPZFunction<std::complex<long double> >;
