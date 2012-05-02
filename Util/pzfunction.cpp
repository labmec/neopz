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

template class TPZFunction<double>;

template class TPZFunction<std::complex<double> >;
