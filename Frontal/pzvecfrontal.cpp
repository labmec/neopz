/**
 * @file
 * @brief Creates vectors of TPZEqnArray and TPZEqnArray* .
 */

#include "pzvec.h"
#include "pzerror.h"

#include "tpzeqnarray.h"

template class TPZVec<TPZEqnArray<float> >;
template class TPZVec<TPZEqnArray<double> >;
template class TPZVec<TPZEqnArray<long double> >;

template class TPZVec<TPZEqnArray<float> *>;
template class TPZVec<TPZEqnArray<double> *>;
template class TPZVec<TPZEqnArray<long double> *>;
