/**
 * @file
 * @brief Creates vectors of TPZEqnArray and TPZEqnArray* .
 */
#include "pzvec.h"
#include "pzerror.h"


#include "tpzeqnarray.h"
template class TPZVec<TPZEqnArray>;
template class TPZVec<TPZEqnArray *>;
