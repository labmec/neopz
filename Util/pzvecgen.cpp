/** 
 * @file 
 * @brief Creates a vectors to integer, REAL, char and TPZString. 
 */

#include "pzvec.h"
#include "pzstring.h" 

template class TPZVec<float>;
template class TPZVec<float * >;
template class TPZVec<double>;
template class TPZVec<double * >;
template class TPZVec<long double>;
template class TPZVec<long double * >;
template class TPZVec<int>;
template class TPZVec<int64_t>;
template class TPZVec<int64_t *>;
template class TPZVec<int *>;
template class TPZVec<char *>;
template class TPZVec<void *>;
template class TPZVec<char>;
template class TPZVec<TPZString>;
