/** 
 * @file 
 * @brief Creates a vectors to integer, REAL, char and TPZString. 
 */

#include "pzvec.h"
#include "pzstring.h" 

template class TPZVec<REAL>;
template class TPZVec<REAL*>;
template class TPZVec<int>;
template class TPZVec<long>;
template class TPZVec<int *>;
template class TPZVec<char *>;
template class TPZVec<void *>;
template class TPZVec<char>;
template class TPZVec<TPZString>;
