/** 
 * @file
 * @brief Creates a free store vector to integers, chars, REAL and TPZString from template TPZManVector.
 */

#include "pzmanvector.h"
#include "pzstring.h" 

template class TPZManVector< int >;
template class TPZManVector< long int >;
template class TPZManVector< int *>;
template class TPZManVector< char *>;
template class TPZManVector< REAL >;
template class TPZManVector< REAL * >;
template class TPZManVector<TPZString>;

