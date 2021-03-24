/** 
 * @file
 * @brief Creates a free store vector to integers, chars, REAL and TPZString from template TPZManVector.
 */

#include "pzmanvector.h"
#include "pzstring.h" 

template class TPZManVector< int >;
template class TPZManVector< int64_t >;
template class TPZManVector< int *>;
template class TPZManVector< char *>;
template class TPZManVector< float >;
template class TPZManVector< float * >;
template class TPZManVector< double >;
template class TPZManVector< double * >;
template class TPZManVector< long double >;
template class TPZManVector< long double * >;
template class TPZManVector<TPZString>;

