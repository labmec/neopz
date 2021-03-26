/** 
 * @file
 * @brief Creates a free store vector to integers, chars, REAL and TPZString from template TPZManVector.
 */

#include "pzmanvector.h"


template class TPZManVector< float,3 >;
template class TPZManVector< double,3 >;
template class TPZManVector< long double,3 >;
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

