/** 
 * @file 
 * @brief Creates a stack vector to integers, REAL and chars. 
 */

#include "pzstack.h"

template class TPZStack<int>;
template class TPZStack<float>;
template class TPZStack<double>;
template class TPZStack<long double>;
template class TPZStack<char *>;
template class TPZStack<int64_t>;
template class TPZStack<char>;

template class TPZStack<std::complex<float> >;
template class TPZStack<std::complex<double> >;

class TPZCompEl;
class TPZGeoEl;

template class TPZStack<TPZCompEl *>;
template class TPZStack<TPZGeoEl *>;
#include "pzcompel.h"
#include "pzgeoelside.h"
template class TPZStack<TPZCompElSide>;
template class TPZStack<TPZGeoElSide>;