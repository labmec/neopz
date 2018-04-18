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
