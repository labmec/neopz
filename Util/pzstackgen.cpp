/** 
 * @file 
 * @brief Creates a stack vector to integers, REAL and chars. 
 */
// $Id: pzstackgen.cpp,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#include "pzstack.h"

template class TPZStack<int>;
template class TPZStack<REAL>;
template class TPZStack<char *>;
template class TPZStack<long int>;
template class TPZStack<char>;

//--| PZ |----------------------------------------------------------------------
