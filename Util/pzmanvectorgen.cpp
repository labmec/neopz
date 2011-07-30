// Generated by Together

/** 
 * @file
 * @brief Creates a free store vector to integers, chars, REAL and TPZString.
 */
// $Id: pzmanvectorgen.cpp,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#include "pzmanvector.h"
#include "pzstring.h" 

template class TPZManVector< int >;
template class TPZManVector< long int >;
template class TPZManVector< int *>;
template class TPZManVector< char *>;
template class TPZManVector< REAL >;
template class TPZManVector< REAL * >;
template class TPZManVector<TPZString>;

//--| PZ |----------------------------------------------------------------------
