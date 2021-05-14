#ifndef TPZMATTYPES_H
#define TPZMATTYPES_H

#include "pzreal.h"
#include <functional>

template <class T>
class TPZVec;
template<class T>
class TPZFMatrix;


//! Alias for forcing function type
template<class TVar>
using ForcingFunctionType =  std::function<void (const TPZVec<REAL> &loc,
                                                 TPZVec<TVar> &result)>;
//! Alias for boundary condition forcing function type
template<class TVar>
using ForcingFunctionBCType =  std::function<void (const TPZVec<REAL> &loc,
                                                   TPZVec<TVar> &rhsVal,
                                                   TPZFMatrix<TVar> &matVal)>;
#endif