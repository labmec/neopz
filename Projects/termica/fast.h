//
// C++ Interface: fast
//
// Description: 
//
//
// Author: Philippe Remy Bernard Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef FASTH
#define FASTH

#include <string>
#include "pzreal.h"
#include "TPZTimer.h"
template<class TVar>
class TPZFYsmpMatrix;

template<class TVar>
class TPZFMatrix;
template<class TVar>
class TPZMatrix;

void RunFast();
// this function will read the matrix
TPZFYsmpMatrix<REAL> *ReadMatrix(const std::string &filename, TPZFMatrix<REAL> &rhs);

void TimeMultiply(TPZMatrix<REAL> *mat, TPZMultiTimer &timer);

void SolveJacobi(TPZMatrix<REAL> *mat, TPZFMatrix<REAL> &rhs, REAL tol, TPZMultiTimer &timer);

void SolveCG(TPZMatrix<REAL> *mat, TPZFMatrix<REAL> &rhs, REAL tol, TPZMultiTimer &timer);

void SolveSSOR(TPZMatrix<REAL> *mat, TPZFMatrix<REAL> &rhs, REAL tol, TPZMultiTimer &timer);

void Compare(TPZMatrix<REAL> *first, TPZMatrix<REAL> *second);

#endif
