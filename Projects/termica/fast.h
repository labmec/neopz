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
#include "pzysmp.h"

void RunFast();
// this function will read the matrix
TPZFYsmpMatrix<STATE> *ReadMatrix(const std::string &filename, TPZFMatrix<STATE> &rhs);

void TimeMultiply(TPZMatrix<STATE> *mat, TPZMultiTimer &timer);

void SolveJacobi(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> &rhs, REAL tol, TPZMultiTimer &timer);

void SolveCG(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> &rhs, REAL tol, TPZMultiTimer &timer);

void SolveSSOR(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> &rhs, REAL tol, TPZMultiTimer &timer);

void Compare(TPZMatrix<REAL> *first, TPZMatrix<REAL> *second);

#endif
