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
class TPZFYsmpMatrix;
class TPZFMatrix;
class TPZMatrix;

void RunFast();
// this function will read the matrix
TPZFYsmpMatrix *ReadMatrix(const std::string &filename, TPZFMatrix &rhs);

void TimeMultiply(TPZMatrix *mat, TPZMultiTimer &timer);

void SolveJacobi(TPZMatrix *mat, TPZFMatrix &rhs, REAL tol, TPZMultiTimer &timer);

void SolveCG(TPZMatrix *mat, TPZFMatrix &rhs, REAL tol, TPZMultiTimer &timer);

void SolveSSOR(TPZMatrix *mat, TPZFMatrix &rhs, REAL tol, TPZMultiTimer &timer);

void Compare(TPZMatrix *first, TPZMatrix *second);

#endif
