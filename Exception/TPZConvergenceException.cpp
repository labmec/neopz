/* 
 * File:   TPZConvergenceException.cpp
 * Author: thiago
 * 
 * Created on 21 de Agosto de 2018, 11:57
 */

#include "TPZConvergenceException.h"

TPZConvergenceException::TPZConvergenceException(REAL tolerance, unsigned int maxIterations, REAL error, unsigned int itNumber, std::string detail) : TPZException(detail),
fTolerance(tolerance), fMaxIterations(maxIterations), fError(error), fItNumber(itNumber) {
}

TPZConvergenceException::~TPZConvergenceException() {
}

