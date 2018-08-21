/* 
 * File:   TPZConvergenceException.h
 * Author: thiago
 *
 * Created on 21 de Agosto de 2018, 11:57
 */

#ifndef TPZCONVERGENCEEXCEPTION_H
#define TPZCONVERGENCEEXCEPTION_H

#include "pzreal.h"
#include "TPZException.h"
#include <string>

class TPZConvergenceException : public TPZException {
public:
    TPZConvergenceException(REAL tolerance, unsigned int maxIterations, REAL error, unsigned int itNumber, std::string detail = "");
    virtual ~TPZConvergenceException();
private:
    REAL fTolerance;
    unsigned int fMaxIterations;
    REAL fError;
    unsigned int fItNumber;
};

#endif /* TPZCONVERGENCEEXCEPTION_H */

