/* 
 * File:   TPZInconsistentStateException.cpp
 * Author: thiago
 * 
 * Created on 21 de Agosto de 2018, 12:13
 */

#include "TPZInconsistentStateException.h"

TPZInconsistentStateException::TPZInconsistentStateException(std::string detail) : TPZException(detail) {
}

TPZInconsistentStateException::~TPZInconsistentStateException() {
}

