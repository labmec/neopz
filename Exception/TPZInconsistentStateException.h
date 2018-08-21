/* 
 * File:   TPZInconsistentStateException.h
 * Author: thiago
 *
 * Created on 21 de Agosto de 2018, 12:13
 */

#ifndef TPZINCONSISTENTSTATEEXCEPTION_H
#define TPZINCONSISTENTSTATEEXCEPTION_H

#include <string>

#include "TPZException.h"

class TPZInconsistentStateException : public TPZException {
public:
    TPZInconsistentStateException(std::string detail);
    virtual ~TPZInconsistentStateException();
private:

};

#endif /* TPZINCONSISTENTSTATEEXCEPTION_H */

