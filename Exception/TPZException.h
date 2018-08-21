/* 
 * File:   TPZException.h
 * Author: thiago
 *
 * Created on 21 de Agosto de 2018, 12:14
 */

#ifndef TPZEXCEPTION_H
#define TPZEXCEPTION_H

#include <string>

class TPZException {
public:
    TPZException(std::string detail);
    virtual std::string GetDetail() const;
    virtual ~TPZException();
private:
    std::string fDetail;
};

#endif /* TPZEXCEPTION_H */

