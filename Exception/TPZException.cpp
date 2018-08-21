/* 
 * File:   TPZException.cpp
 * Author: thiago
 * 
 * Created on 21 de Agosto de 2018, 12:14
 */

#include "TPZException.h"

TPZException::TPZException(std::string detail) : fDetail(detail) {
}

std::string TPZException::GetDetail() const {
    return fDetail;
}

TPZException::~TPZException() {
}

