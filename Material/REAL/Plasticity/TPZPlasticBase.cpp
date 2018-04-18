/* 
 * File:   TPZPlasticBase.cpp
 * Author: quinelato
 * 
 * Created on 30 de Outubro de 2017, 09:32
 */

#include "TPZPlasticBase.h"

int TPZPlasticBase::ClassId() const {
    return Hash("TPZPlasticBase");
}

TPZPlasticBase::~TPZPlasticBase() {

}

int TPZPlasticBase::IntegrationSteps() const {
    return 1;
}
