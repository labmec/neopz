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

void TPZPlasticBase::ResetPlasticStrain() {
    auto state = this->GetState();
    state.m_eps_p.Zero();
    SetState(state);
}
