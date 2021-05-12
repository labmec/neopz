/* 
 * File:   TPZPlasticCriterium.h
 * Author: thiago
 *
 * Created on 28 de Maio de 2018, 18:31
 */

#ifndef TPZPLASTICCRITERION_H
#define TPZPLASTICCRITERION_H

#include "pzreal.h"
#include "TPZSavable.h"
#include "pzvec.h"

class TPZPlasticCriterion : public TPZSavable {
public:

    virtual void YieldFunction(const TPZVec<STATE> &sigma, STATE kprev, TPZVec<STATE> &yield) const = 0;

    virtual int GetNYield() const = 0;
    
    virtual void Print(std::ostream &out) const {
        std::cout << __PRETTY_FUNCTION__ << " Should not be called, please check children classes." << std::endl;
//        DebugStop();
    }
};

#endif /* TPZPLASTICCRITERION_H */

