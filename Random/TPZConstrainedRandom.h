/* 
 * File:   TPZConstrainedRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 13:37
 */

#ifndef TPZCONSTRAINEDRANDOM_H
#define TPZCONSTRAINEDRANDOM_H

#include "TPZRandom.h"

class TPZConstrainedRandom : public TPZRandom {
public:
    TPZConstrainedRandom(REAL begin, REAL end);
    TPZConstrainedRandom(const TPZConstrainedRandom& orig);
    virtual ~TPZConstrainedRandom();
protected :
    REAL begin, end;
};

#endif /* TPZCONSTRAINEDRANDOM_H */

