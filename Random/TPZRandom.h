/* 
 * File:   TPZRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 11:03
 */

#ifndef TPZRANDOM_H
#define TPZRANDOM_H

#include "pzreal.h"

class TPZRandom {
public:
    TPZRandom();
    TPZRandom(const TPZRandom& orig);
    virtual REAL next() = 0;
    virtual REAL pdf(REAL x) = 0;
    virtual ~TPZRandom();
protected :
};

#endif /* TPZRANDOM_H */

