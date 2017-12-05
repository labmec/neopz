/* 
 * File:   TPZConstrainedNormalRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 15:28
 */

#ifndef TPZCONSTRAINEDNORMALRANDOM_H
#define TPZCONSTRAINEDNORMALRANDOM_H

#include "TPZConstrainedRandom.h"
#include "TPZNormalRandom.h"

class TPZConstrainedNormalRandom : public TPZConstrainedRandom, public TPZNormalRandom {
public:
    TPZConstrainedNormalRandom(REAL begin, REAL end, REAL mean, REAL stdev);
    TPZConstrainedNormalRandom(const TPZConstrainedNormalRandom& orig);
    virtual REAL next();
    REAL pdf(REAL x);
    virtual ~TPZConstrainedNormalRandom();

};

#endif /* TPZCONSTRAINEDNORMALRANDOM_H */

