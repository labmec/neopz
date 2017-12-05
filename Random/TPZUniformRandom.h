/* 
 * File:   TPZUniformRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 13:31
 */

#ifndef TPZUNIFORMRANDOM_H
#define TPZUNIFORMRANDOM_H

#include <random>
#include <functional>

#include "TPZConstrainedRandom.h"


class TPZUniformRandom : public TPZConstrainedRandom {
public:
    TPZUniformRandom(REAL begin, REAL end);
    TPZUniformRandom(const TPZUniformRandom& orig);
    REAL next();
    REAL pdf(REAL x);
    virtual ~TPZUniformRandom();
protected:
    std::_Bind_helper<false, std::uniform_real_distribution<REAL>, std::default_random_engine>::type generator;
};

#endif /* TPZUNIFORMRANDOM_H */

