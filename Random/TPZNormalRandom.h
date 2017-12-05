/* 
 * File:   TPZNormalRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 15:25
 */

#ifndef TPZNORMALRANDOM_H
#define TPZNORMALRANDOM_H

#include <random>
#include <functional>
#include "TPZRandom.h"


class TPZNormalRandom : public TPZRandom {
public:
    TPZNormalRandom(REAL mean, REAL stdev);
    TPZNormalRandom(const TPZNormalRandom& orig);
    virtual REAL next();
    REAL cdf(REAL x);
    REAL pdf(REAL x);
    virtual ~TPZNormalRandom();
protected:
    REAL mean, stdev;
private :
    std::_Bind_helper<false, std::normal_distribution<REAL>, std::default_random_engine>::type generator;
};

#endif /* TPZNORMALRANDOM_H */

