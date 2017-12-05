/* 
 * File:   TPZConstrainedNormalRandom.cpp
 * Author: quinelato
 * 
 * Created on 5 de Dezembro de 2017, 15:28
 */

#include "TPZConstrainedNormalRandom.h"

TPZConstrainedNormalRandom::TPZConstrainedNormalRandom(REAL begin, REAL end, REAL mean, REAL stdev) : TPZConstrainedRandom(begin, end), TPZNormalRandom(mean, stdev) {
}

TPZConstrainedNormalRandom::TPZConstrainedNormalRandom(const TPZConstrainedNormalRandom& orig) : TPZConstrainedRandom(orig), TPZNormalRandom(orig) {
}

REAL TPZConstrainedNormalRandom::next() {
    REAL value;
    do {
        value = TPZNormalRandom::next();
    } while (value <= begin || value >= end);
    return value;
}

REAL TPZConstrainedNormalRandom::pdf(REAL x) {
    if (x<begin || x > end) return 0;
    REAL normal_pdf = TPZNormalRandom::pdf(x);
    REAL area = TPZNormalRandom::cdf(end)-TPZNormalRandom::cdf(begin);
    return normal_pdf / area;
}


TPZConstrainedNormalRandom::~TPZConstrainedNormalRandom() {
}

