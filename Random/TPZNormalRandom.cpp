/* 
 * File:   TPZNormalRandom.cpp
 * Author: quinelato
 * 
 * Created on 5 de Dezembro de 2017, 15:25
 */

#include "TPZNormalRandom.h"
#include <random>
#include <functional>

TPZNormalRandom::TPZNormalRandom(REAL mean, REAL stdev) : mean(mean), stdev(stdev), generator(std::bind(std::normal_distribution<REAL>(mean, stdev), std::default_random_engine())) {
    
}

TPZNormalRandom::TPZNormalRandom(const TPZNormalRandom& orig) : mean(orig.mean), stdev(orig.stdev), generator(orig.generator) {
}

REAL TPZNormalRandom::next() {
    return generator();
}

REAL TPZNormalRandom::cdf(REAL x) {
    return 0.5*(1+std::erf((x-mean)/(stdev*M_SQRT2)));
}

REAL TPZNormalRandom::pdf(REAL x) {
    REAL stdev_2 = pow(stdev,2.);
    return 1./(sqrt(2*M_PI*stdev_2))*exp(-pow(x-mean,2)/(2*stdev_2));
}

TPZNormalRandom::~TPZNormalRandom() {
}

