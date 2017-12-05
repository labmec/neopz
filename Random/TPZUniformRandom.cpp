/* 
 * File:   TPZUniformRandom.cpp
 * Author: quinelato
 * 
 * Created on 5 de Dezembro de 2017, 13:31
 */

#include <functional>
#include <random>

#include "TPZUniformRandom.h"

TPZUniformRandom::TPZUniformRandom(REAL begin, REAL end) : TPZConstrainedRandom(begin, end), generator(std::bind(std::uniform_real_distribution<REAL>(begin, end), std::default_random_engine())) {
}

TPZUniformRandom::TPZUniformRandom(const TPZUniformRandom& orig) : TPZConstrainedRandom(orig), generator(orig.generator) {
}

REAL TPZUniformRandom::next(){
    return generator();
}

REAL TPZUniformRandom::pdf(REAL x) {
    return 1./(end-begin);
}

TPZUniformRandom::~TPZUniformRandom() {
}

