/* 
 * File:   TPZConstrainedRandom.cpp
 * Author: quinelato
 * 
 * Created on 5 de Dezembro de 2017, 13:37
 */

#include "TPZConstrainedRandom.h"

TPZConstrainedRandom::TPZConstrainedRandom(REAL begin, REAL end) : TPZRandom(), begin(begin), end(end) {
}

TPZConstrainedRandom::TPZConstrainedRandom(const TPZConstrainedRandom& orig) : TPZRandom(orig), begin(orig.begin), end(orig.end) {
}

TPZConstrainedRandom::~TPZConstrainedRandom() {
}

