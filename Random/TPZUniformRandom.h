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
#include <ctime>

#include "TPZConstrainedRandom.h"

template <typename TVar>
class TPZUniformRandom : public TPZConstrainedRandom<TVar> {
public:
    TPZUniformRandom(TVar begin, TVar end) : TPZConstrainedRandom<TVar>(begin, end), generator(std::bind(std::uniform_real_distribution<TVar>(begin, end), std::default_random_engine(clock()))) {
}
    TPZUniformRandom(const TPZUniformRandom<TVar>& orig): TPZConstrainedRandom<TVar>(orig), generator(orig.generator){
        
    }
    
    virtual TPZRandom<TVar> *clone() override {
        return new TPZUniformRandom<TVar>(*this);
    }
    
    TVar next() override {
        return generator();
    }
    
    TVar pdf(TVar x) override;
    
    virtual ~TPZUniformRandom(){
        
    }
protected:
    std::function<TVar()> generator;
};

template <typename TVar>
TVar TPZUniformRandom<TVar>::pdf(TVar x) {
    return 1./(this->fend - this->fbegin);
}

#endif /* TPZUNIFORMRANDOM_H */

