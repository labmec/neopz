/* 
 * File:   TPZRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 11:03
 */

#ifndef TPZRANDOM_H
#define TPZRANDOM_H

#include "pzreal.h"

template <typename TVar>
class TPZRandom {
public:
    TPZRandom(){
        
    }
    
    TPZRandom(const TPZRandom<TVar>& orig){
        
    }
    
    virtual TPZRandom<TVar> *clone() = 0;
    virtual TVar next() = 0;
    virtual TVar pdf(TVar x) = 0;
    virtual ~TPZRandom(){
        
    }
protected :
};

#endif /* TPZRANDOM_H */

