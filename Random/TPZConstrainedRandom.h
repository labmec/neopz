/* 
 * File:   TPZConstrainedRandom.h
 * Author: quinelato
 *
 * Created on 5 de Dezembro de 2017, 13:37
 */

#ifndef TPZCONSTRAINEDRANDOM_H
#define TPZCONSTRAINEDRANDOM_H

#include "TPZRandom.h"

template <class TVar>
class TPZConstrainedRandom : virtual public TPZRandom<TVar> {
public:
    TPZConstrainedRandom(TVar begin, TVar end) : TPZRandom<TVar>(), fbegin(begin), fend(end){    
    }
    TPZConstrainedRandom(const TPZConstrainedRandom<TVar>& orig) : TPZRandom<TVar>(orig), fbegin(orig.fbegin), fend(orig.fend) {
    }
    
    virtual TVar GetBegin() const {
        return fbegin;
    }

    virtual TVar GetEnd() const {
        return fend;
    }

    virtual ~TPZConstrainedRandom(){
        
    }
protected :
    TVar fbegin, fend;
};

#endif /* TPZCONSTRAINEDRANDOM_H */

