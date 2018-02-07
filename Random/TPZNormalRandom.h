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
#include <ctime>
#include "TPZRandom.h"

template <typename TVar>
class TPZNormalRandom : virtual public TPZRandom<TVar> {
public:
    TPZNormalRandom(TVar mean, TVar stdev) : mean(mean), stdev(stdev), generator(std::bind(std::normal_distribution<TVar>(mean, stdev), std::default_random_engine(clock()))){
        
    }
    
    TPZNormalRandom(const TPZNormalRandom& orig): mean(orig.mean), stdev(orig.stdev), generator(orig.generator){
        
    }
    
    virtual TPZRandom<TVar> *clone(){
        return new TPZNormalRandom<TVar>(*this);
    }
    
    virtual TVar next(){
        return generator();
    }
    
    TVar cdf(TVar x){
        return 0.5*(1+std::erf((x-mean)/(stdev*M_SQRT2)));
    }
    
    TVar pdf(TVar x){
        TVar stdev_2 = pow(stdev,2.);
        return 1./(sqrt(2*M_PI*stdev_2))*exp(-pow(x-mean,2)/(2*stdev_2));
    }
    virtual ~TPZNormalRandom(){
        
    }
protected:
    TVar mean, stdev;
private :
    std::function<TVar()> generator;
};

#endif /* TPZNORMALRANDOM_H */

