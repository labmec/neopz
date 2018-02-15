/* 
 * File:   TPZStochasticSearch.h
 * Author: quinelato
 *
 * Created on 12 de Dezembro de 2017, 15:15
 */

#ifndef TPZSTOCHASTICSEARCH_H
#define TPZSTOCHASTICSEARCH_H

#include <functional>
#include "pzvec.h"
#include "tpzautopointer.h"
#include "TPZRandom.h"

template <typename TVar>
class TPZStochasticSearch {
public:
    TPZStochasticSearch(const unsigned long n_vars) {
        fdistributions.resize(n_vars);
    }
    
    TPZStochasticSearch(const TPZStochasticSearch& orig) : fdistributions(orig.fdistributions) {
    }
    
    void SetDistribution(const unsigned long var, TPZRandom<TVar>& distribution);
    
    TPZVec<TVar> DoSearch(std::function<TVar(TPZVec<TVar>)> objective_function, const unsigned long n_max_iterations, TVar min_relative_error);
    
    virtual ~TPZStochasticSearch(){
        
    }
private:
    TPZVec<TPZAutoPointer<TPZRandom<TVar>>> fdistributions;
};

template <typename TVar>
void TPZStochasticSearch<TVar>::SetDistribution(const unsigned long var, TPZRandom<TVar>& distribution){
#ifdef PZDEBUG
    if (var >= fdistributions.size()){
        DebugStop();
    }
#endif
    TPZAutoPointer<TPZRandom<TVar>> dist(distribution.clone());
    fdistributions[var] = dist;
}

template <typename TVar>
TPZVec<TVar> TPZStochasticSearch<TVar>::DoSearch(std::function<TVar(TPZVec<TVar>)> objective_function, const unsigned long n_max_iterations, TVar min_relative_error) {
    const unsigned long nvars = fdistributions.size();
    TPZVec<TVar> guess(nvars);
    for (unsigned int i = 0; i < nvars; ++i) {
#ifdef PZDEBUG
        if (!fdistributions[i].operator ->()){
            DebugStop();
        }
#endif
        guess[i] = fdistributions[i]->next();
    }

    unsigned long n_iterations = 0;
    TPZVec<TVar> best_guess = guess;
    TVar best_score = objective_function(best_guess);
    TVar relative_error = std::numeric_limits<TVar>::max();
    while (n_iterations < n_max_iterations && relative_error > min_relative_error) {
        for (unsigned int i = 0; i < nvars; ++i) {
            guess[i] = fdistributions[i]->next();
        }
        TVar score = objective_function(guess);
        if (score < best_score){
            relative_error = (best_score-score)/best_score;
            best_score = score;
            best_guess = guess;
        }
        ++n_iterations;
    }
    return best_guess;
}


#endif /* TPZSTOCHASTICSEARCH_H */

