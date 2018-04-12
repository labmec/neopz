/* 
 * File:   TPZFunctionTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 16:57
 */

#ifndef TPZFUNCTIONTRANSLATOR_H
#define TPZFUNCTIONTRANSLATOR_H

#include "TPZChunkTranslator.h"

template<class TVar>
class TPZFunctionTranslator : public TPZChunkTranslator{
public:
    TPZFunctionTranslator();
    TPZFunctionTranslator(const TPZFunctionTranslator<TVar>& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZFunctionTranslator();
private:

};

template<class TVar>
TPZFunctionTranslator<TVar>::TPZFunctionTranslator() {
}

template<class TVar>
TPZFunctionTranslator<TVar>::TPZFunctionTranslator(const TPZFunctionTranslator<TVar>& orig) {
}

template<class TVar>
void TPZFunctionTranslator<TVar>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    
}

template<class TVar>
TPZFunctionTranslator<TVar>::~TPZFunctionTranslator() {
}


#endif /* TPZFUNCTIONTRANSLATOR_H */

