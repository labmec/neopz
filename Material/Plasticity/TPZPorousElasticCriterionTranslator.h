//
//  TPZPorousElasticCriterionTranslator.h
//  pz
//
//  Created by Omar Durán on 9/5/19.
//

#ifndef TPZPorousElasticCriterionTranslator_h
#define TPZPorousElasticCriterionTranslator_h

#include <stdio.h>
#include "TPZChunkTranslator.h"
#include "TPZPlasticStateTranslator.h"
#include "TPZElasticResponseTranslator.h"
#include "TPZPorousElasticResponseTranslator.h"

class TPZPorousElasticCriterionTranslator : public TPZChunkTranslator {
    
public:
    
    TPZPorousElasticCriterionTranslator();
    
    TPZPorousElasticCriterionTranslator(const TPZPorousElasticCriterionTranslator& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZPorousElasticCriterionTranslator();
    
private:
    
    TPZPlasticStateTranslator<STATE> tpzPlasticStateTranslatorSTATE;
    
    TPZElasticResponseTranslator tpzElasticResponseTranslator;
    
    TPZPorousElasticResponseTranslator tpzPorousElasticResponseTranslator;
};

#endif /* TPZPorousElasticCriterionTranslator_h */
