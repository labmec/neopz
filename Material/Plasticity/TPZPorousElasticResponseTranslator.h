//
//  TPZPorousElasticResponseTranslator.h
//  pz
//
//  Created by Omar Durán on 9/5/19.
//

#ifndef TPZPorousElasticResponseTranslator_h
#define TPZPorousElasticResponseTranslator_h

#include <stdio.h>
#include "TPZChunkTranslator.h"

class TPZPorousElasticResponseTranslator : public TPZChunkTranslator {

public:
    
    TPZPorousElasticResponseTranslator();
    
    TPZPorousElasticResponseTranslator(const TPZPorousElasticResponseTranslator & other);
    
    virtual void UpdateAttributes(TPZChunkInTranslation & chunk, const std::map<std::string, uint64_t> & toVersion);
    
    virtual ~TPZPorousElasticResponseTranslator();
    
};

#endif /* TPZPorousElasticResponseTranslator_h */
