//
//  TPZNullMaterialTranslator.h
//  pz
//
//  Created by Omar Durán on 9/26/19.
//

#ifndef TPZNullMaterialTranslator_h
#define TPZNullMaterialTranslator_h

#include <stdio.h>
#include "TPZMaterialTranslator.h"
#include "TPZChunkTranslator.h"

class TPZNullMaterialTranslator : public TPZMaterialTranslator {
    
public:
    
    TPZNullMaterialTranslator();
    
    TPZNullMaterialTranslator(const TPZNullMaterialTranslator & other);
    
    virtual void UpdateAttributes(TPZChunkInTranslation & chunk, const std::map<std::string, uint64_t> & toVersion) override;
    
    virtual ~TPZNullMaterialTranslator();
    
private:
    
    TPZMaterialTranslator parentTranslator;
    
};

#endif /* TPZNullMaterialTranslator_h */
