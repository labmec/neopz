//
//  TPZDiscontinuousGalerkinTranslator.h
//  pz
//
//  Created by Omar Durán on 9/26/19.
//

#ifndef TPZDiscontinuousGalerkinTranslator_h
#define TPZDiscontinuousGalerkinTranslator_h

#include <stdio.h>
#include "TPZMaterialTranslator.h"
#include "TPZChunkInTranslation.h"

class TPZDiscontinuousGalerkinTranslator : public TPZMaterialTranslator {

public:
    
    TPZDiscontinuousGalerkinTranslator();
    
    TPZDiscontinuousGalerkinTranslator(const TPZDiscontinuousGalerkinTranslator & other);
    
    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    
    virtual ~TPZDiscontinuousGalerkinTranslator();
    
private:
    
    TPZMaterialTranslator parentTranslator;
    
};

#endif /* TPZDiscontinuousGalerkinTranslator_h */
