/* 
 * File:   TPZDruckerPragerTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:19
 */

#ifndef TPZDRUCKERPRAGERTRANSLATOR_H
#define TPZDRUCKERPRAGERTRANSLATOR_H

#include "TPZChunkTranslator.h"

class TPZDruckerPragerTranslator : public TPZChunkTranslator {
public:
    TPZDruckerPragerTranslator();
    TPZDruckerPragerTranslator(const TPZDruckerPragerTranslator& orig);
    
    virtual void UpdateAttributes(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion);
    
    virtual ~TPZDruckerPragerTranslator();
private:

};

#endif /* TPZDRUCKERPRAGERTRANSLATOR_H */

