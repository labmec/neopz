/* 
 * File:   TPZMaterialTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 18:51
 */

#ifndef TPZMATERIALTRANSLATOR_H
#define TPZMATERIALTRANSLATOR_H

#include "TPZChunkTranslator.h"

class TPZMaterialTranslator : public TPZChunkTranslator {
public:
    TPZMaterialTranslator();
    TPZMaterialTranslator(const TPZMaterialTranslator& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) override;

    virtual ~TPZMaterialTranslator();
private:

};

#endif /* TPZMATERIALTRANSLATOR_H */

