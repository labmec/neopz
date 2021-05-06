/* 
 * File:   TPZSandlerDimaggioThermoForceATranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 19:45
 */

#ifndef TPZSANDLERDIMAGGIOTHERMOFORCEATRANSLATOR_H
#define TPZSANDLERDIMAGGIOTHERMOFORCEATRANSLATOR_H

#include "TPZChunkTranslator.h"

class TPZSandlerDimaggioThermoForceATranslator : public TPZChunkTranslator {
public:
    TPZSandlerDimaggioThermoForceATranslator();
    TPZSandlerDimaggioThermoForceATranslator(const TPZSandlerDimaggioThermoForceATranslator& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZSandlerDimaggioThermoForceATranslator();
private:

};

#endif /* TPZSANDLERDIMAGGIOTHERMOFORCEATRANSLATOR_H */

