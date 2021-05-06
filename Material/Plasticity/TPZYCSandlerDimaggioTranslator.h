/* 
 * File:   TPZYCSandlerDimaggioTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 19:41
 */

#ifndef TPZYCSANDLERDIMAGGIOTRANSLATOR_H
#define TPZYCSANDLERDIMAGGIOTRANSLATOR_H

#include "TPZChunkTranslator.h"

class TPZYCSandlerDimaggioTranslator : public TPZChunkTranslator {
public:
    TPZYCSandlerDimaggioTranslator();
    TPZYCSandlerDimaggioTranslator(const TPZYCSandlerDimaggioTranslator& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZYCSandlerDimaggioTranslator();
private:

};

#endif /* TPZYCSANDLERDIMAGGIOTRANSLATOR_H */

