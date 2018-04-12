/* 
 * File:   TPZMatWithMemTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 18:43
 */

#ifndef TPZMATWITHMEMTRANSLATOR_H
#define TPZMATWITHMEMTRANSLATOR_H

#include "TPZMaterialTranslator.h"

template <class TMEMTranslator, class TFatherTranslator = TPZMaterialTranslator>
class TPZMatWithMemTranslator : public TFatherTranslator {
public:
    TPZMatWithMemTranslator();
    TPZMatWithMemTranslator(const TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZMatWithMemTranslator();
private:
    TMEMTranslator tMEMTranslator;
    TFatherTranslator parentTranslator;
};

template <class TMEMTranslator, class TFatherTranslator>
TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::TPZMatWithMemTranslator() {
}

template <class TMEMTranslator, class TFatherTranslator>
TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::TPZMatWithMemTranslator(const TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>& orig) {
}

template <class TMEMTranslator, class TFatherTranslator>
void TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    parentTranslator.UpdateStream(chunk, toVersion);
    int updatemem;
    chunk.mOldStream.Read(&updatemem);
    chunk.mNewStream.Write(&updatemem);
    tMEMTranslator.UpdateStream(chunk, toVersion);
    int i, size;
    chunk.mOldStream.Read(&size);
    chunk.mNewStream.Write(&size);
    for (i = 0; i < size; i++) {
        tMEMTranslator.UpdateStream(chunk, toVersion); //fMemory[i]
    }
}

template <class TMEMTranslator, class TFatherTranslator>
TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::~TPZMatWithMemTranslator() {
}



#endif /* TPZMATWITHMEMTRANSLATOR_H */

