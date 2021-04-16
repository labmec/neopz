/* 
 * File:   TPZTensorTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 15:55
 */

#ifndef TPZTENSORTRANSLATOR_H
#define TPZTENSORTRANSLATOR_H

#include "TPZChunkTranslator.h"
#include "pzmanvector.h"
#include "TPZChunkInTranslation.h"

template <typename T>
class TPZTensorTranslator : public TPZChunkTranslator {
public:
    TPZTensorTranslator();
    TPZTensorTranslator(const TPZTensorTranslator<T>& orig);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZTensorTranslator();
private:

};

template <typename T>
TPZTensorTranslator<T>::TPZTensorTranslator() {
}

template <typename T>
TPZTensorTranslator<T>::TPZTensorTranslator(const TPZTensorTranslator<T>& orig) {
}

template <typename T>
void TPZTensorTranslator<T>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    TPZManVector<T, 6> fData;
    chunk.mOldStream.Read(fData);
    chunk.mNewStream.Write(fData);
}

template <typename T>
TPZTensorTranslator<T>::~TPZTensorTranslator() {
}



#endif /* TPZTENSORTRANSLATOR_H */

