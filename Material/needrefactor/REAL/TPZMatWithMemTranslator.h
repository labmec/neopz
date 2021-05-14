/* 
 * File:   TPZMatWithMemTranslator.h
 * Author: thiago
 *
 * Created on 12 de Março de 2018, 18:43
 */

#ifndef TPZMATWITHMEMTRANSLATOR_H
#define TPZMATWITHMEMTRANSLATOR_H

#include "TPZMaterialTranslator.h"
#include "TPZChunkInTranslation.h"
#include "Hash/TPZHash.h"

template <class TMEMTranslator, class TFatherTranslator = TPZMaterialTranslator>
class TPZMatWithMemTranslator : public TFatherTranslator {
public:
    TPZMatWithMemTranslator();
    TPZMatWithMemTranslator(const TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>& orig);

    virtual void UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual void UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

    virtual ~TPZMatWithMemTranslator();
private:
    void UpdateAttributesV3(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);
    void UpdateFromV3(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion);

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
void TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    auto old_version = chunk.mOldVersion["NeoPZ"];
    auto new_version = toVersion.at("NeoPZ");
    switch (old_version) {
        case 1:
        case 2:
            if (new_version <= 3) {
                UpdateAttributesV3(chunk, toVersion);
                break;
            } else {
                DebugStop();
            }
        case 3:
            if (new_version == 3) {
                UpdateAttributesV3(chunk, toVersion);
            } else {
                UpdateFromV3(chunk, toVersion);
            }
            break;
        default:
            UpdateAttributes(chunk, toVersion);
            break;
    }
}

template <class TMEMTranslator, class TFatherTranslator>
void TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::UpdateAttributesV3(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
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
void TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::UpdateFromV3(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    parentTranslator.UpdateStream(chunk, toVersion);
    int updatemem;
    chunk.mOldStream.Read(&updatemem);
    chunk.mNewStream.Write(&updatemem);
    tMEMTranslator.UpdateStream(chunk, toVersion);
    int size;
    chunk.mOldStream.Read(&size);
    auto objId = TPZPersistenceManager::NewChunkInTranslation();
    chunk.mNewStream.Write(&objId);
    int admChunkVectorClassId = Hash("TPZAdmChunkVector") ^ (Hash("TPZChunkVector") ^ tMEMTranslator.GetClassId() << 1 ^ (10 << 2)) << 1;
    TPZAutoPointer<TPZChunkInTranslation> admChunkVectorChunk = new TPZChunkInTranslation(objId, admChunkVectorClassId, chunk.mOldStream, chunk.mOldVersion);
    TPZPersistenceManager::SetChunk(objId, admChunkVectorChunk);
    uint64_t size_64 = size;
    admChunkVectorChunk->mNewStream.Write(&size_64);
    for (int i = 0; i < size; i++) {
        tMEMTranslator.UpdateStream(*admChunkVectorChunk, toVersion); //fMemory[i]
    }
    chunk.mOldStream = admChunkVectorChunk->mOldStream;
    admChunkVectorChunk->mOldStream.clear();
    int compactScheme = 0;
    admChunkVectorChunk->mNewStream.Write(&compactScheme);
    {
        int64_t nel = 0;
        admChunkVectorChunk->mNewStream.Write(&nel);
    }//TPZStack<int> fFree;
    {
        int64_t nel = 0;
        admChunkVectorChunk->mNewStream.Write(&nel);
    } //TPZManVector<int> fNFree;
    chunk.mNewObjIds.push_back(objId);
}

template <class TMEMTranslator, class TFatherTranslator>
void TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    parentTranslator.UpdateStream(chunk, toVersion);
    int updatemem;
    chunk.mOldStream.Read(&updatemem);
    chunk.mNewStream.Write(&updatemem);
    tMEMTranslator.UpdateStream(chunk, toVersion);
    TPZPersistenceManager::TranslateNextPointer(chunk, toVersion); //fMemory[i]
}

template <class TMEMTranslator, class TFatherTranslator>
TPZMatWithMemTranslator<TMEMTranslator, TFatherTranslator>::~TPZMatWithMemTranslator() {
}



#endif /* TPZMATWITHMEMTRANSLATOR_H */

