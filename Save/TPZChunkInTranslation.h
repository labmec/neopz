/* 
 * File:   TPZChunkInTranslator.h
 * Author: quinelato
 *
 * Created on September 18, 2017, 3:54 PM
 */

#ifndef TPZCHUNKINTRANSLATION_H
#define TPZCHUNKINTRANSLATION_H

#include "TPZContBufferedStream.h"


class TPZChunkInTranslation {
public:
    TPZChunkInTranslation(const long int &objId, const int &classId, TPZStream &stream, const size_t &chunkSize, std::map<std::string, long unsigned int> &versionInfo);
    TPZChunkInTranslation(const TPZChunkInTranslation& orig);
    virtual ~TPZChunkInTranslation();
    long int GetObjId() const;
    int GetClassId() const;
private:
    void ReadFromStream(TPZStream &stream, const size_t nBytes);
private:
    TPZContBufferedStream mOldStream;
    TPZContBufferedStream mNewStream;
    
    std::map<std::string, long unsigned int> mOldVersion;
    std::map<std::string, long unsigned int> mNewVersion;
    
    long int mObjId;
    int mClassId;
    
    TPZManVector<long int, 2> mNewObjIds;
    
    friend TPZPersistenceManager;
};

#endif /* TPZCHUNKINTRANSLATION_H */

