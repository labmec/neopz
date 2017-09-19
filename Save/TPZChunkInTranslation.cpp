/* 
 * File:   TPZChunkInTranslator.cpp
 * Author: quinelato
 * 
 * Created on September 18, 2017, 3:54 PM
 */

#include "TPZChunkInTranslation.h"

TPZChunkInTranslation::TPZChunkInTranslation(const long int &objId, const int &classId, TPZStream &stream, const size_t &chunkSize, std::map<std::string, long unsigned int> &versionInfo) :
mObjId(objId),
mClassId(classId),
mNewVersion(versionInfo) {
    this->ReadFromStream(stream, chunkSize);
}

TPZChunkInTranslation::TPZChunkInTranslation(const TPZChunkInTranslation& orig) {
}

TPZChunkInTranslation::~TPZChunkInTranslation() {
}

void TPZChunkInTranslation::ReadFromStream(TPZStream &stream, size_t nBytes) {
    char temp[nBytes];
    stream.Read(temp, nBytes);
    mNewStream.Write(temp, nBytes);
}

long int TPZChunkInTranslation::GetObjId() const {
    return mObjId;
}

int TPZChunkInTranslation::GetClassId() const {
    return mClassId;
}

