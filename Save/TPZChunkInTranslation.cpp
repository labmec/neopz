/* 
 * File:   TPZChunkInTranslator.cpp
 * Author: quinelato
 * 
 * Created on September 18, 2017, 3:54 PM
 */

#include "TPZChunkInTranslation.h"

TPZChunkInTranslation::TPZChunkInTranslation(const int64_t &objId, const int &classId, TPZStream &stream, const size_t &chunkSize, const std::map<std::string, uint64_t> &versionInfo) :
mObjId(objId),
mClassId(classId),
mNewVersion(versionInfo) {
    this->ReadFromStream(stream, chunkSize);
}

TPZChunkInTranslation::TPZChunkInTranslation(const int64_t &objId, const int &classId, const TPZContBufferedStream &oldStream, const std::map<std::string, uint64_t> &versionInfo) :
mOldStream(oldStream),
mObjId(objId),
mClassId(classId),
mOldVersion(versionInfo) {
}

TPZChunkInTranslation::TPZChunkInTranslation(const TPZChunkInTranslation& orig) : mOldStream(orig.mOldStream), mNewStream(orig.mNewStream), mOldVersion(orig.mOldVersion), mNewVersion(orig.mNewVersion), mObjId(orig.mObjId), mClassId(orig.mClassId), mNewObjIds(orig.mNewObjIds) {
}

TPZChunkInTranslation::~TPZChunkInTranslation() {
}

void TPZChunkInTranslation::ReadFromStream(TPZStream &stream, const size_t nBytes) {
    char *temp = new char[nBytes];
    stream.Read(temp, nBytes);
    mNewStream.Write(temp, nBytes);
	delete[] temp;
}

int64_t TPZChunkInTranslation::GetObjId() const {
    return mObjId;
}

int TPZChunkInTranslation::GetClassId() const{
    return mClassId;
}

