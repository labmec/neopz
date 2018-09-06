#include "TPZPersistenceManager.h"
#include <stddef.h>            // for NULL
#include "TPZBFileStream.h"    // for TPZBFileStream
#include "TPZFileStream.h"     // for TPZFileStream
#include "TPZGeneralFStream.h" // for TPZGeneralFStream
#include "TPZSavable.h"       // for TPZSavable
#include "pzerror.h"           // for DebugStop, PZError
#include "pzvec.h"
#include "TPZContBufferedStream.h"             // for TPZVec
#include "TPZChunkInTranslation.h"
#include "TPZChunkTranslator.h"
#include <algorithm>

using namespace TPZPersistenceManagerNS;

TPZGeneralFStream *TPZPersistenceManager::mpStream;

// for READING from file
TPZVec<TPZRestoredInstance> TPZPersistenceManager::mObjVec;
TPZVec<TPZAutoPointer<TPZChunkInTranslation>> TPZPersistenceManager::mChunksVec;
TPZVec<int64_t> TPZPersistenceManager::mMainObjIds;
unsigned int TPZPersistenceManager::mNextMainObjIndex;
std::list<std::map < std::string, uint64_t>> TPZPersistenceManager::mVersionHistory;

// for WRITING to file
std::map<std::string, uint64_t> TPZPersistenceManager::mFileVersionInfo;
TPZContBufferedStream TPZPersistenceManager::mObjectsStream;
TPZVec<const TPZSavable *> TPZPersistenceManager::mPointersToSave;
TPZContBufferedStream TPZPersistenceManager::mCurrentObjectStream;
std::map<const TPZSavable *, int64_t> TPZPersistenceManager::mObjMap;
int64_t TPZPersistenceManager::mNextPointerToSave;

TPZPersistenceManager::TPZPersistenceManager() {
    mpStream = NULL;
    mObjVec.Resize(0);
    mObjMap.clear();
}

/********************************************************************
 *                                                                   *
 *                   WRITE-RELATED METHODS                           *
 *                                                                   *
 ********************************************************************/
void TPZPersistenceManager::OpenWrite(const std::string &fileName,
        const streamType cStreamType) {

    switch (cStreamType) {
        case binary:
        {
            mpStream = new TPZBFileStream();
            break;
        }
        case ascii:
        {
            mpStream = new TPZFileStream();
            break;
        }
        default:
            DebugStop();
    }

    mpStream->OpenWrite(fileName);

    std::string versionString("FileVersion");
    mpStream->Write(versionString.c_str(), 11);

#ifdef PZDEBUG
    if (!mpStream->AmIOpenForWrite()) {
        DebugStop();
    }
    if (mPointersToSave.size() > 0) {
        DebugStop();
    }
#endif    
    mFileVersionInfo.clear();
    mObjectsStream.clear();
    mObjMap.clear();
    mPointersToSave.clear();
    mMainObjIds.clear();
    mNextPointerToSave = 0;
}

void TPZPersistenceManager::WriteToFile(const TPZSavable *obj) {
    auto objId = ScheduleToWrite(obj);
    unsigned int nMainObjIds = mMainObjIds.size();
    mMainObjIds.resize(nMainObjIds + 1);
    mMainObjIds[nMainObjIds] = objId;

#ifdef PZDEBUG    
    uint64_t written = 0x46;
#endif
    for (; mNextPointerToSave < mPointersToSave.size(); ++mNextPointerToSave) {
        // writes obj-id
        mObjectsStream.Write(&mNextPointerToSave, 1);
        auto pointer = mPointersToSave[mNextPointerToSave];
        // writes classId
        auto classId = pointer->ClassId();
        if (classId == -1 || classId == 666) {
            DebugStop();
        }
        mObjectsStream.Write(&classId, 1);
        // writes object data
        mCurrentObjectStream.clear();
        pointer->Write(mCurrentObjectStream, 0);
        unsigned int size = mCurrentObjectStream.Size();
        mObjectsStream.Write(&size, 1);
        mObjectsStream << mCurrentObjectStream;
#ifdef PZDEBUG    
        std::cout << "ObjId: " << mNextPointerToSave
                << " ClassId: " << classId
                << " range: " << written + 1 << " (0x" << std::hex << written + 1 << std::dec << ") - "
                << written + sizeof (int64_t) + 2 * sizeof (int) +size
                << "(0x" << std::hex << written + sizeof (int64_t) + 2 * sizeof (int) +size << ")" << std::dec
                << " size: " << size << std::endl;
        written += sizeof (int64_t) + 2 * sizeof (int) +size;
#endif
    }
}

void TPZPersistenceManager::CloseWrite() {
    mpStream->Write(mFileVersionInfo);

    const uint64_t nObjects = mPointersToSave.size();
    mpStream->Write(&nObjects);

    size_t nObjectBytes = mObjectsStream.Size();
    char *temp = new char[nObjectBytes];
    mObjectsStream.GetDataFromBuffer(temp);
    mpStream->Write(temp, nObjectBytes);
    delete[] temp;

    unsigned int nMainObjects = mMainObjIds.size();
    mpStream->Write(&nMainObjects);
    for (auto mainObjId : mMainObjIds) {
        mpStream->Write(&mainObjId);
    }

    mpStream->CloseWrite();
    mFileVersionInfo.clear();
    mObjectsStream.clear();
    mObjMap.clear();
    mPointersToSave.clear();
    mMainObjIds.clear();
    mNextPointerToSave = 0;
}

int64_t TPZPersistenceManager::ScheduleToWrite(const TPZSavable *obj) {
    if (!obj) return -1;
    auto iMap = mObjMap.find(obj);
    int64_t objId;
    if (iMap == mObjMap.end()) { //object isn't in the map yet
        objId = mObjMap.size();
        if (mPointersToSave.size() != objId) {
            DebugStop();
        }
        mObjMap.insert(std::make_pair(obj, objId));
        mFileVersionInfo.insert(obj->Version());
        mPointersToSave.resize(objId + 1);
        mPointersToSave[objId] = obj;
    } else {
        objId = iMap->second;
    }
    return objId;
}

void TPZPersistenceManager::WritePointer(const TPZSavable *obj, TPZStream *stream) {
    const int64_t objId = ScheduleToWrite(obj);
    stream->Write(&objId);
}

/********************************************************************
 *                                                                   *
 *                    READ-RELATED METHODS                           *
 *                                                                   *
 ********************************************************************/

void TPZPersistenceManager::TranslateNextPointer(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    int64_t objId;
    chunk.mOldStream.Read(&objId);
    chunk.mNewStream.Write(&objId);
}

void TPZPersistenceManager::TranslatePointers(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    uint64_t size;
    chunk.mOldStream.Read(&size);
    chunk.mNewStream.Write(&size);
    for (uint64_t i = 0; i < size; ++i) {
        TranslateNextPointer(chunk, toVersion);
    }
}

unsigned int TPZPersistenceManager::OpenRead(const std::string &fileName,
        const streamType cStreamType) {

    if (TPZSavable::ClassIdMap().size() == 0) {
        //@TODO parallelize
        for (const auto &restoreClass : TPZSavable::RestoreClassSet()) {
#ifdef PZDEBUG
//            std::cout << restoreClass->Restore()->ClassId() << "\t" << typeid(*restoreClass->Restore()).name();
//            if (restoreClass->GetTranslator()){
//                std::cout << "\t" << typeid(*restoreClass->GetTranslator()).name();
//            }
//            std::cout << std::endl;
#endif
            TPZSavable *savable = restoreClass->Restore();
            //@TODO ensure thread-safety
            TPZSavable::RegisterClassId(savable->ClassId(), restoreClass);
            auto objHistory = savable->VersionHistory();
            for (auto versionMap : objHistory) {
                if (std::find(mVersionHistory.begin(), mVersionHistory.end(), versionMap) == mVersionHistory.end()) {
                    //@TODO ensure thread-safety
                    mVersionHistory.push_back(versionMap);
                }
            }
        }

        std::list<std::map < std::string, uint64_t>> mapsToRemove;
        for (auto versionMap : mVersionHistory) {
            for (auto otherMap : mVersionHistory) {
                if (otherMap == versionMap) {
                    continue;
                }
                if (std::search(otherMap.begin(), otherMap.end(), versionMap.begin(), versionMap.end()) != otherMap.end()) {
                    mapsToRemove.push_back(versionMap);
                    break;
                }
            }
        }
        for (auto mapToRemove : mapsToRemove) {
            mVersionHistory.remove(mapToRemove);
        }
    }

    switch (cStreamType) {
        case binary:
        {
            mpStream = new TPZBFileStream();
            break;
        }
        case ascii:
        {
            mpStream = new TPZFileStream();
            break;
        }
        default:
            DebugStop();
    }
    mpStream->OpenRead(fileName);
    std::string versionString("FileVersion");
    char versionRead[12];
    mpStream->Read(versionRead, 11);
    versionRead[11] = '\0'; // terminates c-style string

    if (versionString.compare(versionRead) != 0) { // unversioned file --- not readable
        PZError << "Error reading file - unversioned file" << std::endl;
        DebugStop();
    }

    mFileVersionInfo.clear();
    mpStream->Read(mFileVersionInfo);

    uint64_t nObjects;
    mpStream->Read(&nObjects);
    mChunksVec.Resize(nObjects);

    int64_t objId;
    int classId;
    unsigned int objSize;
    for (int64_t i = 0; i < nObjects; i++) {
        mpStream->Read(&objId);
        mpStream->Read(&classId);
        mpStream->Read(&objSize);
        mChunksVec[i] = new TPZChunkInTranslation(objId, classId, *mpStream, objSize, mFileVersionInfo);
    }

    std::map<std::string, uint64_t> currentVersionInfo = mFileVersionInfo;
    auto versionIt = std::find(mVersionHistory.begin(), mVersionHistory.end(), currentVersionInfo);
    if (versionIt == mVersionHistory.end()) {
        std::cerr << "This file cannot be read by this version of the program." << std::endl;
        DebugStop();
    }
    versionIt++;

    while (versionIt != mVersionHistory.end()) {
        std::map<std::string, uint64_t> nextVersion = *versionIt;
        //@TODO parallelize
        for (int64_t i = 0; i < nObjects; i++) {
            auto chunk = mChunksVec[i];
            chunk->mOldVersion = chunk->mNewVersion;
            chunk->mOldStream = chunk->mNewStream;
            chunk->mNewStream.clear();
        }
        //@TODO parallelize (it's already thread-safe)
        for (int64_t i = 0; i < nObjects; i++) {
            TPZAutoPointer<TPZChunkInTranslation> chunk = mChunksVec[i];
            int classId = chunk->GetClassId();
            if (classId != -1) { // this class still exists in this version
                chunk->mNewVersion = nextVersion;
                TPZChunkTranslator *translator = TPZSavable::ClassIdMap().at(classId)->GetTranslator();
                if (translator) {
                    translator->UpdateStream(*chunk, nextVersion);
                } else {
                    chunk->mNewStream = chunk->mOldStream;
                }
            }
        }
        currentVersionInfo = nextVersion;
        nObjects = mChunksVec.size();
        versionIt++;
    }
    // allocate
    mObjVec.Resize(mChunksVec.size());
    //@TODO parallelize (it's already thread-safe)
    for (unsigned int i = 0; i < mObjVec.size(); ++i) {
        int classId = mChunksVec[i]->GetClassId();
        if (classId != -1) {
            mObjVec[i].ResetReadStatus();
            mObjVec[i].SetInstance(TPZSavable::CreateInstance(classId));
        }
    }

    //@TODO parallelize (with a semaphore to read each chunk)
    for (unsigned int i = 0; i < mObjVec.size(); ++i) {
        int classId = mChunksVec[i]->GetClassId();
        if (classId != -1 && !mObjVec[i].IsAlreadyRead()) {
            mObjVec[i].SetRead();
            mObjVec[i].GetPointerToMyObj()->Read(mChunksVec[i]->mNewStream, NULL);
#ifdef PZDEBUG    
            if (mChunksVec[i]->mNewStream.Size() != 0) {
                DebugStop();
            }
#endif
        }
    }

    unsigned int nMainObjects;
    mpStream->Read(&nMainObjects);
    mMainObjIds.resize(nMainObjects);
    for (unsigned int i = 0; i < nMainObjects; ++i) {
        mpStream->Read(&mMainObjIds[i]);
    }
    mpStream->CloseRead();

    mNextMainObjIndex = 0;
    return nMainObjects;
}

int64_t TPZPersistenceManager::NewChunkInTranslation() {
    auto nChunks = mChunksVec.size();
    mChunksVec.resize(nChunks + 1);
    return nChunks;
}

void TPZPersistenceManager::SetChunk(const int64_t& objId, TPZAutoPointer<TPZChunkInTranslation> chunk) {
    mChunksVec[objId] = chunk;
}


TPZRestoredInstance *TPZPersistenceManager::NewRestoredInstance() {
    auto nObjects = mObjVec.size();
    mObjVec.Resize(nObjects + 1);
    return &mObjVec[nObjects];
}

TPZSavable *TPZPersistenceManager::ReadFromFile() {
    TPZSavable *obj = GetInstance(mMainObjIds[mNextMainObjIndex++]);
    if (mNextMainObjIndex == mMainObjIds.size()) {
        CloseRead();
    }
    return obj;
}

void TPZPersistenceManager::AddInstanceToVec(TPZSavable *obj, const int &cId) {
    mObjVec[cId].SetInstance(obj);
}

TPZSavable *TPZPersistenceManager::GetInstance(const int64_t &objId) {
    if (objId != -1) {
        if (!mObjVec[objId].IsAlreadyRead()) {
            mObjVec[objId].SetRead();
            mObjVec[objId].GetPointerToMyObj()->Read(mChunksVec[objId]->mNewStream, NULL);
            if (mChunksVec[objId]->mNewStream.Size() != 0) {
                DebugStop();
            }
        }
        return mObjVec[objId].GetPointerToMyObj();
    }
    return NULL;
}

TPZSavable *TPZPersistenceManager::GetInstance(TPZStream *stream) {
    int64_t objId;
    stream->Read(&objId);
    return GetInstance(objId);
}

TPZAutoPointer<TPZSavable> TPZPersistenceManager::GetAutoPointer(const int64_t &objId) {
    TPZAutoPointer<TPZSavable> autoPointer;
    if (objId != -1) {
        if (!mObjVec[objId].IsAlreadyRead()) {
            mObjVec[objId].SetRead();
            mObjVec[objId].GetPointerToMyObj()->Read(mChunksVec[objId]->mNewStream, NULL);
            if (mChunksVec[objId]->mNewStream.Size() != 0) {
                DebugStop();
            }
        }
        return mObjVec[objId].GetAutoPointerToMyObj();
    }
    return autoPointer;
}

std::shared_ptr<TPZSavable> TPZPersistenceManager::GetSharedPointer(const int64_t &objId) {
    if (objId != -1) {
        if (!mObjVec[objId].IsAlreadyRead()) {
            mObjVec[objId].SetRead();
            mObjVec[objId].GetPointerToMyObj()->Read(mChunksVec[objId]->mNewStream, NULL);
            if (mChunksVec[objId]->mNewStream.Size() != 0) {
                DebugStop();
            }
        }
        return mObjVec[objId].GetSharedPtrToMyObj();
    }
    std::shared_ptr<TPZSavable> sharedPtr;
    return sharedPtr;
}

TPZAutoPointer<TPZSavable> TPZPersistenceManager::GetAutoPointer(TPZStream *stream) {
    int64_t objId;
    stream->Read(&objId);
    return GetAutoPointer(objId);
}

std::shared_ptr<TPZSavable> TPZPersistenceManager::GetSharedPointer(TPZStream *stream) {
    int64_t objId;
    stream->Read(&objId);
    return GetSharedPointer(objId);
}

void TPZPersistenceManager::CloseRead() {
    mMainObjIds.clear();
    mChunksVec.clear();
    mObjVec.clear();
    mNextMainObjIndex = 0;
}
