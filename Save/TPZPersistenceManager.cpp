#include "TPZPersistenceManager.h"
#include <stddef.h>            // for NULL
#include "TPZBFileStream.h"    // for TPZBFileStream
#include "TPZFileStream.h"     // for TPZFileStream
#include "TPZGeneralFStream.h" // for TPZGeneralFStream
#include "TPZSaveable.h"       // for TPZSaveable
#include "pzerror.h"           // for DebugStop, PZError
#include "pzvec.h"
#include "TPZContBufferedStream.h"             // for TPZVec

using namespace TPZPersistenceManagerNS;

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
    mpStream->Write(&versionString);
}

void TPZPersistenceManager::WriteToFile(const TPZSaveable *obj) {
#ifdef PZDEBUG
    if (!mpStream->AmIOpenForWrite()) {
        DebugStop();
    }
    if (mPointersToSave.size() > 0) {
        DebugStop();
    }
#endif    
    mFileVersionInfo.clear();
    TPZContBufferedStream objectsStream;
    mObjMap.clear();
    mPointersToSave.clear();
    ScheduleToWrite(obj);

    for (long int i = 0; i < mPointersToSave.size(); ++i) {
        // writes obj-id
        objectsStream.Write(&i, 1);
        auto pointer = mPointersToSave[i];
        // writes classId
        auto classId = pointer->ClassId();
        objectsStream.Write(&classId, 1);
        // writes class-id and object
        mCurrentObjectStream.clear();
        pointer->Write(mCurrentObjectStream, false);
        size_t size = mCurrentObjectStream.Size();
        objectsStream.Write(&size, 1);
        objectsStream << mCurrentObjectStream;
    }

    mpStream->Write(mFileVersionInfo);

    const int nObjects = mPointersToSave.size();
    mpStream->Write(&nObjects);

    size_t nObjectBytes = objectsStream.Size();
    char temp[nObjectBytes];
    objectsStream.GetDataFromBuffer(temp);
    mpStream->Write(temp, nObjectBytes);

    mpStream->CloseWrite();
    mFileVersionInfo.clear();
    objectsStream.clear();
    mObjMap.clear();
    mPointersToSave.clear();
}

long int TPZPersistenceManager::ScheduleToWrite(const TPZSaveable *obj) {
    if (!obj) return -1;
    auto iMap = mObjMap.find(obj);
    long int nObjectId;
    if (iMap == mObjMap.end()) { //object isn't in the map yet
        nObjectId = mObjMap.size();
        if (mPointersToSave.size() != nObjectId) {
            DebugStop();
        }
        mObjMap.insert(std::make_pair(obj, nObjectId));
        mFileVersionInfo.insert(obj->Version());
        mPointersToSave.resize(nObjectId + 1);
        mPointersToSave[nObjectId] = obj;
    } else {
        nObjectId = iMap->second;
    }
    return nObjectId;
}

void TPZPersistenceManager::WritePointer(const TPZSaveable *obj) {
    long int nObjectId = ScheduleToWrite(obj);
    mCurrentObjectStream.Write(&nObjectId, 1);
}

/********************************************************************
 *                                                                   *
 *                    READ-RELATED METHODS                           *
 *                                                                   *
 ********************************************************************/

void TPZPersistenceManager::OpenRead(const std::string &fileName,
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
    mpStream->OpenRead(fileName);
    std::string versionString("FileVersion");
    char versionRead[12];
    mpStream->Read(versionRead, 11);
    versionRead[11] = '\0'; // terminates c-style string

    if (versionString.compare(versionRead) != 0) { // unversioned file --- not readable
        PZError << "Error reading file - unversioned file" << std::endl;
        DebugStop();
    }
    // reads how many objects are to be read?

    // update mpStream now! at any point of the program's execution
    // the mpStream is considered to be at NeoPZ's most recent version.
    // the mpStream must NOT leave this ctor without being updated first.
    //TODO: Implement Translators
}

TPZRestoredInstance *TPZPersistenceManager::NewRestoredInstance() {
    auto nObjects = mObjVec.size();
    mObjVec.Resize(nObjects + 1);
    return &mObjVec[nObjects];
}

void TPZPersistenceManager::ReadFromFile(TPZSaveable & dest) {
#ifdef PZDEBUG
    if (!mpStream->AmIOpenForRead()) {
        DebugStop();
    }
#endif
    mFileVersionInfo.clear();
    mpStream->Read(mFileVersionInfo);

    long unsigned int nObjects;
    mpStream->Read(&nObjects);
    mChunksVec.Resize(nObjects);

    long int objId;
    int classId;
    size_t objSize;
    for (long unsigned int i = 0; i < nObjects; i++) {
        mpStream->Read(&objId);
        mpStream->Read(&classId);
        mpStream->Read(&objSize);
        mChunksVec[objId] = new TPZChunkInTranslation(objId, classId, *mpStream, objSize, mFileVersionInfo);
    }

    //    std::map<std::string, long unsigned int> currentVersionInfo = mFileVersionInfo;
    //    std::map<std::string, long unsigned int> nextVersion = ComputeNextVersion(currentVersionInfo);
    //    while (){ //nextVersion is not null
    //        // move newChunk to oldChunk
    //        for (long unsigned int i = 0; i < nObjects; i++) {
    //            //translate
    //        }
    //        currentVersionInfo = nextVersion;
    //        nextVersion = ComputeNextVersion(currentVersionInfo);
    //    }
    // allocate
    mObjVec.Resize(mChunksVec.size());
    for (unsigned int i = 0; i < mObjVec.size(); ++i) {
        int classId = mChunksVec[i]->GetClassId();
        if (classId != -1) {
            mObjVec[i].SetInstance(TPZSaveable::CreateInstance(classId));
        }
        //mObjVec[i].SetInstance();
    }


    // read
    /*TPZSaveable::Restore will read the class_id and
    call the appropriate ::Read(TPZStream*) function,
    where the attributes will be read from file*/
    //TPZSaveable * obj = TPZSaveable::Restore(*mpStream, NULL);

    dest = *(mObjVec[0].GetPointerToMyObj());
}

void TPZPersistenceManager::AddInstanceToVec(TPZSaveable *obj, const int &cId) {
    mObjVec[cId].SetInstance(obj);
}

static TPZSaveable *TPZPersistenceManager::GetInstance(const long int &objId) {
    return (objId == -1) ? NULL : mObjVec[objId].GetPointerToMyObj();
}

static TPZAutoPointer<TPZSaveable>& TPZPersistenceManager::GetAutoPointer(const long int &objId){
    TPZAutoPointer<TPZSaveable> autoPointer;
    if (objId != -1) {
        return mObjVec[objId].GetAutoPointerToMyObj();
    }
    return autoPointer;
}