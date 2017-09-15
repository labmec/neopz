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
    mObjectsStream.clear();
    mObjMap.clear();
    mPointersToSave.clear();
    ScheduleToWrite(obj);

    TPZContBufferedStream thisPointerStream;
    for (long unsigned int i = 0; i < mPointersToSave.size(); ++i) {
        // writes obj-id
        thisPointerStream.Write(&i, 1);
        auto pointer = mPointersToSave[i];
        // writes classId
        auto classId = pointer->ClassId();
        thisPointerStream.Write(&classId, 1);
        // writes class-id and object
        pointer->Write(thisPointerStream, false);
        size_t size = thisPointerStream.Size();
        mObjectsStream.Write(&size, 1);
        mObjectsStream << thisPointerStream;
        thisPointerStream.clear();
    }

    mpStream->Write(mFileVersionInfo);

    const int nObjects = mPointersToSave.size();
    mpStream->Write(&nObjects);

    size_t nObjectBytes = mObjectsStream.Size();
    char temp[nObjectBytes];
    mObjectsStream.GetDataFromBuffer(temp);
    mpStream->Write(temp, nObjectBytes);

    mpStream->CloseWrite();
    mFileVersionInfo.clear();
    mObjectsStream.clear();
    mObjMap.clear();
    mPointersToSave.clear();
}

long unsigned int TPZPersistenceManager::ScheduleToWrite(const TPZSaveable *obj) {
    auto iMap = mObjMap.find(obj);
    long unsigned int nObjectId;
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
    long unsigned int nObjectId = ScheduleToWrite(obj);
    mObjectsStream.Write(&nObjectId, 1);
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

    std::map<std::string, int> fileVersionInfo;
    if (versionString.compare(versionRead) == 0) { // versioned file
        mpStream->Read(fileVersionInfo);
    } else { // unversioned file --- not readable
        PZError << "Error reading file - unversioned file" << std::endl;
        DebugStop();
    }
    // reads how many objects are to be read?

    // update mpStream now! at any point of the program's execution
    // the mpStream is considered to be at NeoPZ's most recent version.
    // the mpStream must NOT leave this ctor without being updated first.
    //TODO: Implement Translators
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
    mpStream->Read(&nObjects, 1);
    mObjVec.Resize(nObjects);

    size_t objSize;
    for (long unsigned int i = 0; i < nObjects; i++) {
        mpStream->Read(&objSize, 1);
        TPZRestoredInstance instance;
        instance.SetFileVersionInfo(mFileVersionInfo);
        instance.ReadFromStream(*mpStream, objSize);
        mObjVec[instance.GetObjId()] = instance; // We should not do this. This causes unnecessary copies.
        
        //TPZSaveable::Restore will read the class_id and
        //call the appropriate ::Read(TPZStream*) function,
        //where the attributes will be read from file
        //TPZSaveable * obj = TPZSaveable::Restore(*mpStream, NULL);
        //AddInstanceToVec(obj, objId);
    }

    dest = *(mObjVec[0].GetPointerToMyObj());
}

TPZSaveable *TPZPersistenceManager::AssignPointers(const int &cId) {
    TPZRestoredInstance &obj = mObjVec[cId];
    if (!obj.HaveMyPtrsBeenAssigned()) {
        obj.SetPtrsAsRestored();
        obj.GetPointerToMyObj()->AssignPointers(obj.MyPointersVec());
    }

    return obj.GetPointerToMyObj();
}

void TPZPersistenceManager::AddInstanceToVec(TPZSaveable *obj, const int &cId) {
    mObjVec[cId].SetInstance(obj);
}