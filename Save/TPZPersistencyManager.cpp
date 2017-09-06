#include "TPZPersistencyManager.h"
#include <stddef.h>            // for NULL
#include "TPZBFileStream.h"    // for TPZBFileStream
#include "TPZFileStream.h"     // for TPZFileStream
#include "TPZGeneralFStream.h" // for TPZGeneralFStream
#include "TPZSaveable.h"       // for TPZSaveable
#include "pzerror.h"           // for DebugStop, PZError
#include "pzvec.h"             // for TPZVec

using namespace TPZPersistencyManagerNS;

TPZPersistencyManager::TPZPersistencyManager() {
    mpStream = NULL;
    mObjVec.Resize(0);
    mObjMap.clear();
}

void TPZPersistencyManager::OpenWrite(const std::string &fileName,
                                      const streamType cStreamType) {

    switch (cStreamType) {
    case binary: {
        mpStream = new TPZBFileStream();
        break;
    }
    case ascii: {
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

void TPZPersistencyManager::PopulateMap(
    TPZSaveable *obj, std::map<std::string, int> &fileVersionInfo) {
}

void TPZPersistencyManager::WriteToFile(TPZSaveable *obj) {
    // write all versionInfo
    std::map<std::string, int> fileVersionInfo;
    PopulateMap(obj, fileVersionInfo);

    mpStream->Write(fileVersionInfo);

    // writes how many objects are to be written?
    // const int mapSize = mObjMap.size();
    // mpStream->Write(mapSize);

    // transverse the map writing everyone to file
    for (std::map<TPZSaveable *, int>::iterator iMap = mObjMap.begin();
         iMap != mObjMap.end(); ++iMap) {
        // writes obj-cId
        mpStream->Write(&(iMap->second));
        // writes class-cId and object
        (iMap->first)->Write(*mpStream, true);
    }
    mpStream->CloseWrite();
    mObjMap.clear(); // map can now be deleted
}

void TPZPersistencyManager::OpenRead(const std::string &fileName,
                                     const streamType cStreamType) {

    switch (cStreamType) {
    case binary: {
        mpStream = new TPZBFileStream();
        break;
    }
    case ascii: {
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
}

TPZSaveable *TPZPersistencyManager::AssignPointers(const int &cId) {
    TPZRestoredInstance &obj = mObjVec[cId];
    if (!obj.HaveIBeenRestored()) {
        obj.SetAsRestored();
        obj.GetPointerToMyObj()->AssignPointers(obj.MyPointersVec());
    }

    return obj.GetPointerToMyObj();
}

void TPZPersistencyManager::AddObjectToVec(TPZSaveable *obj, const int &cId) {
    mObjVec[cId].SetPointerToMyObj(obj);
}
