#include "TPZPersistencyManager.h"
#include <stddef.h>             // for NULL
#include "TPZBFileStream.h"     // for TPZBFileStream
#include "TPZFileStream.h"      // for TPZFileStream
#include "TPZGeneralFStream.h"  // for TPZGeneralFStream
#include "TPZSaveable.h"        // for TPZSaveable
#include "pzerror.h"            // for DebugStop, PZError
#include "pzvec.h"              // for TPZVec

using namespace TPZPersistencyManagerNS;


TPZPersistencyManager::TPZPersistencyManager(){
    stream = NULL;
    objVec.Resize(0);
    objMap.clear();
}

void TPZPersistencyManager::OpenWrite(const std::string &fileName, const streamType sType){
    
    switch (sType) {
        case binary:{
            stream = new TPZBFileStream();
            break;
        }
        case ascii:{
            stream = new TPZFileStream();
            break;
        }
        default:
            DebugStop();
    }
    
    stream->OpenWrite(fileName);
    
    std::string versionString("FileVersion");
    stream->Write(&versionString);
    
}

//TODO: think in a name to combine PopulateMap+WriteToFile
//Having them as separate methods would allow for more than one object (and its pointers and etc.) to be saved. However, now I cannot think in a way to restore more than one object
void TPZPersistencyManager::PopulateMap(TPZSaveable *obj){
    
 
}

void TPZPersistencyManager::WriteToFile(){
    //write all versionInfo
    
    //transverse the map writing everyone to file
    for (std::map<TPZSaveable *,int>::iterator it=objMap.begin(); it!=objMap.end(); ++it){
        stream->Write(&(it->second));
        //(it->first)->Write(stream);
    }
    stream->CloseWrite();
    objMap.clear();//map can now be deleted
}


void TPZPersistencyManager::OpenRead(const std::string &fileName, const streamType sType){
    
    switch (sType) {
        case binary:{
            stream = new TPZBFileStream();
            TPZBFileStream *temp = (TPZBFileStream *)stream;
            temp->OpenWrite(fileName);
            break;
        }
        case ascii:{
            stream = new TPZFileStream();
            TPZFileStream *temp = (TPZFileStream *)stream;
            temp->OpenRead(fileName);
            break;
        }
        default:
            DebugStop();
    }
    
    std::string versionString("FileVersion");
    char versionRead[12];
    stream->Read(versionRead, 11);
    versionRead[11] = '\0';    // terminates c-style string
    
    std::map<std::string, int> fileVersionInfo;
    if (versionString.compare(versionRead) == 0) { // versioned file
        stream->Read(fileVersionInfo);
    } else { // unversioned file --- not readable
        PZError << "Error reading file - unversioned file" << std::endl;
        DebugStop();
    }
    
    //update stream now! at any point of the program's execution
    //the stream is considered to be at NeoPZ's most recent version.
    //the stream must NOT leave this ctor without being updated first.
}


TPZSaveable *
TPZPersistencyManager::AssignPointers(const int &id) {
    TPZRestoredInstance &obj = objVec[id];
    if (!obj.HaveIBeenRestored()) {
        obj.SetAsRestored();
        obj.GetPointerToMyObj()->AssignPointers(obj.MyPointersVec());
    }
    
    return obj.GetPointerToMyObj();
}

void TPZPersistencyManager::AddObject(TPZSaveable *obj, const int &id){
    const int oldSize = objVec.size();
#ifdef PZDEBUG
    if (oldSize!= id) {
        DebugStop();//objects should be added at the same order as they were saved in the first place.
    }
#endif
    objVec.Resize(oldSize + 1);
    objVec[oldSize].SetPointerToMyObj(obj);    
}
