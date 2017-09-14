#ifndef PERSISTENCEMANAGER_H
#define PERSISTENCEMANAGER_H
#include <map>                   // for map
#include <ostream>               // for operator<<, string
#include "TPZRestoredInstance.h" // for TPZRestoredInstance
#include "pzmanvector.h"
#include "TPZContBufferedStream.h"         // for TPZManVector
class TPZSaveable;
class TPZGeneralFStream;

namespace TPZPersistenceManagerNS {
    enum streamType { binary = 1, ascii = 2 };
};

using namespace TPZPersistenceManagerNS;

class TPZPersistenceManager {
    
    /** @brief This static function guarantees that the gFileVersionInfo object is available when needed */
    static std::map<std::string, long unsigned int> &gFileVersionInfo() {
        static std::map<std::string, long unsigned int> gFileVersionInfoObj;
        return gFileVersionInfoObj;
    }
    
    static void RegisterFileVersionInfo(std::string projectName, long unsigned int fileVersion) {
        gFileVersionInfo()[projectName]=fileVersion;
    }
    
    TPZPersistenceManager();
    std::pair<std::string, long unsigned int> getFileVersionInfo(std::string projectName);
    
  protected:
    TPZGeneralFStream *mpStream;
    TPZContBufferedStream mObjectsStream;
    TPZManVector<TPZRestoredInstance, 10> mObjVec; // for READING from file
    
    // for WRITING to file
    std::map<std::string, long unsigned int> mFileVersionInfo;
    TPZManVector<const TPZSaveable *, 10> mPointersToSave;
    std::map<const TPZSaveable *, long unsigned int> mObjMap;          
    // WRITE-RELATED METHODS
  public:
    void OpenWrite(const std::string &fileName, streamType = binary);
    void WriteToFile(const TPZSaveable *);
    void WritePointer(const TPZSaveable *obj);
  protected:
    long unsigned int ScheduleToWrite(const TPZSaveable *obj);

    // READ-RELATED METHODS
  public:
    void OpenRead(const std::string &fileName, streamType = binary);
    void ReadFromFile(TPZSaveable &);
  protected:
    TPZSaveable *AssignPointers(const int &id);
    void AddInstanceToVec(TPZSaveable *, const int &id);
};

#endif // PERSISTENCEMANAGER_H
