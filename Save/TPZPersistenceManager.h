#ifndef PERSISTENCEMANAGER_H
#define PERSISTENCEMANAGER_H
#include <map>                   // for map
#include <ostream>               // for operator<<, string
#include "TPZRestoredInstance.h" // for TPZRestoredInstance
#include "pzmanvector.h"         // for TPZManVector
class TPZSaveable;
class TPZGeneralFStream;

namespace TPZPersistenceManagerNS {
    enum streamType { binary = 1, ascii = 2 };
};

using namespace TPZPersistenceManagerNS;

class TPZPersistenceManager {
    TPZPersistenceManager();
  protected:
    TPZGeneralFStream *mpStream;
    TPZManVector<TPZRestoredInstance, 10> mObjVec; // for READING from file
    
    // for WRITING to file
    TPZManVector<TPZSaveable *, 10> mPointersToSave;
    std::map<TPZSaveable *, int> mObjMap;          
    // WRITE-RELATED METHODS
  public:
    void OpenWrite(const std::string &fileName, streamType = binary);
    void WriteToFile(const TPZSaveable *);
  protected:
    void PopulateMap(const TPZSaveable *obj,
                     std::map<std::string, int> &fileVersionInfo);
    // READ-RELATED METHODS
  public:
    void OpenRead(const std::string &fileName, streamType = binary);
    void ReadFromFile(TPZSaveable &);
  protected:
    TPZSaveable *AssignPointers(const int &id);
    void AddInstanceToVec(TPZSaveable *, const int &id);
};

#endif // PERSISTENCEMANAGER_H
