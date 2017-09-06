#ifndef PERSISTENCYMANAGER_H
#define PERSISTENCYMANAGER_H
#include <map>                   // for map
#include <ostream>               // for operator<<, string
#include "TPZRestoredInstance.h" // for TPZRestoredInstance
#include "pzmanvector.h"         // for TPZManVector
class TPZSaveable;
class TPZGeneralFStream;

namespace TPZPersistencyManagerNS {
    enum streamType { binary = 1, ascii = 2 };
};

using namespace TPZPersistencyManagerNS;

class TPZPersistencyManager {
  protected:
    TPZGeneralFStream *mpStream;
    TPZManVector<TPZRestoredInstance, 10> mObjVec; // for READING from file
    std::map<TPZSaveable *, int> mObjMap;          // for WRITING to file
  public:
    TPZPersistencyManager();

    // WRITE-RELATED METHODS
    void OpenWrite(const std::string &fileName, streamType = binary);
    void WriteToFile(TPZSaveable *);

    // READ-RELATED METHODS
    void OpenRead(const std::string &fileName, streamType = binary);

  protected:
    // WRITE-RELATED METHODS
    void PopulateMap(TPZSaveable *obj,
                     std::map<std::string, int> &fileVersionInfo);
    // READ-RELATED METHODS
    TPZSaveable *AssignPointers(const int &id);
    void AddObjectToVec(TPZSaveable *, const int &id);
};

#endif // PERSISTENCYMANAGER_H
