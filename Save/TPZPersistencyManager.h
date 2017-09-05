#ifndef PERSISTENCYMANAGER_H
#define PERSISTENCYMANAGER_H
#include <map>              // for map
#include <ostream>          // for operator<<, string
#include "TPZRestoredInstance.h"  // for TPZRestoredInstance
#include "pzmanvector.h"    // for TPZManVector
class TPZSaveable;
class TPZGeneralFStream;

namespace TPZPersistencyManagerNS {
    enum streamType {
        binary = 1, ascii = 2
    };
};

using namespace TPZPersistencyManagerNS;

class TPZPersistencyManager {
  protected:
    TPZGeneralFStream *stream;
    TPZManVector<TPZRestoredInstance, 10> objVec; // for READING from file
    std::map<TPZSaveable *, int> objMap; //for WRITING to file
  public:
    TPZPersistencyManager();
    
    //WRITE-RELATED METHODS
    void OpenWrite(const std::string &fileName, streamType = binary);
    void PopulateMap(TPZSaveable *);
    void WriteToFile();
    
    //READ-RELATED METHODS
    void OpenRead(const std::string &fileName, streamType = binary);
    TPZSaveable *AssignPointers(const int &id);
    void AddObject(TPZSaveable *, const int &id);

  protected:
};

#endif//PERSISTENCYMANAGER_H
