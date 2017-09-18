#ifndef PERSISTENCEMANAGER_H
#define PERSISTENCEMANAGER_H
#include <map>                   // for map
#include <ostream>               // for operator<<, string
#include "TPZRestoredInstance.h"

template < class T, int NumExtAlloc>
class TPZManVector;
class TPZSaveable;
class TPZGeneralFStream;
template<class T>
class TPZAutoPointer;
class TPZChunkInTranslation;
class TPZContBufferedStream;

namespace TPZPersistenceManagerNS {
    enum streamType { binary = 1, ascii = 2 };
};

using namespace TPZPersistenceManagerNS;

class TPZPersistenceManager {
    
    TPZPersistenceManager();
    
  protected:
    static TPZGeneralFStream *mpStream;
    
    // for READING from file
    static TPZManVector<TPZRestoredInstance, 10> mObjVec; 
    static TPZManVector<TPZAutoPointer<TPZChunkInTranslation>, 10> mChunksVec; 
    
    // for WRITING to file
    static std::map<std::string, long unsigned int> mFileVersionInfo;
    static TPZManVector<const TPZSaveable *, 10> mPointersToSave;
    static TPZContBufferedStream mCurrentObjectStream;
    static std::map<const TPZSaveable *, long int> mObjMap;          
    // WRITE-RELATED METHODS
  public:
    static void OpenWrite(const std::string &fileName, streamType = binary);
    static void WriteToFile(const TPZSaveable *);
    static void WritePointer(const TPZSaveable *obj);
    static long int ScheduleToWrite(const TPZSaveable *obj);

    // READ-RELATED METHODS
  public:
    static void OpenRead(const std::string &fileName, streamType = binary);
    static TPZRestoredInstance *NewRestoredInstance();
    static void ReadFromFile(TPZSaveable &);
    static TPZSaveable *GetInstance(const long int &objId);
    static TPZAutoPointer<TPZSaveable> GetAutoPointer(const long int &objId);
  protected:
    static void AddInstanceToVec(TPZSaveable *, const int &id);
};

#endif // PERSISTENCEMANAGER_H
