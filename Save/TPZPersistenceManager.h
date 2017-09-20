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
    static TPZManVector<long int,2> mMainObjIds;
    static unsigned int mNextMainObjIndex;
    
    // for WRITING to file
    static std::map<std::string, long unsigned int> mFileVersionInfo;
    static TPZContBufferedStream mObjectsStream;
    static TPZManVector<const TPZSaveable *, 10> mPointersToSave;
    static TPZContBufferedStream mCurrentObjectStream;
    static std::map<const TPZSaveable *, long int> mObjMap;          
    // WRITE-RELATED METHODS
  public:
    static void OpenWrite(const std::string &fileName, streamType = binary);
    static void WriteToFile(const TPZAutoPointer<TPZSaveable>obj);
    static void WriteToFile(const TPZSaveable *);
    static void CloseWrite();
    static void WritePointer(const TPZSaveable *obj);
    static long int ScheduleToWrite(const TPZSaveable *obj, TPZStream *stream);
    
    // READ-RELATED METHODS
  public:
    static unsigned int OpenRead(const std::string &fileName, streamType = binary);
    static TPZRestoredInstance *NewRestoredInstance();
    static TPZSaveable *ReadFromFile();
    static TPZSaveable *GetInstance(const long int &objId, char a);
    static TPZSaveable *GetInstance(TPZStream *stream);
    static TPZAutoPointer<TPZSaveable> GetAutoPointer(const long int &objId, char a);
    static TPZAutoPointer<TPZSaveable> GetAutoPointer(TPZStream *stream);
  protected:
    static void AddInstanceToVec(TPZSaveable *, const int &id);
    
};

#endif // PERSISTENCEMANAGER_H
