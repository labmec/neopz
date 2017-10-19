#ifndef PERSISTENCEMANAGER_H
#define PERSISTENCEMANAGER_H
#include <map>                   // for map
#include <ostream>               // for operator<<, string
#include "TPZRestoredInstance.h"

template < class T, int NumExtAlloc>
class TPZManVector;
class TPZSavable;
class TPZGeneralFStream;
template<class T>
class TPZAutoPointer;
class TPZChunkInTranslation;
class TPZContBufferedStream;

namespace TPZPersistenceManagerNS {
    enum streamType { binary = 1, ascii = 2 };
};

using namespace TPZPersistenceManagerNS;

//TODO: Implement CloseRead()
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
    static TPZManVector<const TPZSavable *, 10> mPointersToSave;
public:
    static TPZContBufferedStream mCurrentObjectStream;
    static std::map<const TPZSavable *, long int> mObjMap;
    static long int mNextPointerToSave;
    // WRITE-RELATED METHODS
  public:
    static void OpenWrite(const std::string &fileName, streamType = binary);
    static void WriteToFile(const TPZSavable *);
    static void CloseWrite();
    static long int ScheduleToWrite(const TPZSavable *obj);
    static void WritePointer(const TPZSavable *obj, TPZStream *stream);
    
    // READ-RELATED METHODS
  public:
    static unsigned int OpenRead(const std::string &fileName, streamType = binary);
    static TPZRestoredInstance *NewRestoredInstance();
    static TPZSavable *ReadFromFile();
    static TPZSavable *GetInstance(const long int &objId);
    static TPZSavable *GetInstance(TPZStream *stream);
    static TPZAutoPointer<TPZSavable> GetAutoPointer(const long int &objId);
    static TPZAutoPointer<TPZSavable> GetAutoPointer(TPZStream *stream);
  protected:
    static void AddInstanceToVec(TPZSavable *, const int &id);
    
};

#endif // PERSISTENCEMANAGER_H
