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
    static TPZVec<TPZRestoredInstance> mObjVec; 
    static TPZVec<TPZAutoPointer<TPZChunkInTranslation>> mChunksVec; 
    static TPZVec<int64_t> mMainObjIds;
    static unsigned int mNextMainObjIndex;
    
    // for WRITING to file
    static std::map<std::string, uint64_t> mFileVersionInfo;
    static TPZContBufferedStream mObjectsStream;
    static TPZVec<const TPZSavable *> mPointersToSave;
public:
    static TPZContBufferedStream mCurrentObjectStream;
    static std::map<const TPZSavable *, int64_t> mObjMap;
    static int64_t mNextPointerToSave;
    // WRITE-RELATED METHODS
  public:
    static void OpenWrite(const std::string &fileName, streamType = binary);
    static void WriteToFile(const TPZSavable *);
    static void CloseWrite();
    static int64_t ScheduleToWrite(const TPZSavable *obj);
    static void WritePointer(const TPZSavable *obj, TPZStream *stream);
    
    // READ-RELATED METHODS
  public:
    static unsigned int OpenRead(const std::string &fileName, streamType = binary);
    static TPZRestoredInstance *NewRestoredInstance();
    static TPZSavable *ReadFromFile();
    static TPZSavable *GetInstance(const int64_t &objId);
    static TPZSavable *GetInstance(TPZStream *stream);
    static TPZAutoPointer<TPZSavable> GetAutoPointer(const int64_t &objId);
    static TPZAutoPointer<TPZSavable> GetAutoPointer(TPZStream *stream);
    static void CloseRead();
  protected:
    static void AddInstanceToVec(TPZSavable *, const int &id);
    
};

#endif // PERSISTENCEMANAGER_H
