#ifndef WORKPOOLH
#define WORKPOOLH

#include <stddef.h>
#include "pzvec.h"

#define WORKPOOLSIZE 32767	// allocation in chunks of 64 kbytes
#define NUMWORKPOOLS 10		// number of chunks in the workpool

class TPZWorkPool {

 public:

  TPZWorkPool();

  ~TPZWorkPool();

  void *NewPointer(size_t Size);
	
  void Release(void *Ptr);

 private:

  size_t fUnitSize;
  long fNumAllocated;
  TPZVec<char *> fStorage;
  TPZVec<char *> fNextFree;
  TPZVec<int> fNumUsed;
  int fLastDeleted;
  int fCurrentSlot;
	
  int FreeSlot();
	
  int FindSlot(void *Ptr);

};

#endif

