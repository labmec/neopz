#include "pzworkpool.h"


TPZWorkPool::TPZWorkPool() : fStorage(NUMWORKPOOLS), fNextFree(NUMWORKPOOLS), fNumUsed(NUMWORKPOOLS) {

  fUnitSize = 0;
  fNumAllocated = 0;
  fLastDeleted = 0;
  for(int i=0; i<NUMWORKPOOLS; i++) {
    fStorage[i] = 0;
    fNextFree[i] = 0;
    fNumUsed[i] = 0;
  }
}

TPZWorkPool::~TPZWorkPool() {
  int nums = fStorage.NElements();
  for(int i = 0; i< nums; i++) {
    delete fStorage[i];
    fStorage[i] = 0;
    fNextFree[i] = 0;
    fNumUsed[i] = 0;
  }
}

void *TPZWorkPool::NewPointer(size_t Size) {
  if(Size != fUnitSize && fUnitSize != 0) {
    return new char[Size];
  }

  if(!fUnitSize) {
    fUnitSize = Size;
    fNumAllocated = WORKPOOLSIZE / Size;
    fCurrentSlot = FreeSlot();
  }

  if(!fNextFree[fCurrentSlot]) fCurrentSlot = FreeSlot();
  char *result = fNextFree[fCurrentSlot];
  if(!result) return NULL;

  char *NewNextFree = *((char **) result);
  fNumUsed[fCurrentSlot]++;
  fNextFree[fCurrentSlot] = NewNextFree;
  return result;
}

void TPZWorkPool::Release(void *Ptr) {
  char **WorkPointer = (char **) Ptr;
  char * Shoot = fStorage[fLastDeleted];
#ifdef __BOORLANDC__
  long offset = (((char huge *)Ptr) - ((char huge *)Shoot))/fUnitSize;
#else	
  long offset = (((char *)Ptr) - ((char *)Shoot))/fUnitSize;
#endif	
  if(offset >= 0 && offset < fNumAllocated) {
    char * &Help = fNextFree[fLastDeleted];
    *WorkPointer = Help;
    Help = (char *) Ptr;
    int &NumUsed = fNumUsed[fLastDeleted];
    NumUsed--;
    if(!NumUsed) {
      delete fStorage[fLastDeleted];
      fStorage[fLastDeleted] = NULL;
      fNextFree[fLastDeleted] = NULL;
    }
    return;
  }

  if(FindSlot(Ptr) >= 0) {
    char * &Help = fNextFree[fLastDeleted];
    *WorkPointer = Help;
    Help = (char *) Ptr;
    int &NumUsed = fNumUsed[fLastDeleted];
    NumUsed--;
    if(!NumUsed) {
      delete fStorage[fLastDeleted];
      fStorage[fLastDeleted] = NULL;
      fNextFree[fLastDeleted] = NULL;
    }
    return;
  }

  //	SysBeep(40);
}

int TPZWorkPool::FindSlot(void *Ptr) {

  long i, LoopSize = fStorage.NElements();
  long offset;
  for (i=0; i<LoopSize; i++){
#ifdef __BOORLANDC__		
    offset = (((char huge *)Ptr) - ((char huge *)fStorage[i]))/fUnitSize;
#else
    offset = (((char *)Ptr) - ((char *)fStorage[i]))/fUnitSize;
#endif
    if(offset >= 0 && offset < fNumAllocated) {
      fLastDeleted = i;
      return i;
    }
  }
  return -1;
}

int TPZWorkPool::FreeSlot(){

  long i,LoopSize;
  LoopSize = fStorage.NElements();
  for(i = 0; i<LoopSize; i++)
    if(fStorage[i] && fNextFree[i]) return i;
  for(i = 0; i< LoopSize; i++)
    if(!fStorage[i]) {
      char *runner = new char[fNumAllocated*fUnitSize];
      fStorage[i] = runner;
      //cout << "TPZWorkPool new storage i = " << i << " address = " << fStorage(i) << "\n";
      //cout.flush();
      void *LastAlloc = runner+fNumAllocated*fUnitSize;
      char *ic;
      for ( ic= runner; ic < (char *) LastAlloc; ic+=fUnitSize) *((void**) ic) = ic+fUnitSize;
      ic -= fUnitSize;
      *((void**)ic)=NULL;
      fNextFree[i] = runner;
      return i;
    }
  fStorage.Resize(LoopSize+10);
  fNextFree.Resize(LoopSize+10);
  fNumUsed.Resize(LoopSize+10);
  int ist;
  for(ist = LoopSize; ist < LoopSize+10; ist++) {
    fStorage[ist] = 0;
    fNextFree[ist] = 0;
    fNumUsed[ist] = 0;
  }
  {
    char *runner = new char[fNumAllocated*fUnitSize];
    fStorage[i] = runner;
    //cout << "TPZWorkPool new storage i = " << i << " address = " << fStorage(i) << "\n";
    //cout.flush();
    void *LastAlloc = runner+fNumAllocated*fUnitSize;
    char *ic;
    for ( ic= runner; ic < (char *) LastAlloc; ic+=fUnitSize) *((void**) ic) = ic+fUnitSize;
    ic -= fUnitSize;
    *((void**)ic) =NULL;
    fNextFree[i] = runner;
    return i;
  }
}

