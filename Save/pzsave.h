// -*- c++ -*-
#ifndef PZSAVEH
#define PZSAVEH

//#include "glib.h"
#include <map>
#include <vector>

#ifdef WIN32
#pragma warning (disable:4786)
#endif

#include "pzvec.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzfilebuffer.h"
#include "pzreal.h"
#include "tpzautopointer.h"

const int TPZSAVEABLEID = -1;


class TPZStream;
class TPZSaveable;

typedef TPZSaveable *(*TPZRestore_t)(TPZStream &,void *);


/// This class defines the interface to save and restore objects from TPZStream objects
class TPZSaveable {

#ifndef ELLIPS
static std::map<int,TPZRestore_t> &Map() {
static std::map<int,TPZRestore_t> gMap;
   return gMap;
}
#endif

public:

virtual ~TPZSaveable()
{
}

virtual int ClassId() const ;

virtual void Write(TPZStream &buf, int withclassid);

virtual void Read(TPZStream &buf, void *context);

template<class T>
static void WriteObjects(TPZStream &buf, TPZVec<T> &vec) 
{
  int c,nc = vec.NElements();
  buf.Write(&nc,1);
  for(c=0; c<nc; c++) vec[c].Write(buf,0);
}

template<class T>
static void WriteObjects(TPZStream &buf, std::vector<T> &vec) 
{
  int c,nc = vec.size();
  buf.Write(&nc,1);
  for(c=0; c<nc; c++) vec[c].Write(buf,0);
}

  template<class T, int EXP>
  static void WriteObjects(TPZStream &buf, TPZChunkVector<T,EXP> &vec)
  {
    int c,nc = vec.NElements();
    buf.Write(&nc,1);
    for(c=0; c<nc; c++) vec[c].Write(buf,0);
  }
	
  template<class T, int EXP>
  static void WriteObjects(TPZStream &buf, TPZAdmChunkVector<T,EXP> &vec)
  {
    int c,nc = vec.NElements();
    buf.Write(&nc,1);
    for(c=0; c<nc; c++) vec[c].Write(buf,0);
    buf.Write(&vec.fCompactScheme,1);
    WriteObjects(buf,vec.fFree);
    WriteObjects(buf,vec.fNFree);
  }
	
      


template<class T>
static void WriteObjectPointers(TPZStream &buf, TPZVec<T *> &vec) 
{
  int c,nc = vec.NElements(),one = -1;
  buf.Write(&nc,1);
  for(c=0; c<nc; c++) 
  {
    if(vec[c]) 
    {
      vec[c]->Write(buf);
    } else {
      buf.Write(&one,1);
    }
  }
}
template<class T>
    static void WriteObjectPointers(TPZStream &buf, std::map<int, TPZAutoPointer<T> > &vec) 
{
  int nc = vec.size(),one = -1;
  buf.Write(&nc,1);
  typedef typename std::map<int, TPZAutoPointer<T> >::iterator vecit_type;
  vecit_type vecit;
  for(vecit=vec.begin(); vecit!= vec.end(); vecit++) 
  {
    int id = vecit->first;
    buf.Write(&id,1);
    if(vecit->second) 
    {
      vecit->second->Write(buf,1);
    } else {
      buf.Write(&one,1);
    }
  }
}

  template<class T, int EXP>
  static void WriteObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec)
  {
    int c,m1=-1,nc = vec.NElements();
    buf.Write(&nc,1);
    for(c=0; c<nc; c++) 
    {
      T *ptr = vec[c];
      if(ptr) ptr->Write(buf);
      else buf.Write(&m1,1);
    }
  }

  template<class T, int EXP>
  static void WriteObjectPointers(TPZStream &buf, TPZAdmChunkVector<T *,EXP> &vec)
  {
    int c,m1=-1,nc = vec.NElements();
    buf.Write(&nc,1);
    for(c=0; c<nc; c++) 
    {
      T *ptr = vec[c];
      if(ptr) ptr->Write(buf,1);
      else buf.Write(&m1,1);
    }
    buf.Write(&vec.fCompactScheme,1);
    WriteObjects(buf,vec.fFree);
    WriteObjects(buf,vec.fNFree);
  }
	
template<class T>
static void ReadObjects(TPZStream &buf, TPZVec<T> &vec, void *context) 
{
  int c,nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  for(c=0; c<nc; c++) 
  {
    vec[c].Read(buf,context);
  }
}

static void ReadObjects(TPZStream &buf, TPZVec<int> &vec) 
{
  int nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  if(nc) buf.Read(&vec[0],nc);
}

template<class T>
static void ReadObjects(TPZStream &buf, std::vector<T> &vec, void *context) 
{
  int c,nc;
  buf.Read(&nc,1);
  vec.resize(nc);
  for(c=0; c<nc; c++) 
  {
    vec[c].Read(buf,context);
  }
}

static void ReadObjects(TPZStream &buf, std::vector<int> &vec) 
{
  int nc;
  buf.Read(&nc,1);
  vec.resize(nc);
  if(nc) buf.Read(&vec[0],nc);
}

static void ReadObjects(TPZStream &buf, std::vector<REAL> &vec) 
{
  int nc;
  buf.Read(&nc,1);
  vec.resize(nc);
  if(nc) buf.Read(&vec[0],nc);
}

static void ReadObjects(TPZStream &buf, TPZVec<REAL> &vec) 
{
  int nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  if(nc) buf.Read(&vec[0],nc);
}

template<int N>
static void ReadObjects(TPZStream &buf, TPZManVector<REAL,N> &vec) 
{
  int nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  if(nc) buf.Read(&vec[0],nc);
}


  template<class T, int EXP>
  static void ReadObjects(TPZStream &buf, TPZChunkVector<T,EXP> &vec, void *context) 
{
  int c,nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  for(c=0; c<nc; c++) 
  {
    vec[c].Read(buf,context);
  }
}

  template<class T, int EXP>
  static void ReadObjects(TPZStream &buf, TPZAdmChunkVector<T,EXP> &vec, void *context)
  {
    int c,nc;
    buf.Read(&nc,1);
    vec.Resize(nc);
    for(c=0; c<nc; c++) vec[c].Read(buf,context);
    buf.Read(&vec.fCompactScheme,1);
    ReadObjects(buf,vec.fFree);
    ReadObjects(buf,vec.fNFree);
  }

template<class T>
static void ReadObjectPointers(TPZStream &buf, TPZVec<T *> &vec, void *context)
{
  int c,nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  for(c=0; c<nc; c++) 
  {
    vec[c] = dynamic_cast<T *>(Restore(buf,context));
  }  
}
 
template<class T>
    static void ReadObjectPointers(TPZStream &buf, std::map<int, TPZAutoPointer<T> > &vec, void *context)
{
  int c,nc;
  buf.Read(&nc,1);
//  vec.Resize(nc);
  for(c=0; c<nc; c++) 
  {
    int id;
    buf.Read(&id,1);
    vec[id] = TPZAutoPointer<T>(dynamic_cast<T *>(Restore(buf,context)));
  }  
}
 
  template<class T, int EXP>
  static void ReadObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec, void *context) 
{
  int c,nc;
  buf.Read(&nc,1);
  vec.Resize(nc);
  for(c=0; c<nc; c++) 
  {
    vec[c] = dynamic_cast<T *>(Restore(buf,context));
  }
}

  template<class T, int EXP>
  static void ReadObjectPointers(TPZStream &buf, TPZAdmChunkVector<T *,EXP> &vec, void *context)
  {
    int c,nc;
    buf.Read(&nc,1);
    vec.Resize(nc);
    for(c=0; c<nc; c++) vec[c] = dynamic_cast<T *>(Restore(buf,context));
    buf.Read(&vec.fCompactScheme,1);
    ReadObjects(buf,vec.fFree);
    ReadObjects(buf,vec.fNFree);
  }

static void WriteObjects(TPZStream &buf, TPZVec<REAL> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
}

static void WriteObjects(TPZStream &buf, std::vector<REAL> &vec) 
{
  int nel = vec.size();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.size());
}

#ifndef ELLIPS
static void WriteObjects(TPZStream &buf, TPZVec<std::TPZFlopCounter> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
}

static void WriteObjects(TPZStream &buf, std::vector<std::TPZFlopCounter> &vec) 
{
  int nel = vec.size();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.size());
}
#endif

static void WriteObjects(TPZStream &buf, TPZVec<int> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
}

static void WriteObjects(TPZStream &buf, std::vector<int> &vec) 
{
  int nel = vec.size();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.size());
}

static void WriteObjects(TPZStream &buf, TPZVec<char> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
}

static void WriteObjects(TPZStream &buf, std::vector<char> &vec) 
{
  int nel = vec.size();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.size());
}


static void Register(int classid, TPZRestore_t fun);

static TPZSaveable *Restore(TPZStream &buf, void *context);

};




template<class T>
TPZSaveable *Restore(TPZStream &buf, void *context) {
  T *ptr = new T;
  ptr->Read(buf,context);
  return ptr;
}

#ifndef ELLIPS
/// this class implements an interface to register a class id and a restore function
/**
A declaration of the type "template class<classname, classid> put in .cpp file does the trick
The static object which is "automatically" created calls the proper interface of the TPZSaveable class
*/
#ifndef BORLAND
template<class T, int N>
class TPZRestoreClass {
public:
TPZRestoreClass()
{
  std::string func_name = __PRETTY_FUNCTION__;
#ifndef WIN32
  std::cout << func_name << std::endl;
#endif
  TPZSaveable::Register(N,Restore);
}

public:
static TPZSaveable *Restore(TPZStream &buf, void *context) {
  T *ptr = new T;
  ptr->Read(buf,context);
  return ptr;
}
private:
   static TPZRestoreClass gRestoreObject;
};

template<>
inline TPZSaveable *TPZRestoreClass<TPZSaveable,-1>::Restore(TPZStream &buf, void *context)
{
  return 0;
}
template<class T, int N>
TPZRestoreClass<T,N> TPZRestoreClass<T,N>::gRestoreObject;

template<>
inline TPZSaveable *Restore<TPZSaveable>(TPZStream &buf, void *context) {
  return 0;
}

template class TPZRestoreClass<TPZSaveable, -1>;
#endif //borland
#endif //ellips

#endif //PZSAVEH


