// -*- c++ -*-
#ifndef PZSAVEH
#define PZSAVEH

#include <map>
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzfilebuffer.h"

const int TPZSAVEABLEID = -1;

using namespace std;

class TPZStream;
class TPZSaveable;

typedef TPZSaveable *(*TPZRestore_t)(TPZStream &,void *);


class TPZSaveable {

static map<int,TPZRestore_t> &Map() {
static map<int,TPZRestore_t> gMap;
   return gMap;
}

public:

virtual int ClassId() const ;

virtual void Write(TPZStream &buf, int withclassid = 1);

virtual void Read(TPZStream &buf, void *context = 0);

template<class T>
static void WriteObjects(TPZStream &buf, TPZVec<T> &vec) 
{
  int c,nc = vec.NElements();
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

  template<class T, int EXP>
  static void WriteObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec)
  {
    int c,m1=-1,nc = vec.NElements();
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

static void WriteObjects(TPZStream &buf, TPZVec<double> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
}

static void WriteObjects(TPZStream &buf, TPZVec<int> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
}


static void WriteObjects(TPZStream &buf, TPZVec<char> &vec) 
{
  int nel = vec.NElements();
  buf.Write(&nel,1);
  if(nel) buf.Write(&vec[0],vec.NElements());
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

template<class T, int N>
class TPZRestoreClass {

TPZRestoreClass()
{
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

template<class T, int N>
TPZRestoreClass<T,N> TPZRestoreClass<T,N>::gRestoreObject;


template<>
inline TPZSaveable *Restore<TPZSaveable>(TPZStream &buf, void *context) {
  return 0;
}

class TPZToto : public TPZSaveable {
  
public:

TPZToto() {
  int i;
  for(i=0; i<10; i++) fNum[i] = i;
  for(i=0; i<20; i++) fDNum[i] = i*30;
  for(i=0; i<20; i++) fStr[i] = 'a'+i;
}

virtual int ClassId() const 
{
  return 1;
}

 virtual void Write(TPZStream &buf, int withclassid = 1);

virtual void Read(TPZStream &buf, void *context);

static int main();

static TPZSaveable *Restore(TPZStream &buf, void *context);

private:
  
  int fNum[10];
  double fDNum[20];
  char fStr[20];

};
#endif //PZSAVEH
