// -*- c++ -*-
#ifndef PZSAVEH
#define PZSAVEH

#include <map>
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzfilebuffer.h"

using namespace std;

class TPZStream;
class TPZSaveable;

typedef TPZSaveable *(*TPZRestore_t)(TPZStream &,void *);


class TPZSaveable {

static map<int,TPZRestore_t> gMap;

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
  int c,nc = vec.NElements();
  buf.Write(&nc,1);
  for(c=0; c<nc; c++) vec[c]->Write(buf);
}

  template<class T, int EXP>
  static void WriteObjectPointers(TPZStream &buf, TPZChunkVector<T *,EXP> &vec)
  {
    int c,nc = vec.NElements();
    for(c=0; c<nc; c++) vec[c]->Write(buf);
  }

  template<class T, int EXP>
  static void WriteObjectPointers(TPZStream &buf, TPZAdmChunkVector<T *,EXP> &vec)
  {
    int c,nc = vec.NElements();
    buf.Write(&nc,1);
    for(c=0; c<nc; c++) vec[c]->Write(buf);
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
  buf.Write(&vec[0],vec.NElements());
}

static void WriteObjects(TPZStream &buf, TPZVec<int> &vec) 
{
  buf.Write(&vec[0],vec.NElements());
}


static void WriteObjects(TPZStream &buf, TPZVec<char> &vec) 
{
  buf.Write(&vec[0],vec.NElements());
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

class TPZToto : public TPZSaveable {
  
public:

TPZToto() {
  int i;
  for(i=0; i<10; i++) fNum[i] = i;
  for(i=0; i<20; i++) fDNum[i] = i*30;
  for(i=0; i<20; i++) fStr[i] = 'a'+i;
}

virtual int ClassId() 
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
