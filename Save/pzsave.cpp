#include "pzsave.h"
#include "pzfilebuffer.h"

#include <iostream>
#include <fstream>
using namespace std;

map<int,TPZRestore_t> TPZSaveable::gMap;

int TPZSaveable::ClassId() const {
  return -1;
}

void TPZSaveable::Write(TPZStream &buf, int withclassid) 
{
  if(withclassid) {
    int var = ClassId();
    buf.Write(&var,1);
  }
}

void TPZSaveable::Read(TPZStream &buf, void *context)
{
}

void TPZSaveable::Register(int classid, TPZRestore_t fun) 
{
  map<int,TPZRestore_t>::iterator it;
  it = gMap.find(classid);
  if(it != gMap.end()) 
  {
    cout << "TPZSaveable::Register duplicate classid " << classid << endl;
    return;
  }
  gMap[classid] = fun;
}

TPZSaveable *TPZSaveable::Restore(TPZStream &buf, void *context) {
  int classid;
  buf.Read(&classid,1);
  map<int,TPZRestore_t>::iterator it;
  it = gMap.find(classid);
  if(it == gMap.end()) 
  {
    cout << "TPZSaveable trying to restore unknown object " << classid << endl;
    return 0;
  }
  TPZRestore_t fun= it->second;
  return (*fun)(buf,context);
}
  
  



void TPZToto::Write(TPZStream &buf, int withclassid)
{
  if(withclassid) TPZSaveable::Write(buf);
  buf.Write(fNum,10);
  buf.Write(fDNum,20);
  buf.Write(fStr,20);
}

void TPZToto::Read(TPZStream &buf, void *context) 
{
  buf.Read(fNum,10);
  buf.Read(fDNum,20);
  buf.Read(fStr,20);
  
}

TPZSaveable *TPZToto::Restore(TPZStream &buf, void *context) {
  TPZToto *t = new TPZToto;
  t->Read(buf,context);
  return t;
}

int TPZToto::main() 
{
  TPZToto test;
  TPZVec<double> tes(10,0.);
  TPZAdmChunkVector<TPZToto> vectes;
  TPZRestore_t fun = ::Restore<TPZToto>;
  TPZSaveable::Register(test.ClassId(),fun);
  {
    TPZFileStream f;
    f.OpenWrite("file.out");
    test.Write(f);
    WriteObjects(f,tes);
    WriteObjects<TPZToto>(f,vectes);
  }
  {
    TPZFileStream f;
    f.OpenRead("file.out");
    //    int classid;
    //    f.Read(&classid,1);
    TPZSaveable *t = TPZSaveable::Restore(f,0);
    TPZToto *toto = (TPZToto *) t;
    cout << toto->fNum[0] << endl;
//    test.Read(f,0);
  }
  return 0;  
  
}  
