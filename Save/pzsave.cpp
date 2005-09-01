#include "pzsave.h"
#include "pzfilebuffer.h"

#include <iostream>
#include <fstream>

#ifdef WIN32
#pragma warning (disable:4786)
#endif


using namespace std;


int TPZSaveable::ClassId() const {
  return -1;
}

void TPZSaveable::Write(TPZStream &buf, int withclassid) 
{
  if(withclassid) { 
    int var = ClassId();
    if(var == -1)
    {
      cout << "TPZSaveable::Write with classid -1 expect trouble\n";
    }
    buf.Write(&var,1);
  }
}

void TPZSaveable::Read(TPZStream &buf, void *context)
{
}

void TPZSaveable::Register(int classid, TPZRestore_t fun) 
{
#ifndef ELLIPS
  map<int,TPZRestore_t>::iterator it;
 it = Map().find(classid);
  if(it != Map().end()) 
  {
    cout << "TPZSaveable::Register duplicate classid " << classid << endl;
    return;
  }
  Map()[classid] = fun;
#endif
}

TPZSaveable *TPZSaveable::Restore(TPZStream &buf, void *context) {
#ifndef ELLIPS
  int classid;
  buf.Read(&classid,1);
  map<int,TPZRestore_t>::iterator it;
  it = Map().find(classid);
  if(it == Map().end()) 
  {
    cout << "TPZSaveable trying to restore unknown object " << classid << endl;
    return 0;
  }
  TPZRestore_t fun= it->second;
  return (*fun)(buf,context);
#else
  return 0;
#endif
}
