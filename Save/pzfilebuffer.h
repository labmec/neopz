#ifndef PZFILEBUFFERH
#define PZFILEBUFFERH
#include <iostream>
#include <fstream>
#include <string>
#include "pzmanvector.h"
using namespace std;



class TPZStream {

 public:

  virtual ~TPZStream() {}

  virtual void Write(int *p, int size=1)=0;

  virtual void Write(double *p, int size=1)=0;

  virtual void Write(const char *p, int size=1)=0;
  
  virtual void Write(string *p, int size=1) = 0;

  void Write(TPZFlopCounter *p, int size=1) 
  {
    int i;
    for(i=0; i<size; i++) Write(&(p[i].fVal),1);
  }

  virtual void Read(int *p, int size=1)=0;

  virtual void Read(double *p, int size=1)=0;

  void Read(TPZFlopCounter *p, int size=1)
  {
    int i;
    for(i=0; i<size; i++)
    {
      Read(&(p[i].fVal),1);
    }
  }

  virtual void Read(char *p, int size=1)=0;

  virtual void Read(string *p, int size=1) = 0;

};

class TPZFileStream : public TPZStream {

  ofstream fo;
  ifstream fi;

 public:

  TPZFileStream(){}

  virtual ~TPZFileStream();

  void OpenWrite(const string &filename) {
    fo.open(filename.c_str());
  }

  void OpenRead(const string &filename) {
    fi.open(filename.c_str());
  }

  virtual void Write(int *p, int size) {
    Writes<int>(p,size);
  }
  
  virtual void Write(double *p, int size) {
    Writes<double>(p,size);
  }

  virtual void Write(const char *p, int size) {
    Writes<char>(p,size);
  }

  virtual void Write(string *p, int size) {
    Writes<string>(p,size);
  }
  
  template<class T>
    void  Writes(const T *p, int size) 
  {
    int c;
    for(c=0; c<size; c++) fo << p[c] << endl;
  }
  
  
  virtual void Read(int *p, int size) {
    Reads<int>(p,size);
  }

  virtual void Read(double *p, int size) {
    Reads<double>(p,size);
  }

  virtual void Read(char *p, int size) {
    Reads<char>(p,size);
  }

  virtual void Read(string *p, int size) {
    int c;
    char buf[2560];
    for(c=0; c<size; c++) 
    {
      fi.getline(buf,2560);
      p[c] = buf;
    }
  }

  template<class T>
    void Reads(T *p, int size) {
    int c;
    char buf[100];
    for(c=0; c<size; c++) fi >> p[c];
    fi.getline(buf,100);
  }

};

#endif
