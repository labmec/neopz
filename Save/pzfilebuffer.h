#ifndef PZFILEBUFFERH
#define PZFILEBUFFERH
#include <iostream>
#include <fstream>
#include <string>
using namespace std;



class TPZStream {

 public:

  virtual ~TPZStream() {}

  virtual void Write(int *p, int size)=0;

  virtual void Write(double *p, int size)=0;

  virtual void Write(char *p, int size)=0;

  virtual void Read(int *p, int size)=0;

  virtual void Read(double *p, int size)=0;

  virtual void Read(char *p, int size)=0;


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

  virtual void Write(char *p, int size) {
    Writes<char>(p,size);
  }

  template<class T>
    void  Writes(T *p, int size) 
  {
    int c;
    for(c=0; c<size; c++) fo << p[c] << ' ';
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

  template<class T>
    void Reads(T *p, int size) {
    int c;
    for(c=0; c<size; c++) fi >> p[c];
  }

};

#endif
