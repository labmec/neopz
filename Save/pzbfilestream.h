//
// C++ Interface: TPZBFileStream
//
// Description: 
//
//
// Author: Thiago M. N. Oliveira <thiago@labmec.fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef STDPZBFILESTREAM_H
#define STDPZBFILESTREAM_H

#include <pzfilebuffer.h>

#include <stdio.h>
/*
extern "C" {
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
}
*/

using namespace std;


/// this class implements the interface to a binary file
/**
@author Thiago M. N. Oliveira
*/
class TPZBFileStream : public TPZStream
{

  FILE *ofd;
  FILE *ifd;

 public:

  TPZBFileStream(){
    ofd=0;
    ifd=0;
  }

  virtual ~TPZBFileStream() {
    if(ofd) fclose(ofd);
    if(ifd) fclose(ifd);
  }

  void OpenWrite(const string &filename) {
    ofd = fopen(filename.c_str(),"wb" );
  }

  void OpenRead(const string &filename) {
    ifd = fopen(filename.c_str(), "rb");
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
    int c;
    for(c=0; c<size; c++) 
    {
      int sz = p[c].size();
      Write(&sz,1);
      Write(p[c].c_str(),p[c].size());
    }
  }
  
  template<class T>
    void  Writes(const T *p, int size) 
  {
    fwrite(p,sizeof(T),size,ofd);
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
  
  virtual void Read(string *p, int size) 
  {
    char buf[1000];
    int c;
    for(c=0; c<size; c++) 
    {
      int sz;
      Read(&sz,1);
      Read(buf,sz);
      buf[sz] = 0;
      p[c] = buf;
    }
  }

  template<class T>
    void Reads(T *p, int size)
  {
    fread(p,sizeof(T),size,ifd);
  }

};



#endif
