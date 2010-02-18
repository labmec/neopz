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

#include "pzfilebuffer.h"

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

  void OpenWrite(const std::string &filename) {
    ofd = fopen(filename.c_str(),"wb" );
  }

  void OpenRead(const std::string &filename) {
    ifd = fopen(filename.c_str(), "rb");
	  if(!ifd)
	  {
		  std::cout << "could not open file " << filename << std::endl;
	  }
  }

  virtual void Write(int *p, int size) {
#ifndef WIN32
    Writes<int>(p,size);
#endif
  }
  
  virtual void Write(double *p, int size) {
#ifndef WIN32
    Writes<double>(p,size);
#endif
  }

  virtual void Write(const char *p, int size) {
#ifndef WIN32
    Writes<char>(p,size);
#endif
  }

  virtual void Write(std::string *p, int size) {
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
#ifndef WIN32
    Reads<int>(p,size);
#endif
  }

  virtual void Read(double *p, int size) {
#ifndef WIN32
    Reads<double>(p,size);
#endif
  }

  virtual void Read(char *p, int size) {
#ifndef WIN32
    Reads<char>(p,size);
#endif
  }
  
  virtual void Read(std::string *p, int size) 
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
	  if(ifd)
	  {
		  fread(p,sizeof(T),size,ifd);
	  }
  }

};



#endif
