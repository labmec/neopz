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

extern "C" {
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
}


using namespace std;



/**
@author Thiago M. N. Oliveira
*/
class TPZBFileStream : public TPZStream
{

  int ofd;
  int ifd;

 public:

  TPZBFileStream(){}

  virtual ~TPZBFileStream() {
    close(ofd);
    close(ifd);
  }

  void OpenWrite(const string &filename) {
    ofd = open(filename.c_str(), O_WRONLY | O_CREAT, 0644);
  }

  void OpenRead(const string &filename) {
    ifd = open(filename.c_str(), O_RDONLY);
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
    int n_bytes = size*sizeof(T);
    char *buff = (char*) calloc (n_bytes+3, sizeof(char));
    strncat(buff, (char*) p, n_bytes);
    write(ofd, buff, n_bytes);
    fsync(ofd);
    free(buff);
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
    void Reads(T *p, int size)
  {
    int n_bytes = size*sizeof(T);
    char *buff = (char*) calloc (n_bytes+3, sizeof(char));
    read(ifd, buff, n_bytes);
    memcpy(p, buff, n_bytes);
    free(buff);
  }

};



#endif
