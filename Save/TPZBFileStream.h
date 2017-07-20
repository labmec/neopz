#ifndef TPZBFILESTREAM_H
#define TPZBFILESTREAM_H

#include <fstream>
#include <string>
#include "TPZStream.h"
/**
 * @ingroup save
 * @brief Implements reading from and writing to a binary file. \ref save
 * "Persistency"
 */
class TPZBFileStream : public TPZStream {
 private:
  std::ifstream fIn;
  std::ofstream fOut;

  template <class T>
  void ReadData(T *p, int howMany);
  
  template <class T>
  void WriteData(const T *p, int howMany);
    
 public:
  TPZBFileStream();
  ~TPZBFileStream();

  void OpenRead(const std::string &fileName);
  void OpenWrite(const std::string &fileName);

  void CloseRead();
  void CloseWrite();
    
    using TPZStream::Read;
    using TPZStream::Write;
  virtual void Read(std::string *p, int howMany) {
    char buf[1000];
    int c;
    for (c = 0; c < howMany; c++) {
      int sz;
      Read(&sz, 1);
      Read(buf, sz);
      buf[sz] = 0;
      p[c] = buf;
    }
  }

  virtual void Write(const std::string *p, int howMany) {
    int c;
    for (c = 0; c < howMany; c++) {
      int sz = p[c].size();
      Write(&sz, 1);
      Write(p[c].c_str(), p[c].size());
    }
  }
};

#endif  // TPZBFILESTREAM_H
