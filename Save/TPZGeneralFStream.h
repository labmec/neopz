#ifndef TPZGENERALFSTREAM_H
#define TPZGENERALFSTREAM_H

#include "TPZStream.h"

class TPZGeneralFStream : public TPZStream{
  public:
    virtual void OpenRead(const std::string &fileName) = 0;
    virtual void OpenWrite(const std::string &fileName) = 0;
    
    virtual bool AmIOpenForRead() = 0;
    virtual bool AmIOpenForWrite() = 0;
    
    virtual void CloseRead() = 0;
    virtual void CloseWrite() = 0;
};
#endif//TPZGENERALFSTREAM_H
