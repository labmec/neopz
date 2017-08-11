#ifndef TPZBUFFEREDSTREAM_H
#define TPZBUFFEREDSTREAM_H
#include "TPZStream.h"
#include <sstream>

class TPZBufferedStream : public TPZStream {
public:
    TPZBufferedStream(TPZStream &fUnderlyingStream);
    TPZBufferedStream(const TPZBufferedStream &readBuffer);
    ~TPZBufferedStream();
    virtual void Read(double *p, const int &size=1);
    virtual void Write(const double *var, const int &size=1);
    TPZBufferedStream &operator=(const TPZBufferedStream &);
    virtual TPZBufferedStream &operator<<(TPZBufferedStream &other);
    virtual TPZBufferedStream &operator<<(const TPZBufferedStream &other);
    void Print();
    void BeginUpdate();
    void EndUpdate(const unsigned long &new_version);
protected:
    virtual void ReadFromBuffer(char *dest, const size_t &nBytes);
    virtual void ConstRead(char *dest, const size_t &nBytes) const;
    virtual void WriteToBuffer(const char *source, const size_t &nBytes);
private:
    char *fBuffer;
    char *fFirst;
    char *fLast;
    size_t fNAllocatedBytes;
    size_t fSize;
    bool fReadFromUnderlyingStream;
    TPZStream &fUnderlyingStream;

    virtual void ConstReadFromBuffer(char *dest, const size_t &nBytes) const;
    static const size_t MIN_SIZE_INCREMENT = size_t(1);
};
#endif // TPZBUFFEREDSTREAM_H
