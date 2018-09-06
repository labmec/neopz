#include "TPZContBufferedStream.h"
#include <iostream>
#include <stdexcept>
#include <string.h>

#ifdef _AUTODIFF
#include "fad.h"
#endif

TPZContBufferedStream::TPZContBufferedStream() {

    fNAllocatedBytes = MIN_SIZE_INCREMENT;
    fBuffer = new char[fNAllocatedBytes];
    fFirst = fBuffer;
    fSize = 0;
    fLast = fBuffer - 1;
}

TPZContBufferedStream::TPZContBufferedStream(const TPZContBufferedStream &other)
: TPZStream(other), fBuffer(NULL) {
    *this = other;
}

TPZContBufferedStream &TPZContBufferedStream::
        operator=(const TPZContBufferedStream &other) {
    fNAllocatedBytes = other.fNAllocatedBytes;
    if (fBuffer)
        delete[] fBuffer;
    fBuffer = new char[fNAllocatedBytes];
    memcpy(fBuffer, other.fBuffer, other.fSize);
    fSize = other.fSize;
    fFirst = fBuffer + (other.fFirst-other.fBuffer);
    fLast = fFirst - 1 + fSize;
    return *this;
}

TPZContBufferedStream::~TPZContBufferedStream() {
    delete[] fBuffer;
}

TPZContBufferedStream &TPZContBufferedStream::
operator<<(TPZContBufferedStream &other) {
    if (&other != this) {
        const unsigned int nBytesOther = other.fSize;
        char *temp = new char[nBytesOther];
        other.ReadFromBuffer(temp, nBytesOther);
        WriteToBuffer(temp, nBytesOther);
        delete[] temp;
    }
    return *this;
}

TPZContBufferedStream &TPZContBufferedStream::
operator<<(const TPZContBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char *temp = new char[nBytesOther];
    other.ConstRead(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
    delete[] temp;
    return *this;
}

void TPZContBufferedStream::ReadFromBuffer(char *dest, const size_t &nBytes) {
    if (nBytes > fSize) {
        std::string msg("TPZContBufferedStream: Cannot read ");
        msg.append(std::to_string(nBytes));
        msg.append(" bytes; there are only ");
        msg.append(std::to_string(fSize));
        msg.append(" available.");
        PZError << msg << std::endl;
        DebugStop();
    }
    memcpy(dest, fFirst, nBytes);
    fFirst += nBytes;
    fSize -= nBytes;
}

void TPZContBufferedStream::ConstRead(char *dest, const size_t &nBytes) const {
    ConstReadFromBuffer(dest, nBytes);
}

void TPZContBufferedStream::ConstReadFromBuffer(char *dest,
        const size_t &nBytes) const {
    if (nBytes > fSize) {
        std::string msg("TPZContBufferedStream: Cannot read ");
        msg.append(std::to_string(nBytes));
        msg.append(" bytes; there are only ");
        msg.append(std::to_string(fSize));
        msg.append(" available.");
        PZError << msg << std::endl;
        DebugStop();
    }
    char *endBuffer = fBuffer + fNAllocatedBytes;
    memcpy(dest, fFirst, nBytes);
}

void TPZContBufferedStream::WriteToBuffer(const char *source,
        const size_t &nBytes) {
    if (fFirst - fBuffer + fSize + nBytes > fNAllocatedBytes) {
        const size_t oldSize = fSize;
        const size_t newAllocatedBytes =
                oldSize * 1.1 + nBytes +
                MIN_SIZE_INCREMENT; // 10% increase + nBytes + MIN_SIZE_INCREMENT
        char *temp = new char[newAllocatedBytes];
        ConstReadFromBuffer(temp, oldSize);
        memcpy(temp + oldSize, source, nBytes);
        delete[] fBuffer;
        fBuffer = temp;
        fNAllocatedBytes = newAllocatedBytes;
        fFirst = fBuffer;
        fSize = oldSize + nBytes;
        fLast = fBuffer - 1 + fSize;
    } else {
        memcpy(fLast + 1, source, nBytes);
        fLast += nBytes;
        fSize += nBytes;
    }
}

void TPZContBufferedStream::Print() {
    std::cout << "fSize=" << fSize << std::endl;
    double *temp = new double[fSize / 8];
    ConstRead(reinterpret_cast<char *> (temp), fSize);
    for (unsigned int i = 0; i < fSize / 8; ++i) {
        std::cout << temp[i] << " ";
    }
    delete[] temp;
    std::cout << std::endl;
}

void TPZContBufferedStream::Write(const int *p, int howMany) {
    WriteData<int>(p, howMany);
}

void TPZContBufferedStream::Write(const unsigned int *p, int howMany) {
    WriteData<unsigned int>(p, howMany);
}

void TPZContBufferedStream::Write(const uint64_t *p, int howMany) {
    WriteData<uint64_t>(p, howMany);
}

void TPZContBufferedStream::Write(const int64_t *p, int howMany) {
    WriteData<int64_t>(p, howMany);
}

void TPZContBufferedStream::Write(const float *p, int howMany) {
    WriteData<float>(p, howMany);
}

void TPZContBufferedStream::Write(const double *p, int howMany) {
    WriteData<double>(p, howMany);
}

void TPZContBufferedStream::Write(const unsigned char *p, int howMany) {
    WriteData<unsigned char>(p, howMany);
}

void TPZContBufferedStream::Write(const char *p, int howMany) {
    WriteData<char>(p, howMany);
}

void TPZContBufferedStream::Write(const std::complex<float> *p, int howMany) {
    WriteData<std::complex<float>>(p, howMany);
}

void TPZContBufferedStream::Write(const std::complex<double> *p, int howMany) {
    WriteData<std::complex<double>>(p, howMany);
}

#ifdef _AUTODIFF

void TPZContBufferedStream::Write(const TFad<1,REAL> *p, int howMany) {
    WriteData<TFad<1,REAL>>(p, howMany);
}

void TPZContBufferedStream::Write(const TFad<6,REAL> *p, int howMany) {
    WriteData<TFad<6,REAL>>(p, howMany);
}

void TPZContBufferedStream::Write(const TFad<8,REAL> *p, int howMany) {
    WriteData<TFad<8,REAL>>(p, howMany);
}

void TPZContBufferedStream::Write(const TFad<9,REAL> *p, int howMany) {
    WriteData<TFad<9,REAL>>(p, howMany);
}

void TPZContBufferedStream::Write(const TFad<10,REAL> *p, int howMany) {
    WriteData<TFad<10,REAL>>(p, howMany);
}

void TPZContBufferedStream::Write(const TFad<14,REAL> *p, int howMany) {
    WriteData<TFad<14,REAL>>(p, howMany);
}

void TPZContBufferedStream::Write(const Fad<float> *p, int howMany) {
    WriteData<Fad<float>>(p, howMany);
}

void TPZContBufferedStream::Write(const Fad<double> *p, int howMany) {
    WriteData<Fad<double>>(p, howMany);
}

#endif

void TPZContBufferedStream::Read(int *p, int howMany) {
    ReadData<int>(p, howMany);
}

void TPZContBufferedStream::Read(unsigned int *p, int howMany) {
    ReadData<unsigned int>(p, howMany);
}

void TPZContBufferedStream::Read(uint64_t *p, int howMany) {
    ReadData<uint64_t>(p, howMany);
}

void TPZContBufferedStream::Read(int64_t *p, int howMany) {
    ReadData<int64_t>(p, howMany);
}

void TPZContBufferedStream::Read(float *p, int howMany) {
    ReadData<float>(p, howMany);
}

void TPZContBufferedStream::Read(double *p, int howMany) {
    ReadData<double>(p, howMany);
}

void TPZContBufferedStream::Read(unsigned char *p, int howMany) {
    ReadData<unsigned char>(p, howMany);
}

void TPZContBufferedStream::Read(char *p, int howMany) {
    ReadData<char>(p, howMany);
}

void TPZContBufferedStream::Read(std::complex<float> *p, int howMany) {
    ReadData<std::complex<float>>(p, howMany);
}

void TPZContBufferedStream::Read(std::complex<double> *p, int howMany) {
    ReadData<std::complex<double>>(p, howMany);
}

#ifdef _AUTODIFF

void TPZContBufferedStream::Read(TFad<1,REAL> *p, int howMany) {
    ReadData<TFad<1,REAL>>(p, howMany);
}

void TPZContBufferedStream::Read(TFad<6,REAL> *p, int howMany) {
    ReadData<TFad<6,REAL>>(p, howMany);
}

void TPZContBufferedStream::Read(TFad<8,REAL> *p, int howMany) {
    ReadData<TFad<8,REAL>>(p, howMany);
}

void TPZContBufferedStream::Read(TFad<9,REAL> *p, int howMany) {
    ReadData<TFad<9,REAL>>(p, howMany);
}

void TPZContBufferedStream::Read(TFad<10,REAL> *p, int howMany) {
    ReadData<TFad<10,REAL>>(p, howMany);
}

void TPZContBufferedStream::Read(TFad<14,REAL> *p, int howMany) {
    ReadData<TFad<14,REAL>>(p, howMany);
}

void TPZContBufferedStream::Read(Fad<float> *p, int howMany) {
    ReadData<Fad<float>>(p, howMany);
}

void TPZContBufferedStream::Read(Fad<double> *p, int howMany) {
    ReadData<Fad<double>>(p, howMany);
}

#endif

void TPZContBufferedStream::GetDataFromBuffer(char *dest) const {
    memcpy(dest, fFirst, fSize);
}

void TPZContBufferedStream::clear() {
    fFirst = fBuffer;
    fSize = 0;
    fLast = fBuffer - 1;
}

size_t TPZContBufferedStream::Size() const {
    return fSize;
}
