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
    fFirst = fBuffer;
    fLast = fBuffer - 1 + fSize;
    return *this;
}

TPZContBufferedStream::~TPZContBufferedStream() {
    delete[] fBuffer;
}

TPZContBufferedStream &TPZContBufferedStream::
operator<<(TPZContBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char temp[nBytesOther];
    other.ReadFromBuffer(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
    return *this;
}

TPZContBufferedStream &TPZContBufferedStream::
operator<<(const TPZContBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char temp[nBytesOther];
    other.ConstRead(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
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
    double temp[fSize / 8];
    ConstRead(reinterpret_cast<char *> (temp), fSize);
    for (unsigned int i = 0; i < fSize / 8; ++i) {
        std::cout << temp[i] << " ";
    }
    std::cout << std::endl;
}

template <class T> void TPZContBufferedStream::ReadData(T *p, int howMany) {
    ReadFromBuffer(reinterpret_cast<char *> (p), howMany * sizeof (T));
}

template <class T>
void TPZContBufferedStream::WriteData(const T *p, int howMany) {
    WriteToBuffer(reinterpret_cast<const char *> (p), howMany * sizeof (T));
}

void TPZContBufferedStream::GetDataFromBuffer(char *dest) const {
    memcpy(dest, fFirst, fSize);
}

void TPZContBufferedStream::clear(){
    fFirst = fBuffer;
    fSize = 0;
    fLast = fBuffer - 1;
}

size_t TPZContBufferedStream::Size() const {
    return fSize;
}
