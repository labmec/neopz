#include "TPZCircBufferedStream.h"
#include <iostream>
#include <stdexcept>
#include <string.h>

#ifdef _AUTODIFF
#include "fad.h"
#endif

TPZCircBufferedStream::TPZCircBufferedStream(){

    fNAllocatedBytes = MIN_SIZE_INCREMENT;
    fBuffer = new char[fNAllocatedBytes];
    fFirst = fBuffer;
    fSize = 0;
    fLast = fBuffer - 1;
}

TPZCircBufferedStream::TPZCircBufferedStream(const TPZStream &other)
    : TPZStream(other), fBuffer(NULL) {
    *this = other;
}

TPZCircBufferedStream &TPZCircBufferedStream::
operator=(const TPZCircBufferedStream &other) {
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

TPZCircBufferedStream::~TPZCircBufferedStream() { delete[] fBuffer; }

TPZCircBufferedStream &TPZCircBufferedStream::operator<<(TPZCircBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char *temp = new char[nBytesOther];
    other.ReadFromBuffer(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
	delete[] temp;
    return *this;
}

TPZCircBufferedStream &TPZCircBufferedStream::
operator<<(const TPZCircBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char *temp = new char[nBytesOther];
    other.ConstRead(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
	delete[] temp;
    return *this;
}

void TPZCircBufferedStream::ReadFromBuffer(char *dest, const size_t &nBytes) {
    if (nBytes > fSize) {
        std::string msg("TPZCircBufferedStream: Cannot read ");
        msg.append(std::to_string(nBytes));
        msg.append(" bytes; there are only ");
        msg.append(std::to_string(fSize));
        msg.append(" available.");
        PZError << msg << std::endl;
        DebugStop();
    }
    char *endBuffer = fBuffer + fNAllocatedBytes;

    if (fFirst + nBytes < endBuffer) {
        // direct reading (we do not need to cycle to the beginning of the
        // buffer to read)
        memcpy(dest, fFirst, nBytes);
        fFirst += nBytes;
    } else {
        // we need to read past the end of the buffer.
        const size_t nBytesRead(endBuffer - fFirst);
        memcpy(dest, fFirst, nBytesRead);
        if (nBytes != nBytesRead) {
            memcpy(dest + nBytesRead, fBuffer, nBytes - nBytesRead);
        }
        fFirst = fBuffer + nBytes - nBytesRead;
    }

    fSize -= nBytes;
}

void TPZCircBufferedStream::ConstRead(char *dest, const size_t &nBytes) const {
    ConstReadFromBuffer(dest, nBytes);
}

void TPZCircBufferedStream::ConstReadFromBuffer(char *dest,
                                            const size_t &nBytes) const {
    if (nBytes > fSize) {
        std::string msg("TPZCircBufferedStream: Cannot read ");
        msg.append(std::to_string(nBytes));
        msg.append(" bytes; there are only ");
        msg.append(std::to_string(fSize));
        msg.append(" available.");
        PZError << msg << std::endl;
        DebugStop();
    }
    char *endBuffer = fBuffer + fNAllocatedBytes;

    if (fFirst + nBytes < endBuffer) {
        // direct reading (we do not need to cycle to the beginning of the
        // buffer to read)
        memcpy(dest, fFirst, nBytes);
    } else {
        // we need to read past the end of the buffer.
        const size_t nBytesRead(endBuffer - fFirst);
        memcpy(dest, fFirst, nBytesRead);
        if (nBytes != nBytesRead) {
            memcpy(dest + nBytesRead, fBuffer, nBytes - nBytesRead);
        }
    }
}

void TPZCircBufferedStream::WriteToBuffer(const char *source,
                                      const size_t &nBytes) {
    if (fSize + nBytes > fNAllocatedBytes) {
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
        char *endBuffer = fBuffer + fNAllocatedBytes;

        if (fLast + nBytes < endBuffer) {
            // direct writing (we do not need to cycle to the beginning of the
            // buffer to write)
            memcpy(fLast + 1, source, nBytes);
            fLast += nBytes;
        } else {
            // we need to write past the end of the buffer.
            const size_t nBytesWritten(endBuffer - fLast - 1);
            memcpy(fLast + 1, source, nBytesWritten);
            if (nBytes != nBytesWritten) {
                memcpy(fBuffer, source + nBytesWritten, nBytes - nBytesWritten);
            }
            fLast = fBuffer - 1 + nBytes - nBytesWritten;
        }

        fSize += nBytes;
    }
}

void TPZCircBufferedStream::Print() {
    std::cout << "fSize=" << fSize << std::endl;
    double *temp = new double[fSize / 8];
    ConstRead(reinterpret_cast<char *>(temp), fSize);
    for (unsigned int i = 0; i < fSize / 8; ++i) {
        std::cout << temp[i] << " ";
    }
	delete[] temp;
    std::cout << std::endl;
}

template <class T> void TPZCircBufferedStream::ReadData(T *p, int howMany) {
    ReadFromBuffer(reinterpret_cast<char *>(p), howMany * sizeof(T));
}

template <class T> void TPZCircBufferedStream::WriteData(const T *p, int howMany) {
    WriteToBuffer(reinterpret_cast<const char *>(p), howMany * sizeof(T));
}
