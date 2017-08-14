#include "TPZBufferedStream.h"
#include <iostream>
#include <stdexcept>
#include <iostream>
#include <string.h>
TPZBufferedStream::TPZBufferedStream(TPZStream &fUnderlyingStream) : TPZStream(fUnderlyingStream.fFromVersion), fUnderlyingStream(fUnderlyingStream) {

    fNAllocatedBytes = MIN_SIZE_INCREMENT;
    fBuffer = new char[fNAllocatedBytes];
    fFirst = fBuffer;
    fSize = 0;
    fLast = fBuffer - 1;
    fReadFromUnderlyingStream = true;
}

TPZBufferedStream::TPZBufferedStream(const TPZBufferedStream &other) : TPZStream(other), fBuffer(NULL), fUnderlyingStream(other.fUnderlyingStream) {
    *this = other;
}

TPZBufferedStream &TPZBufferedStream::operator=(const TPZBufferedStream &other) {
    fUnderlyingStream = other.fUnderlyingStream;
    fNAllocatedBytes = other.fNAllocatedBytes;
    if (fBuffer) delete[] fBuffer;
    fBuffer = new char[fNAllocatedBytes];
    memcpy(fBuffer, other.fBuffer, other.fSize);
    fSize = other.fSize;
    fFirst = fBuffer;
    fLast = fBuffer - 1 + fSize;
    fReadFromUnderlyingStream = other.fReadFromUnderlyingStream;
    return *this;
}

TPZBufferedStream::~TPZBufferedStream() {
    delete []fBuffer;
}

TPZBufferedStream &TPZBufferedStream::operator<<(TPZBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char temp[nBytesOther];
    other.ReadFromBuffer(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
    return *this;
}

TPZBufferedStream &TPZBufferedStream::operator<<(const TPZBufferedStream &other) {
    const unsigned int nBytesOther = other.fSize;
    char temp[nBytesOther];
    other.ConstRead(temp, nBytesOther);
    WriteToBuffer(temp, nBytesOther);
    return *this;
}

void TPZBufferedStream::Read(double *p, const int &size){
    if (fReadFromUnderlyingStream) {
        fUnderlyingStream.Read(p, size);
    } else {
        ReadFromBuffer(reinterpret_cast<char *> (p), size*sizeof (double));
    }
}

void TPZBufferedStream::Write(const double *var, const int &size){
    WriteToBuffer(reinterpret_cast<const char *> (var), size*sizeof (double));
}

void TPZBufferedStream::ReadFromBuffer(char *dest, const size_t &nBytes) {
    if (nBytes > fSize) {
        std::string msg("TPZBufferedStream: Cannot read ");
        msg.append(std::to_string(nBytes));
        msg.append(" bytes; there are only ");
        msg.append(std::to_string(fSize));
        msg.append(" available.");
        throw std::runtime_error(msg);
    }
    char *endBuffer = fBuffer + fNAllocatedBytes;

    if (fFirst + nBytes < endBuffer) {
        // direct reading (we do not need to cycle to the beginning of the buffer to read)
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

void TPZBufferedStream::ConstRead(char *dest, const size_t &nBytes) const {
    if (fReadFromUnderlyingStream) {
        throw std::runtime_error("TPZBufferedStream: We are still reading from the "
                "underlying stream and cannot guarantee that it can be const read!");
    }
    ConstReadFromBuffer(dest, nBytes);
}

void TPZBufferedStream::ConstReadFromBuffer(char *dest, const size_t &nBytes) const {
    if (nBytes > fSize) {
        std::string msg("TPZBufferedStream: Cannot read ");
        msg.append(std::to_string(nBytes));
        msg.append(" bytes; there are only ");
        msg.append(std::to_string(fSize));
        msg.append(" available.");
        throw std::runtime_error(msg);
    }
    char *endBuffer = fBuffer + fNAllocatedBytes;

    if (fFirst + nBytes < endBuffer) {
        // direct reading (we do not need to cycle to the beginning of the buffer to read)
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

void TPZBufferedStream::WriteToBuffer(const char *source, const size_t &nBytes) {
    if (fSize + nBytes > fNAllocatedBytes) {
        const size_t oldSize = fSize;
        const size_t newAllocatedBytes = oldSize * 1.1 + nBytes + MIN_SIZE_INCREMENT;//10% increase + nBytes + 
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
            // direct writing (we do not need to cycle to the beginning of the buffer to write)
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

void TPZBufferedStream::Print() {
    std::cout << "fSize=" << fSize << std::endl;
    double temp[fSize / 8];
    ConstRead(reinterpret_cast<char*> (temp), fSize);
    for (unsigned int i = 0; i < fSize / 8; ++i) {
        std::cout << temp[i] << " ";
    }
    std::cout << std::endl;
}

void TPZBufferedStream::BeginUpdate() {
}

void TPZBufferedStream::EndUpdate(const unsigned long &new_version) {
    fFromVersion = new_version;
    fReadFromUnderlyingStream = false;
}
