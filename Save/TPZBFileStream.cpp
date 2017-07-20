#include "TPZBFileStream.h"
#include "pzerror.h"
#include <exception>
#include <iostream>
#include <string>

TPZBFileStream::TPZBFileStream() {}

TPZBFileStream::~TPZBFileStream() {
    CloseWrite();
    CloseRead();
}

void TPZBFileStream::OpenRead(const std::string &fileName) {
#ifdef PZDEBUG
    if (fIn.is_open()) {
        PZError << "TPZBFileStream: File is already opened" << std::endl;
        DebugStop();
    }
#endif
    fIn.open(fileName.c_str(), std::ifstream::binary);
#ifdef PZDEBUG
    if (!fIn.is_open()) {
        PZError << "TPZBFileStream: Could not open file" << std::endl;
        DebugStop();
    }
#endif

    std::string versionString("FileVersion");
    char versionRead[12];
    fIn.read(versionRead, 11); // reads header of inputfile
    versionRead[11] = '\0';    // terminates c-style string

    if (versionString.compare(versionRead) == 0) { // versioned file
        fIn.read(reinterpret_cast<char *>(&fFromVersion), sizeof(fFromVersion));
    } else { // unversioned file aka V0
        fFromVersion = 0;
        fIn.seekg(0, fIn.beg); // goes back to beginning of file
    }
}
void TPZBFileStream::OpenWrite(const std::string &fileName) {
#ifdef PZDEBUG
    if (fOut.is_open()) {
        PZError << "TPZBFileStream: File is already opened" << std::endl;
        DebugStop();
    }
#endif
    fOut.open(fileName.c_str(), std::ofstream::binary | std::ofstream::trunc);
#ifdef PZDEBUG
    if (!fOut.is_open()) {
        PZError << "TPZBFileStream: Could not open file" << std::endl;
        DebugStop();
    }
#endif
    std::string fileInfo("FileVersion");
    fOut.write(fileInfo.c_str(), fileInfo.length());
    const unsigned long temp = fCurrentVersion;
    fOut.write(reinterpret_cast<const char *>(&temp), sizeof(temp));
}

void TPZBFileStream::CloseRead() {
    if (fIn.is_open()) {
        fIn.close();
    }
}

void TPZBFileStream::CloseWrite() {
    if (fOut.is_open()) {
        fOut.close();
    }
}

/** @brief Reads howMany objects of the class T from pointer location p */
template <class T> void TPZBFileStream::ReadData(T *p, int howMany) {
    fIn.read(reinterpret_cast<char *>(p), howMany * sizeof(T));
#ifdef PZDEBUG
    if (fIn.bad()) {
        PZError << "TBFileStream:Could not read from stream" << std::endl;
        DebugStop();
    }
#endif
}

/** @brief Reads howMany objects of the class T from pointer location p */
template <class T> void TPZBFileStream::WriteData(const T *p, int howMany) {
    fOut.write(reinterpret_cast<const char *>(p), howMany * sizeof(T));
#ifdef PZDEBUG
    if (fOut.bad()) {
        PZError << "TBFileStream:Could not write to stream" << std::endl;
        DebugStop();
    }
#endif
}

template void TPZBFileStream::ReadData<int>(int *p, int howMany);
template void TPZBFileStream::ReadData<unsigned int>(unsigned int *p, int howMany);
template void TPZBFileStream::ReadData<long>(long *p, int howMany);
template void TPZBFileStream::ReadData<float>(float *p, int howMany);
template void TPZBFileStream::ReadData<double>(double *p, int howMany);
template void TPZBFileStream::ReadData<long double>(long double *p, int howMany);
template void TPZBFileStream::ReadData<char>(char *p, int howMany);
template void TPZBFileStream::ReadData<std::complex<float> >(std::complex<float> *p, int howMany);
template void TPZBFileStream::ReadData<std::complex<double> >(std::complex<double> *p, int howMany);
template void TPZBFileStream::ReadData<std::complex<long double> >(std::complex<long double>  *p, int howMany);
#ifdef _AUTODIFF
template void TPZBFileStream::ReadData<Fad <float> >(Fad <float>  *p, int howMany);
template void TPZBFileStream::ReadData<Fad <double> >(Fad <double>  *p, int howMany);
template void TPZBFileStream::ReadData<Fad<long double> >(Fad<long double>  *p, int howMany);
#endif

template void TPZBFileStream::WriteData<int>(const int *p, int howMany);
template void TPZBFileStream::WriteData<unsigned int>(const unsigned int *p, int howMany);
template void TPZBFileStream::WriteData<long>(const long *p, int howMany);
template void TPZBFileStream::WriteData<float>(const float *p, int howMany);
template void TPZBFileStream::WriteData<double>(const double *p, int howMany);
template void TPZBFileStream::WriteData<long double>(const long double *p, int howMany);
template void TPZBFileStream::WriteData<char>(const char *p, int howMany);
template void TPZBFileStream::WriteData<std::complex<float> >(const std::complex<float> *p, int howMany);
template void TPZBFileStream::WriteData<std::complex<double> >(const std::complex<double> *p, int howMany);
template void TPZBFileStream::WriteData<std::complex<long double> >(const std::complex<long double>  *p, int howMany);
#ifdef _AUTODIFF
template void TPZBFileStream::WriteData<Fad <float> >(Fad <float>  *p, int howMany);
template void TPZBFileStream::WriteData<Fad <double> >(Fad <double>  *p, int howMany);
template void TPZBFileStream::WriteData<Fad<long double> >(Fad<long double>  *p, int howMany);
#endif
