#include "TPZBFileStream.h"
#include "fad.h"
#include "tfad.h"

TPZBFileStream::TPZBFileStream() {}

TPZBFileStream::~TPZBFileStream() {
    CloseWrite();
    CloseRead();
}

void TPZBFileStream::OpenRead(const std::string &fileName) {
#ifdef PZDEBUG
    if (AmIOpenForRead()) {
        PZError << "TPZBFileStream: File is already opened" << std::endl;
        DebugStop();
    }
#endif
    fIn.open(fileName.c_str(), std::ifstream::binary);
#ifdef PZDEBUG
    if (!AmIOpenForRead()) {
        PZError << "TPZBFileStream: Could not open file" << std::endl;
        DebugStop();
    }
#endif
}
void TPZBFileStream::OpenWrite(const std::string &fileName) {
#ifdef PZDEBUG
    if (AmIOpenForWrite()) {
        PZError << "TPZBFileStream: File is already opened" << std::endl;
        DebugStop();
    }
#endif
    fOut.open(fileName.c_str(), std::ofstream::binary | std::ofstream::trunc);
#ifdef PZDEBUG
    if (!AmIOpenForWrite()) {
        PZError << "TPZBFileStream: Could not open file" << std::endl;
        DebugStop();
    }
#endif
//    std::string fileInfo("FileVersion");
//    fOut.write(fileInfo.c_str(), fileInfo.length());
//    const uint64_t temp = fCurrentVersion;
//    fOut.write(reinterpret_cast<const char *>(&temp), sizeof(temp));
}

bool TPZBFileStream::AmIOpenForRead(){
    return fIn.is_open();
}

bool TPZBFileStream::AmIOpenForWrite(){
    return fOut.is_open();
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

template
void TPZBFileStream::ReadData<double>(double *p, int howMany);

template
void TPZBFileStream::WriteData<double>(const double *p, int howMany);

template
void TPZBFileStream::WriteData<int>(int const* p, int howMany);

template
void TPZBFileStream::WriteData<long const>(long const* p, int howMany);

template
void TPZBFileStream::ReadData<int>(int* p, int howMany);

template
void TPZBFileStream::ReadData<long>(long* p, int howMany);
