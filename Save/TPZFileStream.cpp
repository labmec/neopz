#include "TPZFileStream.h"
#include "pzerror.h"

#include "fad.h"
#include "tfad.h"

TPZFileStream::TPZFileStream() {
}

TPZFileStream::~TPZFileStream() {
    CloseWrite();
    CloseRead();
}

void TPZFileStream::OpenRead(const std::string &fileName) {
#ifdef PZDEBUG
    if (AmIOpenForRead()) {
        PZError<<"TPZFileStream: File is already opened"<<std::endl;
    }
#endif
    fIn.open(fileName.c_str());
#ifdef PZDEBUG
    if (!AmIOpenForRead()) {
        PZError<<"TPZFileStream: Could not open file"<<std::endl;
    }
#endif
}

void TPZFileStream::OpenWrite(const std::string &fileName) {
#ifdef PZDEBUG
    if (AmIOpenForWrite()) {
        PZError<<"TPZFileStream: File is already opened"<<std::endl;
    }
#endif
    fOut.open(fileName.c_str(), std::ofstream::binary | std::ofstream::trunc);
#ifdef PZDEBUG
    if (!AmIOpenForWrite()) {
        PZError<<"TPZFileStream: Could not open file"<<std::endl;
    }
#endif
}

bool TPZFileStream::AmIOpenForRead(){
    return fIn.is_open();
}

bool TPZFileStream::AmIOpenForWrite(){
    return fOut.is_open();
}

void TPZFileStream::CloseRead() {
    if (fIn.is_open()) {
        fIn.close();
    }
}

void TPZFileStream::CloseWrite() {
    if (fOut.is_open()) {
        fOut.close();
    }
}

template<class T>
    void TPZFileStream::ReadData(T *p, int howMany) {
        int c;
        char buf[100];
        if(!fIn)
        {
            DebugStop();
        }
        if(howMany)
        {
            for(c=0; c<howMany; c++) fIn >> p[c];
            fIn.getline(buf,100);
        }
#ifdef PZDEBUG
        if (fIn.bad()) {
            PZError << "TFileStream:Could not read from stream" << std::endl;
            DebugStop();
        }
#endif
    }

template<class T>
    void  TPZFileStream::WriteData(const T *p, int howMany){
        for(int c=0; c<howMany; c++) fOut << p[c] << std::endl;
#ifdef PZDEBUG
        if (fOut.bad()) {
            PZError << "TFileStream:Could not write to stream" << std::endl;
            DebugStop();
        }
#endif
    }

template void TPZFileStream::ReadData<double>(double *p, int howMany);
template void TPZFileStream::WriteData<double>(double const* p, int howMany);
template void TPZFileStream::WriteData<int>(int const* p, int howMany);
template void TPZFileStream::WriteData<long const>(long const* p, int howMany);
template void TPZFileStream::ReadData<int>(int* p, int howMany);
template void TPZFileStream::ReadData<long>(long* p, int howMany);
