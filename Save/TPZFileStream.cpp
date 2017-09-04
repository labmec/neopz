#include "TPZFileStream.h"
#include "pzerror.h"

TPZFileStream::TPZFileStream() {
}

TPZFileStream::~TPZFileStream() {
    CloseWrite();
    CloseRead();
}

void TPZFileStream::OpenRead(const std::string &fileName) {
#ifdef PZDEBUG
    if (fIn.is_open()) {
        PZError<<"TPZFileStream: File is already opened"<<std::endl;
    }
#endif
    fIn.open(fileName.c_str());
#ifdef PZDEBUG
    if (!fIn.is_open()) {
        PZError<<"TPZFileStream: Could not open file"<<std::endl;
    }
#endif
    fFromVersion = 0;
}

void TPZFileStream::OpenWrite(const std::string &fileName) {
#ifdef PZDEBUG
    if (fOut.is_open()) {
        PZError<<"TPZFileStream: File is already opened"<<std::endl;
    }
#endif
    fOut.open(fileName.c_str(), std::ofstream::binary | std::ofstream::trunc);
#ifdef PZDEBUG
    if (!fOut.is_open()) {
        PZError<<"TPZFileStream: Could not open file"<<std::endl;
    }
#endif
    fOut.precision(15);
    std::string fileInfo("FileVersion");
    fOut.write(fileInfo.c_str(), fileInfo.length());
    const unsigned long temp = fCurrentVersion;
    fOut.write(reinterpret_cast<const char *> (&temp), sizeof (temp));
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
    }

template<class T>
    void  TPZFileStream::WriteData(const T *p, int howMany)
    {
        for(int c=0; c<howMany; c++) fOut << p[c] << std::endl;
    }
