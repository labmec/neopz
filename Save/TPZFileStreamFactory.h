#ifndef TPZFILESTREAMFACTORY_H
#define TPZFILESTREAMFACTORY_H

#include "TPZFileStream.h"
#include "TPZBFileStream.h"

class TPZFileStreamFactory {
public:
    
    /**
     Returns a stream associated with a file to be read.
     It distinguishes from binary/ascii files based on its version.
     Versioned files will be assumed to be binary and vice-versa.
     
     @param fileName Name of input file
     @return pointer to stream associated with fileName file
     */
    static TPZStream *openReadFileStream(const std::string &fileName){
        TPZBFileStream *file = new TPZBFileStream();
        file->OpenRead(fileName);
        if (file->fFromVersion == 0){
            file->CloseRead();
            delete file;
            TPZFileStream *ascii_file = new TPZFileStream();
            ascii_file->OpenRead(fileName);
            return ascii_file;
        } else {
            return file;
        }
    }
    
    /**
     Returns a stream associated with a binary file for writing.
     Binary files are the standard from version 1 onwards.
     
     @param fileName Name of output BINARY file
     @return pointer to stream associated with fileName file.
     */
    static TPZStream *openWriteFileStream(const std::string &fileName){
        TPZBFileStream *file = new TPZBFileStream();
        file->OpenWrite(fileName);
        return file;
    }
private:
    
};

#endif//TPZFILESTREAMFACTORY_H
