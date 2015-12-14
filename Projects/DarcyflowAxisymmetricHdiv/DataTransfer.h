//
//  DataTransfer.h
//  PZ
//
//  Created by Omar on 12/4/15.
//  Class that implements the data transfer from files to XXX
//

#ifndef DataTransfer_hpp
#define DataTransfer_hpp

#include "tpzautopointer.h"
#include <stdio.h>

class DataTransfer{
    
private:
    
    std::string file_name;
    
public:
    
    /** @brief Default constructor $ */
    DataTransfer();
    
    /** @brief Constructor based on a given file name $ */
    DataTransfer(std::string file_name);
    
    /** @brief Default desconstructor $ */
    ~DataTransfer();
    
    /** @brief Copy constructor $ */
    DataTransfer(const DataTransfer& other);

    void ReadFile();
    
};


#endif /* DataTransfer_h */
