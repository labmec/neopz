//
//  TRMGmshReader.hpp
//  PZ
//
//  Created by Omar on 1/15/17.
//
//

#ifndef TRMGmshReader_h
#define TRMGmshReader_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "pzgmesh.h"

class TPZGeoMesh;


struct MaterialDataS {

    int fMatID;
    int fMatID;    
    std::string fMaterial;
    TPZStack<REAL> fProperties;
    
    MaterialDataV() : fMatID(-1), fMaterial(), fProperties()
    {
        
    }
    MaterialDataV(int num) : fMatID(-1), fMaterial(), fProperties()
    {
        
    }
    MaterialDataV(const MaterialDataV &copy) : fMatID(copy.fMatID),
    fMaterial(copy.fMaterial), fProperties(copy.fProperties)
    {
    }
    MaterialDataV &operator=(const MaterialDataV &copy)
    {
        fMatID = copy.fMatID;
        fMaterial = copy.fMaterial;
        fProperties = copy.fProperties;
        return *this;
    }
};


#endif /* TRMGmshReader_h */
