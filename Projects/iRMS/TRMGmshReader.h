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
    TPZStack<std::pair<int ,std::string> >  fMaterial;
    
    MaterialDataS() : fMatID(-1), fMaterial(){
        
    }
    
    MaterialDataS(int num) : fMatID(-1), fMaterial(){
        
    }
    
    MaterialDataS(const MaterialDataS &copy) : fMatID(copy.fMatID),
    fMaterial(copy.fMaterial) {
    }
    
    MaterialDataS &operator=(const MaterialDataS &copy){
        fMatID = copy.fMatID;
        fMaterial = copy.fMaterial;
        return *this;
    }
    
};


#endif /* TRMGmshReader_h */
