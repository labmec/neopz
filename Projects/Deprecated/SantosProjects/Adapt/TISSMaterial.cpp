//
//  TISSMaterial.cpp
//  PZ
//
//  Created by Thiago Dias dos Santos on 8/3/15.
//
//

#include "TISSMaterial.h"


TISSMaterial::TISSMaterial(int matid) : TPZMaterial(matid)
{
    
}

/** @brief Default constructor */
TISSMaterial::TISSMaterial() : TPZMaterial()
{
    
}


TISSMaterial::TISSMaterial(const TISSMaterial &mat) : TPZMaterial(mat)
{
    //nothing here
}

TISSMaterial::~TISSMaterial()
{
    //nothing here
}

void TISSMaterial::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    //nothing here
    DebugStop();
}

void TISSMaterial::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    //nothing here
    DebugStop();
}

void TISSMaterial::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var) );
    
    if(var == 1){
    
        Solout[0] = data.sol[0][0];
        return;
    }
    
    if(var == 2){
        
        Solout[0] = data.sol[0][1];
        return;
    }
    
    if(var == 3){
        
        Solout[0] = data.sol[0][2];
        return;
    }
    
    if(var == 4){
        
        Solout[0] = data.sol[0][3];
        return;
    }
    
    if(var == 5){
        
        Solout[0] = data.sol[0][4];
        return;
    }
    
    if(var == 6){
        
        Solout[0] = data.sol[0][5];
        return;
    }
    
    if(var == 7){
        
        Solout[0] = data.sol[0][6];
        return;
    }
    
    if(var == 8){
        
        Solout[0] = data.sol[0][7];
        return;
    }
 
    
}

/** Returns the variable index associated with the name */
int TISSMaterial::VariableIndex(const std::string &name){
    if(!strcmp("Surface",name.c_str()))           return 1;
    if(!strcmp("Base",name.c_str()))       return 2;
    if(!strcmp("Bed",name.c_str()))      return 3;
    if(!strcmp("Pressure",name.c_str()))      return 4;
    if(!strcmp("Temperature",name.c_str()))      return 5;
    if(!strcmp("Vx",name.c_str()))  return 6;
    if(!strcmp("Vy",name.c_str()))      return 7;
    if(!strcmp("MaskLevelSet",name.c_str()))      return 8;
    
    return TPZMaterial::VariableIndex(name);
}

int TISSMaterial::NSolutionVariables(int var){
    if(var == 1) return 1;
    if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    if(var == 7) return 1;
    if(var == 8) return 1;

    return TPZMaterial::NSolutionVariables(var);
}