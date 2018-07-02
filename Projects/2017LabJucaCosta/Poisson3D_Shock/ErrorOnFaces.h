//
//  ErrorOnFaces.hpp
//  PZ
//
//  Created by labmec on 02/03/18.
//
//

#ifndef ErrorOnFaces_hpp
#define ErrorOnFaces_hpp

#include "pzcmesh.h"
#include "pzanalysis.h"

#include <stdio.h>


bool ComputePressureJumpOnFaces_Hdiv(TPZAnalysis *analysis,int matid,TPZVec<long> &elIndex,TPZVec<int> &sideCoDim1,TPZVec<STATE> &PressureJump);
    
bool ComputeFluxJumpOnFaces_Hdiv(TPZAnalysis *analysis,int matid,TPZVec<long> &elIndex,TPZVec<int> &sideCoDim1,TPZVec<STATE> &PressureJump);
    
#endif /* ErrorOnFaces_hpp */
