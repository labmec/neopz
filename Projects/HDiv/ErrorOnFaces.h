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

/**
 * Faces contains side of codimension 1 however the side is boundary or has one neighboard on the other side.
 * AnotherSideFaces contains NULL if corresponding Faces is boundary or the element and side neighboard on the other side of the faces.
 * Then Faces is a boundary of the computational elements of the higher level.
 **/
bool IdentifyingFaces(TPZCompMesh *cmesh,TPZStack<TPZCompElSide> &Faces, TPZStack<TPZCompElSide> &AnotherSideFaces);

//bool ComputePressureJumpOnFaces_Hdiv(TPZAnalysis *analysis,int matid,TPZVec<long> &elIndex,TPZVec<int> &sideCoDim1,TPZVec<STATE> &PressureJump);
bool ComputePressureJumpOnFaces(TPZCompMesh *cmeshpressure,int matid,TPZStack<TPZCompElSide> &Faces, TPZStack<TPZCompElSide> &AnotherSideFaces,STATE &Error);
bool ComputePressureJumpOnFaces(TPZCompMesh *cmeshpressure,int matid,STATE &Error);

bool ComputeFluxJumpOnFaces_Hdiv(TPZAnalysis *analysis,int matid,TPZVec<long> &elIndex,TPZVec<int> &sideCoDim1,TPZVec<STATE> &PressureJump);
    
#endif /* ErrorOnFaces_hpp */
