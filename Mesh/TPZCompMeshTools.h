//
//  TPZCompMeshTools.h
//  PZ
//
//  Created by Philippe Devloo on 6/4/15.
//
//

#ifndef __PZ__TPZCompMeshTools__
#define __PZ__TPZCompMeshTools__

#include <stdio.h>
#include "pzcmesh.h"

/// class whose methods implement a functionality on a computational mesh
class TPZCompMeshTools
{
public:
    
    void AddHDivPyramidRestraints(TPZCompMesh *cmesh);
    
};

#endif /* defined(__PZ__TPZCompMeshTools__) */
