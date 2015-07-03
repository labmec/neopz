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
#include "pzfunction.h"

/// class whose methods implement a functionality on a computational mesh
class TPZCompMeshTools
{
public:
    
    static void AddHDivPyramidRestraints(TPZCompMesh *cmesh);
    
    static void ExpandHDivPyramidRestraints(TPZCompMesh *cmesh);
    
    static void LoadSolution(TPZCompMesh *cpressure, TPZFunction<STATE> &Forcing);

    
};

#endif /* defined(__PZ__TPZCompMeshTools__) */
