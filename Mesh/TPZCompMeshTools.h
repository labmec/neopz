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

    /// group the elements joining boundary elements and their neighbours
    static void GroupElements(TPZCompMesh *cmesh, std::set<long> elbasis, std::set<long> &grouped);
    
    /// Put the element set into a subcompmesh and make the connects internal
    static void PutinSubmeshes(TPZCompMesh *cmesh, std::set<long> &elindices, long &index, bool KeepOneLagrangian);
    
    /// Put the element set into a subcompmesh and make the connects internal
    static void PutinSubmeshes(TPZCompMesh *cmesh, std::map<long,std::set<long> >&elindices, std::set<long> &indices, bool KeepOneLagrangian);
    
    /// group all embedded elements of the computational mesh
    static void GroupElements(TPZCompMesh *cmesh);
    
    /// created condensed elements for the elements that have internal nodes
    static void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian);

    /// compute the norm of the difference between two meshes
    /// put the computed error in the element solution
    static void ComputeDifferenceNorm(TPZCompMesh *mesh1, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors);

};

#endif /* defined(__PZ__TPZCompMeshTools__) */
