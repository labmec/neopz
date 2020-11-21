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
#include "TPZRenumbering.h"

/**
 * This namespace concentrates general purpose auxiliary methods for computational meshes
 */
namespace TPZCompMeshTools
{
    
    void AddHDivPyramidRestraints(TPZCompMesh *cmesh);
    
    void ExpandHDivPyramidRestraints(TPZCompMesh *cmesh);
    
    void LoadSolution(TPZCompMesh *cpressure, TPZFunction<STATE> &Forcing);

    /// group the elements joining boundary elements and their neighbours
    void GroupElements(TPZCompMesh *cmesh, std::set<int64_t> elbasis, std::set<int64_t> &grouped);
    
    /// Put the element set into a subcompmesh and make the connects internal
    void PutinSubmeshes(TPZCompMesh *cmesh, std::set<int64_t> &elindices, int64_t &index, char KeepOneLagrangian);
    
    /// Put the element set into a subcompmesh and make the connects internal
    // @input elindices
    // elindices.first - identifier of the subdomain
    // elindices.second - indices of the elements that will be transferred to the submesh
    // @output indices
    // indices.first - identifier of the subdomain
    // indices.second - index of the SubCMesh element
    void PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, char KeepOneLagrangian);
    
    /// group elements of the computational mesh
    // for all elements that have dimension of the dimension of cmesh
    // verify if connected elements are included in the current element
    void GroupElements(TPZCompMesh *cmesh);
    
    /// group elements of the computational mesh
    // for all elements in the list
    // group all elements that share a connect with the current element
    // the elements in the list should not share a connect
    // the elements in the mesh cannot share connects with more than one element in the list
    // if either condition isn't met, the program will stop
    void AgglomerateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &ellist, unsigned char lagrangelevel);
    
    /// created condensed elements for the elements that have internal nodes
    void CreatedCondensedElements(TPZCompMesh *cmesh, bool KeepOneLagrangian, bool keepmatrix = false);
    
    /// create a condensed element and do not condense the connect with a given lagrange level or above
    // the method does the same procedure as CreatedCondensedElements, but has different policy for
    // keeping a connect out the condensation loop
    void CondenseElements(TPZCompMesh *cmesh, char LagrangeLevelNotCondensed, bool keepmatrix = false);

    /// ungroup all embedded elements of the computational mesh
    void UnGroupElements(TPZCompMesh *cmesh);
    
    /// uncondensed elements for the elements that have internal nodes
    void UnCondensedElements(TPZCompMesh *cmesh);

    /// compute the norm of the difference between two meshes
    /// put the computed error in the element solution
    void ComputeDifferenceNorm(TPZCompMesh *mesh1, TPZCompMesh *mesh2, TPZVec<STATE> &square_errors);

    /// adjust the polynomial orders of the hdiv elements such that the internal order is higher than the sideorders
    void AdjustFluxPolynomialOrders(TPZCompMesh *fluxmesh, int hdivplusplus);

    /// set the pressure order acording to the order of internal connect of the elements of the fluxmesh
    void SetPressureOrders(TPZCompMesh *fluxmesh, TPZCompMesh *pressuremesh);
    
    void OptimizeBandwidth(TPZCompMesh *cmesh);

    /// Print cmesh solution per geometric element
    void PrintSolutionByGeoElement(TPZCompMesh *cmesh, std::ostream &out);

    /// Print stiffness matrix per geometric element of given material IDs.
    /// If matIDs is empty, all materials are considered.
    void PrintStiffnessMatrixByGeoElement(TPZCompMesh *cmesh, std::ostream &out, std::set<int> matIDs = {});

    void PrintConnectInfoByGeoElement(TPZCompMesh *cmesh, std::ostream &out, std::set<int> matIDs = {},
                                      bool printSeqNumber = true, bool printSolution = true,
                                      bool printLagrangeMult = true);
}

#endif /* defined(__PZ__TPZCompMeshTools__) */
