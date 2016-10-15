//
//  TPZMHMixedMeshControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#ifndef TPZMHMixedMeshControl_hpp
#define TPZMHMixedMeshControl_hpp

#include <stdio.h>

#include "TPZMHMeshControl.h"

class TPZMHMixedMeshControl : public TPZMHMeshControl
{
    
    
    /// computational mesh to contain the pressure elements
    // this mesh is the same as fCMesh if there are no lagrange multipliers assocated with the average pressure
    TPZAutoPointer<TPZCompMesh> fFluxMesh;
    

public:
    
    TPZMHMixedMeshControl() : TPZMHMeshControl()
    {
        
    }
    
    TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices);
    
    TPZMHMixedMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &coarseindices);
    
    
    TPZMHMixedMeshControl(const TPZMHMixedMeshControl &copy) : TPZMHMeshControl(copy)
    {
        fFluxMesh = copy.fFluxMesh;
    }
    
    TPZMHMixedMeshControl &operator=(const TPZMHMixedMeshControl &cp)
    {
        fFluxMesh = cp.fFluxMesh;
        TPZMHMeshControl::operator=(cp);
        return *this;
    }
    
    /// Create all data structures for the computational mesh
    void BuildComputationalMesh(bool usersubstructure);
    
    TPZAutoPointer<TPZCompMesh> FluxMesh()
    {
        return fFluxMesh;
    }
    
    /// Set the flag for creating Lagrange Dofs for the average pressure
    void SetLagrangeAveragePressure(bool flag)
    {
    }
    

    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        meshvec.Resize(2);
        meshvec[0] = fFluxMesh.operator->();
        meshvec[1] = fPressureFineMesh.operator->();
    }

    /// print the data structure
    void Print(std::ostream &out);

protected:
    
    TPZCompMesh *CreateHDivMHMMesh();
    
    TPZCompMesh * CreatePressureMHMMesh();
    
    void DuplicateNeighbouringConnects();

    TPZCompMesh * CreateHDivPressureMHMMesh();

    void HideTheElements();

};

#endif /* TPZMHMixedMeshControl_hpp */
