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
    
    
    TPZMHMixedMeshControl(const TPZMHMixedMeshControl &copy) : TPZMHMeshControl(copy)
    {
        
    }
    
    TPZMHMixedMeshControl &operator=(const TPZMHMixedMeshControl &cp)
    {
        fFluxMesh = cp.fFluxMesh;
        TPZMHMeshControl::operator=(cp);
        return *this;
    }
    
    /// Create all data structures for the computational mesh
    void BuildComputationalMesh(bool usersubstructure);
    
    /// will create 1D elements on the interfaces between the coarse element indices
    void CreateCoarseInterfaces(int matid);
    

    TPZAutoPointer<TPZCompMesh> FluxMesh()
    {
        return fFluxMesh;
    }
    
protected:
    
    TPZCompMesh *CreateHDivMHMMesh(TPZGeoMesh * gmesh, int porder);
    
    TPZCompMesh * CreatePressureMHMMesh(TPZGeoMesh * gmesh, int porder);
    
    void DuplicateNeighbouringConnects(TPZCompMesh * HDivMesh);

    TPZCompMesh * CreateHDivPressureMHMMesh(TPZVec<TPZCompMesh * > & cmeshes);

};

#endif /* TPZMHMixedMeshControl_hpp */
