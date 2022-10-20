

#ifndef PZ_TPZMHMHDIVCREATOR_H
#define PZ_TPZMHMHDIVCREATOR_H

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMHMApproxCreator.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxCreator.h"


class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;

class TPZMHMHDivApproxCreator : public TPZHDivApproxCreator, public TPZMHMApproxCreator {

private:
    int fAvPresLevel;
    int fDistFluxLevel;
    int fPOrderSkeleton;
    int fSkeletonMatId = -12345;

public:

    TPZMHMHDivApproxCreator(TPZGeoMesh *gmesh);
    
    TPZMHMHDivApproxCreator(TPZGeoMesh* gmesh, TPZVec<int64_t>& elPartition);

    ~TPZMHMHDivApproxCreator(){};
    
    TPZMultiphysicsCompMesh * BuildMultiphysicsCMesh();

    void InsertMaterialObjects(TPZAnalyticSolution &analytic);

    void PutinSubstructures(TPZCompMesh &cmesh);

    void CondenseElements(TPZCompMesh &cmesh);

    /// Set skeleton default polynomial order
    void SetPOrderSkeleton(const int ord){
        fPOrderSkeleton = ord;
    }

    /// Get skeleton default polynomial order
    int &GetPOrderSkeleton(){return fPOrderSkeleton;}

    
};

#endif
