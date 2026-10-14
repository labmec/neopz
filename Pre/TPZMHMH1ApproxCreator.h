
#ifndef PZ_TPZMHMH1CREATOR_H
#define PZ_TPZMHMH1CREATOR_H

#include "TPZMHMApproxCreator.h"
#include "TPZH1ApproxCreator.h"
#include "pzerror.h"

class TPZMHMH1ApproxCreator : public TPZMHMApproxCreator, public TPZH1ApproxCreator
{
public:
    TPZMHMH1ApproxCreator(/* args */);

    TPZMHMH1ApproxCreator(TPZVec<int64_t>& elementPartition, TPZGeoMesh * gmesh);
    
    ~TPZMHMH1ApproxCreator();

protected:
    void AddHybridizationGeoElements() override;
    TPZCompMesh *CreateL2Space() override;
    TPZCompMesh *CreateBoundaryHDivSpace() override;
    TPZCompMesh *CreateConstantSpace(const int &lagMult) override;
    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) override;

    virtual void CreateSkeletonElements();
    TPZCompMesh *CriaMalhaTemporaria();
    virtual void CreateWrapElements();
};

#endif
