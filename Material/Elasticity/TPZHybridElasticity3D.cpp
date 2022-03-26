#include "Elasticity/TPZHybridElasticity3D.h"
#include "pzaxestools.h"


void TPZHybridElasticity3D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    if(data.size() != 2)
    {
        std::cout << "Please implement me!\n";
        DebugStop();
    }
    TPZElasticity3D::Contribute(data[1],weight,ek,ef);
}

void TPZHybridElasticity3D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    if(data.size() != 2)
    {
        std::cout << "Please implement me!\n";
        DebugStop();
    }
    TPZElasticity3D::ContributeBC(data[1],weight,ek,ef,bc);
}

void TPZHybridElasticity3D::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity3D::FillDataRequirements(data[1]);
}

void TPZHybridElasticity3D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity3D::FillBoundaryConditionDataRequirements(type, data[1]);
}




void TPZHybridElasticity3D::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data,
                               int var, TPZVec<STATE> &Solout)
{
    TPZElasticity3D::Solution(data[1], var, Solout);
}


TPZMaterial* TPZHybridElasticity3D::NewMaterial() const
{
    return new TPZHybridElasticity3D(*this);
}


int TPZHybridElasticity3D::ClassId() const
{
    return Hash("TPZHybridElasticity3D") ^
        TBase::ClassId() << 1;
}

void TPZHybridElasticity3D::Read(TPZStream &buf, void *context)
{
    TPZElasticity3D::Read(buf,context);
}
    
void TPZHybridElasticity3D::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity3D::Write(buf,withclassid);
}


void TPZHybridElasticity3D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                             TPZVec<REAL> &values) {
    TPZElasticity3D::Errors(data[1],values);
}

