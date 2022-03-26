#include "Elasticity/TPZHybridElasticity2D.h"
#include "pzaxestools.h"


void TPZHybridElasticity2D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    if(data.size() != 2)
    {
        std::cout << "Please implement me!\n";
        DebugStop();
    }
    TPZElasticity2D::Contribute(data[1],weight,ek,ef);
}

void TPZHybridElasticity2D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    if(data.size() != 2)
    {
        std::cout << "Please implement me!\n";
        DebugStop();
    }
    TPZElasticity2D::ContributeBC(data[1],weight,ek,ef,bc);
}

void TPZHybridElasticity2D::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity2D::FillDataRequirements(data[1]);
}

void TPZHybridElasticity2D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity2D::FillBoundaryConditionDataRequirements(type, data[1]);
}




void TPZHybridElasticity2D::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data,
                               int var, TPZVec<STATE> &Solout)
{
    TPZElasticity2D::Solution(data[1], var, Solout);
}


TPZMaterial* TPZHybridElasticity2D::NewMaterial() const
{
    return new TPZHybridElasticity2D(*this);
}


int TPZHybridElasticity2D::ClassId() const
{
    return Hash("TPZHybridElasticity2D") ^
        TBase::ClassId() << 1;
}

void TPZHybridElasticity2D::Read(TPZStream &buf, void *context)
{
    TPZElasticity2D::Read(buf,context);
}
	
void TPZHybridElasticity2D::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity2D::Write(buf,withclassid);
}


void TPZHybridElasticity2D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                             TPZVec<REAL> &values) {
    TPZElasticity2D::Errors(data[1],values);
}
