//
// Created by Gustavo Batistela on 4/28/21.
//

#include "TPZHash.h"
#include "TPZMaterial.h"
#include "TPZMatInterfaceSingleSpace.h"
#include "TPZBndCondT.h"

template<class TVar>
void TPZMatInterfaceSingleSpace<TVar>::ContributeInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft,
    const TPZMaterialDataT<TVar> &dataright,
    REAL weight, TPZFMatrix<TVar> &ef)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeInterface(data,dataleft,dataright,weight,fakeek,ef);
}

template<class TVar>
void TPZMatInterfaceSingleSpace<TVar>::ContributeBCInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft, REAL weight,
    TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(), ef.Rows(), 0.);
    this->ContributeBCInterface(data,dataleft,weight,fakeek,ef,bc);
}




template<class TVar>
int TPZMatInterfaceSingleSpace<TVar>::ClassId() const {
    return Hash("TPZMatInterfaceSingleSpace") ^ ClassIdOrHash<TVar>() << 1;
}

// this method is your chance to verify if the material to which this
// BC interface applies is compatible with this boundary interface
// it is called in the method SetMaterial of class TPZBndCondBase
template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::SetMaterialImpl(TPZMaterial *mat)
{
    auto fMatInterface = dynamic_cast<TPZMatInterfaceSingleSpace<TVar>*>(mat);

    if(!fMatInterface){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Invalid reference for creating BC.\nAborting...\n";
        DebugStop();
    }
}


template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::SolutionInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft,
    const TPZMaterialDataT<TVar> &dataright,
    const int var,
    TPZVec<TVar> &Solout) {}

template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::ContributeInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft,
    const TPZMaterialDataT<TVar> &dataright,
    REAL weight, TPZFMatrix<TVar> &ek,
    TPZFMatrix<TVar> &ef)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}

template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::ContributeInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft,
    const TPZMaterialDataT<TVar> &dataright,
    REAL weight, TPZFMatrix<TVar> &ef)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}

template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::ContributeBCInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft, REAL weight,
    TPZFMatrix<TVar> &ef,
    TPZBndCondT<TVar> &bc)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}

template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::ContributeBCInterface(
    const TPZMaterialDataT<TVar> &data,
    const TPZMaterialDataT<TVar> &dataleft, REAL weight,
    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
    TPZBndCondT<TVar> &bc)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}
    
template<class TVar>
void
TPZMatInterfaceSingleSpaceBC<TVar>::FillDataRequirementsInterface(
    TPZMaterialData &data) const
{}

template class TPZMatInterfaceSingleSpace<STATE>;
template class TPZMatInterfaceSingleSpace<CSTATE>;
template class TPZMatInterfaceSingleSpaceBC<STATE>;
template class TPZMatInterfaceSingleSpaceBC<CSTATE>;
