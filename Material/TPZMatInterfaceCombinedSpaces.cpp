//
// Created by Gustavo Batistela on 4/28/21.
//

#include "TPZHash.h"
#include "TPZMaterial.h"
#include "TPZMatInterfaceCombinedSpaces.h"

template<class TVar>
int TPZMatInterfaceCombinedSpaces<TVar>::ClassId() const {
    return Hash("TPZMatInterfaceCombinedSpaces") ^ ClassIdOrHash<TVar>() << 1;
}


template<class TVar>
void TPZMatInterfaceCombinedSpaces<TVar>::ContributeInterface(
    const TPZMaterialDataT<TVar> &data,
    const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
    const std::map<int, TPZMaterialDataT<TVar>> &dataright, REAL weight,
    TPZFMatrix<TVar> &ef)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(),ef.Rows());
    ContributeInterface(data,dataleft,dataright,weight,fakeek,ef);
}

template<class TVar>
void TPZMatInterfaceCombinedSpaces<TVar>::ContributeBCInterface(
    const TPZMaterialDataT<TVar> &data,
    const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
    REAL weight,
    TPZFMatrix<TVar> &ef,
    TPZBndCondT<TVar> &bc)
{
    TPZFMatrix<TVar> fakeek(ef.Rows(),ef.Rows());
    ContributeBCInterface(data,dataleft,weight,fakeek,ef,bc);
}

template<class TVar>
int TPZMatInterfaceCombinedSpaces<TVar>::GetIntegrationOrder(const TPZVec<int> &porder_left, const TPZVec<int> &porder_right) const
{
    int maxl = 0, maxr = 0;
    for (auto porder: porder_left) {
        maxl = std::max(maxl,porder);
    }
    for (auto porder: porder_right) {
        maxr = std::max(maxr,porder);
    }
    return maxl+maxr;
}

// this method is your chance to verify if the material to which this
// BC interface applies is compatible with this boundary interface
// it is called in the method SetMaterial of class TPZBndCondBase
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
SetMaterialImpl(TPZMaterial *mat)
{
    fMatInterface = dynamic_cast<TPZMatInterfaceCombinedSpaces<TVar>*> (mat);
    if (!fMatInterface) {
      PZError << __PRETTY_FUNCTION__;
      PZError << "\nERROR: Invalid reference for creating BC.\nAborting...\n";
      DebugStop();
    }
}


template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
ContributeInterface(const TPZMaterialDataT<TVar> &data,
                    const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                    const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                    REAL weight,
                    TPZFMatrix<TVar> &ek,
                    TPZFMatrix<TVar> &ef)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}
    
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
ContributeInterface(const TPZMaterialDataT<TVar> &data,
                    const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                    const std::map<int, TPZMaterialDataT<TVar>> &dataright,
                    REAL weight,
                    TPZFMatrix<TVar> &ef)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}
//throws
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                      const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                      REAL weight,
                      TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}
//throws
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
ContributeBCInterface(const TPZMaterialDataT<TVar> &data,
                      const std::map<int, TPZMaterialDataT<TVar>> &dataleft,
                      REAL weight,
                      TPZFMatrix<TVar> &ef,
                      TPZBndCondT<TVar> &bc)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}
    
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
FillDataRequirementsInterface(TPZMaterialDataT<TVar> &data,
                              std::map<int, TPZMaterialDataT<TVar>> &datavec_left,
                              std::map<int, TPZMaterialDataT<TVar>> &datavec_right)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}
    
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
SolutionInterface(const TPZMaterialDataT<TVar> &data,
                  const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                  const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                  int var, TPZVec<TVar> &Solout)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}

    
template<class TVar>
void TPZMatInterfaceCombinedSpacesBC<TVar>::
SolutionInterface(const TPZMaterialDataT<TVar> &data,
                  const std::map<int, TPZMaterialDataT<TVar>> &dataleftvec,
                  const std::map<int, TPZMaterialDataT<TVar>> &datarightvec,
                  int var, TPZVec<TVar> &Solout,
                  TPZCompEl *left,TPZCompEl *right)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<< "should not be called! Aborting...\n";
    DebugStop();
}


template class TPZMatInterfaceCombinedSpaces<STATE>; 
template class TPZMatInterfaceCombinedSpaces<CSTATE>;
template class TPZMatInterfaceCombinedSpacesBC<STATE>; 
template class TPZMatInterfaceCombinedSpacesBC<CSTATE>;
