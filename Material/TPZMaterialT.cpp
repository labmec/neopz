#include "TPZMaterialT.h"
#include "TPZBndCond.h"
#include "Hash/TPZHash.h"
template<class TVar>
TPZBndCondT<TVar>* TPZMaterialT<TVar>::CreateBC(TPZMaterial *reference,
                                                int id, int type,
                                                const TPZFMatrix<TVar> &val1,
                                                const TPZVec<TVar> &val2)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" should NOT be called.";
    PZError<<"All materials should derive from TPZMatBase.";
    PZError<<"Read the documentation for more info.";
    PZError<<"Aborting...";
    DebugStop();
    return nullptr;
}

template<class TVar>
void TPZMaterialT<TVar>::SetForcingFunction(ForcingFunctionType<TVar> f, int pOrder){
    auto *bnd = dynamic_cast<TPZBndCond*>(this);
    if(bnd){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<" ERROR: use SetForcingFunctionBC for boundary materials.\n";
        PZError<<"Aborting...\n";
        DebugStop();
    }
    fForcingFunction = f;
    fForcingFunctionPOrder = pOrder;
}

template<class TVar>
int TPZMaterialT<TVar>::ClassId() const
{
    return Hash("TPZMaterialT") ^
        ClassIdOrHash<TVar>() << 1 ^
        TPZMaterial::ClassId() << 2;
}
template class TPZMaterialT<STATE>;
template class TPZMaterialT<CSTATE>;