#include "TPZMatErrorSingleSpace.h"
#include "Hash/TPZHash.h"

template<class TVar>
int TPZMatErrorSingleSpace<TVar>::ClassId() const
{
    return Hash("TPZMatErrorSingleSpace") ^
        TPZMatError<TVar>::ClassId();
}

template<class TVar>
void TPZMatErrorSingleSpace<TVar>::Errors(const TPZMaterialDataT<TVar> &data,
                                          TPZVec<REAL> &errors)
{
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
    return Errors(data.x, data.sol[0],
                  data.dsol[0], data.axes,
                  errors );
}
template class TPZMatErrorSingleSpace<STATE>;
template class TPZMatErrorSingleSpace<CSTATE>;
template class TPZMatErrorSingleSpaceBC<STATE>;
template class TPZMatErrorSingleSpaceBC<CSTATE>;