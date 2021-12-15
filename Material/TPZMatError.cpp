#include "TPZMatError.h"
#include "Hash/TPZHash.h"

template<class TVar>
void TPZMatError<TVar>::SetExactSol(ExactSolType<TVar> f, int pOrder){
    fExactSol = f;
    fExactPOrder = pOrder;
}

//template<class TVar>
//int TPZMatError<TVar>::ClassId() const{
//    return Hash("TPZMatError") ^
//        ClassIdOrHash<TVar>() << 1;
//}

template class TPZMatError<STATE>;
template class TPZMatError<CSTATE>;
