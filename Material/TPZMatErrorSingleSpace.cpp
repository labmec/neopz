#include "TPZMatErrorSingleSpace.h"
#include "Hash/TPZHash.h"

//template<class TVar>
//int TPZMatErrorSingleSpace<TVar>::ClassId() const
//{
//    return Hash("TPZMatErrorSingleSpace") ^
//        TPZMatError<TVar>::ClassId();
//}

template class TPZMatErrorSingleSpace<STATE>;
template class TPZMatErrorSingleSpace<CSTATE>;
template class TPZMatErrorSingleSpaceBC<STATE>;
template class TPZMatErrorSingleSpaceBC<CSTATE>;
