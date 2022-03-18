#include "TPZMatErrorCombinedSpaces.h"
#include "Hash/TPZHash.h"

//template<class TVar>
//int TPZMatErrorCombinedSpaces<TVar>::ClassId() const
//{
//    return Hash("TPZMatErrorCombinedSpaces") ^
//        TPZMatError<TVar>::ClassId();
//}

template class TPZMatErrorCombinedSpaces<STATE>;
template class TPZMatErrorCombinedSpaces<CSTATE>;
template class TPZMatErrorCombinedSpacesBC<STATE>;
template class TPZMatErrorCombinedSpacesBC<CSTATE>;
