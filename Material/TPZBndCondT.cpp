#include "TPZBndCondT.h"
#include "TPZMaterial.h"
#include "TPZMaterialDataT.h"
#include "Hash/TPZHash.h"

#include <iostream>

template<class TVar>
TPZBndCondT<TVar>::TPZBndCondT(int type,
                               const TPZFMatrix<TVar> &val1,
                               const TPZVec<TVar> &val2) :
    TPZRegisterClassId(&TPZBndCond::ClassId),
    TPZBndCond(type), fBCVal1(val1), fBCVal2(val2)
{

}

//! Prints data associated with the material.
template<class TVar>
void TPZBndCondT<TVar>::Print(std::ostream &out) const {
    TPZBndCond::Print(out);
    out << "Has forcing function " << HasForcingFunctionBC() << std::endl;
    out << "Forcing function p order " << fForcingFunctionBCPOrder << std::endl;
    fBCVal1.Print("Val1", out);
    out << "Val2 " << fBCVal2 << std::endl;
}



template<class TVar>
int TPZBndCondT<TVar>::ClassId() const{
    return Hash("TPZBndCondT") ^
        ClassIdOrHash<TVar>() << 1 ^
        TPZBndCond::ClassId() << 2;
}

template<class TVar>
void TPZBndCondT<TVar>::Read(TPZStream& buf, void* context){
    TPZBndCond::Read(buf,context);
}

template<class TVar>
void TPZBndCondT<TVar>::Write(TPZStream& buf, int withclassid) const{
    TPZBndCond::Write(buf,withclassid);
}


template class TPZBndCondT<STATE>;
template class TPZBndCondT<CSTATE>;
