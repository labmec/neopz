// $Id: TPZTensor.cpp,v 1.2 2008-03-08 03:12:52 erick Exp $

#include "TPZTensor.h"

template<class T>
STATE TPZTensor<T>::TPZDecomposed::gEigval[3] = {T(0.)};

template<>
void TPZTensor<STATE>::Read(TPZStream& buf, void* context) {
    buf.Read(fData);
}

template<>
void TPZTensor<STATE>::Write(TPZStream& buf, int withclassid) const {
    buf.Write(fData);
}

template class TPZTensor<STATE>;