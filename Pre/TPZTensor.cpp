// $Id: TPZTensor.cpp,v 1.2 2008-03-08 03:12:52 erick Exp $

#include "TPZTensor.h"

template<class T>
STATE TPZTensor<T>::TPZDecomposed::gEigval[3] = {0.};

template class TPZTensor<STATE>;

#include "tfad.h"
template class TPZTensor<TFad<6,double>>;