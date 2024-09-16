#include "TPZLinearEigenSolver.h"
#include "TPZYSMPPardiso.h"
#include "TPZSYSMPPardiso.h"

#ifdef USING_MKL
template<class TVar>
TPZPardisoSolver<TVar> *
TPZLinearEigenSolver<TVar>::GetPardisoControl(TPZAutoPointer<TPZMatrix<TVar>> mat){
  auto sym = TPZAutoPointerDynamicCast<TPZSYsmpMatrixPardiso<TVar>>(mat);
  if(sym){
    return  &(sym->GetPardisoControl());
  }
  auto nsym = TPZAutoPointerDynamicCast<TPZFYsmpMatrixPardiso<TVar>>(mat);
  if(nsym){
    return  &(nsym->GetPardisoControl());
  }
  return nullptr;
}
#else
template<class TVar>
TPZPardisoSolver<TVar> *
TPZLinearEigenSolver<TVar>::GetPardisoControl(TPZAutoPointer<TPZMatrix<TVar>> mat){
  std::cout<<__PRETTY_FUNCTION__
           <<"\nNeoPZ was not configured with MKL!"
           <<std::endl;
  return nullptr;
}
#endif

template class TPZLinearEigenSolver<float>;
template class TPZLinearEigenSolver<double>;
template class TPZLinearEigenSolver<long double>;

template class TPZLinearEigenSolver<std::complex<float>>;
template class TPZLinearEigenSolver<std::complex<double>>;
template class TPZLinearEigenSolver<std::complex<long double>>;