/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrix methods.
 */

#ifdef USING_MKL
#include "TPZSYSMPPardiso.h"
#include "pzfmatrix.h"
// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrixPardiso<TVar>::TPZSYsmpMatrixPardiso() : TPZRegisterClassId(&TPZSYsmpMatrixPardiso::ClassId),
TPZSYsmpMatrix<TVar>() {

}

template<class TVar>
void TPZSYsmpMatrixPardiso<TVar>::CopyFrom(const TPZMatrix<TVar> *  mat)
{                                                           
  auto *from = dynamic_cast<const TPZSYsmpMatrixPardiso<TVar> *>(mat);
  if (from) {
    *this = *from;
  }
  else
  {
    auto *from = dynamic_cast<const TPZSYsmpMatrix<TVar> *>(mat);
    if (from && from->IsDecomposed() == ENoDecompose) {
      *this = *from;
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Called with incompatible type\n.";
    PZError<<"Aborting...\n";
    DebugStop();
  }
}

template<class TVar>
int TPZSYsmpMatrixPardiso<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrixPardiso") ^ TPZSYsmpMatrix<TVar>::ClassId() << 1;
}


template<class TVar>
void TPZSYsmpMatrixPardiso<TVar>::SetIsDecomposed(DecomposeType val){
  TPZBaseMatrix::SetIsDecomposed(val);
  if(val){fPardisoControl.fDecomposed = true;}
}

template<class TVar>
int TPZSYsmpMatrixPardiso<TVar>::Decompose(const DecomposeType dt)
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }

  if(!fPardisoControl.HasCustomSettings()){
    typename TPZPardisoSolver<TVar>::MStructure str =
      TPZPardisoSolver<TVar>::MStructure::ESymmetric;
    typename TPZPardisoSolver<TVar>::MSystemType sysType =
      TPZPardisoSolver<TVar>::MSystemType::ESymmetric;
    typename TPZPardisoSolver<TVar>::MProperty prop =
      this->IsDefPositive() ?
      TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
      TPZPardisoSolver<TVar>::MProperty::EIndefinite;
    fPardisoControl.SetStructure(str);
    fPardisoControl.SetMatrixType(sysType,prop);
  }
  fPardisoControl.Decompose(this);
  this->SetIsDecomposed(dt);
  return 0;
}

template<class TVar>
int TPZSYsmpMatrixPardiso<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt)
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }
  if(this->fDecomposed == ENoDecompose) this->Decompose(dt);
  const TPZSYsmpMatrixPardiso<TVar>* this_ct = const_cast<const TPZSYsmpMatrixPardiso<TVar>*>(this);
  this_ct->SolveDirect(F,dt);
  return 0;
}

template<class TVar>
int TPZSYsmpMatrixPardiso<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }
  if(!this->fDecomposed){
    this->Error(__PRETTY_FUNCTION__,"matrix should've been decomposed already");
  }
  TPZFMatrix<TVar> x(F);
  fPardisoControl.Solve(this,x,F);
  return 0;
}


template class TPZSYsmpMatrixPardiso<double>;
template class TPZSYsmpMatrixPardiso<float>;
template class TPZSYsmpMatrixPardiso<long double>;
template class TPZSYsmpMatrixPardiso<std::complex<float>>;
template class TPZSYsmpMatrixPardiso<std::complex<double>>;
template class TPZSYsmpMatrixPardiso<std::complex<long double>>;

#endif
