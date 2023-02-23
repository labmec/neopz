/**
 * @file
 * @brief Contains the implementation of the TPZFYsmpMatrix methods.
 */

#ifdef USING_MKL
#include "TPZYSMPPardiso.h"
#include "pzfmatrix.h"
// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZFYsmpMatrixPardiso<TVar>::TPZFYsmpMatrixPardiso() : TPZRegisterClassId(&TPZFYsmpMatrixPardiso::ClassId),
TPZFYsmpMatrix<TVar>() {

}

template<class TVar>
void TPZFYsmpMatrixPardiso<TVar>::CopyFrom(const TPZMatrix<TVar> *  mat)
{                                                           
  auto *from = dynamic_cast<const TPZFYsmpMatrixPardiso<TVar> *>(mat);
  if (from) {
    *this = *from;
  }
  else
  {
    auto *from = dynamic_cast<const TPZFYsmpMatrix<TVar> *>(mat);
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
int TPZFYsmpMatrixPardiso<TVar>::ClassId() const{
    return Hash("TPZFYsmpMatrixPardiso") ^ TPZFYsmpMatrix<TVar>::ClassId() << 1;
}


template<class TVar>
void TPZFYsmpMatrixPardiso<TVar>::SetIsDecomposed(DecomposeType val){
  TPZBaseMatrix::SetIsDecomposed(val);
  if(val){fPardisoControl.fDecomposed = true;}
}

template<class TVar>
int TPZFYsmpMatrixPardiso<TVar>::Decompose(const DecomposeType dt)
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }

  if(!fPardisoControl.HasCustomSettings()){
    typename TPZPardisoSolver<TVar>::MStructure str = this->IsSymmetric() ? TPZPardisoSolver<TVar>::MStructure::ESymmetric : TPZPardisoSolver<TVar>::MStructure::ENonSymmetric;
    typename TPZPardisoSolver<TVar>::MSystemType sysType = this->IsSymmetric() ? TPZPardisoSolver<TVar>::MSystemType::ESymmetric : TPZPardisoSolver<TVar>::MSystemType::ENonSymmetric;
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
int TPZFYsmpMatrixPardiso<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt)
{
  if(this->fDecomposed && this->fDecomposed != dt){
    this->Error(__PRETTY_FUNCTION__,"matrix is already decomposed with other scheme");
  }
  if(this->fDecomposed == ENoDecompose) this->Decompose(dt);
  const TPZFYsmpMatrixPardiso<TVar>* this_ct = const_cast<const TPZFYsmpMatrixPardiso<TVar>*>(this);
  this_ct->SolveDirect(F,dt);
  return 0;
}

template<class TVar>
int TPZFYsmpMatrixPardiso<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const
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


template class TPZFYsmpMatrixPardiso<double>;
template class TPZFYsmpMatrixPardiso<float>;
template class TPZFYsmpMatrixPardiso<long double>;
template class TPZFYsmpMatrixPardiso<std::complex<float>>;
template class TPZFYsmpMatrixPardiso<std::complex<double>>;
template class TPZFYsmpMatrixPardiso<std::complex<long double>>;

#endif
