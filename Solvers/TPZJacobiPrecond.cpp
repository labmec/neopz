#include "TPZJacobiPrecond.h"
#include "TPZFMatrixRef.h"

template<class TVar>
TPZJacobiPrecond<TVar>::TPZJacobiPrecond(TPZAutoPointer<TPZMatrix<TVar> >  refmat)
{
  if(!refmat){
    PZError<<__PRETTY_FUNCTION__
           <<"\nInvalid matrix as input! nullptr"<<std::endl;
    DebugStop();
  }
  this->fReferenceMatrix = refmat;
  this->SetMatrix(new TPZFMatrix<TVar>());
  UpdateFrom(refmat);
}


template<class TVar>
void TPZJacobiPrecond<TVar>::UpdateFrom(TPZAutoPointer<TPZBaseMatrix> mat)
{
  if(mat && mat == this->fReferenceMatrix){
    const auto sz = mat->Rows();
    if(sz!= this->Matrix()->Rows()){this->Matrix()->Resize(sz,1);}
    auto ms = this->Matrix()->Storage().Elem();
    auto matval = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(mat);
    for(int i = 0; i < sz; i++){
      ms[i] = 1./(matval->GetVal(i,i));
    }
  }
}

template<class TVar>
void TPZJacobiPrecond<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual)
{
  const auto nr = F.Rows();

#ifdef PZDEBUG
  if(nr != this->Matrix()->Rows()){
    PZError<<__PRETTY_FUNCTION__
           <<"\nIncompatible sizes!\n"
           <<"\tmatrix diagonal size: "<<this->Matrix()->Rows()<<'\n'
           <<"\trhs size: "<<nr<<'\n'
           <<"Aborting..."<<std::endl;
    DebugStop();
  }
#endif
  const auto nc = F.Cols();
  if(result.Rows() != nr || result.Cols() != nc){
    result.Resize(nr,nc);
  }
  auto ms = this->Matrix()->Storage().Elem();
  for(int ic = 0; ic < nc; ic++){
    for(int ir = 0; ir < nr; ir++){
      result.PutVal(ir,ic, F.GetVal(ir,ic)*ms[ir]);
    }
  }
}

template class TPZJacobiPrecond<STATE>;
template class TPZJacobiPrecond<CSTATE>;