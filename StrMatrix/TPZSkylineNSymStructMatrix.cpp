
#include "TPZSkylineNSymStructMatrix.h"
#include "pzskylnsymmat.h"
#include "pzcmesh.h"

template<class TVar, class TPar>
TPZSkylineNSymStructMatrix<TVar,TPar>::TPZSkylineNSymStructMatrix(TPZCompMesh *cmesh)
  : TPZSkylineStructMatrix<TVar,TPar>(cmesh)
{
  ///nothing here
}

template<class TVar, class TPar>
TPZSkylineNSymStructMatrix<TVar,TPar>::TPZSkylineNSymStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh)
  : TPZSkylineStructMatrix<TVar,TPar>(cmesh)
{
  ///nothing here
}

template<class TVar, class TPar>
TPZMatrix<TVar> * TPZSkylineNSymStructMatrix<TVar,TPar>::ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline)
{
  return new TPZSkylNSymMatrix<TVar>(neq,skyline);
}

template<class TVar, class TPar>
TPZStructMatrix * TPZSkylineNSymStructMatrix<TVar,TPar>::Clone()
{
  return new TPZSkylineNSymStructMatrix<TVar,TPar>(*this);
}

template<class TVar, class TPar>
int TPZSkylineNSymStructMatrix<TVar,TPar>::ClassId() const{
  return Hash("TPZSkylineNSymStructMatrix") ^
    TPZSkylineStructMatrix<TVar,TPar>::ClassId() << 1;
}

#include "pzstrmatrixot.h"
#include "pzstrmatrixflowtbb.h"

template class TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixOR<STATE>>;
template class TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixOT<STATE>>;
template class TPZSkylineNSymStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>>;
template class TPZSkylineNSymStructMatrix<CSTATE,TPZStructMatrixOR<CSTATE>>;
template class TPZSkylineNSymStructMatrix<CSTATE,TPZStructMatrixOT<CSTATE>>;
template class TPZSkylineNSymStructMatrix<CSTATE,TPZStructMatrixTBBFlow<CSTATE>>;