
#include "TPZSkylineNSymStructMatrix.h"
#include "pzskylnsymmat.h"
#include "pzcmesh.h"

TPZSkylineNSymStructMatrix::TPZSkylineNSymStructMatrix(TPZCompMesh *cmesh)
                           : TPZSkylineStructMatrix(cmesh)
{
  ///nothing here
}

TPZSkylineNSymStructMatrix::TPZSkylineNSymStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh)
                            : TPZSkylineStructMatrix(cmesh)
{
  ///nothing here
}


TPZSkylineNSymStructMatrix::TPZSkylineNSymStructMatrix(const TPZSkylineStructMatrix &cp):TPZSkylineStructMatrix(cp)
{
  ///nothing here
}

TPZSkylineNSymStructMatrix::~TPZSkylineNSymStructMatrix()
{
  ///nothing here
}

TPZMatrix<STATE> * TPZSkylineNSymStructMatrix::ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline)
{
  return new TPZSkylNSymMatrix<STATE>(neq,skyline);
}

