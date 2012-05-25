
#include "TPZSkylineNSymStructMatrix.h"
#include "pzskylnsymmat.h"

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

TPZMatrix<STATE> * TPZSkylineNSymStructMatrix::ReallyCreate(int neq, const TPZVec<int> &skyline)
{
  return new TPZSkylNSymMatrix<STATE>(neq,skyline);
}

