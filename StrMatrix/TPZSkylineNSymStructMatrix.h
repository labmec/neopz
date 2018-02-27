//---------------------------------------------------------------------------

#ifndef TPZSkylineNSymStructMatrixH
#define TPZSkylineNSymStructMatrixH

#include "pzskylstrmatrix.h"

/**
 * Implements Non symmetrical SkyLine Structural Matrices
 * ingroup structural
 */
class TPZSkylineNSymStructMatrix : public TPZSkylineStructMatrix {

protected:

  /** Returns the skyline matrix object */
  virtual TPZMatrix<STATE> * ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline);

public:

  TPZSkylineNSymStructMatrix(TPZCompMesh *cmesh);

  TPZSkylineNSymStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);

  TPZSkylineNSymStructMatrix(const TPZSkylineStructMatrix &cp);

  ~TPZSkylineNSymStructMatrix();

  virtual TPZStructMatrix * Clone()
  {
    return new TPZSkylineNSymStructMatrix(*this);
  }

    

};

#endif
