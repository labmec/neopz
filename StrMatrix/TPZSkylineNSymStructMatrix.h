//---------------------------------------------------------------------------

#ifndef TPZSkylineNSymStructMatrixH
#define TPZSkylineNSymStructMatrixH

#include "pzskylstrmatrix.h"

/**
 * Implements a skyline structural matrix using TPZSkylNSymMatrix as a storage format.
 * ingroup structural
 */
template<class TVar=STATE, class TPar=TPZStructMatrixOR<TVar>>
class TPZSkylineNSymStructMatrix : public TPZSkylineStructMatrix<TVar,TPar> {

protected:

  /** Returns the skyline matrix object */
  TPZMatrix<TVar> * ReallyCreate(int64_t neq, const TPZVec<int64_t> &skyline) override;

public:

  TPZSkylineNSymStructMatrix(TPZCompMesh *cmesh);

  TPZSkylineNSymStructMatrix(TPZAutoPointer<TPZCompMesh> cmesh);

  TPZStructMatrix * Clone() override;

  int ClassId() const override;
};

#endif
