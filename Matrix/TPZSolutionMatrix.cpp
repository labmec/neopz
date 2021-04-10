#include "TPZSolutionMatrix.h"
#include "TPZStream.h"

TPZSolutionMatrix::TPZSolutionMatrix(int nrows, int ncols,
                                                  bool is_complex)
    : fIsComplex(is_complex) {
  if (!fIsComplex) {
    fBaseMatrix = &fRealMatrix;
  }
  else{DebugStop();}
  // else{ fBaseMatrix = &fComplexMatrix;}
  fBaseMatrix->Resize(nrows, ncols);
}

TPZSolutionMatrix::TPZSolutionMatrix(TPZFMatrix<STATE> &sol)
    : fIsComplex(false), fRealMatrix(sol), fBaseMatrix(&fRealMatrix)
{

}

// TPZSolutionMatrix::TPZSolutionMatrix(TPZFMatrix<CSTATE> &sol)
//     : fIsComplex(true), fComplexMatrix(sol), fBaseMatrix(&fComplexMatrix)
// {

// }

TPZSolutionMatrix::TPZSolutionMatrix(const TPZSolutionMatrix &cp) :
    fIsComplex(cp.fIsComplex), fRealMatrix(cp.fRealMatrix)
    // , fComplexMatrix(cp.fComplexMatrix)
{
    if (!fIsComplex) {
    fBaseMatrix = &fRealMatrix;
  }
    else{DebugStop();}
  // else{ fBaseMatrix = &fComplexMatrix;}
}

TPZSolutionMatrix&
TPZSolutionMatrix::operator=(const TPZSolutionMatrix &cp)
{
    fIsComplex = cp.fIsComplex;
    fRealMatrix = cp.fRealMatrix;
    // fComplexMatrix = copy.fComplexMatrix;
    if (!fIsComplex) {
    fBaseMatrix = &fRealMatrix;
  }
    else{DebugStop();}
  // else{ fBaseMatrix = &fComplexMatrix;}
    return *this;
}

void TPZSolutionMatrix::Read(TPZStream &buf, void *context) {
  buf.Read(fIsComplex);
  if (!fIsComplex)
    return fRealMatrix.Read(buf, context);
  else{DebugStop();}
  // else return fComplexMatrix.Read(buf,context);
}
void TPZSolutionMatrix::Write(TPZStream &buf, int withclassid) const {
  buf.Write(fIsComplex);
  if (!fIsComplex)
    return fRealMatrix.Write(buf, withclassid);
  else{DebugStop();}
  // else return fComplexMatrix.Write(buf,withclassid);
}