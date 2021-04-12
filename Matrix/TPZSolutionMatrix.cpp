#include "TPZSolutionMatrix.h"
#include "TPZStream.h"

TPZSolutionMatrix::TPZSolutionMatrix()
    : fSolType(EUndefined) {
    fBaseMatrix = nullptr;
}

TPZSolutionMatrix::TPZSolutionMatrix(bool is_complex) {
  fSolType = is_complex ? EComplex : EReal;
  if(fSolType == EReal)
    {
      fBaseMatrix = &fRealMatrix;
    }
  else{DebugStop();}
  // else{ fBaseMatrix = &fComplexMatrix;}
}
TPZSolutionMatrix::TPZSolutionMatrix(int nrows, int ncols,
                                                  bool is_complex)
    : fSolType(is_complex ? EComplex : EReal) {
    if(fSolType == EReal) {
    fBaseMatrix = &fRealMatrix;
  }
  else{DebugStop();}
  // else{ fBaseMatrix = &fComplexMatrix;}
  fBaseMatrix->Resize(nrows, ncols);
}

TPZSolutionMatrix::TPZSolutionMatrix(const TPZFMatrix<STATE> &sol)
    : fSolType(EReal), fRealMatrix(sol), fBaseMatrix(&fRealMatrix)
{

}

// TPZSolutionMatrix::TPZSolutionMatrix(const TPZFMatrix<CSTATE> &sol)
//     : fSolType(EComplex), fComplexMatrix(sol), fBaseMatrix(&fComplexMatrix)
// {

// }

TPZSolutionMatrix::TPZSolutionMatrix(const TPZSolutionMatrix &cp) :
    fSolType(cp.fSolType), fRealMatrix(cp.fRealMatrix)
    // , fComplexMatrix(cp.fComplexMatrix)
{
    if(fSolType == EReal) {
    fBaseMatrix = &fRealMatrix;
  }
    else{DebugStop();}
  // else{ fBaseMatrix = &fComplexMatrix;}
}

TPZSolutionMatrix&
TPZSolutionMatrix::operator=(const TPZSolutionMatrix &cp)
{
    if(fSolType == EUndefined){
        fSolType = cp.fSolType;
    }
    if(fSolType == cp.fSolType){
        fRealMatrix = cp.fRealMatrix;
        // fComplexMatrix = copy.fComplexMatrix;
        if(fSolType == EReal) {
          fBaseMatrix = &fRealMatrix;
        }
        else{DebugStop();}
        // else{ fBaseMatrix = &fComplexMatrix;}
        return *this;
  } else{DebugStop();}
    return *this;
}

template<class TVar>
TPZSolutionMatrix&
TPZSolutionMatrix::operator=(const TPZFMatrix<TVar> &mat)
{
    if constexpr (std::is_same<TVar,STATE>::value
                  // ||std::is_same<TVar,CSTATE>::value
                  ){
        if constexpr(std::is_same<TVar,STATE>::value){
            if(fSolType == EReal || fSolType == EUndefined){
                fSolType = EReal;
                fRealMatrix = mat;
                fBaseMatrix = &fRealMatrix;
                return *this;
            }
            // else if(fSolType == EComplex || fSolType == EUndefined){
            //     fSolType = EComplex;
            //     fRealMatrix = mat;
            //     fBaseMatrix = &fComplexMatrix;
            //     return *this;
            // }
        }
    }
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" called with incompatible type\n";
    DebugStop();
    return *this;
}

TPZSolutionMatrix::operator TPZBaseMatrix &(){
    return *fBaseMatrix;
}

TPZSolutionMatrix::operator TPZFMatrix<STATE> &(){
    if(fSolType != EReal){
        PZError<<__PRETTY_FUNCTION__;
        PZError << " called with incompatible type\n";
        DebugStop();
    }
    return fRealMatrix;
}

TPZSolutionMatrix::operator const TPZFMatrix<STATE> &()const{
    if(fSolType != EReal){
        PZError<<__PRETTY_FUNCTION__;
        PZError << " called with incompatible type\n";
        DebugStop();
    }
    return fRealMatrix;
}

// TPZSolutionMatrix::operator TPZFMatrix<CSTATE> &(){
//     if(fSolType != EComplex){
//         PZError<<__PRETTY_FUNCTION__;
//         PZError << " called with incompatible type\n";
//         DebugStop();
//     }
//     return fComplexMatrix;
// }

void TPZSolutionMatrix::Read(TPZStream &buf, void *context) {
  bool isComplex;
  buf.Read(isComplex);
  if(isComplex) fSolType = EComplex;
  else fSolType = EReal;
  
  if (fSolType == EReal)
    return fRealMatrix.Read(buf, context);
  else{DebugStop();}
  // else return fComplexMatrix.Read(buf,context);
}
void TPZSolutionMatrix::Write(TPZStream &buf, int withclassid) const {
  if(fSolType == EUndefined){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" cannot write TPZSolutionMatrix without defining its type!\n";
    DebugStop();
  }
  const bool isComplex = fSolType == EComplex ? true : false;
  buf.Write(isComplex);
  if (!isComplex)
    return fRealMatrix.Write(buf, withclassid);
  else{DebugStop();}
  // else return fComplexMatrix.Write(buf,withclassid);
}

template
TPZSolutionMatrix &TPZSolutionMatrix::operator=<STATE>(const TPZFMatrix<STATE> &mat);