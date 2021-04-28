#include "TPZSolutionMatrix.h"
#include "TPZStream.h"

TPZSolutionMatrix::TPZSolutionMatrix() :
  fSolType(EUndefined), fBaseMatrix(nullptr){
}

TPZSolutionMatrix::TPZSolutionMatrix(bool is_complex) {
  fSolType = is_complex ? EComplex : EReal;
  if(fSolType == EReal)
    {
      fBaseMatrix = &fRealMatrix;
    }
  else{ fBaseMatrix = &fComplexMatrix;}
}
TPZSolutionMatrix::TPZSolutionMatrix(int nrows, int ncols,
                                     bool is_complex)
  : fSolType(is_complex ? EComplex : EReal) {
  if(fSolType == EReal)
    {
      fBaseMatrix = &fRealMatrix;
    }
    else{ fBaseMatrix = &fComplexMatrix;}
  fBaseMatrix->Resize(nrows, ncols);
}

TPZSolutionMatrix::TPZSolutionMatrix(const TPZFMatrix<STATE> &sol)
    : fSolType(EReal), fRealMatrix(sol),
      fBaseMatrix(&fRealMatrix)
{

}

TPZSolutionMatrix::TPZSolutionMatrix(const TPZFMatrix<CSTATE> &sol)
    : fSolType(EComplex), fComplexMatrix(sol),
      fBaseMatrix(&fComplexMatrix)
{

}

TPZSolutionMatrix::TPZSolutionMatrix(const TPZSolutionMatrix &cp) :
    fSolType(cp.fSolType), fRealMatrix(cp.fRealMatrix)
    , fComplexMatrix(cp.fComplexMatrix)
{
  if(fSolType == EReal)
    {
      fBaseMatrix = &fRealMatrix;
    }
  else{ fBaseMatrix = &fComplexMatrix;}
}

TPZSolutionMatrix&
TPZSolutionMatrix::operator=(const TPZSolutionMatrix &cp)
{
    if(fSolType == EUndefined){
        fSolType = cp.fSolType;
    }
    if(fSolType == cp.fSolType){
        fRealMatrix = cp.fRealMatrix;
        fComplexMatrix = cp.fComplexMatrix;
        if(fSolType == EReal) {
          fBaseMatrix = &fRealMatrix;
        }
        else{ fBaseMatrix = &fComplexMatrix;}
        return *this;
  } else{DebugStop();}
    return *this;
}

template<class TVar>
TPZSolutionMatrix&
TPZSolutionMatrix::operator=(const TPZFMatrix<TVar> &mat)
{
  if constexpr(std::is_same<TVar,STATE>::value){
    if(fSolType == EReal || fSolType == EUndefined){
      fSolType = EReal;
      fRealMatrix = mat;
      fBaseMatrix = &fRealMatrix;
      return *this;
    }
  }
  else if constexpr(std::is_same<TVar,CSTATE>::value){
    if(fSolType == EComplex || fSolType == EUndefined){
      fSolType = EComplex;
      fComplexMatrix = mat;
      fBaseMatrix = &fComplexMatrix;
      return *this;
    }
  }
  PZError<<__PRETTY_FUNCTION__;
  PZError<<" called with incompatible type\n";
  DebugStop();
  return *this;
}

template <class TVar>
TPZSolutionMatrix &TPZSolutionMatrix::operator+=(const TPZFMatrix<TVar> &mat){
  if (fSolType == EUndefined){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" called in undefined solution matrix.\n";
    PZError<<"Aborting...\n";
    DebugStop();
  }
  if constexpr(std::is_same<TVar,STATE>::value){
    if(fSolType == EReal){
      fRealMatrix += mat;
      return *this;
    }
  }
  else if constexpr(std::is_same<TVar,CSTATE>::value){
    if(fSolType == EComplex){
      fComplexMatrix += mat;              
      return *this;
    }
  }
    
  PZError<<__PRETTY_FUNCTION__;
  PZError<<" called with incompatible type\n";
  DebugStop();
  return *this;
}

TPZSolutionMatrix &TPZSolutionMatrix::operator+=(const TPZSolutionMatrix &sol){
  if(fSolType == EReal && sol.fSolType == EReal){
    fRealMatrix += sol.fRealMatrix;
    return *this;
  }
  else if(fSolType == EComplex && sol.fSolType == EComplex){
    fComplexMatrix += sol.fComplexMatrix;
    return *this;
  }
  PZError << __PRETTY_FUNCTION__;
  PZError << " called with incompatible type\n";
  DebugStop();
  return *this;
}

TPZSolutionMatrix::operator TPZBaseMatrix &(){
  if(!fBaseMatrix){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" trying to access undefined matrix\n";
    PZError<<"Aborting...";
    DebugStop();
  }
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

TPZSolutionMatrix::operator TPZFMatrix<CSTATE> &(){
    if(fSolType != EComplex){
        PZError<<__PRETTY_FUNCTION__;
        PZError << " called with incompatible type\n";
        DebugStop();
    }
    return fComplexMatrix;
}

TPZSolutionMatrix::operator const TPZFMatrix<CSTATE> &() const{
    if(fSolType != EComplex){
        PZError<<__PRETTY_FUNCTION__;
        PZError << " called with incompatible type\n";
        DebugStop();
    }
    return fComplexMatrix;
}


void TPZSolutionMatrix::ExpandAndSetSol(const TPZSolutionMatrix & sol, const int64_t nrows){
  const auto origRows = this->Rows();
  const auto solRows = sol.Rows();
  const auto solCols = sol.Cols();
#ifdef PZ_DEBUG
  if(origRows < solRows){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<" called with original solution smaller than it should be.\n";
    PZError<<"   original #rows: "<<origRows<<'\n';
    PZError<<"   new      #rows: "<<solRows<<'\n';
    PZerror<<"Aborting...\n";
    DebugStop();
  }
#endif
  this->Resize(origRows,solCols);
  if(fSolType == EReal){
    const TPZFMatrix<STATE> &solst = sol;
    for(auto j=0;j<solCols;j++)
    {
        for(auto i=0;i<solRows;i++)
        {
          //let us skip boundary checking
          fRealMatrix.PutVal(i,j,solst.GetVal(i, j));
        }
    }
  }
  else if(fSolType == EComplex){
    const TPZFMatrix<CSTATE> &solst = sol;
    for(auto j=0;j<solCols;j++)
    {
        for(auto i=0;i<solRows;i++)
        {
          fComplexMatrix.PutVal(i,j,solst.GetVal(i, j));
        }
    }
  }
  else{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<": solution does not have a defined type! Aborting...\n";
    DebugStop();
  }
}

void TPZSolutionMatrix::Read(TPZStream &buf, void *context) {
  bool isComplex;
  buf.Read(isComplex);
  if(isComplex) fSolType = EComplex;
  else fSolType = EReal;
  
  if (fSolType == EReal)
    return fRealMatrix.Read(buf, context);
  else return fComplexMatrix.Read(buf,context);
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
  else return fComplexMatrix.Write(buf,withclassid);
}




template
TPZSolutionMatrix &TPZSolutionMatrix::operator=<STATE>(const TPZFMatrix<STATE> &mat);
template
TPZSolutionMatrix &TPZSolutionMatrix::operator+=<STATE>(const TPZFMatrix<STATE> &mat);

template
TPZSolutionMatrix &TPZSolutionMatrix::operator=<CSTATE>(const TPZFMatrix<CSTATE> &mat);
template
TPZSolutionMatrix &TPZSolutionMatrix::operator+=<CSTATE>(const TPZFMatrix<CSTATE> &mat);