#include "TPZSpectralTransform.h"
#include "Hash/TPZHash.h"
#include "pzvec.h"
#include "TPZStream.h"
#include "pzmatrix.h"
#include "TPZFMatrixRef.h"
#include "TPZYSMPMatrix.h"
#include "TPZSYSMPMatrix.h"
#include <complex>


template<class TVar>
int TPZSpectralTransform<TVar>::ClassId() const{
  return Hash("TPZSpectralTransform") ^
    ClassIdOrHash<TVar>() << 1;
}

template<class TVar>
void
TPZSTShiftOrigin<TVar>::CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A, TPZAutoPointer<TPZMatrix<TVar>>B) const
{
  const auto sp= B->GetSymmetry();
  const bool use_lu = sp == SymProp::Herm || (!std::is_same_v<TVar,RTVar> && sp == SymProp::Sym);
  if (use_lu) B->Decompose(ELU);
  else B->Decompose(ELDLt);
  //b-1 * shiftedA will be computed at the arnoldi iteration
  const auto &shift = Shift();
  const auto nRows = A->Rows();
  for(int i = 0; i < nRows; i++) A->PutVal(i,i,A->GetVal(i,i)-shift);
}

template<class TVar>
void
TPZSTShiftOrigin<TVar>::CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A) const
{
  const auto nRows = A->Rows();
  //calculating A-sigmaB
  const auto &shift = Shift();
  for(int i = 0; i < nRows; i++) A->PutVal(i,i,A->GetVal(i,i)-shift);
}

template<class TVar>
void TPZSTShiftOrigin<TVar>::TransformEigenvalues(TPZVec<CTVar> &w) const
{
  for(auto &mappedw : w) mappedw += fShift;
}

template<class TVar>
int TPZSTShiftOrigin<TVar>::ClassId() const
{
  return Hash("TPZShiftOrigin") ^
    TPZSpectralTransform<TVar>::ClassId() << 1;
}

template<class TVar>
void TPZSTShiftOrigin<TVar>::Write(TPZStream &buf, int withclassid) const
{
  TPZSpectralTransform<TVar>::Write(buf,withclassid);
  buf.Write(&fShift);
}

template<class TVar>
void
TPZSTShiftOrigin<TVar>::Read(TPZStream &buf, void *context)
{
  TPZSpectralTransform<TVar>::Read(buf,context);
  buf.Read(&fShift);
}

template<class TVar>
void
TPZSTShiftAndInvert<TVar>::CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A, TPZAutoPointer<TPZMatrix<TVar>>B) const
{
  const auto &shift = this->Shift();
  if(A->Storage().Rows() != B->Storage().Rows()){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Matrices have uncompatible storage formats.\nAborting...\n";
    DebugStop();
  }
  A->Storage() -= B->Storage() * shift;

  const auto sp= A->GetSymmetry();
  const bool use_lu = sp == SymProp::Herm || (!std::is_same_v<TVar,RTVar> && sp == SymProp::Sym);
  if (use_lu) A->Decompose(ELU);
  else A->Decompose(ELDLt);
}

template<class TVar>
void
TPZSTShiftAndInvert<TVar>::CalcMatrix(TPZAutoPointer<TPZMatrix<TVar>>A) const
{
  const auto &shift = this->Shift();
  const auto nRows = A->Rows();
  for(int i = 0; i < nRows; i++) A->PutVal(i,i,A->GetVal(i,i)-shift);
  const auto sp= A->GetSymmetry();
  const bool use_lu = sp == SymProp::Herm || (!std::is_same_v<TVar,RTVar> && sp == SymProp::Sym);
  if (use_lu) A->Decompose(ELU);
  else A->Decompose(ELDLt);
}

template<class TVar>
void
TPZSTShiftAndInvert<TVar>::TransformEigenvalues(TPZVec<CTVar> &w) const
{
  const auto &s = this->fShift;
  for(auto &mappedw : w) mappedw = ((CTVar)1.0)/mappedw + s;
}

template<class TVar>
int
TPZSTShiftAndInvert<TVar>::ClassId() const
{
  return Hash("TPZSTShiftAndInvert") ^
    TPZSTShiftOrigin<TVar>::ClassId() << 1;
}


#define INSTANTIATE_TEMPLATES(TCLASS)           \
  template class TCLASS<float>;                 \
  template class TCLASS<double>;                \
  template class TCLASS<std::complex<float>>;   \
  template class TCLASS<std::complex<double>>;  \


INSTANTIATE_TEMPLATES(TPZSTShiftOrigin)
INSTANTIATE_TEMPLATES(TPZSTShiftAndInvert)
#undef INSTANTIATE_TEMPLATES
