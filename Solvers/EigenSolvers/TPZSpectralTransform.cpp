#include "TPZSpectralTransform.h"
#include "Hash/TPZHash.h"
#include "pzvec.h"
#include "TPZStream.h"
#include "pzmatrix.h"

#include <complex>


template<class TVar>
int TPZSpectralTransform<TVar>::ClassId() const{
  return Hash("TPZSpectralTransform") ^
    ClassIdOrHash<TVar>() << 1;
}

template<class TVar>
TPZMatrix<TVar> *
TPZSTShiftOrigin<TVar>::CalcMatrix(TPZMatrix<TVar> &A, TPZMatrix<TVar> &B) const
{
  if (B.IsSymmetric()) B.Decompose_LDLt();
  else B.Decompose_LU();
  TPZMatrix<TVar> * shiftedMat = A.NewMatrix();
  //b-1 * shiftedA will be computed at the arnoldi iteration
  const auto &shift = Shift();
  const auto nRows = A.Rows();
  for(int i = 0; i < nRows; i++) shiftedMat->PutVal(i,i,shiftedMat->GetVal(i,i)-shift);
  return shiftedMat;
}

template<class TVar>
TPZMatrix<TVar> *
TPZSTShiftOrigin<TVar>::CalcMatrix(TPZMatrix<TVar> &A) const
{
  const auto nRows = A.Rows();
  TPZMatrix<TVar> * shiftedMat = &A;
  //calculating A-sigmaB
  const auto &shift = Shift();
  for(int i = 0; i < nRows; i++) shiftedMat->PutVal(i,i,shiftedMat->GetVal(i,i)-shift);
  return shiftedMat;
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
TPZMatrix<TVar> *
TPZSTShiftAndInvert<TVar>::CalcMatrix(TPZMatrix<TVar> &A, TPZMatrix<TVar> &B) const
{
  TPZMatrix<TVar> * shiftedMat = A.NewMatrix();
  const auto &shift = this->Shift();
  B*=shift;
  A.Subtract(B,*shiftedMat);
  if (shiftedMat->IsSymmetric()) shiftedMat->Decompose_LDLt();
  else shiftedMat->Decompose_LU();
  return shiftedMat;
}

template<class TVar>
TPZMatrix<TVar> *
TPZSTShiftAndInvert<TVar>::CalcMatrix(TPZMatrix<TVar> &A) const
{
  TPZMatrix<TVar> * shiftedMat = A.NewMatrix();
  const auto &shift = this->Shift();
  const auto nRows = A.Rows();
  for(int i = 0; i < nRows; i++) shiftedMat->PutVal(i,i,A.GetVal(i,i)-shift);
  if (shiftedMat->IsSymmetric()) shiftedMat->Decompose_LDLt();
  else shiftedMat->Decompose_LU();
  return shiftedMat;
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