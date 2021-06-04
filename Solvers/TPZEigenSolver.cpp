#include "TPZEigenSolver.h"
#include "TPZKrylovEigenSolver.h"
#include "TPZSpectralTransform.h"
#include "Hash/TPZHash.h"
#include <numeric>

template <class TVar>
void TPZEigenSolver<TVar>::SetAsGeneralised(bool isGeneralised) {
  fIsGeneralised = isGeneralised;
}

template <class TVar> int TPZEigenSolver<TVar>::ClassId() const {
  return Hash("TPZEigenSolver") ^ ClassIdOrHash<TVar>() << 1 ^
         TPZSolver::ClassId() << 2;
}

template<class TVar>
void TPZEigenSolver<TVar>::ResetMatrix()
{
  TPZAutoPointer<TPZMatrix<TVar>> newA, newB;
  fMatrixA = newA;
  fMatrixB = newB;
}

template<class TVar>
void TPZEigenSolver<TVar>::SortEigenvalues(TPZVec<CTVar> &w, TPZVec<int> &indices)
{
  const CTVar target = Target();
  
  const auto eigOrder = EigenSorting();
  auto sortFunc = [eigOrder,target](const CTVar a, const CTVar b) {
    switch (eigOrder) {
    case TPZEigenSort::AbsAscending:
      return fabs(a) < fabs(b);
    case TPZEigenSort::AbsDescending:
      return fabs(a) > fabs(b);
    case TPZEigenSort::RealAscending:
      return a.real() < b.real();
    case TPZEigenSort::RealDescending:
      return a.real() > b.real();
    case TPZEigenSort::ImagAscending:
      return a.imag() < b.imag();
    case TPZEigenSort::ImagDescending:
      return a.imag() > b.imag();
    case TPZEigenSort::TargetRealPart:
      return fabs(a.real() - target.real()) < fabs(b.real() - target.real());
    case TPZEigenSort::TargetImagPart:
      return fabs(a.imag() - target.imag()) < fabs(b.imag() - target.imag());
    case TPZEigenSort::TargetMagnitude:
      return fabs(fabs(a) - fabs(target)) < fabs(fabs(b) - fabs(target));
    }
    unreachable();
  };
  // sorting eigenvalues
  indices.Resize(w.size());
  std::iota(indices.begin(), indices.end(), 0); // Initializing
  std::stable_sort(
      indices.begin(), indices.end(),
      [&w, &sortFunc](int i, int j) { return sortFunc(w[i], w[j]); });
  std::stable_sort(w.begin(), w.end(),
                   [&sortFunc](auto i, auto j) { return sortFunc(i, j); });

  w.Resize(NEigenpairs());
}

template class TPZEigenSolver<float>;
template class TPZEigenSolver<double>;
template class TPZEigenSolver<long double>;

template class TPZEigenSolver<std::complex<float>>;
template class TPZEigenSolver<std::complex<double>>;
template class TPZEigenSolver<std::complex<long double>>;