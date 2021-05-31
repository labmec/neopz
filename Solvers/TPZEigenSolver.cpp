#include "TPZEigenSolver.h"
#include "Hash/TPZHash.h"
#include <numeric>

template<class TVar>
void TPZEigenSolver<TVar>::SetAsGeneralised(bool isGeneralised)
{
    fIsGeneralised = isGeneralised;
}

template<class TVar>
int TPZEigenSolver<TVar>::ClassId() const
{
    return Hash("TPZEigenSolver") ^
        ClassIdOrHash<TVar>() << 1 ^
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
    const auto eigOrder = EigenSorting();
    auto sortFunc = [eigOrder](const CTVar a, const CTVar b) {
      switch (eigOrder) {
      case TPZEigenSort::EAbsAscending:
        return fabs(a) < fabs(b);
      case TPZEigenSort::EAbsDescending:
        return fabs(a) > fabs(b);
      case TPZEigenSort::ERealAscending:
        return a.real() < b.real();
      case TPZEigenSort::ERealDescending:
        return a.real() > b.real();
      case TPZEigenSort::EImagAscending:
        return a.imag() < b.imag();
      case TPZEigenSort::EImagDescending:
        return a.imag() > b.imag();
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