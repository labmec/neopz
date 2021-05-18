#include "TPZEigenSolver.h"
#include "Hash/TPZHash.h"

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

template class TPZEigenSolver<float>;
template class TPZEigenSolver<double>;
template class TPZEigenSolver<long double>;

template class TPZEigenSolver<std::complex<float>>;
template class TPZEigenSolver<std::complex<double>>;
template class TPZEigenSolver<std::complex<long double>>;