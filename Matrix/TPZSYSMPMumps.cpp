/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrixMumps methods.
 */

#ifdef USING_MUMPS
#include "TPZSYSMPMumps.h"
#include "pzfmatrix.h"

#include <dmumps_c.h>

// Note: MUMPS does not have optimized sparse matrix-vector multiplication
// like MKL's sparse BLAS. We'll use the base class implementation for MultAdd.

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************

template <class TVar>
TPZSYsmpMatrixMumps<TVar>::TPZSYsmpMatrixMumps() : TPZRegisterClassId(&TPZSYsmpMatrixMumps::ClassId),
                                                   TPZSYsmpMatrix<TVar>()
{
}

template <class TVar>
TPZSYsmpMatrixMumps<TVar>::TPZSYsmpMatrixMumps(const TPZSYsmpMatrixMumps<TVar> &cp) : TPZRegisterClassId(&TPZSYsmpMatrixMumps::ClassId),
                                                                                      TPZSYsmpMatrix<TVar>(cp),
                                                                                      fCOOValid(false)
{
    // Note: fMumpsControl is NOT copied - each matrix gets a fresh MUMPS instance
    // The copied matrix will need to be decomposed again if needed
    // COO format will be recomputed when needed
}

template <class TVar>
TPZSYsmpMatrixMumps<TVar>::TPZSYsmpMatrixMumps(TPZSYsmpMatrixMumps<TVar> &&cp) : TPZRegisterClassId(&TPZSYsmpMatrixMumps::ClassId),
                                                                                 TPZSYsmpMatrix<TVar>(std::move(cp)),
                                                                                 fMumpsControl(std::move(cp.fMumpsControl)),
                                                                                 fIRN1Based(std::move(cp.fIRN1Based)),
                                                                                 fJCN1Based(std::move(cp.fJCN1Based)),
                                                                                 fCOOValid(cp.fCOOValid)
{
    cp.fCOOValid = false;
}

template <class TVar>
TPZSYsmpMatrixMumps<TVar> &TPZSYsmpMatrixMumps<TVar>::operator=(const TPZSYsmpMatrixMumps<TVar> &copy)
{
    if (this != &copy)
    {
        TPZSYsmpMatrix<TVar>::operator=(copy);
        // Note: fMumpsControl is NOT copied - each matrix gets a fresh MUMPS instance
        fMumpsControl = TPZMumpsSolver<TVar>();
        // Reset decomposition state since we have a new MUMPS instance
        this->SetIsDecomposed(ENoDecompose);
        fCOOValid = false;
    }
    return *this;
}

template <class TVar>
TPZSYsmpMatrixMumps<TVar> &TPZSYsmpMatrixMumps<TVar>::operator=(TPZSYsmpMatrixMumps<TVar> &&copy)
{
    if (this != &copy)
    {
        TPZSYsmpMatrix<TVar>::operator=(std::move(copy));
        fMumpsControl = std::move(copy.fMumpsControl);
        fIRN1Based = std::move(copy.fIRN1Based);
        fJCN1Based = std::move(copy.fJCN1Based);
        fCOOValid = copy.fCOOValid;
        copy.fCOOValid = false;
    }
    return *this;
}

template <class TVar>
void TPZSYsmpMatrixMumps<TVar>::CopyFrom(const TPZMatrix<TVar> *mat)
{
    auto *from = dynamic_cast<const TPZSYsmpMatrixMumps<TVar> *>(mat);
    if (from)
    {
        *this = *from;
        return;
    }

    auto *fromBase = dynamic_cast<const TPZSYsmpMatrix<TVar> *>(mat);
    if (fromBase && fromBase->IsDecomposed() == ENoDecompose)
    {
        *this = *fromBase;
        return;
    }

    PZError << __PRETTY_FUNCTION__;
    PZError << "\nERROR: Called with incompatible type\n";
    PZError << "Aborting...\n";
    DebugStop();
}

template <class TVar>
int TPZSYsmpMatrixMumps<TVar>::ClassId() const
{
    return Hash("TPZSYsmpMatrixMumps") ^ TPZSYsmpMatrix<TVar>::ClassId() << 1;
}

template <class TVar>
void TPZSYsmpMatrixMumps<TVar>::MultAdd(const TPZFMatrix<TVar> &x,
                                        const TPZFMatrix<TVar> &y,
                                        TPZFMatrix<TVar> &z,
                                        const TVar alpha, const TVar beta,
                                        const int opt) const
{
    // MUMPS does not provide sparse matrix-vector multiplication utilities
    // like MKL's sparse BLAS. We have two options:
    // 1. Use the base class implementation (TPZSYsmpMatrix::MultAdd)
    // 2. Implement our own CSR matrix-vector product

    // For simplicity and correctness, we'll use the base class implementation
    // which already handles symmetric matrices correctly
    TPZSYsmpMatrix<TVar>::MultAdd(x, y, z, alpha, beta, opt);

    // Note: If performance is critical, you could implement an optimized
    // CSR symmetric matrix-vector multiplication here using OpenMP
}

template<class TVar>
void TPZSYsmpMatrixMumps<TVar>::SetData(const TPZVec<int64_t> &IA, const TPZVec<int64_t> &JA, const TPZVec<TVar> &A) {
    TPZSYsmpMatrix<TVar>::SetData(IA, JA, A);
    fCOOValid = false;
}

template<class TVar>
void TPZSYsmpMatrixMumps<TVar>::SetData(TPZVec<int64_t> &&IA, TPZVec<int64_t> &&JA, TPZVec<TVar> &&A) {
    TPZSYsmpMatrix<TVar>::SetData(std::move(IA), std::move(JA), std::move(A));
    fCOOValid = false;
}

template <class TVar>
void TPZSYsmpMatrixMumps<TVar>::SetIsDecomposed(DecomposeType val)
{
    TPZBaseMatrix::SetIsDecomposed(val);
    if (val == ENoDecompose)
    {
        // Matrix may have been modified, invalidate COO format
        fCOOValid = false;
    }
    if (val)
    {
        fMumpsControl.fDecomposed = true;
    }
}

template <class TVar>
int TPZSYsmpMatrixMumps<TVar>::Decompose(const DecomposeType dt)
{
    // Check if already decomposed with a different scheme
    if (this->fDecomposed && this->fDecomposed != dt)
    {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix is already decomposed with other scheme");
    }

    // Set up MUMPS parameters based on matrix properties
    if (!fMumpsControl.HasCustomSettings())
    {
        const auto sysType = this->GetSymmetry();
        typename TPZMumpsSolver<TVar>::MProperty prop =
            this->IsDefPositive() ? TPZMumpsSolver<TVar>::MProperty::EPositiveDefinite : TPZMumpsSolver<TVar>::MProperty::EIndefinite;
        fMumpsControl.SetMatrixType(sysType, prop);
    }

    if (!fCOOValid) {
       UpdateCOOFormat();
    }

    // Perform decomposition
    fMumpsControl.Decompose(this);
    fIRN1Based.Resize(0);
    fJCN1Based.Resize(0);
    fCOOValid = false;

    this->SetIsDecomposed(dt);

    return 0;
}

template <class TVar>
int TPZSYsmpMatrixMumps<TVar>::SolveDirect(TPZFMatrix<TVar> &F,
                                           const DecomposeType dt)
{
    // Non-const version

    // Check if already decomposed with a different scheme
    if (this->fDecomposed && this->fDecomposed != dt)
    {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix is already decomposed with other scheme");
    }

    // Decompose if not already done
    if (this->fDecomposed == ENoDecompose)
    {
        this->Decompose(dt);
    }

    // Call const version
    const TPZSYsmpMatrixMumps<TVar> *this_ct =
        const_cast<const TPZSYsmpMatrixMumps<TVar> *>(this);
    this_ct->SolveDirect(F, dt);

    return 0;
}

template <class TVar>
int TPZSYsmpMatrixMumps<TVar>::SolveDirect(TPZFMatrix<TVar> &F,
                                           const DecomposeType dt) const
{
    // Const version - actual solve

    // Check if already decomposed with a different scheme
    if (this->fDecomposed && this->fDecomposed != dt)
    {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix is already decomposed with other scheme");
    }

    // Check if decomposed at all
    if (!this->fDecomposed)
    {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix should've been decomposed already");
    }

    // Create temporary for solution (MUMPS overwrites RHS)
    TPZFMatrix<TVar> x(F);

    // Solve using MUMPS
    fMumpsControl.Solve(this, x, F);

    return 0;
}

// Explicit template instantiations
template class TPZSYsmpMatrixMumps<double>;
template class TPZSYsmpMatrixMumps<float>;
template class TPZSYsmpMatrixMumps<std::complex<float>>;
template class TPZSYsmpMatrixMumps<std::complex<double>>;
template class TPZSYsmpMatrixMumps<long double>;
template class TPZSYsmpMatrixMumps<std::complex<long double>>;

template <class TVar>
void TPZSYsmpMatrixMumps<TVar>::UpdateCOOFormat()
{
    const long long n = this->Rows();
    const long long nnz = this->fA.size();

    fIRN1Based.Resize(nnz);
    fJCN1Based.Resize(nnz);

    // Convert CSR to COO format
    long long k = 0;
    for (long long i = 0; i < n; i++)
    {
        const long long row_start = this->fIA[i];
        const long long row_end = this->fIA[i + 1];
        for (long long pos = row_start; pos < row_end; pos++)
        {
            fIRN1Based[k] = static_cast<MUMPS_INT>(i + 1);              // MUMPS uses 1-based indexing
            fJCN1Based[k] = static_cast<MUMPS_INT>(this->fJA[pos] + 1); // MUMPS uses 1-based indexing
            k++;
        }
    }

    if (k != nnz)
    {
        std::cerr << "ERROR: CSR to COO conversion mismatch: k=" << k << " nnz=" << nnz << std::endl;
        DebugStop();
    }

    fCOOValid = true;
}

// Explicit template instantiations for all supported types
template void TPZSYsmpMatrixMumps<double>::UpdateCOOFormat();
template void TPZSYsmpMatrixMumps<float>::UpdateCOOFormat();
template void TPZSYsmpMatrixMumps<std::complex<float>>::UpdateCOOFormat();
template void TPZSYsmpMatrixMumps<std::complex<double>>::UpdateCOOFormat();

// Note: long double instantiations for completeness, even though MUMPS
// doesn't natively support long double (uses double precision internally)
template void TPZSYsmpMatrixMumps<long double>::UpdateCOOFormat();
template void TPZSYsmpMatrixMumps<std::complex<long double>>::UpdateCOOFormat();

#endif // USING_MUMPS