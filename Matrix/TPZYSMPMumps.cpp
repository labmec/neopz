/**
 * @file
 * @brief Contains the implementation of the TPZFYsmpMatrixMumps methods.
 */

#ifdef USING_MUMPS
#include "TPZYSMPMumps.h"
#include "pzfmatrix.h"

#include <dmumps_c.h>

// Note: For non-symmetric matrices, MUMPS doesn't provide optimized
// sparse matrix-vector multiplication like MKL's sparse BLAS.

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZFYsmpMatrixMumps<TVar>::TPZFYsmpMatrixMumps() : 
    TPZRegisterClassId(&TPZFYsmpMatrixMumps::ClassId),
    TPZFYsmpMatrix<TVar>() {
}

template<class TVar>
TPZFYsmpMatrixMumps<TVar>::TPZFYsmpMatrixMumps(const TPZFYsmpMatrixMumps<TVar> &cp) :
    TPZRegisterClassId(&TPZFYsmpMatrixMumps::ClassId),
    TPZFYsmpMatrix<TVar>(cp),
    fCOOValid(false) {
    // Note: fMumpsControl is NOT copied - each matrix gets a fresh MUMPS instance
}

template<class TVar>
TPZFYsmpMatrixMumps<TVar>::TPZFYsmpMatrixMumps(TPZFYsmpMatrixMumps<TVar> &&cp) :
    TPZRegisterClassId(&TPZFYsmpMatrixMumps::ClassId),
    TPZFYsmpMatrix<TVar>(std::move(cp)),
    fMumpsControl(std::move(cp.fMumpsControl)),
    fIRN1Based(std::move(cp.fIRN1Based)),
    fJCN1Based(std::move(cp.fJCN1Based)),
    fCOOValid(cp.fCOOValid) {
    cp.fCOOValid = false;
}

template<class TVar>
TPZFYsmpMatrixMumps<TVar>& TPZFYsmpMatrixMumps<TVar>::operator=(const TPZFYsmpMatrixMumps<TVar> &copy) {
    if (this != &copy) {
        TPZFYsmpMatrix<TVar>::operator=(copy);
        // Note: fMumpsControl is NOT copied - each matrix gets a fresh MUMPS instance
        fMumpsControl = TPZMumpsSolver<TVar>();
        // Reset decomposition state
        this->SetIsDecomposed(ENoDecompose);
        fCOOValid = false;
    }
    return *this;
}

template<class TVar>
TPZFYsmpMatrixMumps<TVar>& TPZFYsmpMatrixMumps<TVar>::operator=(TPZFYsmpMatrixMumps<TVar> &&copy) {
    if (this != &copy) {
        TPZFYsmpMatrix<TVar>::operator=(std::move(copy));
        fMumpsControl = std::move(copy.fMumpsControl);
        fIRN1Based = std::move(copy.fIRN1Based);
        fJCN1Based = std::move(copy.fJCN1Based);
        fCOOValid = copy.fCOOValid;
        copy.fCOOValid = false;
    }
    return *this;
}

template<class TVar>
void TPZFYsmpMatrixMumps<TVar>::CopyFrom(const TPZMatrix<TVar> *mat) {
    auto *from = dynamic_cast<const TPZFYsmpMatrixMumps<TVar> *>(mat);
    if (from) {
        *this = *from;
        return;
    }
    
    auto *fromBase = dynamic_cast<const TPZFYsmpMatrix<TVar> *>(mat);
    if (fromBase && fromBase->IsDecomposed() == ENoDecompose) {
        *this = *fromBase;
        return;
    }
    
    PZError << __PRETTY_FUNCTION__;
    PZError << "\nERROR: Called with incompatible type\n";
    PZError << "Aborting...\n";
    DebugStop();
}

template<class TVar>
int TPZFYsmpMatrixMumps<TVar>::ClassId() const {
    return Hash("TPZFYsmpMatrixMumps") ^ TPZFYsmpMatrix<TVar>::ClassId() << 1;
}

template<class TVar>
void TPZFYsmpMatrixMumps<TVar>::MultAdd(const TPZFMatrix<TVar> &x,
                                         const TPZFMatrix<TVar> &y,
                                         TPZFMatrix<TVar> &z,
                                         const TVar alpha, const TVar beta,
                                         const int opt) const {
    // computes z = beta * y + alpha * opt(this)*x
    // z and x cannot share storage
    this->MultAddChecks(x, y, z, alpha, beta, opt);

    // MUMPS does not provide sparse matrix-vector multiplication utilities
    // like MKL's sparse BLAS. We have two options:
    // 1. Use the base class implementation (TPZFYsmpMatrix::MultAdd)
    // 2. Implement our own CSR matrix-vector product
    
    // For simplicity and correctness, we'll use the base class implementation
    TPZFYsmpMatrix<TVar>::MultAdd(x, y, z, alpha, beta, opt);
    
    // Note: If performance is critical, you could implement an optimized
    // CSR non-symmetric matrix-vector multiplication here using OpenMP.
    // The original MKL code was:
    // - Using mkl_sparse_X_create_csr to create sparse matrix handle
    // - Using mkl_sparse_X_mv for matrix-vector multiplication
    // - Supporting transpose, conjugate transpose operations
    // 
    // A custom OpenMP implementation would look like:
    // #pragma omp parallel for
    // for (int i = 0; i < nrows; i++) {
    //     TVar sum = 0;
    //     for (int j = ia[i]; j < ia[i+1]; j++) {
    //         sum += a[j] * x[ja[j]];
    //     }
    //     z[i] = alpha * sum + beta * y[i];
    // }
}

template<class TVar>
void TPZFYsmpMatrixMumps<TVar>::SetData(TPZVec<int64_t> &IA, TPZVec<int64_t> &JA, TPZVec<TVar> &A) {
    TPZFYsmpMatrix<TVar>::SetData(IA, JA, A);
    fCOOValid = false;
}

template<class TVar>
void TPZFYsmpMatrixMumps<TVar>::SetData(TPZVec<int64_t> &&IA, TPZVec<int64_t> &&JA, TPZVec<TVar> &&A) {
    TPZFYsmpMatrix<TVar>::SetData(std::move(IA), std::move(JA), std::move(A));
    fCOOValid = false;
}

template<class TVar>
void TPZFYsmpMatrixMumps<TVar>::SetIsDecomposed(DecomposeType val) {
    TPZBaseMatrix::SetIsDecomposed(val);
    if (val == ENoDecompose) {
        // Matrix may have been modified, invalidate COO format
        fCOOValid = false;
    }
    if (val) {
        fMumpsControl.fDecomposed = true;
    }
}

template<class TVar>
int TPZFYsmpMatrixMumps<TVar>::Decompose(const DecomposeType dt) {
    // Check if already decomposed with a different scheme
    if (this->fDecomposed && this->fDecomposed != dt) {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix is already decomposed with other scheme");
    }

    // Set up MUMPS parameters based on matrix properties
    if (!fMumpsControl.HasCustomSettings()) {
        const auto sysType = this->GetSymmetry();
        typename TPZMumpsSolver<TVar>::MProperty prop =
            this->IsDefPositive() ?
            TPZMumpsSolver<TVar>::MProperty::EPositiveDefinite :
            TPZMumpsSolver<TVar>::MProperty::EIndefinite;
        fMumpsControl.SetMatrixType(sysType, prop);
    }
    
    if (!fCOOValid) {
       UpdateCOOFormat();
    }
    
    // Perform decomposition using MUMPS. fIRN1Based/fJCN1Based remain valid
    // and cached: they will be reused if Decompose is called again with the
    // same sparsity pattern (e.g., reassembly with new values). They are
    // invalidated by SetIsDecomposed(ENoDecompose) and by any structural change.
    fMumpsControl.Decompose(this);

    this->SetIsDecomposed(dt);
    
    return 0;
}

template<class TVar>
int TPZFYsmpMatrixMumps<TVar>::SolveDirect(TPZFMatrix<TVar> &F,
                                            const DecomposeType dt) {
    // Non-const version
    
    // Check if already decomposed with a different scheme
    if (this->fDecomposed && this->fDecomposed != dt) {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix is already decomposed with other scheme");
    }
    
    // Decompose if not already done
    if (this->fDecomposed == ENoDecompose) {
        this->Decompose(dt);
    }
    
    // Call const version
    const TPZFYsmpMatrixMumps<TVar> *this_ct = 
        const_cast<const TPZFYsmpMatrixMumps<TVar>*>(this);
    this_ct->SolveDirect(F, dt);
    
    return 0;
}

template<class TVar>
int TPZFYsmpMatrixMumps<TVar>::SolveDirect(TPZFMatrix<TVar> &F,
                                            const DecomposeType dt) const {
    // Const version - actual solve
    
    // Check if already decomposed with a different scheme
    if (this->fDecomposed && this->fDecomposed != dt) {
        this->Error(__PRETTY_FUNCTION__,
                    "matrix is already decomposed with other scheme");
    }
    
    // Check if decomposed at all
    if (!this->fDecomposed) {
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
#ifdef MUMPS_HAVE_SINGLE
template class TPZFYsmpMatrixMumps<float>;
#endif
#ifdef MUMPS_HAVE_DOUBLE
template class TPZFYsmpMatrixMumps<double>;
template class TPZFYsmpMatrixMumps<long double>;
#endif
#ifdef MUMPS_HAVE_COMPLEX
template class TPZFYsmpMatrixMumps<std::complex<float>>;
#endif
#ifdef MUMPS_HAVE_COMPLEX16
template class TPZFYsmpMatrixMumps<std::complex<double>>;
template class TPZFYsmpMatrixMumps<std::complex<long double>>;
#endif

template<class TVar>
void TPZFYsmpMatrixMumps<TVar>::UpdateCOOFormat() {
    const long long n = this->Rows();
    const long long nnz = this->fA.size();
    
    fIRN1Based.Resize(nnz);
    fJCN1Based.Resize(nnz);
    
    // Convert CSR to COO format
    long long k = 0;
    for (long long i = 0; i < n; i++) {
        const long long row_start = this->fIA[i];
        const long long row_end = this->fIA[i + 1];
        for (long long pos = row_start; pos < row_end; pos++) {
            fIRN1Based[k] = static_cast<MUMPS_INT>(i + 1);  // MUMPS uses 1-based indexing
            fJCN1Based[k] = static_cast<MUMPS_INT>(this->fJA[pos] + 1);  // MUMPS uses 1-based indexing
            k++;
        }
    }
    
    if (k != nnz) {
        std::cerr << "ERROR: CSR to COO conversion mismatch: k=" << k << " nnz=" << nnz << std::endl;
        DebugStop();
    }
    
    fCOOValid = true;
}

#ifdef MUMPS_HAVE_DOUBLE
template void TPZFYsmpMatrixMumps<double>::UpdateCOOFormat();
template void TPZFYsmpMatrixMumps<long double>::UpdateCOOFormat();
#endif
#ifdef MUMPS_HAVE_SINGLE
template void TPZFYsmpMatrixMumps<float>::UpdateCOOFormat();
#endif
#ifdef MUMPS_HAVE_COMPLEX
template void TPZFYsmpMatrixMumps<std::complex<float>>::UpdateCOOFormat();
#endif
#ifdef MUMPS_HAVE_COMPLEX16
template void TPZFYsmpMatrixMumps<std::complex<double>>::UpdateCOOFormat();
template void TPZFYsmpMatrixMumps<std::complex<long double>>::UpdateCOOFormat();
#endif

#endif