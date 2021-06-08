//
//  TPZPardisoSolver.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#include "TPZPardisoSolver.h"

#include "pzsysmp.h"
#include "pzysmp.h"
#include "pzlog.h"


#ifdef USING_MKL

#include "mkl_pardiso.h"
#else
#define NOMKL                                                   \
    PZError<<"The class TPZPardisoSolver should not be used ";  \
    PZError<<"if NeoPZ was configured with USING_MKL=OFF\n";    \
    PZError<<"Aborting..."<<std::endl;                          \
    DebugStop();
#endif 
/*auxiliary functions*/
void Error_check(int error);

template<class TVar>
int DataType([[maybe_unused]]TVar a);

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.pardisocontrol");
#endif

/// empty constructor (non symetric and LU decomposition
template<class TVar>
TPZPardisoSolver<TVar>::TPZPardisoSolver() :
    TPZMatrixSolver<TVar>(),fPardisoControl(),fParam(64,0)
{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
}


template<class TVar>
TPZPardisoSolver<TVar>::TPZPardisoSolver(MSystemType systemtype,
                                         MStructure structure,
                                         MProperty prop) :
    fSystemType(systemtype), fStructure(structure), fProperty(prop),
    fPardisoControl(), fParam(64,0)
{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
}


template<class TVar>
TPZPardisoSolver<TVar>::~TPZPardisoSolver()
{
#ifdef USING_MKL
    long long phase = -1;
    long long n=1;
    long long av,bv,xv;
    void *a= &av,*b = &bv, *x = &xv;
    long long ia,ja,perm,nrhs = 1;
    long long Error = 0;
    if(fPardisoInitialized)
        pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType,
                    &phase, &n, a, &ia, &ja, &perm,
                    &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    
    if (Error) {
        DebugStop();
    }
#endif
    //we should NOT delete fSymmetricSystem and fNonSymmetricSystem
}
template<class TVar>
void TPZPardisoSolver<TVar>::SetMatrix(TPZAutoPointer<TPZBaseMatrix> refmat)
{
    auto *symSystem =
        dynamic_cast<TPZMatrix<TVar>*>(refmat.operator->());
    auto *nSymSystem =
        dynamic_cast<TPZMatrix<TVar>*>(refmat.operator->());
#ifdef PZDEBUG
    if(!symSystem && !nSymSystem){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"This solver is only compatible with sparse matrices.\nAborting...\n";
        DebugStop();
    }
#endif
    
    fDecomposed = refmat->IsDecomposed();
    /*the following variables could have been initialized by the user*/
    if (fStructure == MStructure::ENonInitialized)
        fStructure = symSystem ? MStructure::ESymmetric : MStructure::ENonSymmetric;
    
    if(fSystemType == MSystemType::ENonInitialized ||
       fProperty == MProperty::ENonInitialized){
        const MProperty prop = refmat->IsDefPositive() ?
            MProperty::EPositiveDefinite : MProperty::EIndefinite;
        const MSystemType sym = refmat->IsSymmetric() ?
            MSystemType::ESymmetric : MSystemType::ENonSymmetric;
        SetMatrixType(sym,prop);
    }
    TPZMatrixSolver<TVar>::SetMatrix(refmat);
}


template<class TVar>
void TPZPardisoSolver<TVar>::Solve(const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol,
                                   TPZFMatrix<TVar> *res)
{
    if(!this->Matrix()) {
		PZError << __PRETTY_FUNCTION__;
        PZError<< " called without a matrix pointer\n";
		DebugStop();
	}
    if(!fDecomposed) Decompose();
    sol = rhs;
    Solve(this->Matrix().operator->(),rhs,sol);

    if(res) res->Redim(rhs.Rows(),rhs.Cols());
}
template<class TVar>
void TPZPardisoSolver<TVar>::Decompose()
{
    auto *tmpSym =
        dynamic_cast<TPZSYsmpMatrix<TVar>*>(this->Matrix().operator->());
    auto *tmpNSym =
        dynamic_cast<TPZFYsmpMatrix<TVar>*>(this->Matrix().operator->());
    if(tmpSym){
        Decompose(tmpSym);
    }else if(tmpNSym){
        Decompose(tmpNSym);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"This solver is only compatible with sparse matrices.\nAborting...\n";
        DebugStop();
    }
}

template<class TVar>
void TPZPardisoSolver<TVar>::Decompose(TPZMatrix<TVar> *mat)
{
#ifndef USING_MKL
    NOMKL
#else
    auto *symSystem =
        dynamic_cast<TPZSYsmpMatrix<TVar>*>(mat);
    auto *nSymSystem =
        dynamic_cast<TPZFYsmpMatrix<TVar>*>(mat);

    long long n=0;
    TVar bval = 0., xval = 0.;
    TVar *a,*b = &bval, *x = &xval;
    long long *ia,*ja;
    if (symSystem) {
        if (symSystem->Rows()==0) {
            return;
        }
        a = &(symSystem->fA[0]);
        ia = (long long *) &(symSystem->fIA[0]);
        ja = (long long *) &(symSystem->fJA[0]);
        n = symSystem->Rows();
    }
    if (nSymSystem) {
        a = &(nSymSystem->fA[0]);
        ia = (long long *) &(nSymSystem->fIA[0]);
        ja = (long long *) &(nSymSystem->fJA[0]);
        n = nSymSystem->Rows();
        
    }

    long long *perm = 0,nrhs = 0;
    long long Error = 0;
    nrhs = 0;    
    
    /// analyse and factor the equations
    long long phase = 12;
    fPermutation.resize(n);
    for (long long i=0; i<n; i++) {
        fPermutation[i] = i;
    }
    perm = &fPermutation[0];

    /*@orlandini: most values were taken from eigen PardisoSupport.
      There are other values being set at:
      TPZPardisoSolver<TVar>::MatrixType()
     */

    
    //fParam[0] No default values
    fParam[0] = 1;
    //fParam[1]  use Metis for the ordering
    fParam[1] = 2;
    /*fParam[3]  Preconditioned CGS/CG. 
       0 = // No iterative-direct algorithm
       10*L+K
       L = stoppping criterion: 10^-L
       K = 
          0: The factorization is always computed as required by phase
          1: CGS iteration replaces the computation of LU. 
             The preconditioner is LU that was computed at a previous step
             (the first step or last step with a failure) in a sequence of
             solutions needed for identical sparsity patterns.
          2: CGS iteration for symmetric positive definite matrices
             Replaces the computation of LLt. The preconditioner is LLT
             that was computed at a previous step
             (the first step or last step with a failure)
             in a sequence of solutions needed for identical sparsity patterns. 
     */
    //fParam[4]  No user fill-in reducing permutation
    if constexpr (!is_complex<TVar>::value){
        fParam[3] = fSystemType == MSystemType::ESymmetric ? 10*6+2 : 10*6+1;
        if(fProperty == MProperty::EIndefinite) fParam[4] =1;
    }else{
        fParam[3] = 0;
        fParam[4] = 0;
    }
    //fParam[9]  Perturb the pivot elements with 1E-fParam[9]
    fParam[9] = 13;
    //fParam[26] Whether to check matrix data
    fParam[26] = 1;
    //fParam[59]  Do not use OOC
    fParam[59] = 0;
    
    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm,
                &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    if (Error) {
        Error_check(int(Error));
        std::cout << __PRETTY_FUNCTION__ << " error code " << Error << std::endl;
        DebugStop();
    }
#ifdef PZDEBUG
    std::cout << "Pardiso:: decomposition complete. \n";
#endif
    fDecomposed = true;
#endif
}

/// Use the decomposed matrix to invert the system of equations
template<class TVar>
void TPZPardisoSolver<TVar>::Solve(const TPZMatrix<TVar> *mat,
                                   const TPZFMatrix<TVar> &rhs,
                                   TPZFMatrix<TVar> &sol) const
{
#ifndef USING_MKL
    NOMKL
#else
#ifdef PZDEBUG
    if(!fDecomposed){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nError: Matrix has not been decomposed.\nAborting..."<<std::endl;
        DebugStop();
    }
#endif
    auto *symSystem =
        dynamic_cast<const TPZSYsmpMatrix<TVar>*>(mat);
    auto *nSymSystem =
        dynamic_cast<const TPZFYsmpMatrix<TVar>*>(mat);
    long long n=0;
    TVar *a,*b, *x;
    long long *ia,*ja;
    if (symSystem) {
        if(symSystem->Rows() == 0)
        {
            return;
        }
        a = &(symSystem->fA[0]);
        ia = (long long *) &(symSystem->fIA[0]);
        ja = (long long *) &(symSystem->fJA[0]);
        n = symSystem->Rows();
    }
    if (nSymSystem) {
        a = &(nSymSystem->fA[0]);
        ia = (long long *) &(nSymSystem->fIA[0]);
        ja = (long long *) &(nSymSystem->fJA[0]);
        n = nSymSystem->Rows();
    }
    
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "The pardiso control vector is\n";
        sout << fParam << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    long long *perm,nrhs;
    long long Error = 0;
    nrhs = rhs.Cols();
    n = rhs.Rows();
    b = &rhs.g(0,0);
    x = &sol.g(0,0);
    perm = &fPermutation[0];
    /// forward and backward substitution
    long long phase = 33;
    
    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm,
                &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    
    if(fParam[19]>150){
        std::cout << "Pardiso:: Number of iterations " << fParam[19] << " > 150, calling numerical factorization... " << std::endl;
        phase = 23;
        pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm,
                    &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    }
    
    int rest = abs(fParam[19]%10); // CG/CGS error report
    if(fParam[19] <= 0){
        switch (rest) {
            case 1:{
                std::cout << "Pardiso:: fluctuations of the residuum are too large. " << std::endl;
            }
                break;
                
            case 2:{
                std::cout << "Pardiso:: Slow convergence - Main matrix and matrix for preconditioner differ a lot. " << std::endl;
            }
                break;
            case 3:{
                std::cout <<"Pardiso:: stopping criterion is not reached at max_it_cgs. " << std::endl;
            }
                break;
            case 4:{
                std::cout << "Pardiso:: perturbed pivots caused iterative refinement. " << std::endl;
            }
                break;
                
            case 5:{
                std::cout << "Pardiso:: factorization is too fast for this matrix. It is better to use the factorization method with iparm[3] = 0 " << std::endl;
                fParam[3] = 0;
            }
                break;
            case 6:{
                std::cout << "Pardiso:: There is not a diagnostig. " << std::endl;
            }
                break;
            default:
                break;
        }
        
    }
    
    
    if (Error<0) {
        Error_check(int(Error));
        std::cout << "Pardiso:: Calling a numerical factorization. \n";
        phase = 23;
        pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm,
                    &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    }
    
    if (Error) {
        Error_check(int(Error));
        DebugStop();
    }
    
#ifdef PZDEBUG
    std::cout << "Pardiso:: linear solve complete. \n";
#endif
    
#endif//USING_MKL
}

template<class TVar>
void TPZPardisoSolver<TVar>::SetMessageLevel(int lvl){
    lvl == 0 ? fMessageLevel = 0 : fMessageLevel = 1;
}

template<class TVar>
TPZPardisoSolver<TVar> *TPZPardisoSolver<TVar>::Clone() const{
    return new TPZPardisoSolver<TVar>(*this);
}

template<class TVar>
void TPZPardisoSolver<TVar>::SetMatrixType(MSystemType systemtype, MProperty prop)
{
    if(fSystemType != MSystemType::ENonInitialized){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR:\n";
        PZError<<"this function should not be called on an initialized instance.\n";
        PZError<<"Aborting..."<<std::endl;
        DebugStop();
    }
    fSystemType = systemtype;
    fProperty = prop;
    fMatrixType = MatrixType();
}

template<class TVar>
long long TPZPardisoSolver<TVar>::MatrixType()
{
    if (fStructure == MStructure::ENonSymmetric){
        if(fSystemType == MSystemType::ESymmetric){
            if constexpr (is_complex<TVar>::value){
                fMatrixType = 3;
            }else{
                fMatrixType = 1;
            }
        }else{
            if constexpr (is_complex<TVar>::value){
                fMatrixType = 13;
            }else{
                fMatrixType = 11;
            }
        }
    }else{
        if(fSystemType != MSystemType::ESymmetric){
            if constexpr (is_complex<TVar>::value){
                fMatrixType = 3;
            }else{
                fMatrixType = 1;
            }
        }
        else if(fProperty == MProperty::EPositiveDefinite){
            if constexpr (is_complex<TVar>::value){
                fMatrixType = 4;
            }else{
                fMatrixType = 2;
            }
        }else{
            if constexpr (is_complex<TVar>::value){
                fMatrixType = -4;
            }else{
                fMatrixType = -2;
            }
        }
    }
    
    int param[64] = {0};
    int matrixtype = fMatrixType;
#ifdef USING_MKL
    pardisoinit(fHandle,&matrixtype,param);
    fPardisoInitialized = true;
#endif
    for (int i=0; i<64; i++) {
        fParam[i] = param[i];
    }
    //fParam[10]  Use nonsymmetric permutation and scaling MPS
    fParam[10] = fSystemType == MSystemType::ESymmetric ? 0 : 1;
    //fParam[12]  Maximum weighted matching algorithm is switched-off (default for symmetric).
    fParam[12] = fSystemType == MSystemType::ESymmetric ? 0 : 1;
    //fParam[27] float or double
    fParam[27] = DataType((TVar)0);
    //fParam[34]  zero-based indexing
    fParam[34] = 1;
    //Use CSR
    fParam[36] = 0;
    return fMatrixType;
}


template<class TVar>
int DataType([[maybe_unused]]TVar a)
{
    DebugStop();
	return 0;
}


template<>
int DataType([[maybe_unused]]double a)
{
    return 0;
}

template<>
int DataType([[maybe_unused]]float a)
{
    return 1;
}

template<>
int DataType([[maybe_unused]]std::complex<double> a)
{
    return 0;
}

template<>
int DataType([[maybe_unused]]std::complex<float> a)
{
    return 1;
}

void Error_check(int error) {
    
    switch (error) {
        case -1:
            std::cout << "Pardiso:: Input inconsistent." << std::endl;
            break;
        case -2:
            std::cout << "Pardiso:: Not enough memory." << std::endl;
            break;
        case -3:
            std::cout << "Pardiso:: Reordering problem." << std::endl;
            break;
        case -4:
            std::cout << "Pardiso:: Zero pivot, numerical fact. or iterative refinement problem. " << std::endl;
            break;
        case -5:
            std::cout << "Pardiso:: Unclassified (internal) error. " << std::endl;
            break;
        case -6:
            std::cout << "Pardiso:: Preordering failed (matrix types 11, 13 only). " << std::endl;
            break;
        case -7:
            std::cout << "Pardiso:: Diagonal matrix problem. " << std::endl;
            break;
        case -8:
            std::cout << "Pardiso:: 32-bit integer overflow problem. " << std::endl;
            break;
        default:
            std::cout << "Pardiso:: There is not a explanation. " << std::endl;
            break;
    }
    
}

template class TPZPardisoSolver<double>;
template class TPZPardisoSolver<long double>;
template class TPZPardisoSolver<float>;
template class TPZPardisoSolver<std::complex<float>>;
template class TPZPardisoSolver<std::complex<double>>;
template class TPZPardisoSolver<std::complex<long double>>;
