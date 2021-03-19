//
//  TPZPardisoControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#include "TPZPardisoControl.h"

#ifdef USING_MKL

#include "mkl_pardiso.h"
#include "pzsysmp.h"
#include "pzysmp.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.pardisocontrol");
#endif
//#define Release_Memory_Q

/// empty constructor (non symetric and LU decomposition
template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl() : fSystemType(ENonSymmetric),
fStructure(EStructureSymmetric), fProperty(EIndefinite), fPardisoControl(), fHandle(0),
fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0),
fNonSymmetricSystem(0), fSymmetricSystem(0)
{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
//    fMatrixType = MatrixType();
}


template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl(MSystemType systemtype, MProperty prop) : fSystemType(systemtype),
        fStructure(EStructureSymmetric), fProperty(prop), fPardisoControl(), fHandle(0),
        fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0),
        fNonSymmetricSystem(0), fSymmetricSystem(0)

{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
//    fMatrixType = MatrixType();
}

/// change the matrix type
// this method should only be called if the pardiso control is zero (non initialized)
template<class TVar>
void TPZPardisoControl<TVar>::SetMatrixType(MSystemType systemtype, MProperty prop)
{
    fSystemType = systemtype;
    fProperty = prop;
    fMatrixType = MatrixType();
}

//fSystemType(systemtype),
//fStructure(EStructureSymmetric), fProperty(prop), fPardisoControl(), fHandle(0),
//fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0)

template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl(const TPZPardisoControl<TVar> &copy) : fSystemType(copy.fSystemType),
fStructure(copy.fStructure), fProperty(copy.fProperty), fPardisoControl(copy.fPardisoControl), fHandle(copy.fHandle),
fParam(copy.fParam), fMax_num_factors(copy.fMax_num_factors), fMatrix_num(copy.fMatrix_num), fMessageLevel(copy.fMessageLevel),
fError(copy.fError), fPermutation(copy.fPermutation), fMatrixType(copy.fMatrixType),
fSymmetricSystem(copy.fSymmetricSystem), fNonSymmetricSystem(copy.fNonSymmetricSystem)
{
    
}

template<class TVar>
TPZPardisoControl<TVar> &TPZPardisoControl<TVar>::operator=(const TPZPardisoControl &copy)
{
    fSystemType = copy.fSystemType;
    fStructure = copy.fStructure;
    fProperty = copy.fProperty;
    fPardisoControl = copy.fPardisoControl;
    fHandle = copy.fHandle;
    fParam = copy.fParam;
    fMax_num_factors = copy.fMax_num_factors;
    fMatrix_num = copy.fMatrix_num;
    fMessageLevel = copy.fMessageLevel;
    fError = copy.fError;
    fPermutation = copy.fPermutation;
    fMatrixType = copy.fMatrixType;
    fSymmetricSystem = copy.fSymmetricSystem;
    fNonSymmetricSystem = copy.fNonSymmetricSystem;
    return *this;
}


//enum MSystemType {ESymmetric, EHermitian, EnonSymmetric};
//
//enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
//
//enum MProperty {EPositiveDefinite, EIndefinite};

template<class TVar>
int DataType(TVar a)
{
    DebugStop();
	return 0;
}


template<>
int DataType(double a)
{
    return 0;
}

template<>
int DataType(float a)
{
    return 1;
}


template<class TVar>
long long TPZPardisoControl<TVar>::MatrixType()
{
    // should not happen
    if (fStructure == EStructureNonSymmetric) {
        DebugStop();
    }
    if (fSystemType == ESymmetric && fProperty == EIndefinite) {
        if(fMatrixType != 0 && fMatrixType != -2)
        {
            DebugStop();
        }
        else if(fMatrixType == -2)
        {
            return -2;
        }
        fMatrixType = -2;
    }
    if (fSystemType == ESymmetric && fProperty == EPositiveDefinite) {
        if(fMatrixType != 0 && fMatrixType != 2)
        {
            DebugStop();
        }
        else if(fMatrixType == 2)
        {
           return 2;
        }
        fMatrixType = 2;
    }
    if (fSystemType == ENonSymmetric && fStructure == EStructureSymmetric) {
        fMatrixType = 11;
    }
    if (fSystemType == ENonSymmetric && fProperty == EPositiveDefinite) {
        DebugStop();
    }
    
//    for (long long i=0; i<64; i++) {
//        long long val = fHandle[i];
//        if (val) {
//            DebugStop();
//        }
//    }
    
    int param[64] = {0};
    int matrixtype = fMatrixType;
    pardisoinit(fHandle,&matrixtype,param);
    for (int i=0; i<64; i++) {
        fParam[i] = param[i];
    }

    fParam[34] = 1;
    return fMatrixType;
}


template<class TVar>
void TPZPardisoControl<TVar>::Decompose()
{
    long long n=0;
    TVar bval = 0., xval = 0.;
    TVar *a,*b = &bval, *x = &xval;
    long long *ia,*ja;
    if (fSymmetricSystem) {
        if (fSymmetricSystem->Rows()==0) {
            return;
        }
        a = &(fSymmetricSystem->fA[0]);
        ia = (long long *) &(fSymmetricSystem->fIA[0]);
        ja = (long long *) &(fSymmetricSystem->fJA[0]);
        n = fSymmetricSystem->Rows();
    }
    if (fNonSymmetricSystem) {
        a = &(fNonSymmetricSystem->fA[0]);
        ia = (long long *) &(fNonSymmetricSystem->fIA[0]);
        ja = (long long *) &(fNonSymmetricSystem->fJA[0]);
        n = fNonSymmetricSystem->Rows();
        
    }

    long long *perm = 0,nrhs = 0;
    long long Error = 0;
    nrhs = 0;
    fPermutation.resize(n);
    perm = &fPermutation[0];
    fParam[34] = 1;
    // Do not use OOC
    fParam[59] = 0;
    /// analyse and factor the equations
    long long phase = 12;
    fPermutation.resize(n);
    for (long long i=0; i<n; i++) {
        fPermutation[i] = i;
    }
    perm = &fPermutation[0];
    
    /// analyse and factor the equations
    // LU preconditioned CGS (10*L+K) where K={1:CGS,2:CG} and L=10^-L stopping threshold
    if (fProperty == EIndefinite) {
        fParam[4] = 1;
        if(fSystemType == ESymmetric){ // The factorization is always computed as required by phase.
            fParam[3 ] = 10*6+2;
        }else{ // CGS iteration replaces the computation of LU. The preconditioner is LU that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
            fParam[3 ] = 10*6+1;
            fParam[10] = 1;
            fParam[12] = 1;
        }
    }else{
        
        if(fSystemType == ESymmetric){ // CGS iteration for symmetric positive definite matrices replaces the computation of LLT. The preconditioner is LLT that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
            fParam[3 ] = 10*6+2;
        }else{
            fParam[3 ] = 10*6+1;
            fParam[10] = 1;
            fParam[12] = 1;
        }
    }
    
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
}

/// Use the decomposed matrix to invert the system of equations
template<class TVar>
void TPZPardisoControl<TVar>::Solve(TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const
{
    long long n=0;
    TVar *a,*b, *x;
    long long *ia,*ja;
    if (fSymmetricSystem) {
        if(fSymmetricSystem->Rows() == 0)
        {
            return;
        }
        a = &(fSymmetricSystem->fA[0]);
        ia = (long long *) &(fSymmetricSystem->fIA[0]);
        ja = (long long *) &(fSymmetricSystem->fJA[0]);
        n = fSymmetricSystem->Rows();
    }
    if (fNonSymmetricSystem) {
        a = &(fNonSymmetricSystem->fA[0]);
        ia = (long long *) &(fNonSymmetricSystem->fIA[0]);
        ja = (long long *) &(fNonSymmetricSystem->fJA[0]);
        n = fNonSymmetricSystem->Rows();
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
    b = &rhs(0,0);
    x = &sol(0,0);
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
    
#ifdef Release_Memory_Q
    phase = -1;
    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm, &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    if (fSymmetricSystem) {
        fSymmetricSystem->SetIsDecomposed(0);
    }
    if (fNonSymmetricSystem) {
        fNonSymmetricSystem->SetIsDecomposed(0);
    }
#ifdef PZDEBUG
    std::cout << "Pardiso:: release memory complete. \n";
#endif
#endif
}

template<class TVar>
TPZPardisoControl<TVar>::~TPZPardisoControl()
{
    long long phase = -1;
    long long n=1;
    long long av,bv,xv;
    void *a= &av,*b = &bv, *x = &xv;
    long long ia,ja,perm,nrhs = 1;
    long long Error = 0;
    
    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, &ia, &ja, &perm,
                &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    
    if (Error) {
        DebugStop();
    }
    
}

template<class TVar>
void TPZPardisoControl<TVar>::Error_check(int error) const {
    
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


template class TPZPardisoControl<double>;
template class TPZPardisoControl<long double>;
template class TPZPardisoControl<float>;



#endif
