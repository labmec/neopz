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

#define ISM_new

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
    
//    void pardiso (_MKL_DSS_HANDLE_t pt, const MKL_INT *maxfct, const MKL_INT *mnum, const
//                  MKL_INT *mtype, const MKL_INT *phase, const MKL_INT *n, const void *a, const MKL_INT
//                  *ia, const MKL_INT *ja, MKL_INT *perm, const MKL_INT *nrhs, MKL_INT *iparm, const
//                  MKL_INT *msglvl, void *b, void *x, MKL_INT *error);
    
    for (long long i=0; i<64; i++) {
        long long val = fHandle[i];
        if (val) {
            DebugStop();     
        }
    }
//    void pardiso_64( _MKL_DSS_HANDLE_t,       const long long int *, const long long int *, const long long int *,
//                    const long long int *, const long long int *, const void *,          const long long int *,
//                    const long long int *, long long int *, const long long int *, long long int *,
//                    const long long int *, void *,                void *,                long long int * );

    int param[64] = {0};
    int matrixtype = fMatrixType;
    pardisoinit(fHandle,&matrixtype,param);
    for (int i=0; i<64; i++) {
        fParam[i] = param[i];
    }
//    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, &a, &ia, &ja, &perm,
//                &nrhs, &fParam[0], &fMessageLevel, &b, &x, &Error);
    
    TVar toto = 0;
    fParam[27] = ::DataType(toto);
    /// establish that the datastructures are zero based
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
//    for (int i=0; i<n+1; i++) {
//        std::cout << ia[i] << ' ';
//    }
//    std::cout << std::endl;
//    for (int i=0; i<ia[n]; i++) {
//        std::cout << ja[i] << ' ' << a[i] << "| ";
//    }
//    std::cout << std::endl;
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
    
#ifndef ISM_new
    fParam[4 ] = 1; // user permutation PERM
#else
    
    /// analyse and factor the equations
    // LU preconditioned CGS (10*L+K) where K={1:CGS,2:CG} and L=10^-L stopping threshold
    if (fProperty == EIndefinite) {
        if(fSystemType == ESymmetric){ // The factorization is always computed as required by phase.
            fParam[3 ] = 10*6+2;
            fParam[10] = 1;
            fParam[12] = 1;
        }else{ // CGS iteration replaces the computation of LU. The preconditioner is LU that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
            fParam[3 ] = 10*6+1;
            fParam[10] = 1;
            fParam[12] = 1;
        }
    }else{

        if(fSystemType == ESymmetric){ // CGS iteration for symmetric positive definite matrices replaces the computation of LLT. The preconditioner is LLT that was computed at a previous step (the first step or last step with a failure) in a sequence of solutions needed for identical sparsity patterns.
            fParam[3 ] = 10*6+2;
            fParam[10] = 0;
            fParam[12] = 0;
        }else{
            fParam[3 ] = 10*6+1;
            fParam[10] = 1;
            fParam[12] = 1;
        }
    }
#endif

    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm,
                &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    if (Error) {
        std::cout << __PRETTY_FUNCTION__ << " error code " << Error << std::endl;
        DebugStop();
    }
    std::cout << "Pardiso:: decomposition complete. \n";
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
    }
    if (fNonSymmetricSystem) {
        a = &(fNonSymmetricSystem->fA[0]);
        ia = (long long *) &(fNonSymmetricSystem->fIA[0]);
        ja = (long long *) &(fNonSymmetricSystem->fJA[0]);
        
    }

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
    
//    std::cout << "Norm RHS " << Norm(rhs) << std::endl;
//    std::cout << "Norm sol " << Norm(sol) << std::endl;
//    rhs.Print("rhs");
//    sol.Print("sol");
    if (Error) {
        DebugStop();
    }
    std::cout << "Pardiso:: linear solve complete. \n";
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
    
    //    void pardiso_64( _MKL_DSS_HANDLE_t,       const long long int *, const long long int *, const long long int *,
    //                    const long long int *, const long long int *, const void *,          const long long int *,
    //                    const long long int *, long long int *, const long long int *, long long int *,
    //                    const long long int *, void *,                void *,                long long int * );
    
    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, &ia, &ja, &perm,
                &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);
    
    if (Error) {
        DebugStop();
    }
    
}


template class TPZPardisoControl<double>;
template class TPZPardisoControl<long double>;
template class TPZPardisoControl<float>;



#endif
