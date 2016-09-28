//
//  TPZPardisoControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#include "TPZPardisoControl.h"
#ifdef USING_MKL

#include "pzsysmp.h"
#include "pzysmp.h"

/// empty constructor (non symetric and LU decomposition
template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl() : fSystemType(ENonSymmetric),
fStructure(EStructureSymmetric), fProperty(EIndefinite), fPardisoControl(), fHandle(0),
fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0),
fNonSymmetricSystem(0), fSymmetricSystem(0)
{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
    fMatrixType = MatrixType();
}


template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl(MSystemType systemtype, MProperty prop) : fSystemType(systemtype),
        fStructure(EStructureSymmetric), fProperty(prop), fPardisoControl(), fHandle(0),
        fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0),
        fNonSymmetricSystem(0), fSymmetricSystem(0)

{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
    fMatrixType = MatrixType();
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
        fMatrixType = -2;
    }
    if (fSystemType == ESymmetric && fProperty == EPositiveDefinite) {
        fMatrixType = 2;
    }
    if (fSystemType == ENonSymmetric && fStructure == EStructureSymmetric) {
        fMatrixType = 1;
    }
    if (fSystemType == ENonSymmetric && fProperty == EPositiveDefinite) {
        DebugStop();
    }
    
//    void pardiso (_MKL_DSS_HANDLE_t pt, const MKL_INT *maxfct, const MKL_INT *mnum, const
//                  MKL_INT *mtype, const MKL_INT *phase, const MKL_INT *n, const void *a, const MKL_INT
//                  *ia, const MKL_INT *ja, MKL_INT *perm, const MKL_INT *nrhs, MKL_INT *iparm, const
//                  MKL_INT *msglvl, void *b, void *x, MKL_INT *error);
    
    long long phase = 0;
    long long n=1;
    long long av,bv,xv;
    void *a= &av,*b = &bv, *x = &xv;
    long long ia,ja,perm,nrhs = 1;
    long long Error = 0;
    
    for (long i=0; i<64; i++) {
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
    
    TVar toto;
    fParam[27] = ::DataType(toto);
    /// establish that the datastructures are zero based
    fParam[34] = 1;

    if (Error) {
        DebugStop();
    }
    return fMatrixType;
}


template<class TVar>
void TPZPardisoControl<TVar>::Decompose()
{
    if(sizeof(long long) != sizeof(long))
    {
        DebugStop();
    }
    long long n=0;
    TVar bval = 0., xval = 0.;
    TVar *a,*b = &bval, *x = &xval;
    long long *ia,*ja;
    if (fSymmetricSystem) {
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
    for (int i=0; i<n; i++) {
        bool hasdiag = false;
        if (ia[i+1] == ia[i]+1 && ja[ia[i]] == i ) {
            hasdiag = true;
        }
        for (int j = ia[i]; j < ia[i+1]-1; j++) {
            if (ja[j] == i || ja[j+1] == i) {
                hasdiag = true;
            }
            if (ja[j] >= ja[j+1]) {
                DebugStop();
            }
        }
        if (hasdiag == false) {
            std::cout << "i = " << i << std::endl;
            for (int k=ia[i]; k<ia[i+1]; k++) {
                std::cout << ja[k] << " ";
            }
            std::cout << std::endl;
            DebugStop();
        }
    }
//    std::cout << std::endl;
    long long *perm = 0,nrhs = 0;
    long long Error = 0;
    nrhs = 0;
    fPermutation.resize(n);
    
    for(long i = 0; i < n; i++){
        fPermutation[i] = i;
    }
    
    perm = &fPermutation[0];
    fParam[34] = 1;
    
    
    
    /// analyse and factor the equations
// https://github.com/boyle/meagre-crowd/blob/master/src/solver_pardiso.c
    /// testing RTFM = READ THE FUCKING MANUAL
    if (fProperty == EIndefinite && fSystemType == ESymmetric) {
        fParam[4] = 1;
        
//        fParam[10] = 10;
//        fParam[5] = 1;
        
//        fParam[12] = 0;

        // Note: other values unused
//        fParam[1 ] = 0; // use default values (2, 4..64), 3 MUST be set .. this will overwrite the following config with defaults (it mostly here for documentation)
//        fParam[2 ] = 2; // fill-in reducing ordering (0: min-degree, 2: METIS)
//        fParam[3 ] = 1; // number of processors: must match OMP_NUM_THREADS TODO  -- NOTE this is an *upper-limit* on the number of processors...
//        fParam[4 ] = 0; // LU preconditioned CGS (10*L+K) where K=1:CGS,2:CG L=10^-L stopping threshold
//        fParam[5 ] = 1; // user permutation PERM
//        fParam[6 ] = 0; // overwrite b with x "write solution on x"
        
//        fParam[5] = 1;
//        fParam[10] = 8;
//        fParam[11] = 1;
//        fParam[12] = 1;
//        fParam[13] = 1;
        
//        p->IPARM( 10 ) = 8;
//        p->IPARM( 11 ) = 1; // scaling vectors (MTYPE=11,13, if(IPARM(13) then also -2,-4,6)
//        // use IPARM(11) = IPARM(13) = 1 for highly indefinite matrices: requires numerical values of A in analysis
//        p->IPARM( 12 ) = 0; // solve transpose matrix A^t x = b TODO support this
//        // improved accuracy, 1: normal matchines, 2: advanced matchings (highly indefinite matrices)
//        if (( p->MTYPE == REAL_UNSYM ) || ( p->MTYPE == COMPLEX_UNSYM ) )
//            p->IPARM( 13 ) = 1;
    }
    long long phase = 12;
    
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    std::cout  << " fParam = "<< fParam[14] << std::endl;
    
    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, a, ia, ja, perm,
                &nrhs, &fParam[0], &fMessageLevel, b, x, &Error);

    std::cout  << " fParam = "<< fParam[14] << std::endl;
    std::cout << "Pardiso:: Done." << std::endl;
    if (Error) {
        DebugStop();
    }
}

/// Use the decomposed matrix to invert the system of equations
template<class TVar>
void TPZPardisoControl<TVar>::Solve(TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &sol) const
{
    long long n=0;
    TVar *a,*b, *x;
    long long *ia,*ja;
    if (fSymmetricSystem) {
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
    
    double toto;
    //    void pardiso_64( _MKL_DSS_HANDLE_t,       const long long int *, const long long int *, const long long int *,
    //                    const long long int *, const long long int *, const void *,          const long long int *,
    //                    const long long int *, long long int *, const long long int *, long long int *,
    //                    const long long int *, void *,                void *,                long long int * );
    
    int matrixtype = fMatrixType;
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
