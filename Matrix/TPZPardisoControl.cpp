//
//  TPZPardisoControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#include "TPZPardisoControl.h"

template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl(MSystemType systemtype, MProperty prop) : fSystemType(systemtype),
        fStructure(EStructureSymmetric), fProperty(prop), fPardisoControl(64,0), fHandle(0),
        fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0)
{
    fHandle = &fPardisoControl[0];
    fMatrixType = MatrixType();
}

//enum MSystemType {ESymmetric, EHermitian, EnonSymmetric};
//
//enum MStructure {EStructureSymmetric, EStructureNonSymmetric};
//
//enum MProperty {EPositiveDefinite, EIndefinite};


template<>
long long TPZPardisoControl<double>::MatrixType()
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
    if (fSystemType == EnonSymmetric && fStructure == EStructureSymmetric) {
        fMatrixType = 1;
    }
    if (fSystemType == EnonSymmetric && fProperty == EPositiveDefinite) {
        DebugStop();
    }
    
//    void pardiso (_MKL_DSS_HANDLE_t pt, const MKL_INT *maxfct, const MKL_INT *mnum, const
//                  MKL_INT *mtype, const MKL_INT *phase, const MKL_INT *n, const void *a, const MKL_INT
//                  *ia, const MKL_INT *ja, MKL_INT *perm, const MKL_INT *nrhs, MKL_INT *iparm, const
//                  MKL_INT *msglvl, void *b, void *x, MKL_INT *error);
    
    long long phase = 1;
    long long n=0;
    long long av,bv,xv;
    void *a= &av,*b = &bv, *x = &xv;
    long long ia,ja,perm,nrhs;
    long long Error = 0;
    
//    void pardiso_64( _MKL_DSS_HANDLE_t,       const long long int *, const long long int *, const long long int *,
//                    const long long int *, const long long int *, const void *,          const long long int *,
//                    const long long int *, long long int *, const long long int *, long long int *,
//                    const long long int *, void *,                void *,                long long int * );

    pardiso_64 (fHandle,  &fMax_num_factors, &fMatrix_num, &fMatrixType, &phase, &n, &a, &ia, &ja, &perm,
                &nrhs, &fParam[0], &fMessageLevel, &b, &x, &Error);
    
    
}

template<class TVar>
long long TPZPardisoControl<TVar>::MatrixType()
{
    DebugStop();
}