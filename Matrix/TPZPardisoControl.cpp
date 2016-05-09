//
//  TPZPardisoControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/5/16.
//
//

#include "TPZPardisoControl.h"
#ifdef USING_MKL

/// empty constructor (non symetric and LU decomposition
template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl() : fSystemType(ENonSymmetric),
fStructure(EStructureNonSymmetric), fProperty(EIndefinite), fPardisoControl(), fHandle(0),
fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0)
{
    fPardisoControl = new TPZManVector<long long,64>(64,0);
    fHandle = &fPardisoControl.operator->()->operator[](0);
    fMatrixType = MatrixType();
}


template<class TVar>
TPZPardisoControl<TVar>::TPZPardisoControl(MSystemType systemtype, MProperty prop) : fSystemType(systemtype),
        fStructure(EStructureSymmetric), fProperty(prop), fPardisoControl(), fHandle(0),
        fParam(64,0), fMax_num_factors(1), fMatrix_num(1), fMessageLevel(0), fError(0), fPermutation(), fMatrixType(0)
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
fError(copy.fError), fPermutation(copy.fPermutation), fMatrixType(copy.fMatrixType)
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
    
    long long phase = 1;
    long long n=0;
    long long av,bv,xv;
    void *a= &av,*b = &bv, *x = &xv;
    long long ia,ja,perm,nrhs;
    long long Error = 0;
    
    for (long i=0; i<64; i++) {
        if (fHandle[i]) {
            DebugStop();     
        }
    }
    double toto;
    fParam[27] = ::DataType(toto);
    fParam[34] = 1;
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


template class TPZPardisoControl<double>;
#endif
