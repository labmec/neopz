//
// Created by natalia on 14/05/19.
//
#include "TPZIrregularBlocksMatrix.h"

#ifdef USING_MKL
#include <mkl.h>
#endif

TPZIrregularBlocksMatrix::TPZIrregularBlocksMatrix() : TPZMatrix<REAL>(), fBlocksInfo() {
    this->Resize(0,0);
    fBlocksInfo.fNumBlocks = -1;
    fBlocksInfo.fStorage.resize(0);
    fBlocksInfo.fRowSizes.resize(0);
    fBlocksInfo.fColSizes.resize(0);
    fBlocksInfo.fMatrixPosition.resize(0);
    fBlocksInfo.fRowFirstIndex.resize(0);
    fBlocksInfo.fColFirstIndex.resize(0);
    fBlocksInfo.fRowPtr.resize(0);
    fBlocksInfo.fColInd.resize(0);
}

TPZIrregularBlocksMatrix::TPZIrregularBlocksMatrix(const int64_t rows,const int64_t cols) : TPZMatrix<REAL>(rows,cols), fBlocksInfo() {
    fBlocksInfo.fNumBlocks = -1;
    fBlocksInfo.fStorage.resize(0);
    fBlocksInfo.fRowSizes.resize(0);
    fBlocksInfo.fColSizes.resize(0);
    fBlocksInfo.fMatrixPosition.resize(0);
    fBlocksInfo.fRowFirstIndex.resize(0);
    fBlocksInfo.fColFirstIndex.resize(0);
    fBlocksInfo.fRowPtr.resize(0);
    fBlocksInfo.fColInd.resize(0);
}

TPZIrregularBlocksMatrix::~TPZIrregularBlocksMatrix() {}

TPZIrregularBlocksMatrix::TPZIrregularBlocksMatrix(const TPZIrregularBlocksMatrix &copy) {
    TPZMatrix<REAL>::operator=(copy);
    fBlocksInfo = copy.fBlocksInfo;
}

TPZIrregularBlocksMatrix &TPZIrregularBlocksMatrix::operator=(const TPZIrregularBlocksMatrix &copy) {
    if(&copy == this){
        return *this;
    }
    TPZMatrix<REAL>::operator=(copy);
    fBlocksInfo = copy.fBlocksInfo;

    return *this;
}

void TPZIrregularBlocksMatrix::Multiply(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &res, int opt) {
    char trans;
    char matdescra[] = {'G',' ',' ','C'};
    REAL alpha, beta;

    if (opt == false) {
        trans = 'N';
        alpha = 1.;
    }
    else if (opt == true) {
        trans = 'T';
        alpha = -1.;
    }

    const int m = fRow;
    const int n = 1;
    const int k = fCol;
    beta = 0.;
//    mkl_dcsrmm (&trans, &m, &n, &k, &alpha, matdescra, &fStorage[0], &fColInd[0], &fRowPtr[0], &fRowPtr[1], &A(0,0), &n, &beta, &res(0,0), &n);
    mkl_dcsrmv(&trans, &m, &k, &alpha, matdescra , &fBlocksInfo.fStorage[0], &fBlocksInfo.fColInd[0], &fBlocksInfo.fRowPtr[0], &fBlocksInfo.fRowPtr[1], A, &beta, &res(0,0));
}