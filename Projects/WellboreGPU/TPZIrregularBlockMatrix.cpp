//
// Created by natalia on 14/05/19.
//
#include "pzintel.h"
#include "TPZIrregularBlockMatrix.h"

#ifdef USING_MKL
#include <mkl.h>
#endif

TPZIrregularBlockMatrix::TPZIrregularBlockMatrix() : fNumBlocks(-1), fStorage(0), fRowSizes(0), fColSizes(0),
                                                     fMatrixPosition(0), fRowFirstIndex(0), fColFirstIndex(0),
                                                     fRow(-1), fCol(-1), fRowPtr(0), fColInd(0){}

TPZIrregularBlockMatrix::~TPZIrregularBlockMatrix() {
}

TPZIrregularBlockMatrix::TPZIrregularBlockMatrix(const TPZIrregularBlockMatrix &copy) {
    fNumBlocks = copy.fNumBlocks;
    fStorage = copy.fStorage;
    fRowSizes = copy.fRowSizes;
    fColSizes = copy.fColSizes;
    fMatrixPosition = copy.fMatrixPosition;
    fRowFirstIndex = copy.fRowFirstIndex;
    fColFirstIndex = copy.fColFirstIndex;
    fRow = copy.fRow;
    fCol = copy.fCol;
    fRowPtr = copy.fRowPtr;
    fColInd = copy.fColInd;
}

TPZIrregularBlockMatrix &TPZIrregularBlockMatrix::operator=(const TPZIrregularBlockMatrix &copy) {
    if(&copy == this){
        return *this;
    }

    fNumBlocks = copy.fNumBlocks;
    fStorage = copy.fStorage;
    fRowSizes = copy.fRowSizes;
    fColSizes = copy.fColSizes;
    fMatrixPosition = copy.fMatrixPosition;
    fRowFirstIndex = copy.fRowFirstIndex;
    fColFirstIndex = copy.fColFirstIndex;
    fRow = copy.fRow;
    fCol = copy.fCol;
    fRowPtr = copy.fRowPtr;
    fColInd = copy.fColInd;

    return *this;
}

void TPZIrregularBlockMatrix::Multiply(REAL *A, REAL *res, int opt) {
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
    mkl_dcsrmv(&trans, &m, &k, &alpha, matdescra , &fStorage[0], &fColInd[0], &fRowPtr[0], &fRowPtr[1], A, &beta, res);
}