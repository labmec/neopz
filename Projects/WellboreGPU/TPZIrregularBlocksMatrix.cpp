//
// Created by natalia on 14/05/19.
//
#include "TPZIrregularBlocksMatrix.h"

#include "MatMul.h"

#ifdef USING_MKL
#include <mkl.h>
#endif

TPZIrregularBlocksMatrix::TPZIrregularBlocksMatrix() : TPZMatrix<REAL>(), fBlocksInfo() {
#ifdef USING_CUDA
    dBlocksInfo.dStorage.resize(0);
    #ifdef USING_SPARSE
    dBlocksInfo.dRowPtr.resize(0);
    dBlocksInfo.dColInd.resize(0);
    dBlocksInfo.dRowRowPtr.resize(0);
    dBlocksInfo.dRowRowInd.resize(0);
    dBlocksInfo.dColColPtr.resize(0);
    #else
    dBlocksInfo.dRowSizes.resize(0);
    dBlocksInfo.dColSizes.resize(0);
    dBlocksInfo.dRowFirstIndex.resize(0);
    dBlocksInfo.dColFirstIndex.resize(0);
    dBlocksInfo.dMatrixPosition.resize(0);
    dBlocksInfo.dRowRowPosition.resize(0);
    dBlocksInfo.dColColPosition.resize(0);
    #endif

    fCudaCalls = new TPZCudaCalls();
#endif
    this->Resize(0,0);
    fBlocksInfo.fNumBlocks = -1;
    fBlocksInfo.fStorage.resize(0);
    fBlocksInfo.fRowSizes.resize(0);
    fBlocksInfo.fColSizes.resize(0);
    fBlocksInfo.fMatrixPosition.resize(0);
    fBlocksInfo.fRowFirstIndex.resize(0);
    fBlocksInfo.fColFirstIndex.resize(0);
    fBlocksInfo.fRowRowPosition.resize(0);
    fBlocksInfo.fColColPosition.resize(0);

    fBlocksInfo.fColColPtr.resize(0);
    fBlocksInfo.fRowPtr.resize(0);
    fBlocksInfo.fColInd.resize(0);
    fBlocksInfo.fRowRowPtr.resize(0);
    fBlocksInfo.fRowRowInd.resize(0);
}

TPZIrregularBlocksMatrix::TPZIrregularBlocksMatrix(const int64_t rows,const int64_t cols) : TPZMatrix<REAL>(rows,cols), fBlocksInfo() {
#ifdef USING_CUDA
    dBlocksInfo.dStorage.resize(0);
    #ifdef USING_SPARSE
    dBlocksInfo.dColColPtr.resize(0);
    dBlocksInfo.dRowPtr.resize(0);
    dBlocksInfo.dColInd.resize(0);
    dBlocksInfo.dRowRowPtr.resize(0);
    dBlocksInfo.dRowRowInd.resize(0);
    #else
    dBlocksInfo.dRowSizes.resize(0);
    dBlocksInfo.dColSizes.resize(0);
    dBlocksInfo.dRowFirstIndex.resize(0);
    dBlocksInfo.dColFirstIndex.resize(0);
    dBlocksInfo.dMatrixPosition.resize(0);
    dBlocksInfo.dRowRowPosition.resize(0);
    dBlocksInfo.dColColPosition.resize(0);
    #endif

    fCudaCalls = new TPZCudaCalls();
#endif
    fBlocksInfo.fNumBlocks = -1;
    fBlocksInfo.fStorage.resize(0);
    fBlocksInfo.fRowSizes.resize(0);
    fBlocksInfo.fColSizes.resize(0);
    fBlocksInfo.fMatrixPosition.resize(0);
    fBlocksInfo.fRowFirstIndex.resize(0);
    fBlocksInfo.fColFirstIndex.resize(0);
    fBlocksInfo.fRowRowPosition.resize(0);
    fBlocksInfo.fColColPosition.resize(0);

    fBlocksInfo.fColColPtr.resize(0);
    fBlocksInfo.fRowPtr.resize(0);
    fBlocksInfo.fColInd.resize(0);
    fBlocksInfo.fRowRowPtr.resize(0);
    fBlocksInfo.fRowRowInd.resize(0);
    fBlocksInfo.fColColPtr.resize(0);
}

TPZIrregularBlocksMatrix::~TPZIrregularBlocksMatrix() {
}

TPZIrregularBlocksMatrix::TPZIrregularBlocksMatrix(const TPZIrregularBlocksMatrix &copy) {
    TPZMatrix<REAL>::operator=(copy);
    fBlocksInfo = copy.fBlocksInfo;

#ifdef USING_CUDA
    dBlocksInfo = copy.dBlocksInfo;
    fCudaCalls = copy.fCudaCalls;
#endif
}

TPZIrregularBlocksMatrix &TPZIrregularBlocksMatrix::operator=(const TPZIrregularBlocksMatrix &copy) {
    if(&copy == this){
        return *this;
    }
    TPZMatrix<REAL>::operator=(copy);
    fBlocksInfo = copy.fBlocksInfo;

#ifdef USING_CUDA
    dBlocksInfo = copy.dBlocksInfo;
    fCudaCalls = copy.fCudaCalls;
#endif

    return *this;
}

void TPZIrregularBlocksMatrix::MultiplyVector(REAL *A, REAL *res, int opt) {
    int nblocks = fBlocksInfo.fNumBlocks;

    TPZVec<int> one(nblocks);
    one.Fill(1);

#ifdef USING_CUDA
    int rows = this->Rows();
    int cols = this->Cols();

    #ifdef USING_SPARSE
    if(opt == 0) {
        fCudaCalls->SpMV(0, rows, cols, dBlocksInfo.dStorage.getSize(), 1., dBlocksInfo.dStorage.getData(), dBlocksInfo.dRowPtr.getData(), dBlocksInfo.dColInd.getData(), A, res); 
    } else {
        fCudaCalls->SpMV(1, rows, cols, dBlocksInfo.dStorage.getSize(), -1., dBlocksInfo.dStorage.getData(), dBlocksInfo.dRowPtr.getData(), dBlocksInfo.dColInd.getData(), A, res); 
    }
    #else
    TPZVecGPU<int> dOne(nblocks);
    dOne.set(&one[0], nblocks);

    if(opt == 0) {
        fCudaCalls->Multiply(opt, dBlocksInfo.dRowSizes.getData(), dOne.getData(), dBlocksInfo.dColSizes.getData(), dBlocksInfo.dStorage.getData(), dBlocksInfo.dMatrixPosition.getData(), A, dBlocksInfo.dColFirstIndex.getData(), res, dBlocksInfo.dRowFirstIndex.getData(), 1., nblocks); 
    } else {
        fCudaCalls->Multiply(opt, dBlocksInfo.dColSizes.getData(), dOne.getData(), dBlocksInfo.dRowSizes.getData(), dBlocksInfo.dStorage.getData(), dBlocksInfo.dMatrixPosition.getData(), A, dBlocksInfo.dRowFirstIndex.getData(), res, dBlocksInfo.dColFirstIndex.getData(), -1., nblocks); 
    }
    #endif 
 #else
    if(opt == 0) {
        MatrixMultiplication(opt, &fBlocksInfo.fRowSizes[0], &one[0], &fBlocksInfo.fColSizes[0], &fBlocksInfo.fStorage[0], &fBlocksInfo.fMatrixPosition[0], A, &fBlocksInfo.fColFirstIndex[0], res, &fBlocksInfo.fRowFirstIndex[0], 1., nblocks);
    } else {
        MatrixMultiplication(opt, &fBlocksInfo.fColSizes[0], &one[0], &fBlocksInfo.fRowSizes[0], &fBlocksInfo.fStorage[0], &fBlocksInfo.fMatrixPosition[0], A, &fBlocksInfo.fRowFirstIndex[0], res, &fBlocksInfo.fColFirstIndex[0], -1., nblocks);
    }

#endif
}

void TPZIrregularBlocksMatrix::KMatrix(REAL *A, REAL *res) {
    int nblocks = fBlocksInfo.fNumBlocks;

    int rows = this->Rows();
    int cols = this->Cols();

#ifdef USING_CUDA
    TPZVecGPU<REAL> aux(fBlocksInfo.fMatrixPosition[nblocks]);

    #ifdef USING_SPARSE
    fCudaCalls->SpMSpM(0, rows, rows, cols, fBlocksInfo.fRowRowPosition[nblocks], A, dBlocksInfo.dRowRowPtr.getData(), dBlocksInfo.dRowRowInd.getData(), 
        dBlocksInfo.dStorage.getSize(), dBlocksInfo.dStorage.getData(), dBlocksInfo.dRowPtr.getData(), dBlocksInfo.dColInd.getData(), 
        dBlocksInfo.dStorage.getSize(), aux.getData(), dBlocksInfo.dRowPtr.getData()); 

    fCudaCalls->SpMSpM(1, cols, cols, rows, dBlocksInfo.dStorage.getSize(), dBlocksInfo.dStorage.getData(), dBlocksInfo.dRowPtr.getData(), dBlocksInfo.dColInd.getData(), 
        aux.getSize(), aux.getData(), dBlocksInfo.dRowPtr.getData(), dBlocksInfo.dColInd.getData(), fBlocksInfo.fColColPosition[nblocks], res, dBlocksInfo.dColColPtr.getData()); 

    #else
    fCudaCalls->Multiply(0, dBlocksInfo.dRowSizes.getData(), dBlocksInfo.dColSizes.getData(), dBlocksInfo.dRowSizes.getData(), A, dBlocksInfo.dRowRowPosition.getData(), dBlocksInfo.dStorage.getData(), dBlocksInfo.dMatrixPosition.getData(), aux.getData(), dBlocksInfo.dMatrixPosition.getData(), 1., nblocks);
    fCudaCalls->Multiply(1, dBlocksInfo.dColSizes.getData(), dBlocksInfo.dColSizes.getData(), dBlocksInfo.dRowSizes.getData(), dBlocksInfo.dStorage.getData(), dBlocksInfo.dMatrixPosition.getData(), aux.getData(), dBlocksInfo.dMatrixPosition.getData(), res, dBlocksInfo.dColColPosition.getData(), 1., nblocks);
    #endif 
#else
    TPZVec<REAL> aux(fBlocksInfo.fMatrixPosition[nblocks]);

    MatrixMultiplication(0, &fBlocksInfo.fRowSizes[0], &fBlocksInfo.fColSizes[0], &fBlocksInfo.fRowSizes[0], A, &fBlocksInfo.fRowRowPosition[0], &fBlocksInfo.fStorage[0], &fBlocksInfo.fMatrixPosition[0], &aux[0], &fBlocksInfo.fMatrixPosition[0], 1., nblocks); /// aux -> outvar /// A ->Dep, Storage -> eps:  C*eps or Dep*eps
    MatrixMultiplication(1, &fBlocksInfo.fColSizes[0], &fBlocksInfo.fColSizes[0], &fBlocksInfo.fRowSizes[0], &fBlocksInfo.fStorage[0], &fBlocksInfo.fMatrixPosition[0], &aux[0], &fBlocksInfo.fMatrixPosition[0], res, &fBlocksInfo.fColColPosition[0], 1., nblocks); /// res -> outvar /// eps_t*C*eps
#endif
}



#ifdef USING_CUDA
void TPZIrregularBlocksMatrix::TransferDataToGPU() {
    dBlocksInfo.dStorage.resize(fBlocksInfo.fStorage.size());
    dBlocksInfo.dStorage.set(&fBlocksInfo.fStorage[0], fBlocksInfo.fStorage.size());

    #ifdef USING_SPARSE
    dBlocksInfo.dRowPtr.resize(fBlocksInfo.fRowPtr.size());
    dBlocksInfo.dRowPtr.set(&fBlocksInfo.fRowPtr[0], fBlocksInfo.fRowPtr.size());

    dBlocksInfo.dColInd.resize(fBlocksInfo.fColInd.size());
    dBlocksInfo.dColInd.set(&fBlocksInfo.fColInd[0], fBlocksInfo.fColInd.size());

    dBlocksInfo.dRowRowPtr.resize(fBlocksInfo.fRowRowPtr.size());
    dBlocksInfo.dRowRowPtr.set(&fBlocksInfo.fRowRowPtr[0], fBlocksInfo.fRowRowPtr.size());

    dBlocksInfo.dRowRowInd.resize(fBlocksInfo.fRowRowInd.size());
    dBlocksInfo.dRowRowInd.set(&fBlocksInfo.fRowRowInd[0], fBlocksInfo.fRowRowInd.size());

    dBlocksInfo.dColColPtr.resize(fBlocksInfo.fColColPtr.size());
    dBlocksInfo.dColColPtr.set(&fBlocksInfo.fColColPtr[0], fBlocksInfo.fColColPtr.size());
    #else
    dBlocksInfo.dRowSizes.resize(fBlocksInfo.fRowSizes.size());
    dBlocksInfo.dRowSizes.set(&fBlocksInfo.fRowSizes[0], fBlocksInfo.fRowSizes.size());

    dBlocksInfo.dColSizes.resize(fBlocksInfo.fColSizes.size());
    dBlocksInfo.dColSizes.set(&fBlocksInfo.fColSizes[0], fBlocksInfo.fColSizes.size());

    dBlocksInfo.dMatrixPosition.resize(fBlocksInfo.fMatrixPosition.size());
    dBlocksInfo.dMatrixPosition.set(&fBlocksInfo.fMatrixPosition[0], fBlocksInfo.fMatrixPosition.size());

    dBlocksInfo.dRowFirstIndex.resize(fBlocksInfo.fRowFirstIndex.size());
    dBlocksInfo.dRowFirstIndex.set(&fBlocksInfo.fRowFirstIndex[0], fBlocksInfo.fRowFirstIndex.size());

    dBlocksInfo.dColFirstIndex.resize(fBlocksInfo.fColFirstIndex.size());
    dBlocksInfo.dColFirstIndex.set(&fBlocksInfo.fColFirstIndex[0], fBlocksInfo.fColFirstIndex.size());

    dBlocksInfo.dRowRowPosition.resize(fBlocksInfo.fRowRowPosition.size());
    dBlocksInfo.dRowRowPosition.set(&fBlocksInfo.fRowRowPosition[0], fBlocksInfo.fRowRowPosition.size());

    dBlocksInfo.dColColPosition.resize(fBlocksInfo.fColColPosition.size());
    dBlocksInfo.dColColPosition.set(&fBlocksInfo.fColColPosition[0], fBlocksInfo.fColColPosition.size());
    #endif
}
#endif

void TPZIrregularBlocksMatrix::CSRVectors() {
    int64_t nblocks = fBlocksInfo.fNumBlocks;
    int64_t rows = this->Rows();
    int64_t cols = this->Cols();

    fBlocksInfo.fRowPtr.resize(rows + 1);
    fBlocksInfo.fColInd.resize(fBlocksInfo.fMatrixPosition[nblocks]);

    fBlocksInfo.fRowRowPtr.resize(rows + 1);
    fBlocksInfo.fRowRowInd.resize(fBlocksInfo.fRowRowPosition[nblocks]);

    fBlocksInfo.fColColPtr.resize(cols + 1);

    for (int iel = 0; iel < nblocks; ++iel) {
        for (int irow = 0; irow < fBlocksInfo.fRowSizes[iel]; ++irow) {
            fBlocksInfo.fRowPtr[irow + fBlocksInfo.fRowFirstIndex[iel]] = fBlocksInfo.fMatrixPosition[iel] + irow * fBlocksInfo.fColSizes[iel];

            for (int icol = 0; icol < fBlocksInfo.fColSizes[iel]; ++icol) {
                fBlocksInfo.fColInd[icol + fBlocksInfo.fMatrixPosition[iel] + irow * fBlocksInfo.fColSizes[iel]] = icol + fBlocksInfo.fColFirstIndex[iel];
            }
        }
        for (int irow = 0; irow < fBlocksInfo.fRowSizes[iel]; ++irow) {
            fBlocksInfo.fRowRowPtr[irow + fBlocksInfo.fRowFirstIndex[iel]] = fBlocksInfo.fRowRowPosition[iel] + irow * fBlocksInfo.fRowSizes[iel];

            for (int icol = 0; icol < fBlocksInfo.fRowSizes[iel]; ++icol) {
                fBlocksInfo.fRowRowInd[icol + fBlocksInfo.fRowRowPosition[iel] + irow * fBlocksInfo.fRowSizes[iel]] = icol + fBlocksInfo.fRowFirstIndex[iel];
            }
        }

        for (int irow = 0; irow < fBlocksInfo.fColSizes[iel]; ++irow) {
            fBlocksInfo.fColColPtr[irow + fBlocksInfo.fColFirstIndex[iel]] = fBlocksInfo.fColColPosition[iel] + irow * fBlocksInfo.fColSizes[iel];
        }
    }
    fBlocksInfo.fRowPtr[rows] = fBlocksInfo.fMatrixPosition[nblocks];
    fBlocksInfo.fRowRowPtr[rows] = fBlocksInfo.fRowRowPosition[nblocks];
    fBlocksInfo.fColColPtr[cols] = fBlocksInfo.fColColPosition[nblocks];
}
