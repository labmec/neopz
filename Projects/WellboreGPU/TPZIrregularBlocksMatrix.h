//
// Created by natalia on 14/05/19.
//
#ifndef INTPOINTSFEM_TPZIRREGULARBLOCKSMATRIX_H
#define INTPOINTSFEM_TPZIRREGULARBLOCKSMATRIX_H
#include "pzmatrix.h"

#ifdef USING_CUDA
#include "TPZVecGPU.h"
#include "TPZCudaCalls.h"
#endif

class TPZIrregularBlocksMatrix : public TPZMatrix<REAL> {

public:
    
    /** @brief Irregular blocks information */
    struct IrregularBlocks {
        int64_t fNumBlocks; //number of blocks
        TPZVec<REAL> fStorage; // blocks values
        TPZVec<int> fRowSizes; // blocks row sizes
        TPZVec<int> fColSizes; // blocks columns sizes
        TPZVec<int> fMatrixPosition; // blocks start position in fStorage vector
        TPZVec<int> fRowFirstIndex; // blocks first row index
        TPZVec<int> fColFirstIndex; // blocks first column index

        TPZVec<int> fRowPtr; // vector of the start of every row and the end of the last row plus one (this is for CSR format)
        TPZVec<int> fColInd; // vector of column indices for each non-zero element of the matrix (this is for CSR format)
    };


#ifdef USING_CUDA
    struct IrregularBlocksDev {
    TPZVecGPU<REAL> dStorage; // blocks values
    TPZVecGPU<int> dRowSizes; // blocks row sizes
    TPZVecGPU<int> dColSizes; // blocks columns sizes
    TPZVecGPU<int> dMatrixPosition; // blocks start position in fStorage vector
    TPZVecGPU<int> dRowFirstIndex; // blocks first row index
    TPZVecGPU<int> dColFirstIndex; // blocks first column index

    TPZVecGPU<int> dRowPtr; // vector of the start of every row and the end of the last row plus one (this is for CSR format)
    TPZVecGPU<int> dColInd; // vector of column indices for each non-zero element of the matrix (this is for CSR format)
    };
#endif

    TPZIrregularBlocksMatrix();

    TPZIrregularBlocksMatrix(const int64_t rows, const int64_t cols);

    ~TPZIrregularBlocksMatrix();

    virtual TPZMatrix<REAL> * Clone() const {
        return new TPZIrregularBlocksMatrix(*this);
    }

    TPZIrregularBlocksMatrix(const TPZIrregularBlocksMatrix &copy);

    TPZIrregularBlocksMatrix &operator=(const TPZIrregularBlocksMatrix &copy);

    void MultiplyVector(REAL *A, REAL *res, int opt);

    void CSRVectors();

    void SetBlocks(struct IrregularBlocks & blocks) {
        fBlocksInfo = blocks;
#ifdef USING_SPARSE
        CSRVectors();
#endif
    }

    struct IrregularBlocks & Blocks() {
        return fBlocksInfo;
    }

#ifdef USING_CUDA    
    struct IrregularBlocksDev & BlocksDev() {
        return dBlocksInfo;
    }    

    void TransferDataToGPU();
#endif

private:
    struct IrregularBlocks fBlocksInfo;

#ifdef USING_CUDA
    struct IrregularBlocksDev dBlocksInfo;
    TPZCudaCalls *fCudaCalls;
#endif
};


#endif //INTPOINTSFEM_TPZIRREGULARBLOCKSMATRIX_H
