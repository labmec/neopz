//
// Created by natalia on 14/05/19.
//
#ifndef INTPOINTSFEM_TPZIRREGULARBLOCKSMATRIX_H
#define INTPOINTSFEM_TPZIRREGULARBLOCKSMATRIX_H
#include "pzmatrix.h"

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

    /** @brief Default constructor */
    TPZIrregularBlocksMatrix();

    /**
     @brief Constructor with initialization parameters
     @param rows number of rows
     @param columns Number of columns
     */
    TPZIrregularBlocksMatrix(const int64_t rows, const int64_t cols);

    /** @brief Default destructor */
    ~TPZIrregularBlocksMatrix();

    /** @brief Clone */
    virtual TPZMatrix<REAL> * Clone() const {
        return new TPZIrregularBlocksMatrix(*this);
    }

    /** @brief Creates a TPZIrregularBlockMatrix with copy constructor
     * @param copy : original TPZIrregularBlockMatrix
     */
    TPZIrregularBlocksMatrix(const TPZIrregularBlocksMatrix &copy);

    /** @brief operator= */
    TPZIrregularBlocksMatrix &operator=(const TPZIrregularBlocksMatrix &copy);

    /** @brief Performs the following operation: res = alpha * this * A + beta * res
     * @param A : matrix that will be multipied
     * @param res : result of the multiplication
     * @param alpha : scalar parameter
     * @param beta : scalar parameter
     * @param opt : indicates if transpose or not
     */
    void Multiply(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &res, int opt);

    /** @brief Set method */
    void SetBlocks(struct IrregularBlocks blocks) {
        fBlocksInfo = blocks;
    }

    /** @brief Access method */
    struct IrregularBlocks Blocks() {
        return fBlocksInfo;
    }

protected:
    struct IrregularBlocks fBlocksInfo;
};


#endif //INTPOINTSFEM_TPZIRREGULARBLOCKSMATRIX_H