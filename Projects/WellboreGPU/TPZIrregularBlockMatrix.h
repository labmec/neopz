//
// Created by natalia on 14/05/19.
//
#include "pzcmesh.h"
#include "TPZPlasticStepPV.h"
#include "TPZYCMohrCoulombPV.h"
#include "TPZElastoPlasticMem.h"
#include "TPZMatElastoPlastic2D.h"

#ifndef INTPOINTSFEM_TPZIRREGULARBLOCKMATRIX_H
#define INTPOINTSFEM_TPZIRREGULARBLOCKMATRIX_H


class TPZIrregularBlockMatrix {
public:
    /** @brief Default constructor */
    TPZIrregularBlockMatrix();

    /** @brief Default destructor */
    ~TPZIrregularBlockMatrix();

    /** @brief Creates a irregular block matrix with copy constructor
     * @param copy : original irregular block matrix
     */
    TPZIrregularBlockMatrix(const TPZIrregularBlockMatrix &copy);

    /** @brief operator= */
    TPZIrregularBlockMatrix &operator=(const TPZIrregularBlockMatrix &copy);

    /** @brief Performs the following operation: res = alpha * BMatrix * A + beta * res
     * @param A : matrix that will be multipied
     * @param res : result of the multiplication
     * @param alpha : scalar parameter
     * @param beta : scalar parameter
     * @param opt : indicates if transpose or not
     */
    void Multiply(REAL *A, REAL *res, int opt);

    /** @brief Access methods */
    int64_t NumBlocks() {
        return fNumBlocks;
    }

    void SetNumBlocks(int64_t nb) {
        fNumBlocks = nb;
    }

    TPZVec<REAL> Storage() {
        return fStorage;
    }

    void SetStorage(TPZVec<REAL> storage) {
        fStorage = storage;
    }

    TPZVec<int> RowSizes() {
        return fRowSizes;
    }

    void SetRowSizes(TPZVec<int> rowsizes) {
        fRowSizes = rowsizes;
    }

    TPZVec<int> ColSizes() {
        return fColSizes;
    }

    void SetColSizes(TPZVec<int> colsizes) {
        fColSizes = colsizes;
    }

    TPZVec<int> MatrixPosition() {
        return fMatrixPosition;
    }

    void SetMatrixPosition(TPZVec<int> matrixpos) {
        fMatrixPosition = matrixpos;
    }

    TPZVec<int> RowFirstIndex() {
        return fRowFirstIndex;
    }

    void SetRowFirstIndex(TPZVec<int> rowfirstindex) {
        fRowFirstIndex = rowfirstindex;
    }

    TPZVec<int> ColFirstIndex() {
        return fColFirstIndex;
    }

    void SetColFirstIndex(TPZVec<int> colfirstindex) {
        fColFirstIndex = colfirstindex;
    }

    int64_t Rows() {
        return fRow;
    }

    void SetRows(int64_t rows) {
        fRow = rows;
    }

    int64_t Cols() {
        return fCol;
    }

    void SetCols(int64_t cols) {
        fCol = cols;
    }

    TPZVec<int> RowPtr() {
        return fRowPtr;
    }

    void SetRowPtr(TPZVec<int> rowptr) {
        fRowPtr = rowptr;
    }

    TPZVec<int> ColInd() {
        return fColInd;
    }
    void SetColInd(TPZVec<int> colind) {
        fColInd = colind;
    }

protected:
    /** @brief Number of blocks */
    int64_t fNumBlocks;

    /** @brief Vector of all matrices values */
    TPZVec<REAL> fStorage;

    /** @brief Vector of number of rows of each matrix */
    TPZVec<int> fRowSizes;

    /** @brief Vector of number of columns of each matrix */
    TPZVec<int> fColSizes;

    /** @brief Vector of matrix position in fStorage */
    TPZVec<int> fMatrixPosition;

    /** @brief Vector of first row index of each matrix */
    TPZVec<int> fRowFirstIndex;

    /** @brief Vector of first column index of each matrix */
    TPZVec<int> fColFirstIndex;

    /** @brief Number of rows of the irregular block matrix */
    int64_t fRow;

    /** @brief Number of columns of the irregular block matrix */
    int64_t fCol;

    /** @brief Vector of the start of every row and the end of the last row plus one (this is for CSR format) */
    TPZVec<int> fRowPtr;

    /** @brief Vector of column indices for each non-zero element of the matrix (this is for CSR format)*/
    TPZVec<int> fColInd;

};


#endif //INTPOINTSFEM_TPZIRREGULARBLOCKMATRIX_H