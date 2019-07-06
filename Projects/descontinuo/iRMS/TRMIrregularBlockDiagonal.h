//
//  TRMIrregularBlockDiagonal.h
//  PZ
//
//  Created by Omar Durán on 5/26/16.
//  @brief A irregular block diagonal matrix, also called a diagonal block matrix, when is a square diagonal matrix in which the diagonal elements are square
//  matrices of any size (possibly even 1×1), and the off-diagonal elements are 0. A irregular block diagonal matrix is therefore a non square block matrix.
//

#ifndef TRMIrregularBlockDiagonal_h
#define TRMIrregularBlockDiagonal_h

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzlog.h"


template <class TVar>
class TRMIrregularBlockDiagonal : public TPZMatrix<TVar>
{
    
public:
    /** @brief Simple constructor */
    TRMIrregularBlockDiagonal ();
    /**
     * @brief Constructor with initialization parameters
     * @param blocksizes Size of blocks on Block Diagonal matrix
     * @param glob Global matrix which will be blocked
     */
    TRMIrregularBlockDiagonal (TPZVec< std::pair<int64_t, int64_t> > &blocksizes, const TPZFMatrix<TVar> &glob );
    /**
     * @brief Constructor with initialization parameters
     * @param blocksizes Size of blocks on Block Diagonal matrix
     */
    TRMIrregularBlockDiagonal (TPZVec< std::pair<int64_t, int64_t> > &blocksizes);
    /** @brief Copy constructor */
    TRMIrregularBlockDiagonal (const TRMIrregularBlockDiagonal & );

    /** @brief Copy constructor */
    TRMIrregularBlockDiagonal &operator=(const TRMIrregularBlockDiagonal & copy)
    {
        DebugStop();
        return *this;
    }

    CLONEDEF(TRMIrregularBlockDiagonal)
    /** @brief Simple destructor */
    ~TRMIrregularBlockDiagonal();
    
    int    Put(const int64_t row,const int64_t col,const TVar& value );
    const TVar &Get(const int64_t row,const int64_t col ) const;
    
    TVar &operator()(const int64_t row, const int64_t col);
    virtual TVar &s(const int64_t row, const int64_t col);
    
    /** @brief This method don't make verification if the element exist. It is fast than Put */
    int    PutVal(const int64_t row,const int64_t col,const TVar& value );
    /** @brief This method don't make verification if the element exist. It is fast than Get */
    const  TVar &GetVal(const int64_t row,const int64_t col ) const;
    
    /** @brief Computes z = alpha * opt(this)*x + beta * y */
    /** @note z and x cannot overlap in memory */
    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const ;
    
    int64_t Dim() const     { return this->Rows(); }
    
    /** @brief Zeroes all the elements of the matrix. */
    int Zero();
    
    /**
     * @brief Return the choosen block size
     * @param blockid - block index
     */
    std::pair<int64_t, int64_t> GetSizeofBlock(int64_t blockid) {return fBlockSize[blockid];}
    
    void Transpose(TPZMatrix<TVar> *const T) const;
    virtual int Decompose_LU();
    virtual int Decompose_LU(std::list<int64_t> &singular);
    
    /** @brief Makes the backward and forward substitutions whether the matrix was LU decomposed */
    virtual int Substitution( TPZFMatrix<TVar> * B ) const;
    
    /** @brief Updates the values of the matrix based on the values of the matrix */
    virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
    
    /** @brief This method checks the working of the class */
    static int main();
    
    /** Fill the matrix with random values (non singular matrix) */
    void AutoFill(int64_t dim, int64_t dimj, int symmetric);
    
private:
    
    /** @brief Clean data matrix. Zeroes number of columns and rows. */
    int Clear();
    
public:
    /**
     * @brief Initializes current matrix based on blocksize
     * @param blocksize Used to initialize current matrix
     */
    void Initialize(TPZVec< std::pair<int64_t, int64_t> > &blocksize);
    
    /**
     * @brief Partial initialization of the current matrix based on the number of blocks
     * @param nblock of blocks
     */
    void Initialize(int64_t nblock);
    
    /**
     * @brief Adds a block to current matrix
     * @param i Adds in ith position
     * @param block Block to be added
     */
    void AddBlock(int64_t i, TPZFMatrix<TVar> &block);
    /**
     * @brief Sets a block in the current matrix
     * @param i Adds in ith position
     * @param block Block to be added
     */
    void SetBlock(int64_t i, TPZFMatrix<TVar> &block);
    
    /**
     * @brief Sets block dimensions in the current matrix
     * @param i Adds in ith position
     * @param block_size seize to be fixed
     */
    void SetBlockSize(int64_t i, std::pair<int64_t, int64_t> &block_size);
    
    /**
     * @brief Gets a block from current matrix
     * @param i Returns teh ith block
     * @param block Contains returned block
     */
    void GetBlock(int64_t i, TPZFMatrix<TVar> &block);
    
    /**
     @brief Builds a block from matrix
     @param matrix Matrix to build from
     */
    void BuildFromMatrix(TPZMatrix<TVar> &matrix);
    /**
     * @brief Prints current matrix data
     * @param message Message to be printed
     * @param out Output device
     * @param format Output format to print
     */
    virtual void Print(const char *message, std::ostream &out = std::cout, const MatrixOutputFormat format =EFormatted) const;
    
    int64_t NumberofBlocks() {return fBlockSize.NElements();}
    
protected:
    /** @brief Stores matrix data */
    TPZVec<TVar> fStorage;
    
    /** @brief Stores blocks data */
    TPZVec<int64_t> fBlockPos;
    
    /** @brief Stores block isizes and jside data */
    TPZVec< std::pair<int64_t, int64_t> > fBlockSize;
};

template<class TVar>
inline TVar &TRMIrregularBlockDiagonal<TVar>::s(const int64_t row, const int64_t col) {
    // verificando se o elemento a inserir esta dentro da matriz
    return this->operator()(row,col);
}


#endif /* TRMIrregularBlockDiagonal_h */
