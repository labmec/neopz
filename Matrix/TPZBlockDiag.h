//
//  TPZBlockDiag.h
//  PZ
//
//  Created by Omar on 5/1/16.
//
//

#ifndef TPZBlockDiag_h
#define TPZBlockDiag_h


/**
 * @file
 * @brief Contains TPZBlockDiagonal class which defines block diagonal matrices.
 */


#include <stdio.h>
#include <iostream>
#include "pzmatrix.h"
#include "pzvec.h"

/**
 * @brief Defines block diagonal matrices. \ref matrix "Matrix"
 * @ingroup matrix
 * @author Philippe Devloo.
 * @since 01/2001
 */
template <class TVar>
class TPZBlockDiag : public TPZMatrix<TVar>
{
    
public:
    /** @brief default constructor */
    TPZBlockDiag();

    /**
     * @brief Constructor with initialization parameters
     * @param blocksizes Size of blocks on Block Diagonal matrix
     * @param glob Global matrix which will be blocked
     */
    TPZBlockDiag (const TPZVec< std::pair<long, long> > &blocksizes, const TPZFMatrix<TVar> &glob );

    /**
     * @brief Constructor with initialization parameters
     * @param blocksizes Size of blocks on Block Diagonal matrix
     */
    TPZBlockDiag(const TPZVec<std::pair<long, long> > &blocksizes);

    /** @brief Copy constructor */
    TPZBlockDiag (const TPZBlockDiag & );
    
    CLONEDEF(TPZBlockDiag)
    /** @brief Simple destructor */
    ~TPZBlockDiag();
    
    int    Put(const long row,const long col,const TVar& value );
    const TVar &Get(const long row,const long col ) const;
    
    TVar &operator()(const long row, const long col);
    virtual TVar &s(const long row, const long col);
    
    /** @brief This method don't make verification if the element exist. It is fast than Put */
    int    PutVal(const long row,const long col,const TVar& value );
    /** @brief This method don't make verification if the element exist. It is fast than Get */
    const  TVar &GetVal(const long row,const long col ) const;
    
    /** @brief Computes z = alpha * opt(this)*x + beta * y */
    /** @note z and x cannot overlap in memory */
    void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 ) const ;
    
    long Dim_i() const     { return this->Rows(); }
    
    long Dim_j() const     { return this->Cols(); }
    
    /** @brief Zeroes all the elements of the matrix. */
    int Zero();
    
    /**
     * @brief Return the choosen block size
     * @param blockid - block index
     */
    int GetSizeofBlock(long blockid) {return fBlockSize[blockid];}
    
    void Transpose(TPZMatrix<TVar> *const T) const;
    virtual int Decompose_LU();
    virtual int Decompose_LU(std::list<long> &singular);
    
    /** @brief Makes the backward and forward substitutions whether the matrix was LU decomposed */
    virtual int Substitution( TPZFMatrix<TVar> * B ) const;
    
    /** @brief Updates the values of the matrix based on the values of the matrix */
    virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);
    
    /** @brief This method checks the working of the class */
    static int main();
    
    /** Fill the matrix with random values (non singular matrix) */
    void AutoFill();
    
private:
    
    /** @brief Clean data matrix. Zeroes number of columns and rows. */
    int Clear();
    
public:
    /**
     * @brief Initializes current matrix based on blocksize
     * @param blocksize Used to initialize current matrix
     */
    void Initialize(const TPZVec< std::pair<long, long> > &blocksize);
    /**
     * @brief Adds a block to current matrix
     * @param i Adds in ith position
     * @param block Block to be added
     */
    void AddBlock(long i, TPZFMatrix<TVar> &block);
    /**
     * @brief Sets a block in the current matrix
     * @param i Adds in ith position
     * @param block Block to be added
     */
    void SetBlock(long i, TPZFMatrix<TVar> &block);
    
    /**
     * @brief Gets a block from current matrix
     * @param i Returns teh ith block
     * @param block Contains returned block
     */
    void GetBlock(long i, TPZFMatrix<TVar> &block);
    
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
    
    long NumberofBlocks() {return fBlockSize.NElements();}
    
protected:
    /** @brief Stores matrix data */
    TPZVec<TVar> fStorage;
    /** @brief Stores blocks data */
    TPZVec<long> fBlockPos;
    /** @brief Stores block sizes data */
    TPZVec< std::pair<long, long> > fBlockSize;
};

template<class TVar>
inline TVar &TPZBlockDiag<TVar>::s(const long row, const long col) {
    // verificando se o elemento a inserir esta dentro da matriz
    return this->operator()(row,col);
}




#endif /* TPZBlockDiag_h */
