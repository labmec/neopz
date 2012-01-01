/**
 * @file
 * @brief Contains TPZBlockDiagonal class which defines block diagonal matrices.
 */
// Author: Philippe Devloo.
//
// File:   pzblockdiag.h
//
// Class:  TPZBlockDiagonal
//
// Obs.:   Implements a block diagonal matrix
//
// Version : 1 - 2001
//


#ifndef _TBLOCKDIAG_
#define _TBLOCKDIAG_


#include <iostream>
#include "pzmatrix.h"
#include "pzvec.h"

/**
 @brief Defines block diagonal matrices. \ref matrix "Matrix"
 @ingroup matrix
 */
class TPZBlockDiagonal : public TPZMatrix
{
	
public:
	/** @brief Simple constructor */
	TPZBlockDiagonal ();
	/**
     @brief Constructor with initialization parameters
     @param blocksizes Size of blocks on Block Diagonal matrix
     @param glob Global matrix which will be blocked
	 */
	TPZBlockDiagonal (const TPZVec<int> &blocksizes, const TPZFMatrix &glob );
	/**
     @brief Constructor with initialization parameters
     @param blocksizes Size of blocks on Block Diagonal matrix
	 */
	TPZBlockDiagonal (const TPZVec<int> &blocksizes);
	/** @brief Copy constructor */
	TPZBlockDiagonal (const TPZBlockDiagonal & );
	
	CLONEDEF(TPZBlockDiagonal)
	/** @brief Simple destructor */
	~TPZBlockDiagonal();
	
	int    Put(const int row,const int col,const REAL& value );
	const REAL &Get(const int row,const int col ) const;
	
	REAL &operator()(const int row, const int col);
	virtual REAL &s(const int row, const int col);
	//estos metodos nao verificam a existencia do elemento
	//sao mas rapidos que Put e Get
	int    PutVal(const int row,const int col,const REAL& value );
	const REAL &GetVal(const int row,const int col ) const;
	
	void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
				 const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 ) const ;
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	
	// Peforms the product (*this)T x D x (*this).
	//  TPZBlockDiagonal  InnerProd(TPZBlockDiagonal &D );
	
	
	int Dim() const     { return Rows(); }
	
	// Zera os elementos da matriz
	int Zero();
	
	/**
	 * @brief Return the choosen block size
	 * @param blockid - block index
	 */
	int GetSizeofBlock(int blockid) {return fBlockSize[blockid];}
	
	void Transpose(TPZMatrix *const T) const;
	virtual int Decompose_LU();
	virtual int Decompose_LU(std::list<int> &singular);
	
	// Faz o Backward e Forward substitutions para a matriz
	// decomposta com LU
	/// Makes the backward and forward substitutions whether the matrix was LU decomposed
	virtual int Substitution( TPZFMatrix * B ) const;
	
	/** @brief Updates the values of the matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix> mat);
	
	/**
	 * method which checks the working of the class
	 */
	static int main();
    
    /** Fill the matrix with random values (non singular matrix) */
    void AutoFill();

	
private:
	
	
	//static int Error(const char *msg1,const char *msg2="" );
	/** @brief Clean data matrix. Zeroes number of columns and rows. */
	int Clear();
public:
	/**
     @brief Initializes current matrix based on blocksize
     @param blocksize Used to initialize current matrix
	 */
	void Initialize(const TPZVec<int> &blocksize);
	/**
     * @brief Adds a block to current matrix
     * @param i Adds in ith position
     * @param block Block to be added
	 */
	void AddBlock(int i, TPZFMatrix &block);
	/**
     * @brief Sets a block in the current matrix
     * @param i Adds in ith position
     * @param block Block to be added
	 */
	void SetBlock(int i, TPZFMatrix &block);
	
	/**
     * @brief Gets a block from current matrix
     * @param i Returns teh ith block
     * @param block Contains returned block
	 */
	void GetBlock(int i, TPZFMatrix &block);
	
	/**
     @brief Builds a block from matrix
     @param matrix Matrix to build from
	 */
	void BuildFromMatrix(TPZMatrix &matrix);
	/**
	 * @brief Prints current matrix data
     * @param message Message to be printed
     * @param out Output device
	 * @param format Output format to print
	 */
	virtual void Print(const char *message, std::ostream &out = std::cout, const MatrixOutputFormat format =EFormatted) const;
	
	int NumberofBlocks() {return fBlockSize.NElements();}
	
protected:
	/**
     @brief Stores matrix data
	 */
	TPZVec<REAL> fStorage;
	/**
     @brief Stores blocks data
	 */
	TPZVec<int> fBlockPos;
	/**
     @brief Stores block sizes data
	 */
	TPZVec<int> fBlockSize;
};


inline REAL &TPZBlockDiagonal::s(const int row, const int col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}


#endif

