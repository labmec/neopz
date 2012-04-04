/**
 * @file pzadmchunk.h
 * @brief Declarates the TPZBlock<REAL>class which implements block matrices.
 */

#ifndef _TBLOCKHH_
#define _TBLOCKHH_

#include "pzmatrix.h"
#include "pzmanvector.h"
#include "pzreal.h"
#include "pzsave.h"

/** @brief Id of the diagonal block matrix */
/** @ingroup matrixutility */
const int TPZBLOCKID = 101;

/**
 * @brief Implements block matrices. \ref matrixutility "Matrix utility"
 * @author Misael Luis Santana Mandujano
 * @since 12/1994
 * @ingroup matrixutility
 */
template<class TVar>
class TPZBlock : public TPZSaveable
{
public: 
	/**
	 * @brief For each elements on matrix a size 1 block is created
	 * @param matrix_to_represent Indicates which matrix is to be represented
	 * @param num_of_blocks Indicates number of blocks
	 * @param initial_blocks_size Indicates initial block size, default value is 1
	 */
	TPZBlock(TPZMatrix<TVar> *const matrix_to_represent = 0,const int num_of_blocks = 0,
			 const int initial_blocks_size = 1 );
	
	/**
	 * @brief Copy constructor
	 * @param bl New object is created based on bl
	 */
	TPZBlock(const TPZBlock<TVar> &bl);
	
	/** @brief Simple Destrutor */
	virtual ~TPZBlock();
	
	/**
     * @brief Changes pointer to other
     * @param other New matrix to be pointed to
	 */
	virtual void SetMatrix(TPZMatrix<TVar> *const other);
	
	/** @brief Returns a pointer to current matrix */
	TPZMatrix<TVar> *Matrix(){ return fpMatrix;}
	
	/**
     * @brief Sets number of blocks on diagonal matrix
     * @param num_of_blocks Number of blocks
	 */
	int SetNBlocks(const int num_of_blocks );
	
	/**
     * @brief Modifies existing block dimensions or creates a new block with given index
     * @param index Given index to be redimensioned or created
     * @param dim New dimension
     * @param pos New position
	 */
	int Set(const int index,const int dim,const int pos = -1 );
	
	/**
     * @brief Computes blocks sequence
     * @param dimensions Contains blocks sequence
	 */
	int SetAll( TPZVec<int> & dimensions );

	/**
     * @brief Resequences blocks positioning
     * @param start Starting position
	 */
	int Resequence(const int start=0);
	
	/**
     * @brief Removes a block
     * @param index Index of the block to be removed
	 */
	int Remove(const int index );
	
	/**
     * @brief Verifies if blocks are sequential and does not overcome matrix size
	 */
	int Verify() const;
	
	TVar & operator()(const int block_row,const int block_col,const int r,const int c ) const;
	
	/** @brief Gets a element from matrix verifying */
	const TVar & Get(const int block_row,const int block_col,const int r,const int c ) const;
	/** @brief Puts a element to matrix verifying */
	int Put(const int block_row,const int block_col,const int r,const int c,const TVar& value );
	
	/** @brief Gets a element from a matrix verifying the existence */
	const TVar & Get(const int block_row,const int r,const int c ) const;
	/** @brief Puts a element to matrix verifying the existence */
	int Put(const int block_row,const int r,const int c,const TVar& value );
	
	/** @brief Gets a element from matrix but not verify the existence */
	const TVar & GetVal(const int bRow,const int bCol,const int r,const int c ) const;
	/** @brief Puts a element to matrix but not verify the existence */
	int PutVal(const int bRow,const int bCol,const int r,const int c,const TVar& value );
	
	/**
     * @brief Puts a block on current matrix
     * @param block_row Contains block row
     * @param block_col Contains block column
     * @param block Block to be inserted
	 */
	int PutBlock(const int block_row,const int block_col,const TPZFMatrix<TVar> & block );
	/**
     * @brief Gets a block on current matrix
     * @param block_row Contains block row
     * @param block_col Contains block column
     * @param block Block to be inserted
	 */
	int GetBlock(const int block_row,const int block_col, TPZFMatrix<TVar> *const block ) const;
	/**
     * @brief Adds a block on current matrix
     * @param block_row Contains block row
     * @param block_col Contains block column
     * @param block Block to be inserted
	 */
	int AddBlock(const int block_row,const int block_col, const TPZFMatrix<TVar> & block );
		
	/**
     * @brief Inserts a block (block_row , block_col) on current matrix target
     * @param block_row Contains block row
     * @param block_col Contains block column
     * @param target Block to be inserted
     * @param row Starting row position
     * @param col Starting column position
	 */
	int InsertBlock(const int block_row,const int block_col,
					const int row,const int col, TPZMatrix<TVar> &target) const;
	
	TPZBlock<TVar>&operator=(const TPZBlock<TVar>& ); 
	 
	/**
     * @brief Prints a matrix block
     * @param block_row Contains block row
     * @param block_col Contains block column
     * @param title Title on printed output device
     * @param out Output device
	 */
	int  PrintBlock(const int block_row,const int block_col,const char *title = "",TPZostream &out = std::cout ) const;
	
	/// Prints all the blocks of the matrix
	void Print(const char *title = "",TPZostream &out = std::cout,TPZMatrix<TVar> *mat=NULL);
	
	void PrintSolution(const char *title, TPZostream &out);
	
	/** @brief Returns the max number of blocks on diagonal */
	int MaxBlockSize() const {return fBlock.NElements();}
	/** @brief Returns number of blocks on diagonal */
	int NBlocks() const {return fBlock.NElements();}
	
	/**
     * @brief Returns block dimension
     * @param block_diagonal Inquired block_diagonal
	 */
	int Size(const int block_diagonal) const { return fBlock[block_diagonal].dim; }
	
	/**
     * @brief Returns the position of first element block dependent on matrix diagonal
     * @param block_diagonal Inquired block_diagonal
	 */
	int Position(const int block_diagonal) const { return fBlock[block_diagonal].pos;}
	
	/** @brief Returns matrix dimension pointed by block */
	int Dim() const {return fBlock.NElements() ? fBlock[fBlock.NElements()-1].pos+fBlock[fBlock.NElements()-1].dim : 0; }
	
	/** @brief returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);

private:
	
	/**
     * @struct TNode
     * @brief Defines a node
	 */
	struct TNode
	{
		int pos; /**< Position of node */
		int dim; /**< Dimension of node */
		
		TNode() {
			pos=0;
			dim=0;
		}
		void Read(TPZStream &buf, void *context);
		void Write(TPZStream &buf, void *context);
	};
	
	/** @brief Nodes vector */
	TPZManVector<TNode>    fBlock;
	/**  @brief Pointer to TPZMatrix */
	TPZMatrix<TVar> *fpMatrix;
	static REAL gZero;//zero
	
};

#endif
