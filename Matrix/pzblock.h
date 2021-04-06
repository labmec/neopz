/**
 * @file pzadmchunk.h
 * @brief Declarates the TPZBlock<REAL>class which implements block matrices.
 */

#ifndef _TBLOCKHH_
#define _TBLOCKHH_

#include "pzmatrix.h"
#include "pzmanvector.h"
#include "pzreal.h"
#include "TPZSavable.h"


/**
 * @brief Implements block matrices. \ref matrixutility "Matrix utility"
 * @author Misael Luis Santana Mandujano
 * @since 12/1994
 * @ingroup matrixutility
 */
template<class TVar>
class TPZBlock : public TPZSavable
{
public:
    TPZBlock() : TPZRegisterClassId(&TPZBlock::ClassId), fBlock(), fpMatrix(0)
    {
        
    }
	/**
	 * @brief For each elements on matrix a size 1 block is created
	 * @param matrix_to_represent Indicates which matrix is to be represented
	 * @param num_of_blocks Indicates number of blocks
	 * @param initial_blocks_size Indicates initial block size, default value is 1
	 */
	TPZBlock(TPZMatrix<TVar> *const matrix_to_represent,const int num_of_blocks = 0,
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
	void SetMatrix(TPZMatrix<TVar> *const other)
    {
        fpMatrix = other;
    }
	
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
	
	TPZBlock<TVar>&operator=(const TPZBlock<TVar>& ); 		
	
    void PrintStructure(std::ostream &out = std::cout);
	
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

        int ClassId() const override;
	
        /** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;

private:
	
	/**
     * @struct TNode
     * @brief Defines a node
	 */
	class TNode : public TPZSavable {
            public :
                
		int pos; /**< Position of node */
		int dim; /**< Dimension of node */
		
		TNode() {
			pos=0;
			dim=0;
		}
                
        int ClassId() const override{
            return Hash("TNode") ^ ClassIdOrHash<TPZBlock<TVar>>() << 1;
        }
                
		void Read(TPZStream &buf, void *context) override;
		void Write(TPZStream &buf, int withclassid) const override;
	};
	
	/** @brief Nodes vector */
	TPZManVector<TNode>    fBlock;
	/**  @brief Pointer to TPZMatrix */
	TPZMatrix<TVar> *fpMatrix;
	static REAL gZero;//zero
	
};

template<class TVar>
int TPZBlock<TVar>::ClassId() const {
    return Hash("TPZBlock") ^ ClassIdOrHash<TVar>() << 1;
}

#endif
