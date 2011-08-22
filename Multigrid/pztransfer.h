/**
 * @file
 * @brief Contains the TPZTransfer class which implements a rectangular sparse block matrix.
 */
#ifndef TRANSFERH
#define TRANSFERH

#include "pzblock.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzmatrix.h"

/**
 * @brief Implements rectangular matrix which extends a solution vector of the coarse mesh to a solution vector in the fine mesh. \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Implements a rectangular sparse block matrix it is assumed that the data is entered one row at a time \n
 * the matrix structure cannot be modified after being defined
 */
class TPZTransfer : public TPZMatrix {
	
	public :
	
    TPZTransfer();
	
	/** @brief The sparse matrix blocks are defined by row, col */
	TPZTransfer(TPZBlock &row, TPZBlock &col,int nvar, int nrowblocks, int ncolblocks);
	
	TPZTransfer(const TPZTransfer &cp) : TPZMatrix(cp),
	fNStateVar(cp.fNStateVar), fRowBlock(cp.fRowBlock),
	fColBlock(cp.fColBlock),fColPosition(cp.fColPosition),
	fNumberofColumnBlocks(cp.fNumberofColumnBlocks),
	fColumnBlockNumber(cp.fColumnBlockNumber),
	fColumnBlockLastUsed(cp.fColumnBlockLastUsed),
	fDoubleValues(cp.fDoubleValues),
	fDoubleValLastUsed(cp.fDoubleValLastUsed)
	{
	}
	
	CLONEDEF(TPZTransfer)
	
	//TPZMatrix : EFormatted, EInputFormat, EMathematicaInput
	virtual void Print(const char *name = NULL, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const;
	
	/** @brief Identifies the number of equations per shapefunction*/
	void SetNStateVariables(int statevar) { fNStateVar = statevar; }
	
	/** 
	 * @brief This operation will reset the matrix to zero with no rows defined
	 * @param row starting row of the block
	 * @param col starting column of the block
	 * @param nvar indicates the number of state variables of the problem
	 * @param nrowblocks number of rows of the block
	 * @param ncolblocks number of columns of the block
     * @note the stride of the matrix will be initialized by nvar
	 */
	void SetBlocks(TPZBlock &row, TPZBlock &col, int nvar, int nrowblocks, int ncolblocks);
	
	/** @brief Returns 1 if the row is defined (i.e. has column entries)*/
	int HasRowDefinition(int row);
	
	/** @brief Will specify the sparsity pattern of row*/
	void AddBlockNumbers(int row, TPZVec<int> &colnumbers);
	
	/** @brief Sets the row,col block equal to matrix mat
	 * if row col was not specified by AddBlockNumbers, \n
	 * an error will be issued and exit
	 */
	void SetBlockMatrix(int row, int col, TPZFMatrix &mat);
	
	/** @brief Multiplies the transfer matrix and puts the result in z*/
	void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
				 const REAL alpha, REAL beta, int opt = 0, int stride = 1) const ;
	
	/**
	 * @brief Will transfer the solution, taking into acount there may be more than
	 * one state variable
	 */
	void TransferSolution(const TPZFMatrix &coarsesol, TPZFMatrix &finesol);
	
	/**
	 * @brief Will transfer the residual, taking into acount there may be more than
	 * one state variable
	 */
	void TransferResidual(const TPZFMatrix &fine, TPZFMatrix &coarse);
	
	void Multiply(const TPZFMatrix &A, TPZFMatrix&B, int opt,
				  int stride) const;
	
private:
	
	/** @brief Increases the storage allocated
     int fColPosition to include numcol more values*/
	void ExpandColumnVectorEntries(int numcol);
	
	/** @brief Increases the storage space available in the fDoubleValues vector to
     include numval entries*/
	void ExpandDoubleValueEntries(int numval);
	
	/** @brief Number of variables associated with each shape function*/
	int fNStateVar;
	/** @brief Block sizes of the rows*/
	TPZBlock fRowBlock;
	/** @brief Block sizes of the columns*/
	TPZBlock fColBlock;
	/** @brief Vector indicating the starting column block for each row*/
	TPZVec<int> fColPosition;
	/** @brief Vector indicating the number of column blocks associated with each row*/
	TPZVec<int> fNumberofColumnBlocks;
	/** @brief Vector indicating the starting point of each column block*/
	TPZManVector<int> fColumnBlockPosition;
	/** @brief Vector indicating the number of the column corresponding to the block*/
	TPZManVector<int> fColumnBlockNumber;
	/** @brief Indicates the next free position*/
	int fColumnBlockLastUsed;
	/** @brief Storage space for the matrix blocks*/
	TPZManVector<REAL> fDoubleValues;
	/** @brief Indicates the next free position of fDoubleValues*/
	int fDoubleValLastUsed;
	
};

#endif

