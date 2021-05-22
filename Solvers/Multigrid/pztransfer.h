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
template<class TVar>
class TPZTransfer : public TPZMatrix<TVar> {
	
public :
    int ClassId() const override;

	/** @brief Default constructor */
    TPZTransfer();
	
	/** @brief The sparse matrix blocks are defined by row, col */
	//TPZTransfer(TPZBlock &row, TPZBlock &col,int nvar, int nrowblocks, int ncolblocks);
	TPZTransfer(TPZBlock &row, TPZBlock &col,int nvar, int nrowblocks, int ncolblocks);
    
    /// copy constructor
	TPZTransfer(const TPZTransfer &cp) : TPZRegisterClassId(&TPZTransfer::ClassId),
    TPZMatrix<TVar>(cp),
	fNTVarVar(cp.fNTVarVar), fRowBlock(cp.fRowBlock),
	fColBlock(cp.fColBlock),fColPosition(cp.fColPosition),
	fNumberofColumnBlocks(cp.fNumberofColumnBlocks),
	fColumnBlockNumber(cp.fColumnBlockNumber),
	fColumnBlockLastUsed(cp.fColumnBlockLastUsed),
	fDoubleValues(cp.fDoubleValues),
	fDoubleValLastUsed(cp.fDoubleValLastUsed)
	{
	}
	inline TPZTransfer*NewMatrix() const override {return new TPZTransfer{};}
	virtual TPZMatrix<TVar> *Clone() const  override { return new TPZTransfer(*this); }

	/** @brief Creates a copy from another TPZTransfer*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override        
  {                                                           
    auto *from = dynamic_cast<const TPZTransfer<TVar> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }
	
	//TPZMatrix<REAL> : EFormatted, EInputFormat, EMathematicaInput
	virtual void Print(const char *name = NULL, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const override;
	
	/** @brief Identifies the number of equations per shapefunction*/
	void SetNTVarVariables(int TVarvar) { fNTVarVar = TVarvar; }
	
	/** 
	 * @brief This operation will reset the matrix to zero with no rows defined
	 * @param row starting row of the block
	 * @param col starting column of the block
	 * @param nvar indicates the number of TVar variables of the problem
	 * @param nrowblocks number of rows of the block
	 * @param ncolblocks number of columns of the block
	 */
	//void SetBlocks(TPZBlock &row, TPZBlock &col, int nvar, int nrowblocks, int ncolblocks);
	void SetBlocks(TPZBlock &row, TPZBlock &col, int nvar, int nrowblocks, int ncolblocks);
	
	/** @brief Returns 1 if the row is defined (i.e. has column entries)*/
	int HasRowDefinition(int row);
	
	/** @brief Will specify the sparsity pattern of row*/
	void AddBlockNumbers(int row, TPZVec<int> &colnumbers);
	//it previously had no implementation, it would return zero anyway
	const TVar GetVal(int64_t,int64_t) const override{return (TVar)0;}
	
	/** @brief Sets the row,col block equal to matrix mat
	 * if row col was not specified by AddBlockNumbers, \n
	 * an error will be issued and exit
	 */
	void SetBlockMatrix(int row, int col, TPZFMatrix<TVar> &mat);
	
	/** @brief Multiplies the transfer matrix and puts the result in z*/
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha,const TVar beta, const int opt = 0) const  override;
	
    /** @brief Multiplies the transfer matrix and puts the result in z*/
    void MultAddScalar(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                 const TVar alpha,const TVar beta, const int opt = 0) const;
	/**
	 * @brief Will transfer the solution, taking into acount there may be more than
	 * one TVar variable
	 */
	void TransferSolution(const TPZFMatrix<TVar> &coarsesol, TPZFMatrix<TVar> &finesol);
	
	/**
	 * @brief Will transfer the residual, taking into acount there may be more than
	 * one TVar variable
	 */
	void TransferResidual(const TPZFMatrix<TVar> &fine, TPZFMatrix<TVar> &coarse);
	
//	void Multiply(const TPZFMatrix<TVar> &A, TPZFMatrix<TVar> &B, int opt) const override;
	
    void MultiplyScalar(const TPZFMatrix<TVar> &A, TPZFMatrix<TVar> &B, int opt) const;
protected:
	inline TVar *&Elem() override
  {
    return fDoubleValues.begin();
  }
  inline const TVar *Elem() const override
  {
    return fDoubleValues.begin();
  }
  inline int64_t Size() const override
  {
    return fDoubleValues.size();
  }
private:
	
	/** @brief Increases the storage allocated
     int fColPosition to include numcol more values*/
	void ExpandColumnVectorEntries(int numcol);
	
	/** @brief Increases the storage space available in the fDoubleValues vector to
     include numval entries*/
	void ExpandDoubleValueEntries(int numval);
	
	/** @brief Number of variables associated with each shape function*/
	int fNTVarVar;
	/** @brief Block sizes of the rows*/
	//TPZBlock fRowBlock;
	TPZBlock fRowBlock;
	/** @brief Block sizes of the columns*/
	//TPZBlock fColBlock;
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
	TPZManVector<TVar> fDoubleValues;
	/** @brief Indicates the next free position of fDoubleValues*/
	int fDoubleValLastUsed;
	
};

template<class TVar>
int TPZTransfer<TVar>::ClassId() const{
    return Hash("TPZTransfer") ^ TPZMatrix<TVar>::ClassId() << 1;
}

#endif

