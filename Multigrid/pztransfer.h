#ifndef TRANSFERH
#define TRANSFERH

#include "pzblock.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzmatrix.h"

/**this class implements a rectangular sparse block matrix
it is assumed that the data is entered one row at a time
the matrix structure cannot be modified after being defined*/
class TPZTransfer : public TPZMatrix {

  public :

    TPZTransfer();

  /**the sparse matrix blocks are defined by row, col*/
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

  /**identifies the number of equations per shapefunction*/
  void SetNStateVariables(int statevar) { fNStateVar = statevar; }

  /**this operation will reset the matrix to zero with no rows defined
     nvar indicates the number of state variables of the problem
     the stride of the matrix will be initialized by nvar*/
  void SetBlocks(TPZBlock &row, TPZBlock &col, int nvar, int nrowblocks, int ncolblocks);

  /**returns 1 if the row is defined (i.e. has column entries)*/
  int HasRowDefinition(int row);

  /**will specify the sparsity pattern of row*/
  void AddBlockNumbers(int row, TPZVec<int> &colnumbers);

  /**sets the row,col block equal to matrix mat
     if row col was not specified by AddBlockNumbers, an error
     will be issued and exit*/
  void SetBlockMatrix(int row, int col, TPZFMatrix &mat);

  /**multiplies the transfer matrix and puts the result in z*/
  void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
	       const REAL alpha, REAL beta, int opt = 0, int stride = 1) const ;

  /**
   * Will transfer the solution, taking into acount there may be more than
   * one state variable
   */
  void TransferSolution(const TPZFMatrix &coarsesol, TPZFMatrix &finesol);

  /**
   * Will transfer the residual, taking into acount there may be more than
   * one state variable
   */
  void TransferResidual(const TPZFMatrix &fine, TPZFMatrix &coarse);

  void Multiply(const TPZFMatrix &A, TPZFMatrix&B, int opt,
		 int stride) const;

 private:

  /**increases the storage allocated
     int fColPosition to include numcol more values*/
  void ExpandColumnVectorEntries(int numcol);

  /**increases the storage space available in the fDoubleValues vector to
     include numval entries*/
  void ExpandDoubleValueEntries(int numval);

  /**number of variables associated with each shape function*/
  int fNStateVar;
  /**block sizes of the rows*/
  TPZBlock fRowBlock;
  /**block sizes of the columns*/
  TPZBlock fColBlock;
  /**vector indicating the starting column block for each row*/
  TPZVec<int> fColPosition;
  /**vector indicating the number of column blocks associated with each row*/
  TPZVec<int> fNumberofColumnBlocks;
  /**vector indicating the starting point of each column block*/
  TPZManVector<int> fColumnBlockPosition;
  /**vector indicating the number of the column corresponding to the block*/
  TPZManVector<int> fColumnBlockNumber;
  /**indicates the next free position*/
  int fColumnBlockLastUsed;
  /**storage space for the matrix blocks*/
  TPZManVector<REAL> fDoubleValues;
  /**indicates the next free position of fDoubleValues*/
  int fDoubleValLastUsed;

};

#endif

