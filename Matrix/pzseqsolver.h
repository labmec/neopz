/* Generated by Together */

#ifndef TPZSEQUENCESOLVER_H
#define TPZSEQUENCESOLVER_H
#include "pzsolve.h"
#include "pzstack.h"

class TPZFMatrix;

/**
   Defines sequence solvers
   @ingroup solver
*/
class TPZSequenceSolver : public TPZMatrixSolver {
public:
  /**
     Constructor with initialization parameter
     @param refmat Sets reference matrix to NILL
  */
  TPZSequenceSolver(TPZMatrix *refmat = 0);
  /**
     Copy constructor
     @param copy Model object to be copied from
  */
  TPZSequenceSolver(const TPZSequenceSolver & copy);
  
  void Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual = 0);

  /**
  This method will reinitialize the solver object, including the solution procedure
  */  
  void ResetSolver();
  
  /**
  This method will reset the matrix associated with the solver
  This is useful when the matrix needs to be recomputed in a non linear problem
  */
  virtual void ResetMatrix();
  
  /**
  * Updates the values of the preconditioner based on the values of the matrix
  */
  virtual void UpdateFrom(TPZMatrix *mat);
  /**
  This method gives a preconditioner to share a matrix with the referring solver object
  */
//  virtual void SetMatrix(TPZMatrixSolver *solver);
  
  void AppendSolver(TPZMatrixSolver & solve);
  
  virtual TPZSolver * Clone() const;

 private:    
  TPZStack < TPZMatrixSolver * > fSolvers;
};
#endif //TPZSEQUENCESOLVER_H
