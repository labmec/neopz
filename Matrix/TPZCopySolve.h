/* Generated by Together */

#ifndef TPZCOPYSOLVE_H
#define TPZCOPYSOLVE_H
#include "pzsolve.h"
class TPZFMatrix;
class TPZCopySolve : public TPZMatrixSolver {
public:

  TPZCopySolve(TPZMatrix *mat) : TPZMatrixSolver(mat)
  {
  }
  TPZCopySolve(const TPZCopySolve &other) : TPZMatrixSolver(other)
  {
  }
    /**
     *      Solves the system of linear equations stored in current matrix
     *      As this class implements only a copy operation, it just copies u to F;
     *      @param F contains Force vector
     *      @param result contains the solution
     */
    void Solve(const TPZFMatrix &F, TPZFMatrix &result, TPZFMatrix *residual)
      {
        result=F;
      }  

    /**
     *      Clones the current object returning a pointer of type TPZSolver
      */
    TPZSolver *Clone() const 
    {
      return new TPZCopySolve(*this);
    }


};
#endif //TPZCOPYSOLVE_H
