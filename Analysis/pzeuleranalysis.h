#ifndef PZEULERANALYSIS_H
#define PZEULERANALYSIS_H

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzflowcmesh.h"
#include "pzadmchunk.h"
#include "pzmaterial.h"
#include "pzconslaw.h"
#include "pzstrmatrix.h"

#include <iostream>

class TPZAnalysis :: TPZEulerAnalysis
{
#ifdef NOTDEFINED

variables in the parent class
    protected:
  /**
   * Geometric Mesh
   */
  TPZGeoMesh *fGeoMesh;
  /**
   * Computational mesh
   */
  TPZCompMesh *fCompMesh;
  /**
   * Graphical mesh
   */
  TPZGraphMesh *fGraphMesh[3];
  /**
   * Load vector ???
   */
  TPZFMatrix fRhs;
  /**
   * Solution vector
   */
  TPZFMatrix fSolution;
  /**
   * Type of solver to be applied
   */
  TPZMatrixSolver *fSolver;
  /**
   * Scalar variables names - to post process
   */
  TPZVec<char *> fScalarNames[3];
  /**
   * Vector variables names - to post process
   */
  TPZVec<char *> fVectorNames[3];
  /**
   * Step ???
   */
  int fStep;
  /**
   * Time ?? (time step ??)
   */
  REAL fTime;

  /**
   * Structural matrix
   */
  TPZStructMatrix * fStructMatrix;

  struct TTablePostProcess {
    TPZVec<int> fGeoElId;
    TPZVec<TPZCompEl *> fCompElPtr;
    int fDimension;
    TPZVec<REAL> fLocations;
    TPZVec<char *> fVariableNames;
    ostream *fOutfile;
    TTablePostProcess();
    ~TTablePostProcess();
  } fTable;
#endif

public:

   TPZEulerAnalysis();

   TPZEulerAnalysis(TPZFlowCompMesh *mesh,std::ostream &out = cout)

   ~TPZEulerAnalysis();

   /**
    * see declaration in the base class.
    */
   virtual void Run(ostream &out = cout);

   /**
    * Sets the solution vector to be the one
    * representing the Advanced State (Wn+1)
    * or the advanced implicit solution.
    *
    */
   void SetAdvancedState();

   /**
    * Sets the solution vector to be the one
    * representing the Current State (Wn)
    * or the last explicit solution.
    *
    */
   void SetLastState();

   /**
    * Adds deltaSol to the last solution and stores it as the current
    * Solution, preparing it to contribute again.
    *
    */
   void UpdateSolution(TPZFMatrix & deltaSol);

   /**
    * After a call to UpdateSolution, this method
    * shifts the history in one solution, so that
    * the current solution in the last iteration
    * becomes the last solution in the newest
    * iteration.
    * This method does not reset the current solution.
    * To reset it, UpdateSolution with a deltaSol
    * containing any valid value (zeros, for instance)
    * must be called.
    */
   void UpdateHistory();

   /**
    *Assemble the stiffness matrix
    **/
   virtual  void Assemble();

   /**
    * Solves an assembled stiffness matrix.
    * Not to be misunderstood as a global
    * solver, because it performs only part
    * of one linear iteration.
    */
   void Solve(REAL & res);

   /**
    * Implements the Newton's method.
    * Returns 1 if the desired tolerance was
    * reached, 0 if it exceeded the maximal
    * number of iterations.
    *
    * @param epsilon [in/out] expected tolerance
    * in the Newton's method. Returns the value
    * of the achieved tolerance.
    *
    * @param numIter [in/out] maximal number of
    * iterations of the Newton's method. Returns
    * the number of iterations needed by the method.
    */
   int RunNewton(REAL & epsilon, int & numIter);

   /**
    * Main Analysis routine
    */
   virtual void Run();

protected:

   /**
    * Stores a pointer to the computational
    * flow mesh. The inherited fCompMesh is also
    * updated with valid values to keep compatibility
    * with inherited methods, but in this class
    * fFlowCompMesh is used whenever a method exclusive
    * from this class is needed.
    */
   TPZFlowCompMesh * fFlowCompMesh;

   /**
    * Matrix to hols a second solution vector.
    * In the Euler residual, this second vector
    * implements the history of solutions for the
    * Euler residual.
    */
   TPZFMatrix fSolution2;

   /**
    * These two pointers to solutions switches
    * among the stored vectors, so that the
    * history can be made without unnecessary
    * copy operations.
    */
   TPZFMatrix * fpCurrSol, * fpLastSol, * fpSolution;

};

#endif
