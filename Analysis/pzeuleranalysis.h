//$Id: pzeuleranalysis.h,v 1.5 2003-10-21 18:10:58 erick Exp $

#ifndef PZEULERANALYSIS_H
#define PZEULERANALYSIS_H

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzflowcmesh.h"
#include "pzadmchunk.h"
#include "pzmaterial.h"
#include "pzconslaw.h"
#include "pzstrmatrix.h"
#include "pzsolve.h"

#include <iostream>

class TPZEulerAnalysis : public TPZAnalysis
{

public:

   TPZEulerAnalysis();

   TPZEulerAnalysis(TPZFlowCompMesh *mesh,std::ostream &out = cout);

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
    * Informs the Analysis class the time at which
    * the current solution in the computational
    * mesh belongs, so that the materials can
    * choose whether to contribute implicitly
    * or explicitly.
    */
   void SetContributionTime(TPZContributeTime time);

   /**
    * Adds deltaSol to the last solution and stores it as the current
    * Solution, preparing it to contribute again.
    * Also updates the CompMesh soution.
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
    * Buffers the assemblage of the rhs with
    * respect to the last state.
    * Sets the current solution in the mesh
    * when finished.
    */
   void BufferLastStateAssemble();

   /**
    * Assembles the stiffness matrix
    * BufferLastStateAssemble or UpdateHistory
    * must be called first.
    **/
   virtual void Assemble();

   /**
    * Assembles the right hand side vector.
    * BufferLastStateAssemble or UpdateHistory
    * must be called first.
    */
   virtual void Assemble(TPZFMatrix & rhs);

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
    * Defines the time step for each material in the mesh
    */
   void ComputeTimeStep();

   /**
    * Settings for the linear system solver
    */
   void SetLinSysCriteria(REAL epsilon, int maxIter);

   /**
    * Settings for the Newton's method.
    */
   void SetNewtonCriteria(REAL epsilon, int maxIter);

   /**
    * Settings for the time integration method.
    */
   void SetTimeIntCriteria(REAL epsilon, int maxIter);

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
    * Vector to hold the contribution of last state
    * to the rhs.
    */

   TPZFMatrix fRhsLast;

   /**
    * These two pointers to solutions switches
    * among the stored vectors, so that the
    * history can be made without unnecessary
    * copy operations.
    */
   TPZFMatrix * fpCurrSol, * fpLastSol, * fpSolution;

   /**
    * Stop criteria for the Newton and time
    * integration loops.
    *
    */
   REAL fLinSysEps;
   int fLinSysMaxIter;

   REAL fNewtonEps;
   int fNewtonMaxIter;

   REAL fTimeIntEps;
   int fTimeIntMaxIter;
};

#endif
