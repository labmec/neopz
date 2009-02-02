//$Id: pzeuleranalysis.h,v 1.20 2009-02-02 10:20:33 phil Exp $

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
#include "pzdxmesh.h"
#include "pzstepsolver.h"
#include "pzblockdiag.h"
#include "pzsave.h"

#include <iostream>

/// This class implements an analysis procedure for computing the steady state solution of a compressible Euler flow simulation
/**
This class implements several tricks to obtain the steady state solution as fast as possible
It evoluates the CFL of the simulation according to the residual of the system of equations and according to the number of iterations of the Newton method at each step
@author Erick Raggio Slis dos Santos
*/
class TPZEulerAnalysis : public TPZAnalysis
{

public:

   TPZEulerAnalysis();

   TPZEulerAnalysis(TPZFlowCompMesh *mesh,std::ostream &out = std::cout);

   ~TPZEulerAnalysis();


   /**
    * Writes the computational mesh onto disk
    */
   void WriteCMesh( const char * str);

   /**
    * see declaration in the base class.
    */
   virtual void Run(std::ostream &out, const std::string & dxout, int dxRes);

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
    * Indicates whether the CFL is to evolute or not
    */
   void SetEvolCFL(int EvolCFL);


   /**
    * Adds deltaSol to the last solution and stores it as the current
    * Solution, preparing it to contribute again.
    * Also updates the CompMesh soution.
    *
    * @param deltaSol [in] result of previous call to Solve.
    * @param epsilon [out] norm of resultant fRhs.
    *
    */
   void UpdateSolAndRhs(TPZFMatrix & deltaSol, REAL & epsilon);

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
    * Evaluates the flux part of the residual for
    * convergence check.
    */
   REAL EvaluateFluxEpsilon();

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
   virtual void AssembleRhs();

   /**
    * Solves an assembled stiffness matrix.
    * Not to be misunderstood as a global
    * solver, because it performs only part
    * of one linear iteration.
    * @param res [out] residual of the linear system solution
    * @param residual [in] pointer to a temporary matrix storage
    * @param delSol [out] returns the delta Solution vector.
    * @return 1 if everything fine, 0 otherwise (no better solution found)
    */
   int Solve(REAL & res, TPZFMatrix * residual, TPZFMatrix & delSol);

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
   REAL ComputeTimeStep();

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

   /**
    * Prepares the DX graph mesh
    */
   TPZDXGraphMesh * PrepareDXMesh(const std::string &dxout, int dxRes);

   /**
    * Informs a block diagonal to be used as preconditioning
    */
   void SetBlockDiagonalPrecond(TPZBlockDiagonal * blockDiag);
   
   /**
   This method will search for the $sol0+\alpha dir$ solution which minimizes the residual
   @param [in/out] residual initial residual (of sol0) and on exiting final residual
   @param [in/out] sol0 solution vector
   @param [in] dir search direction
   @return 1 if a direction was found, 0 otherwise
   */
   int LineSearch(REAL &residual, TPZFMatrix &sol0, TPZFMatrix &dir);
    void CompareRhs();

   /**
    * Evaluates the CFL control based on the newest residual norm
    * (epsilon) and last residual norm (lastEpsilon). Scales are
    * also applied to timeStep
    *
    * @param lastEpsilon [in/out] norm of the last flux vector
    * it is set to epsilon after CFL update.
    * @param epsilon [in] norm of the current flux vector.
    * @param epsilon_Newton [in] residual of the nonlinear invertion.
    * @param timeStep [in/out] parameter to accumulate the same
    * scales.
    *
    */
   void CFLControl(REAL & lastEpsilon, REAL & epsilon, REAL & epsilon_Newton, REAL & timeStep);
    void SetGMResFront(REAL tol, int numiter, int numvectors);
    void SetFrontalSolver();
    void SetGMResBlock(REAL tol, int numiter, int numvec);

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
    * Vector to hold the contribution of last state
    * to the rhs.
    */

   TPZFMatrix fRhsLast;

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


   
   /**
    * Indicates whether the CFL is to evolute or not.
    */
   int fEvolCFL;

   /**
    * Preconditioner
    */
   TPZBlockDiagonal * fpBlockDiag;
   
   /**
    * Total number of newton iterations during this run
    */
   int fTotalNewton;
   
   /**
    * Indication if a frontal matrix is being used as a preconditioner
    */
   int fHasFrontalPreconditioner;
   
};

#endif
