//$Id: pzeuleranalysis.cpp,v 1.7 2003-10-24 19:03:10 erick Exp $

#include "pzeuleranalysis.h"
#include "pzerror.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"

TPZEulerAnalysis::TPZEulerAnalysis():
TPZAnalysis(), fFlowCompMesh(NULL),
fSolution2(), fRhsLast(), fpCurrSol(NULL),
fpLastSol(NULL), fpSolution(NULL),
fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100)
{

}

TPZEulerAnalysis::TPZEulerAnalysis(TPZFlowCompMesh *mesh, std::ostream &out):
TPZAnalysis(mesh, out), fFlowCompMesh(mesh),
fSolution2(), fRhsLast(), fpCurrSol(NULL), fpLastSol(NULL),
fpSolution(NULL), fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100)
{

}

TPZEulerAnalysis::~TPZEulerAnalysis()
{
}

void TPZEulerAnalysis::SetAdvancedState()
{
   fpSolution = fpCurrSol;
   SetContributionTime(Advanced_CT);
   fCompMesh->LoadSolution(*fpSolution);
}

void TPZEulerAnalysis::SetLastState()
{
   fpSolution = fpLastSol;
   SetContributionTime(Last_CT);
   fCompMesh->LoadSolution(*fpSolution);
}

void TPZEulerAnalysis::SetContributionTime(TPZContributeTime time)
{
   fFlowCompMesh->SetContributionTime(time);
}

void TPZEulerAnalysis::UpdateSolution(TPZFMatrix & deltaSol)
{
   (*fpCurrSol) += deltaSol;

   fCompMesh->LoadSolution(*fpCurrSol);
}

void TPZEulerAnalysis::UpdateHistory()
{
   TPZFMatrix * pBuff;

   // switching the current and last solution storages.
   pBuff = fpCurrSol;
   fpCurrSol = fpLastSol;
   fpLastSol = pBuff;

#ifdef RESTART_ZEROED
   //Zeroeing the newest iterative solution
   fpCurrSol.Zero();
#else
   // copying the laststate to the current state, as
   // a first approximation to Wn+1
   // Notice that this manner of computing the history
   // did not needed the pointer swap. Anyway, this
   // structure will be kept to allow a Zero as the
   // first guess of initial solution without a great
   // effort.
   (*fpCurrSol) = (*fpLastSol);
#endif
   // The Current Solution does not have any valid data
   // UpdateSolution must be called afterwards.

   // buffering the contribution to the RHS
   BufferLastStateAssemble();
}

void TPZEulerAnalysis::BufferLastStateAssemble()
{
   fRhsLast.Zero();
   fRhsLast.Redim(fCompMesh->NEquations(),1);
   SetLastState();
   fFlowCompMesh->Assemble(fRhsLast);
   SetAdvancedState();
}

void TPZEulerAnalysis::Assemble()
{
   if(!fCompMesh)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Computational Mesh\n";
      return;
   }

   if(!fStructMatrix)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
      return;
   }

   if(!fSolver)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Solver\n";
      return;
   }


   // retrieving the matrix stored in the Solver,
   // attempting to reuse it.
   TPZMatrix * pTangentMatrix = fSolver->Matrix();

   // verifies if the matrix really exists
   if(!pTangentMatrix)
   {
      // creates a matrix prepared for contribution
      pTangentMatrix = fStructMatrix->Create();

      // attaches the matrix to a solver.
      fSolver->SetMatrix(pTangentMatrix);
   }else{

      pTangentMatrix->Zero();
   }

   // redimensions and zeroes Rhs
   fRhs.Redim(fCompMesh->NEquations(),1);

   // Contributing referring to the advanced state
   // (n+1 index)
   fStructMatrix->Assemble(*pTangentMatrix, fRhs);

   // Contributing referring to the last state (n index)
   fRhs+=/*.Add(fRhsLast, fRhsLast)*/ fRhsLast;
}


void TPZEulerAnalysis::AssembleRhs()
{
   if(!fCompMesh) return;

   // redimensions and zeroes Rhs
   fRhs.Redim(fCompMesh->NEquations(),1);

   // Contributing referring to the advanced state
   // (n+1 index)
   fFlowCompMesh->Assemble(fRhs);

   // Contributing referring to the last state (n index)
   fRhs+=/*.Add(fRhsLast, */fRhsLast/*)*/;
}

void TPZEulerAnalysis::Solve(REAL & res) {
   int numeq = fCompMesh->NEquations();
   if(fRhs.Rows() != numeq ) return;

   TPZFMatrix residual(fRhs);
   TPZFMatrix delu(numeq,1);
   if(fpSolution->Rows() != numeq) {
     fpSolution->Redim(numeq,1);
   } else {
     fSolver->Matrix()->Residual(*fpSolution,fRhs,residual);
   }
   //      REAL normres  = Norm(residual);
   //	cout << "TPZAnalysis::Solve residual : " << normres << " neq " << numeq << endl;
   fSolver->Solve(residual, delu);
   //fSolution += delu;
   UpdateSolution(delu);

   //fCompMesh->LoadSolution(*fpSolution);
   res = Norm(residual);
}

int TPZEulerAnalysis::RunNewton(REAL & epsilon, int & numIter)
{
   int i = 0;
   REAL res;// residual of linear invertion.

   Assemble(); // assembles the linear system.
   // It is assumed that BufferLastStateAssemble(); was already called.
   epsilon = fNewtonEps * 2.;// ensuring the loop will be
   // performed at least once.

   while(i < fNewtonMaxIter && epsilon > fNewtonEps)
   {
      //Solves the linearized system and updates the solution.
      Solve(res);

      // Linearizes the system with the newest iterative solution
      Assemble();

      epsilon = Norm(fRhs);
      cout << "\nEpsilon:" << epsilon;
      i++;
   }

   // updates the history of state variable vectors
   UpdateHistory();

   numIter = i;

   if(epsilon > fNewtonEps)return 0;
   return 1;
}

void TPZEulerAnalysis::Run(ostream &out)
{
   // this analysis loop encloses several calls
   // to Newton's linearizations, updating the
   // time step.

   out << "\nBeginning time integration";

   REAL epsilon_Newton/*, epsilon_Global*/;
   int numIter_Newton/*, numIter_Global*/;
   REAL epsilon;

   int i = 0;
   epsilon = 2. * fTimeIntEps;

   // initializing the solution pointers
   fpCurrSol = & fSolution;
   fpLastSol = & fSolution2;

   // copying the solution from the mesh into the sol vector.
   fpLastSol->operator=(fFlowCompMesh->Solution());

#ifdef RESTART_ZEROED
   fpCurrSol->Zero();
#else
   // the history must be rebuilt -> a current state does not exist yet.
   fpCurrSol->operator=(*fpLastSol); // deltaState = 0;
#endif
/*
// The lines below initialize all the variables with Ones.
   for(int j = 0; j < fpLastSol->Rows(); j++)
   fpLastSol->operator()(j,0) = 1.;
   fFlowCompMesh->LoadSolution(*fpLastSol);
*/
   // Buffers the contribution of the last state
   BufferLastStateAssemble();

   // evaluates the time step based on the solution
   // (last Sol, since Curr sol equals it)
   ComputeTimeStep();

   while(i < fTimeIntMaxIter && epsilon > fTimeIntEps)
   {
      // Solves the nonlinear system, updates the solution,
      // history and last state assemble buffer.
      RunNewton(epsilon_Newton, numIter_Newton);

      // resetting the forcing function for iterations greater than the 1st
      fFlowCompMesh->SetFlowforcingFunction(NULL);

      // Computing the time step, verifying the convergency
      // using the newest time step.
      ComputeTimeStep();

      AssembleRhs(); // computing the residual only
      epsilon = Norm(fRhs);

      out << "\niter:" << i
          << " eps=" << epsilon
	  << " |NewtonEps=" << epsilon_Newton
	  << " nIter=" << numIter_Newton;
   }

   out.flush();

   fSolver->ResetMatrix(); // deletes the memory allocated
   // for the storage of the tangent matrix.
}

void TPZEulerAnalysis::ComputeTimeStep()
{
  fFlowCompMesh->ComputeTimeStep();
}

void TPZEulerAnalysis::SetLinSysCriteria(REAL epsilon, int maxIter)
{
   fLinSysEps     = epsilon;
   fLinSysMaxIter = maxIter;
}

void TPZEulerAnalysis::SetNewtonCriteria(REAL epsilon, int maxIter)
{
   fNewtonEps     = epsilon;
   fNewtonMaxIter = maxIter;
}

void TPZEulerAnalysis::SetTimeIntCriteria(REAL epsilon, int maxIter)
{
   fTimeIntEps     = epsilon;
   fTimeIntMaxIter = maxIter;
}

