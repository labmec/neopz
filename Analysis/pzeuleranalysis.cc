#include "pzeuleranalysis.h"
#include "pzerror.h"
#include "TPZCompElDisc.h"

TPZEulerAnalysis::TPZEulerAnalysis():
TPZAnalysis(), fFlowCompMesh(NULL),
fSolution2(), fpCurrSol(NULL),
fpLastSol(NULL), fpSolution(NULL),
fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100)
{

}

TPZEulerAnalysis::TPZEulerAnalysis(TPZFlowCompMesh *mesh, std::ostream &out):
TPZAnalysis(mesh, out), fFlowCompMesh(mesh),
 fSolution2(), fpCurrSol(NULL), fpLastSol(NULL),
fpSolution(NULL), fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100)
{

}

TPZEulerAnalysis::~TPZEulerAnalysis()
{

}

void RunNewton(REAL & epsilon, int & numIter)
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
   SetContributionTime(Current_CT);
   fCompMesh->LoadSolution(*fpSolution);
}

void TPZEulerAnalysis::SetContributionTime(TPZContributeTime time)
{
   fFlowCompMesh->SetContributionTime(time);
/*   int i, NumFluid;
   NumFluid = fFluidMaterial.NElements();
   for(i = 0; i < NumFluid; i++)
   {
      fFluidMaterial[i]->SetContributionTime(time);
   }*/
}

void TPZEulerAnalysis::UpdateSolution(TPZFMatrix & deltaSol)
{
   (*fpCurrSol) = (*fpLastSol);
   (*fpCurrSol) += deltaSol;
}

void TPZEulerAnalysis::UpdateHistory()
{
   TPZFMatrix * pBuff;

   // switching the current and last solution storages.
   pBuff = fpCurrSol;
   fpCurrSol = fpLastSol;
   fpLastSol = pBuff;

   // The Current Solution does not have any valid data
   // UpdateSolution must be called afterwards.
}

void TPZEulerAnalysis::Assemble()
{
   if(!fCompMesh || !fStructMatrix || !fSolver) return;

   // redimensions and zeroes Rhs
   fRhs.Redim(fCompMesh->NEquations(),1);

   // resets the solver matrix
   fSolver->SetMatrix(0);

   // creates a matrix prepared for contribution
   TPZMatrix * pTangentMatrix = fStructMatrix->Create();

   // Contributing referring to the last state (n index)
   SetLastState();
   fStructMatrix->Assemble(*pTangentMatrix, fRhs);

   // Contributing referring to the advanced state
   // (n+1 index)
   SetAdvancedState();
   fStructMatrix->Assemble(*pTangentMatrix, fRhs);
   //fSolver->SetMatrix(fStructMatrix->CreateAssemble(fRhs));

   // attaches the matrix to a solver.
   fSolver->SetMatrix(pTangentMatrix);
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

   //fCompMesh->LoadSolution(fSolution);
   res = Norm(residual); // pzfmatrix.h
}

int TPZEulerAnalysis::RunNewton(REAL & epsilon, int & numIter)
{
   int i = 0;
   REAL res;// residual of linear invertion.

// the lines below are to be performed by the external loop
//   fpLastSol->Zero();
//   fpCurrSol->Zero();

   Assemble(); // assembles the linear system.
   epsilon = fNewtonEps * 2.;// ensuring the loop will be
   // performed at least once.

   while(i < fNewtonMaxIter && epsilon > fNewtonEps)
   {
      //Solves the linearized system;
      Solve(res);

      // Linearizes the system with the newest iterative solution
      Assemble();

      epsilon = Norm(fRhs);
      i++;
   }

   fSolver->ResetMatrix(); // deletes the memory allocated
   // for the storage of the tangent matrix.

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

   REAL epsilon_Newton, epsilon_Global;
   int numIter_Newton, numIter_Global;
   REAL epsilon;

   int i;
   epsilon_Global = 2. * fTimeIntEps;

   ComputeTimeStep();

   while(i < fTimeIntMaxIter && epsilon_Global > fTimeIntEps)
   {

      RunNewton(epsilon_Newton, numIter_Newton);

      // Computing the time step, verifying the convergency
      // using the newest time step.
      ComputeTimeStep();

      fFlowCompMesh->Assemble(fRhs); // computing the residual only
      epsilon = Norm(fRhs);

      out << "\niter:" << i
          << " eps=" << epsilon
	  << " |NewtonEps=" << epsilon_Newton
	  << " nIter=" << numIter_Newton;

   }

   out.flush();

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

