#include "pzeuleranalysis.h"

TPZEulerAnalysis::TPZEulerAnalysis():TPZAnalysis(), fFlowCompMesh(NULL),
fSolution2(), fpCurrSol(NULL), fpLastSol(NULL), fpSolution(NULL)
{

}

TPZEulerAnalysis::TPZEulerAnalysis(TPZFlowCompMesh *mesh, std::ostream &out = cout):
TPZAnalysis(mesh, out):fFlowCompMesh(mesh), fSolution2(), fpCurrSol(NULL), fpLastSol(NULL),
fpSolution(NULL)
{

}

TPZEulerAnalysis::~TPZEulerAnalysis()
{

}

virtual void TPZEulerAnalysis::Run(ostream &out = cout)
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
   int i, NumFluid;
   NumFluid = fFluidMaterial.NElements();
   for(i = 0; i < NumFluid; i++)
   {
      fFluidMaterial[i]->SetContributionTime(time);
   }
}

void TPZEulerAnalysis::UpdateSolution(TPZFMatrix & deltaSol)
{
   (*fpCurrSol) = (*pfLastSol);
   fpCurrSol->Add(1.0, deltaSol);
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
   TPZMatrix * pTangentMatrix = StructMatrix->Create();

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
   int i = 0, maxIter = numIterl;
   REAL res; // norm of residual of linear invertion.
   REAL resNorm;// norm of residual of Euler

   *fpLastSol.Zero();
   *fpCurrSol.Zero();

   Assemble(); // assembles the linear system.

   while(i < maxIter && normRes > epsilon)
   {
      //Solves the linearized system;
      Solve(res);

      // Linearizes the system with the newest iterative solution
      Assemble();

      normRes = Norm(fRhs);
      i++;
   }

   fSolver->ResetMatrix(); // deletes the memory allocated
   // for the storage of the tangent matrix.

   epsilon = normRes;

   if(i==maxIter)return 0;
   maxIter = i;
   return 1;
}

void TPZEulerAnalysis::Run()
{
   REAL epsilon_Newton, epsilon_Global;
   int numIter_Newton, numIter_Global;

???// Como calcular o resíduo para poder avaliar convergência quando deltaT é atualizado??? sem muito esforço computacional??

}

void TPZEulerAnalysis::SetTimeStep()
{
  TPZFlowCompMesh *fm  = dynamic_cast<TPZFlowCompMesh *>(fCompMesh);

  if(!fm){
     PZErr << "\nInvalid CompMesh type -> TPZFlowCompMesh expected\n";
  }

  REAL maxveloc = fm->MaxVelocityOfMesh();
  REAL deltax = fm->LesserEdgeOfMesh();// deltax

  TPZCompElDisc *disc;
  int degree = disc->gDegree;
  TPZConservationLaw2 *law = dynamic_cast<TPZConservationLaw2 *>(mat);
  law->SetDeltaTime(maxveloc,deltax,degree);
}
