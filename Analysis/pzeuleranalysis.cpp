//$Id: pzeuleranalysis.cpp,v 1.31 2004-05-21 13:32:35 erick Exp $

#include "pzeuleranalysis.h"
#include "pzerror.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "TPZParFrontStructMatrix.h"

TPZEulerAnalysis::TPZEulerAnalysis():
TPZAnalysis(), fFlowCompMesh(NULL),
fSolution2(), fRhsLast(), fpCurrSol(NULL),
fpLastSol(NULL), fpSolution(NULL),
fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100),
fEvolCFL(0), fpBlockDiag(NULL)
{

}

TPZEulerAnalysis::TPZEulerAnalysis(TPZFlowCompMesh *mesh, std::ostream &out):
TPZAnalysis(mesh, out), fFlowCompMesh(mesh),
fSolution2(), fRhsLast(), fpCurrSol(NULL), fpLastSol(NULL),
fpSolution(NULL), fLinSysEps(1e-10), fLinSysMaxIter(20),
fNewtonEps(1e-9),  fNewtonMaxIter(10),
fTimeIntEps(1e-8), fTimeIntMaxIter(100),
fEvolCFL(0), fpBlockDiag(NULL)
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

}

void TPZEulerAnalysis::BufferLastStateAssemble()
{
   fRhsLast.Zero();
   fRhsLast.Redim(fCompMesh->NEquations(),1);
   SetLastState();
   fFlowCompMesh->Assemble(fRhsLast);
   SetAdvancedState();
}

REAL TPZEulerAnalysis::EvaluateFluxEpsilon()
{
// enabling the flux evaluation type
   fFlowCompMesh->SetResidualType(Flux_RT);
   TPZFMatrix Flux;
   Flux.Redim(fCompMesh->NEquations(),1);
   Flux.Zero();
// contribution of last state
   SetLastState();
   fFlowCompMesh->Assemble(Flux);

// constribution of adv state
   SetAdvancedState();
   fFlowCompMesh->Assemble(Flux);

// reseting the residual type to residual
   fFlowCompMesh->SetResidualType(Residual_RT);

   return Norm(Flux);
}

#ifdef NOTDEFINED
void TPZEulerAnalysis::Assemble()
{
   if(!fCompMesh)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Computational Mesh\n";
      return;
      exit(-1);
   }

   if(!fStructMatrix)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
      exit(-1);
      return;
   }

   if(!fSolver)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Solver\n";
      exit(-1);
      return;
   }

   // redimensions and zeroes Rhs
   //fRhs.Redim(fCompMesh->NEquations(),1);
   fRhs = fRhsLast;

   TPZMatrix * pTangentMatrix = fSolver->Matrix();

   if(!pTangentMatrix || dynamic_cast<TPZParFrontStructMatrix <TPZFrontNonSym> *>(fStructMatrix))
   {
      pTangentMatrix = fStructMatrix->CreateAssemble(fRhs);
      fSolver->SetMatrix(pTangentMatrix);
   }
   else
   {
      if(!pTangentMatrix)
      {
         PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
         exit(-1);
         return;

      }

      pTangentMatrix->Zero();

      // Contributing referring to the advanced state
      // (n+1 index)
      fStructMatrix->Assemble(*pTangentMatrix, fRhs);
   }

   if(fpBlockDiag)
   {
      fpBlockDiag->Zero();
      //fpBlockDiag->SetIsDecomposed(0); // Zero already makes it
      fpBlockDiag->BuildFromMatrix(*pTangentMatrix);
   }

   // Contributing referring to the last state (n index)
   //fRhs+=/*.Add(fRhsLast, fRhsLast)*/ fRhsLast;

/*
ofstream Mout("Matriz.out");
ofstream Vout("Vetor.out");

   pTangentMatrix->Print("Matrix", Mout);//EMathematicaInput);

   fRhs.Print("Rhs", Vout);

   Mout.close();
   Vout.close();*/
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

#endif


void TPZEulerAnalysis::Assemble()
{
   if(!fCompMesh)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Computational Mesh\n";
      return;
      exit(-1);
   }

   if(!fStructMatrix)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
      exit(-1);
      return;
   }

   if(!fSolver)
   {
      PZError << "TPZEulerAnalysis::Assemble Error: No Solver\n";
      exit(-1);
      return;
   }

   // redimensions and zeroes Rhs
   //fRhs.Redim(fCompMesh->NEquations(),1);
   fRhs = fRhsLast;

   TPZMatrix * pTangentMatrix = fSolver->Matrix();

   if(!pTangentMatrix || dynamic_cast<TPZParFrontStructMatrix <TPZFrontNonSym> *>(fStructMatrix))
   {
      pTangentMatrix = fStructMatrix->CreateAssemble(fRhs);
      fSolver->SetMatrix(pTangentMatrix);
   }
   else
   {
      if(!pTangentMatrix)
      {
         PZError << "TPZEulerAnalysis::Assemble Error: No Structural Matrix\n";
         exit(-1);
         return;

      }

      pTangentMatrix->Zero();

      // Contributing referring to the advanced state
      // (n+1 index)
      fStructMatrix->Assemble(*pTangentMatrix, fRhs);
   }

   if(fpBlockDiag)
   {
      fpBlockDiag->Zero();
      //fpBlockDiag->SetIsDecomposed(0); // Zero already makes it
      fpBlockDiag->BuildFromMatrix(*pTangentMatrix);
   }

   // Contributing referring to the last state (n index)
   //fRhs+=/*.Add(fRhsLast, fRhsLast)*/ fRhsLast;

/*
ofstream Mout("Matriz.out");
ofstream Vout("Vetor.out");

   pTangentMatrix->Print("Matrix", Mout);//EMathematicaInput);

   fRhs.Print("Rhs", Vout);

   Mout.close();
   Vout.close();*/
}


void TPZEulerAnalysis::AssembleRhs()
{
   if(!fCompMesh) return;

   // Contributing referring to the last state (n index)
   fRhs = fRhsLast;

   // Contributing referring to the advanced state
   // (n+1 index)
   fFlowCompMesh->Assemble(fRhs);

/*
      if(!fCompMesh) return;

   // redimensions and zeroes Rhs
   fRhs.Redim(fCompMesh->NEquations(),1);

   // Contributing referring to the advanced state
   // (n+1 index)
   fFlowCompMesh->Assemble(fRhs);

   // Contributing referring to the last state (n index)
   fRhs+=fRhsLast;*/
}

//ofstream eulerout("Matrizes.out");

void TPZEulerAnalysis::Solve(REAL & res, TPZFMatrix * residual) {
   int numeq = fCompMesh->NEquations();
   if(fRhs.Rows() != numeq ) return;

   TPZFMatrix rhs(fRhs);

   TPZFMatrix delu(numeq,1);

   fSolver->Solve(rhs, delu);

/*
   fSolver->Matrix()->Print("Matriz", eulerout, EMathematicaInput);
   delu.Print("delu", eulerout, EMathematicaInput);
   fRhs.Print("Rhs", eulerout, EMathematicaInput);
   eulerout.flush();
*/
   UpdateSolution(delu);

   if(residual)
   {   // verifying the inversion of the linear system
      residual->Redim(numeq,1);
      fSolver->Matrix()->Residual(delu, fRhs, *residual);
      res = Norm(*residual);

      if(res > fLinSysEps)
      cout << "\n\tLinear system invertion did not achieve expected tolerance:" << res;
      cout.flush();
   }
}

int TPZEulerAnalysis::RunNewton(REAL & epsilon, int & numIter)
{

/*
   int i = 0;
   REAL res;// residual of linear invertion.

   Assemble(); // assembles the linear system.
   // It is assumed that BufferLastStateAssemble(); was already called.
   epsilon = fNewtonEps * 2.;// ensuring the loop will be
   // performed at least once.

   TPZFMatrix residual;

   while(i < fNewtonMaxIter && epsilon > fNewtonEps)
   {
      //Solves the linearized system and updates the solution.
      Solve(res, &residual);
//      cout << "\n      LinEpsilon:" << res;
      // Linearizes the system with the newest iterative solution
      Assemble();

      epsilon = Norm(fRhs);
      cout << "\n   NonLinEpsilon:" << epsilon;
      i++;
   }

   numIter = i;

   if(epsilon > fNewtonEps)return 0;
   return 1;
*/

   int i = 0;
   REAL res;// residual of linear invertion.

//   Assemble(); // assembles the linear system.
   // It is assumed that BufferLastStateAssemble(); was already called.
   epsilon = fNewtonEps * 2.;// ensuring the loop will be
   // performed at least once.

   TPZFMatrix residual;

   while(i < fNewtonMaxIter && epsilon > fNewtonEps)
   {
      // Linearizes the system with the newest iterative solution
      Assemble();
      //Solves the linearized system and updates the solution.
      Solve(res, &residual);
      //cout << "\n      LinEpsilon:" << res;

      AssembleRhs();

      epsilon = Norm(fRhs);
      cout << "\n   NonLinEpsilon:" << epsilon;
      i++;
   }

   numIter = i;

   if(epsilon > fNewtonEps)return 0;
   return 1;
}

TPZDXGraphMesh * TPZEulerAnalysis::PrepareDXMesh(ofstream &dxout, int dxRes)
{
  TPZVec<char *> scalar(4),vector(0);
  scalar[0] = "density";
  scalar[1] = "pressure";
  scalar[2] = "normvelocity";
  scalar[3] = "Mach";

  TPZMaterial * mat = fFlowCompMesh->GetFlowMaterial(0);
  int dim = mat->Dimension();
  //ResetReference(Mesh());//retira referências para criar graph consistente
  TPZDXGraphMesh * graph = new TPZDXGraphMesh (Mesh(),dim,mat,scalar,vector);
  //SetReference(Mesh());//recupera as referências retiradas
  //ofstream *dxout = new ofstream("ConsLaw.dx");
  //cout << "\nDX output file : ConsLaw.dx\n";
  graph->SetOutFile(dxout);
  //int resolution = 2;
  graph->SetResolution(dxRes);
  graph->DrawMesh(dim);

  return graph;
}

void TPZEulerAnalysis::Run(ostream &out, ofstream & dxout, int dxRes)
{
   // this analysis loop encloses several calls
   // to Newton's linearizations, updating the
   // time step.

   out << "\nBeginning time integration";

   clock_t startTime, endTime;
   startTime = clock();//time(&startTime);

   TPZDXGraphMesh * graph = PrepareDXMesh(dxout, dxRes);
   int numIterDX = 10000000;
   int outputStep = 0;

   REAL epsilon_Newton/*, epsilon_Global*/;
   int numIter_Newton/*, numIter_Global*/;
   REAL epsilon, lastEpsilon=0;

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

   // evaluates the time step based on the solution
   // (last Sol, since Curr sol equals it)
   ComputeTimeStep();

   // Buffers the contribution of the last state
   BufferLastStateAssemble();

   out << "\niter\teps(dU/dt)\tNewtonEps=\tnNewtonIter";


   while(i < fTimeIntMaxIter && epsilon > fTimeIntEps)
   {


      if(i%numIterDX==0)
      {
         graph->DrawSolution(outputStep, i);
	 graph->Out()->flush();
	 outputStep++;
      }

      // Solves the nonlinear system, updates the solution,
      // history and last state assemble buffer.
      RunNewton(epsilon_Newton, numIter_Newton);

      // resetting the forcing function for iterations greater than the 1st
      fFlowCompMesh->SetFlowforcingFunction(NULL);

      // updates the history of state variable vectors
      UpdateHistory();

      // buffering the contribution to the RHS
      BufferLastStateAssemble();

      // Computing the time step, verifying the convergency
      // using the newest time step.
      ComputeTimeStep();

      epsilon = EvaluateFluxEpsilon();

      // CFL control based on Flux reduction
      if((lastEpsilon>0.&&epsilon>0.) && fEvolCFL >= 1)
      {
         fFlowCompMesh->ScaleCFL(lastEpsilon/epsilon);
      }

      // CFL control based on Nonlinear invertion
      if(epsilon_Newton < fNewtonEps && fEvolCFL >= 2)
      {
         fFlowCompMesh->ScaleCFL(2.);
      }
      if(epsilon_Newton > fNewtonEps && fEvolCFL >= 3)
      {
         fFlowCompMesh->ScaleCFL(.5);
      }

      lastEpsilon = epsilon;

      out << "\n" << i
          << "\t" << epsilon
	  << "\t" << epsilon_Newton
	  << "\t" << numIter_Newton;

      cout << "\niter:" << i
           << "\t eps(dU/dt)=" << epsilon
	   << "\t |NewtonEps=" << epsilon_Newton
	   << "\t nIter=" << numIter_Newton;

      i++;
   }

   fSolver->ResetMatrix(); // deletes the memory allocated
   // for the storage of the tangent matrix.

   //fpCurrSol->Print("Solution", cout);

   graph->DrawSolution(outputStep, i);
   graph->Out()->flush();
   outputStep++;

   //graph->Close();
   delete graph;

   endTime = clock();//time(&endTime);
   const double CPS = CLOCKS_PER_SEC;

   out  << "\nElapsed time in seconds: " << ((double)(endTime - startTime))/CPS << endl;

   cout << "\nElapsed time in seconds: " << ((double)(endTime - startTime))/CPS << endl;

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

void TPZEulerAnalysis::SetEvolCFL(int EvolCFL)
{
   fEvolCFL = EvolCFL;
}

void TPZEulerAnalysis::SetBlockDiagonalPrecond(TPZBlockDiagonal * blockDiag)
{
   fpBlockDiag = blockDiag;
}

void TPZEulerAnalysis::WriteCMesh( const char * str)
{
   TPZFileStream fstr;
   fstr.OpenWrite(str);
   fGeoMesh->Write(fstr,1);
   fFlowCompMesh->Write(fstr,1);
}
