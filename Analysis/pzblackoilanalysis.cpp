//$Id: pzblackoilanalysis.cpp,v 1.1 2008-11-12 12:46:03 fortiago Exp $

#include "pzblackoilanalysis.h"
#include "pzblackoil2p3d.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"

using namespace std;

TPZBlackOilAnalysis::TPZBlackOilAnalysis(TPZCompMesh *mesh, double TimeStep, std::ostream &out):TPZNonLinearAnalysis(mesh,out){
  this->fTimeStep = TimeStep;
  this->SetConvergence(0, 0.);
  this->SetNewtonConvergence(0, 0.);
  this->SetInitialSolutionAsZero();
  this->SetSaveFrequency(0,0);
}

TPZBlackOilAnalysis::~TPZBlackOilAnalysis(){
///nothing to be done here
}

void TPZBlackOilAnalysis::SetInitialSolution(TPZFMatrix & InitialSol){
  const int nrows = this->Mesh()->Solution().Rows();
  const int ncols = this->Mesh()->Solution().Cols();
  if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
    PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
  }
  else{
    this->fSolution = InitialSol;
  }
}

void TPZBlackOilAnalysis::SetInitialSolutionAsZero(){
  TPZFMatrix & MeshSol = this->Mesh()->Solution();
  this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
  this->fSolution.Zero();
}

void TPZBlackOilAnalysis::Run(std::ostream &out, bool linesearch){

  this->SetAllMaterialsDeltaT();

  TPZFMatrix prevsol, laststate, lastsol;
  for(this->fCurrentStep = 0; this->fCurrentStep < this->fNIter; this->fCurrentStep++){

    ///Computing residual of last state solution
    this->SetLastState();
    this->Assemble();
    laststate = this->fRhs;
    prevsol = fSolution;
    lastsol = fSolution;
    ///Newton's method
    this->SetCurrentState();
    REAL error = this->fNewtonTol * 2. + 1.;    
    int iter = 0;
    while(error > this->fNewtonTol && iter < this->fNewtonMaxIter){

      fSolution.Redim(0,0);
      this->Assemble();
      this->fRhs += laststate;
      this->Solve();

      if (linesearch){
        TPZFMatrix nextSol;
        REAL LineSearchTol = 1e-3 * Norm(fSolution);
        const int niter = 10;
        this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
        fSolution = nextSol;
      }
      else{
        fSolution += prevsol;
      }

      prevsol -= fSolution;
      REAL norm = Norm(prevsol);
      out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << endl;

      prevsol = fSolution;
      TPZAnalysis::LoadSolution();

      error = norm;
      iter++;
   }///Newton's iterations

   prevsol = fSolution;

   if (this->fSaveFrequency){
    if (!(this->fCurrentStep % fSaveFrequency)){
      this->PostProcess(this->fDXResolution);
    }
   }

   prevsol -= lastsol;
   REAL steadynorm = Norm(prevsol);
   std::cout << "*********** Steady state error at iteration " << this->fCurrentStep << " = " << steadynorm << "\n\n";
   if (!fForceAllSteps){
    if (steadynorm < this->fSteadyTol){
      std::cout << "Steady state solution achieved\n\n";
      this->fNIter = this->fCurrentStep;
      break;
    }
   }
   std::cout.flush();

  }///time step iterations

}///method


void TPZBlackOilAnalysis::SetLastState(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second.operator->());
    if (blackoilmat){
      blackoilmat->SetLastState();
    }
  }
}


void TPZBlackOilAnalysis::SetCurrentState(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second.operator->());
    if (blackoilmat){
      blackoilmat->SetCurrentState();
    }
  }
}

void TPZBlackOilAnalysis::SetAllMaterialsDeltaT(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZBlackOil2P3D * blackoilmat = dynamic_cast< TPZBlackOil2P3D *>(matit->second.operator->());
    if (blackoilmat){
      blackoilmat->SetTimeStep(this->TimeStep());
    }
  }
}
 

void TPZBlackOilAnalysis::PostProcess(int resolution, int dimension){
    REAL T = this->fCurrentStep * this->TimeStep();
    this->fTime = T;
    TPZAnalysis::PostProcess(resolution, dimension);
}//method


void TPZBlackOilAnalysis::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
    REAL T = this->fCurrentStep * this->TimeStep();
    out << "\nSOLUTION #" << this->fCurrentStep << " AT TIME = " << T << std::endl;
    TPZAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
}//method

void TPZBlackOilAnalysis::Assemble(){
  if(!fCompMesh || !fStructMatrix || !fSolver){
    cout << "TPZBlackOilAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
         << " fStructMatrix " << (bool) fStructMatrix << " fSolver " << (bool) fSolver << " at file " 
         << __FILE__ << " line " << __LINE__ << endl;
    return;
  }

  int sz = fCompMesh->NEquations();
  fRhs.Redim(sz,1);

  ///
  bool exist = false;
  if(fSolver->Matrix()) if (fSolver->Matrix()->Rows()==sz) exist = true;
  if (exist){
    fSolver->Matrix()->Zero();
    fStructMatrix->Assemble(fSolver->Matrix(),fRhs);
  }
  else{
    TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs);
    fSolver->SetMatrix(mat);
  }
  fSolver->UpdateFrom(fSolver->Matrix());
}

void TPZBlackOilAnalysis::SetConvergence(int niter, REAL eps, bool ForceAllSteps){
  this->fNIter = niter;
  this->fSteadyTol = eps;
  this->fForceAllSteps = ForceAllSteps;
}

void TPZBlackOilAnalysis::SetSaveFrequency(int SaveFrequency, int resolution){
  this->fSaveFrequency = SaveFrequency;
  this->fDXResolution = resolution;
}

void TPZBlackOilAnalysis::SetNewtonConvergence(int niter, REAL eps){
  this->fNewtonMaxIter = niter;
  this->fNewtonTol = eps;
}

REAL & TPZBlackOilAnalysis::TimeStep(){
  return this->fTimeStep;
}

