// -*- c++ -*-

//$Id: pztransientanalysis.cpp,v 1.8 2009-05-06 20:07:01 fortiago Exp $

#include "pztransientanalysis.h"
#include "pztransientmat.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"
#include "checkconv.h"

using namespace std;

template<class TRANSIENTCLASS>
double TPZTransientAnalysis<TRANSIENTCLASS>::gTime = 0.;

template<class TRANSIENTCLASS>
TPZTransientAnalysis<TRANSIENTCLASS>::TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear, std::ostream &out):/*TPZAnalysis*/TPZNonLinearAnalysis(mesh,out), fSavedSolutionVec(){
  this->fTimeStep = 0.;
  this->fCurrentIter = 0;
  this->SetConvergence(0, 0.);
  this->SetNewtonConvergence(0, 0.);
  this->SetInitialSolutionAsZero();
  this->fIsLinearProblem = IsLinear;
  this->SetSaveFrequency(0,0);
  this->fSaveSolutionVecFrequency = 0; 
}

template<class TRANSIENTCLASS>
TPZTransientAnalysis<TRANSIENTCLASS>::~TPZTransientAnalysis(){
  fSavedSolutionVec.clear();
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetInitialSolution(TPZFMatrix & InitialSol){
  const int nrows = this->Mesh()->Solution().Rows();
  const int ncols = this->Mesh()->Solution().Cols();
  if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
    PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
  }
  else{
    this->fSolution = InitialSol;
  }
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetInitialSolutionAsZero(){
  TPZFMatrix & MeshSol = this->Mesh()->Solution();
  this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
  this->fSolution.Zero();
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::Run(std::ostream &out, bool FromBegining, bool linesearch){

  this->SetImplicit();

  if (FromBegining){
    this->fCurrentIter = 0;
    this->fSavedSolutionVec.clear();
  }

  {
  TPZVec<REAL> coefs(1,1.);
  TPZFMatrix cpSol(fSolution);
  TPZFMatrix range(fCompMesh->NEquations(),1,1.);
  this->SetLastState();
  CheckConvergence(*this,cpSol,range,coefs);
  this->SetCurrentState();
  CheckConvergence(*this,cpSol,range,coefs);
  }

  this->LoadSolution(fSolution);
  
  TPZTransientAnalysis::gTime =this->GetCurrentIter() * this->TimeStep();
//   this->PostProcess(this->fDXResolution);
  
  this->SetAllMaterialsDeltaT();
  
  if (this->fIsLinearProblem){
    this->ComputeLinearTangentMatrix();  
  }

  TPZFMatrix prevsol, laststate, lastsol;
  for( ; this->fCurrentIter < this->fNIter;){
  
    this->fCurrentIter++;
    TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;

    //Computing residual of last state solution
//     this->fSolution = prevsol;
    this->SetLastState();
    this->Assemble();
    laststate = this->fRhs;
    prevsol = fSolution;
    lastsol = fSolution;
    //Newton's method
    this->SetCurrentState();        
    REAL error = this->fNewtonTol * 2. + 1.;    
    int iter = 0;
    while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {

      fSolution.Redim(0,0);
      this->Assemble();
      this->fRhs += laststate;

/*    for(int i = 0; i < this->Solver().Matrix()->Rows(); i++)
      for(int j = 0; j < this->Solver().Matrix()->Cols(); j++)
        this->Solver().Matrix()->operator ( )(i,j) *= 1e6;
    this->Rhs().operator*=( 1e6 );*/

    /*{ofstream file("new.nb");
    this->Solver().Matrix()->Print("rigidezNew = ", file, EMathematicaInput);
    this->Rhs().Print("rhsNew = ", file);
    file.flush();}*/

      this->Solve();

      if (linesearch){
        TPZFMatrix nextSol;
        REAL LineSearchTol = 1e-3 * Norm(fSolution);
        const int niter = 100;
        this->LineSearch(prevsol, fSolution, nextSol, LineSearchTol, niter);
        fSolution = nextSol;
      }
      else{
        fSolution += prevsol;
      }

      prevsol -= fSolution;
      REAL norm = Norm(prevsol);
      out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << std::endl;

      prevsol = fSolution;
      TPZAnalysis::LoadSolution();

      error = norm;
      iter++;
   }//Newton's iterations
   
   prevsol = fSolution;
   
   if (this->fSaveFrequency){
    if (!(this->fCurrentIter % fSaveFrequency)){
      this->PostProcess(this->fDXResolution);
    }
   }
   this->SaveCurrentSolutionVec();
   
   prevsol -= lastsol;
   REAL steadynorm = Norm(prevsol);
   std::cout << "*********** Steady state error at iteration " << this->fCurrentIter << " = " << steadynorm << "\n\n";
   if (!fForceAllSteps){
    if (steadynorm < this->fSteadyTol){
      std::cout << "Steady state solution achieved\n\n";
      this->fNIter = fCurrentIter;
      break;
    }
   }
   std::cout.flush();   
   
   
  }//time iterations
  

  
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetImplicit(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetImplicit();
    }
  }
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetExplicit(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetExplicit();
    }
  }
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetLastState(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetLastState();
    }
  }
}
  
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetCurrentState(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetCurrentState();
    }
  }
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetMassMatrix(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetMassMatrix();
    }
  }
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetFluxOnly(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetFluxOnly();
    }
  }
}
  
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetAllMaterialsDeltaT(){
  TPZCompMesh * mesh = this->Mesh();
  std::map<int, TPZAutoPointer<TPZMaterial> >::iterator matit;
  for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
  {
    if(!matit->second) continue;
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(matit->second.operator->());
    if (trans){
      trans->SetTimeStep(this->TimeStep());
    }
  }
}
 
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::PostProcess(int resolution, int dimension){
    REAL T = this->GetCurrentIter() * this->TimeStep();
    this->fTime = T;
    TPZAnalysis::PostProcess(resolution, dimension);
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
    REAL T = this->GetCurrentIter() * this->TimeStep();
    this->gTime = T;
    out << "\nSOLUTION #" << this->GetCurrentIter() << " AT TIME = " << T << std::endl;
    TPZAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
}//method

  #include "TPZParFrontStructMatrix.h"
  #include "pzstepsolver.h"
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::Assemble(){
  if(!fCompMesh || !fStructMatrix || !fSolver){
    cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
         << " fStructMatrix " << (bool) fStructMatrix << " fSolver " << (bool) fSolver << " at file " 
         << __FILE__ << " line " << __LINE__ << endl;
    return;
  }
  
  int sz = fCompMesh->NEquations();
  fRhs.Redim(sz,1);
  
  ///pelo amor de Deus
//   delete fSolver->Matrix();
//   fSolver->ResetMatrix();
//   delete this->fStructMatrix;
//   fStructMatrix = NULL;
//   delete fSolver;
//   fSolver = NULL;
//   TPZCompMesh * cmesh = this->Mesh();
//   TPZParFrontStructMatrix<TPZFrontNonSym>  /*TPZFStructMatrix*/ StrMat(cmesh);
//   this->SetStructuralMatrix(StrMat);
//   TPZStepSolver step;
//   step.SetDirect(ELU);
//   this->SetSolver(step);
  ///
  bool exist = false;
  if(fSolver->Matrix()) if (fSolver->Matrix()->Rows()==sz) exist = true;
  if (exist){
    if (fIsLinearProblem){
      TPZStructMatrix::Assemble(fRhs, *Mesh());
    }
    else{
      fSolver->Matrix()->Zero();
      fStructMatrix->Assemble(fSolver->Matrix(),fRhs);
    }
  }
  else{
    if (this->fIsLinearProblem){
      std::cout << __PRETTY_FUNCTION__ << " @ " << __LINE__ << " Error! StrMatrix must be created using" 
                                                            << " methodTPZTransientAnalysis::ComputeLinearTangentMatrix()"
                                                            << " when (this->fIsLinearProblem == true)\n";
    }
    TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs);
    fSolver->SetMatrix(mat);
  }
  fSolver->UpdateFrom(fSolver->Matrix());
}

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeLinearTangentMatrix(){
  if (!fIsLinearProblem) return;      
  this->SetCurrentState();
  const int sz = this->Mesh()->NEquations();
  fRhs.Redim(sz,1);
  TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs);
  fSolver->SetMatrix(mat);
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeMassMatrix(){
  this->SetMassMatrix();
  const int sz = this->Mesh()->NEquations();
  fRhs.Redim(sz,1);
  TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs);
  fSolver->SetMatrix(mat);
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::ComputeFluxOnly(){
  if(!fCompMesh || !fStructMatrix || !fSolver){
    cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
         << " fStructMatrix " << (bool) fStructMatrix << " fSolver " << (bool) fSolver << " at file " 
         << __FILE__ << " line " << __LINE__ << endl;
    return;
  }
  
  this->SetFluxOnly();  
  int sz = fCompMesh->NEquations();
  fRhs.Redim(sz,1);
  if(fSolver->Matrix() && fSolver->Matrix()->Rows()==sz){
    TPZStructMatrix::Assemble(fRhs, *Mesh());
  }//if
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::RunExplicit(std::ostream &out, bool FromBegining){

  this->SetExplicit();

  if (FromBegining){
    this->fCurrentIter = 0;
    this->fSavedSolutionVec.clear();
  }
  TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;
  this->PostProcess(this->fDXResolution);

  this->SetAllMaterialsDeltaT();

  TPZFMatrix prevsol;
  for( this->fCurrentIter++ ; this->fCurrentIter < this->fNIter; this->fCurrentIter++){

    this->ComputeMassMatrix();

    TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;

    this->SetFluxOnly();

    //Computing residual of last state solution
    prevsol = fSolution;
    TPZAnalysis::LoadSolution();
    this->ComputeFluxOnly();

    /*{ofstream file("new.txt");
    this->Solver().Matrix()->Print("rigidezNew = ", file);
    this->Rhs().Print("rhsNew = ", file);
    file.flush();}
    exit(-1);*/

    this->Solve();
    //now fSolution = deltaSol
    fSolution += prevsol;
    
    TPZAnalysis::LoadSolution();
    if (this->fSaveFrequency){
      if (!(this->fCurrentIter % fSaveFrequency)){
        this->PostProcess(this->fDXResolution);
      }
    }
    this->SaveCurrentSolutionVec();
   
    prevsol -= fSolution;
    REAL steadynorm = Norm(prevsol);
    std::cout << "*********** Steady state error at iteration " << (this->fCurrentIter) << " = " << steadynorm << "\n\n";
    if (!fForceAllSteps){
      if (steadynorm < this->fSteadyTol){
        std::cout << "Steady state solution achieved\n\n";
        this->fNIter = fCurrentIter;
        break;
      }
    }
    std::cout.flush();   
   
  }//time iterations
  
}//method

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SetSaveSolution(int SaveFrequency){
  this->fSaveSolutionVecFrequency = SaveFrequency;
}

template<class TRANSIENTCLASS>
std::list< std::pair<TPZFMatrix, REAL> > & TPZTransientAnalysis<TRANSIENTCLASS>::GetSavedSolutions(){
  return this->fSavedSolutionVec;
}

#include <sstream>
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::SaveCurrentSolutionVec(){
  if(!this->fSaveSolutionVecFrequency) return;
  if(this->fCurrentIter % this->fSaveSolutionVecFrequency == 0){
    std::pair< TPZFMatrix, REAL > mypair;
    mypair.first = this->Solution();
    mypair.second = TPZTransientAnalysis::gTime;
    this->fSavedSolutionVec.push_back(mypair);

   ofstream file("currentsol.txt");
   stringstream mess; mess << "sol( " << TPZTransientAnalysis::gTime << " ) = ";
   this->Solution().Print(mess.str().c_str(), file);

  }
}

#include "pztransrichardsmat.h"

template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::RunDualTimeStep(REAL PseudoTimeStep, std::ostream &out, bool FromBegining){

  this->SetImplicit();

  if (FromBegining){
    this->fCurrentIter = 0;
    this->fSavedSolutionVec.clear();
  }

  this->LoadSolution(fSolution);

  TPZTransientAnalysis::gTime =this->GetCurrentIter() * this->TimeStep();
//   this->PostProcess(this->fDXResolution);

  const REAL TimeStep = this->TimeStep();
  this->SetAllMaterialsDeltaT();

  /// Auxiliar vectors. They store multiplier coefficients
  TPZFMatrix ActualUn, ActualUnPlus, DTSUn, DTSUnPlus;

  /// Auxiliar vector. It stores integral values
  TPZFMatrix laststateRes, /*IntegDTSUn,*/ Residual;

  for( ; this->fCurrentIter < this->fNIter ; ){

    this->TimeStep() = TimeStep;
    this->SetAllMaterialsDeltaT();

    /// Storing last actual solution
    ActualUn = this->fSolution;

    this->fCurrentIter++;
    TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;

    /// Computing residual of last state solution
    this->SetLastState();
    this->Assemble();
    laststateRes = this->fRhs;

    /// Dual Time Step
    DTSUn = ActualUn;
    DTSUnPlus = ActualUn;
    REAL error = this->fNewtonTol * 2. + 1.;
    int iter = 0;
    while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {

      /// Computing residual
      this->TimeStep() = TimeStep;
      this->SetAllMaterialsDeltaT();
      this->SetCurrentState();
      fSolution.Redim(0,0);
      this->Assemble();
      this->fRhs += laststateRes;
      Residual = this->fRhs;

      /// Integrating dual time step vector DTSUn
      this->TimeStep() = PseudoTimeStep;
      this->SetAllMaterialsDeltaT();

/*      this->fSolution = DTSUn;
      TPZAnalysis::LoadSolution();
      this->SetLastState();
      this->fRhs.Zero();
      this->Assemble();
      IntegDTSUn = this->fRhs;
*/
      /// Now, we perform DTSUnPlus/PseudoTimeStep = DTSUn/PseudoTimeStep + Residual
        /// DTSUnPlus matrix
      TPZTransientRichardsMaterial::gPureMassMatrix = true;
      this->SetMassMatrix();
      this->fSolver->Matrix()->Zero();      
      this->Assemble();
      TPZTransientRichardsMaterial::gPureMassMatrix = false;

        /// Computing DTSUnPlus
      this->fRhs = Residual;
      this->Solve();
      DTSUnPlus = this->fSolution;
      DTSUnPlus += DTSUn;

      /// Updating values
      this->fSolution = DTSUnPlus;
      TPZAnalysis::LoadSolution();

      /// checking DTS convergence
      DTSUn -= DTSUnPlus;
      REAL norm = Norm(DTSUn);
      REAL normRes = Norm(Residual);
      out << "\tIteracao n : " << (iter+1) << " : norma de delta U |deltaU|: " << norm 
          << " Norma do Rhs: " << normRes << std::endl;
      error = norm;
      iter++;

      /// Updating values
      DTSUn = DTSUnPlus;

   }///DualTimeStep iterations

   /// Getting DTS result
   ActualUnPlus = this->fSolution;

   /// Saving results
   if (this->fSaveFrequency){
    if (!(this->fCurrentIter % fSaveFrequency)){
      this->PostProcess(this->fDXResolution);
    }
   }
   this->SaveCurrentSolutionVec();

   /// Checking steady state convergence
   ActualUn -= ActualUnPlus;
   REAL steadynorm = Norm(ActualUn);
   std::cout << "*********** Steady state error at iteration " << this->fCurrentIter << " = " << steadynorm << "\n\n";
   if (!fForceAllSteps){
    if (steadynorm < this->fSteadyTol){
      std::cout << "Steady state solution achieved\n\n";
      this->fNIter = fCurrentIter;
      break;
    }
   }
   std::cout.flush();


  }/// time iterations

}/// method

/*
template<class TRANSIENTCLASS>
void TPZTransientAnalysis<TRANSIENTCLASS>::RunDualTimeStep(REAL PseudoTimeStep, std::ostream &out, bool FromBegining){

  this->SetImplicit();

  if (FromBegining){
    this->fCurrentIter = 0;
    this->fSavedSolutionVec.clear();
  }

  this->LoadSolution(fSolution);

  TPZTransientAnalysis::gTime =this->GetCurrentIter() * this->TimeStep();
//   this->PostProcess(this->fDXResolution);

  this->SetAllMaterialsDeltaT();

  /// Auxiliar vectors. They store multiplier coefficients
  TPZFMatrix ActualUn, ActualUnPlus, DTSUn, DTSUnPlus;

  /// Auxiliar vector. It stores integral values
  TPZFMatrix laststateRes, Residual;

  for( ; this->fCurrentIter < this->fNIter ; ){

    /// Storing last actual solution
    ActualUn = this->fSolution;

    this->fCurrentIter++;
    TPZTransientAnalysis::gTime = this->TimeStep() * this->fCurrentIter;

    /// Computing residual of last state solution
    this->SetLastState();
    this->Assemble();
    laststateRes = this->fRhs;

    /// Dual Time Step
    DTSUn = ActualUn;
    DTSUnPlus = ActualUn;
    REAL error = this->fNewtonTol * 2. + 1.;
    int iter = 0;
    while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {

      /// Computing residual
      this->SetCurrentState();
      fSolution.Redim(0,0);
      this->Assemble();
      this->fRhs += laststateRes;
      Residual = this->fRhs;

      /// Now, we perform DTSUnPlus = DTSUn - PseudoTimeStep * Residual
      DTSUnPlus = Residual;
      DTSUnPlus *= -1.;
      DTSUnPlus *= PseudoTimeStep;
      DTSUnPlus += DTSUn;

      /// Updating values
      this->fSolution = DTSUnPlus;
      TPZAnalysis::LoadSolution();
      DTSUn = DTSUnPlus;

      /// checking DTS convergence
      REAL norm = Norm(this->fRhs);
      out << "\tIteracao n : " << (iter+1) << " : norma do residuo |RHS|: " << norm << std::endl;
      error = norm;
      iter++;

   }///DualTimeStep iterations

   /// Getting DTS result
   ActualUnPlus = this->fSolution;

   /// Saving results
   if (this->fSaveFrequency){
    if (!(this->fCurrentIter % fSaveFrequency)){
      this->PostProcess(this->fDXResolution);
    }
   }
   this->SaveCurrentSolutionVec();

   /// Checking steady state convergence
   ActualUn -= ActualUnPlus;
   REAL steadynorm = Norm(ActualUn);
   std::cout << "*********** Steady state error at iteration " << this->fCurrentIter << " = " << steadynorm << "\n\n";
   if (!fForceAllSteps){
    if (steadynorm < this->fSteadyTol){
      std::cout << "Steady state solution achieved\n\n";
      this->fNIter = fCurrentIter;
      break;
    }
   }
   std::cout.flush();


  }/// time iterations

}/// method
*/
///instantiations
#include "pzpoisson3d.h"
template class TPZTransientAnalysis< TPZMatPoisson3d >;

#include "pznonlinearpoisson3d.h"
template class TPZTransientAnalysis< TPZNonLinearPoisson3d >;

#include "pzburger.h"
template class TPZTransientAnalysis< TPZBurger >;

#include "pzrichardsequation.h"
template class TPZTransientAnalysis< TPZRichardsEquation >;
