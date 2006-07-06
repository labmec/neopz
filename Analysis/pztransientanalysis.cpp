// -*- c++ -*-

//$Id: pztransientanalysis.cpp,v 1.2 2006-07-06 15:59:09 tiago Exp $

#include "pztransientanalysis.h"
#include "pztransientmat.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"

#include "pznonlinearpoisson3d.h"
#define TRANSIENTCLASS TPZNonLinearPoisson3d

using namespace std;

double TPZTransientAnalysis::gTime = 0.;

TPZTransientAnalysis::TPZTransientAnalysis(TPZCompMesh *mesh, bool IsLinear, std::ostream &out):TPZAnalysis(mesh,out),fLastState(){

  this->fTimeStep = 0.;
  this->fCurrentIter = 0;
  this->SetConvergence(0, 0.);
  this->SetNewtonConvergence(0, 0.);
  this->SetInitialSolutionAsZero();
  this->fIsLinearProblem = IsLinear;
}

TPZTransientAnalysis::~TPZTransientAnalysis(){

}

void TPZTransientAnalysis::SetInitialSolution(TPZFMatrix & InitialSol){
  const int nrows = this->Mesh()->Solution().Rows();
  const int ncols = this->Mesh()->Solution().Cols();
  if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
    PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
  }
  else{
    this->fAllSolutions[0] = InitialSol;
  }
}

void TPZTransientAnalysis::SetInitialSolutionAsZero(){
  TPZFMatrix & MeshSol = this->Mesh()->Solution();
  this->GetSolution(0).Redim( MeshSol.Rows(), MeshSol.Cols() );
  this->GetSolution(0).Zero();
}

void TPZTransientAnalysis::Run(std::ostream &out){

  TPZTransientAnalysis::gTime = 0.;

  for(int i = 1; i < this->fAllSolutions.NElements(); i++){
    this->GetSolution(i).Zero();
  }//i
  this->SetAllMaterialsDeltaT();
  
  if (this->fIsLinearProblem){
    this->fSolution = this->GetSolution(0);
    this->ComputeLinearTangentMatrix();  
  }  

  TPZFMatrix prevsol;
  for(this->fCurrentIter = 0; this->fCurrentIter < this->fNIter; this->fCurrentIter++){
  
    TPZTransientAnalysis::gTime = this->TimeStep() * (this->fCurrentIter+1);

    //Computing residual of last state solution
    prevsol = this->GetSolution(this->fCurrentIter);
    this->fSolution = prevsol;
    this->SetLastState();
    this->Assemble();
    this->fLastState = this->fRhs;
    
    //Newton's method
    this->SetCurrentState();        
    REAL error = this->fNewtonTol * 2. + 1.;    
    int iter = 0;
    while(error > this->fNewtonTol && iter < this->fNewtonMaxIter) {

      fSolution.Redim(0,0);
      this->Assemble();
      this->fRhs += this->fLastState;
      this->Solve();

      REAL norm = Norm(fSolution);
      out << "Iteracao n : " << (iter+1) << " : norma da solucao |Delta(Un)|: " << norm << std::endl;

      fSolution += prevsol;
      prevsol = fSolution;
      TPZAnalysis::LoadSolution();
//       if(norm < this->fNewtonTol) {
//          out << "Tolerancia atingida na iteracao : " << (iter+1) << endl;
//          out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
// 
//       } else
//       if( (norm - error) > 1.e-9 ) {
//          out << "\nDivergent Method or Stalled\n";
//       }
      error = norm;
      iter++;
   }//Newton's iterations
   
   this->GetSolution(fCurrentIter+1) = fSolution;
   
   fSolution -= this->GetSolution(fCurrentIter);
   REAL steadynorm = Norm(fSolution);
   std::cout << "*********** Steady state error at iteration " << (this->fCurrentIter+1) << " = " << steadynorm << "\n\n";
   if (steadynorm < this->fSteadyTol){
    std::cout << "Steady state solution achieved\n\n";
    this->fNIter = fCurrentIter+1;
    break;
   }
   
   
   
  }//time iterations
  

  
}//method

void TPZTransientAnalysis::SetLastState(){
  TPZCompMesh * mesh = this->Mesh();
  const int nmat = mesh->MaterialVec().NElements();  
  for(int i = 0; i < nmat; i++){
    TPZMaterial * mat = mesh->MaterialVec()[i];
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(mat);
    if (trans){
      trans->SetLastState();
    }
  }
}
  
void TPZTransientAnalysis::SetCurrentState(){
  TPZCompMesh * mesh = this->Mesh();
  const int nmat = mesh->MaterialVec().NElements();  
  for(int i = 0; i < nmat; i++){
    TPZMaterial * mat = mesh->MaterialVec()[i];
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(mat);
    if (trans){
      trans->SetCurrentState();
    }
  }
}
  
void TPZTransientAnalysis::SetAllMaterialsDeltaT(){
  TPZCompMesh * mesh = this->Mesh();
  const int nmat = mesh->MaterialVec().NElements();  
  for(int i = 0; i < nmat; i++){
    TPZMaterial * mat = mesh->MaterialVec()[i];
    TPZTransientMaterial< TRANSIENTCLASS > * trans = dynamic_cast<TPZTransientMaterial< TRANSIENTCLASS > *>(mat);
    if (trans){
      trans->SetTimeStep( this->TimeStep() );
    }
  }
}
 
void TPZTransientAnalysis::PostProcess(int resolution, int dimension){
  for(int it = 0; it <= this->fNIter; it++){
    REAL T = it * this->TimeStep();
    this->fTime = T;
    this->fSolution = this->GetSolution(it);
    TPZAnalysis::LoadSolution();
    TPZAnalysis::PostProcess(resolution, dimension);
  }//it
}//method
  
void TPZTransientAnalysis::PostProcess(TPZVec<REAL> &loc, std::ostream &out){
  for(int it = 0; it <= this->fNIter; it++){
    REAL T = it * this->TimeStep();
    this->gTime = T;
    this->fSolution = this->GetSolution(it);
    TPZAnalysis::LoadSolution();
    out << "\nSOLUTION #" << it << " AT TIME = " << T << std::endl;
    TPZAnalysis::PostProcess(loc, out);
    out << "\n***************************************\n" << std::endl;
  }//it
}//method

void TPZTransientAnalysis::Assemble(){
  if(!fCompMesh || !fStructMatrix || !fSolver){
    cout << "TPZTransientAnalysis::Assemble lacking definition for Assemble fCompMesh "<< (void *) fCompMesh 
         << " fStructMatrix " << (void *) fStructMatrix << " fSolver " << (void *) fSolver << " at file " 
         << __FILE__ << " line " << __LINE__ << endl;
    return;
  }
  
  int sz = fCompMesh->NEquations();
  fRhs.Redim(sz,1);
  if(fSolver->Matrix() && fSolver->Matrix()->Rows()==sz){
    if (fIsLinearProblem){
      this->Mesh()->Assemble(fRhs);
    }
    else{
      fSolver->Matrix()->Zero();
      fStructMatrix->Assemble(*(fSolver->Matrix()),fRhs);
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

void TPZTransientAnalysis::ComputeLinearTangentMatrix(){
  if (!fIsLinearProblem) return;      
  this->SetCurrentState();
  const int sz = this->Mesh()->NEquations();
  fRhs.Redim(sz,1);
  TPZMatrix *mat = fStructMatrix->CreateAssemble(fRhs);
  fSolver->SetMatrix(mat);
}//method
