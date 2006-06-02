// -*- c++ -*-

//$Id: pztransientanalysis.cpp,v 1.1 2006-06-02 17:03:27 tiago Exp $

#include "pztransientanalysis.h"
#include "pztransientmat.h"
#include "pzpoisson3d.h"

TPZTransientAnalysis::TPZTransientAnalysis(TPZCompMesh *mesh,std::ostream &out):TPZAnalysis(mesh,out){

  this->fTimeStep = 0.;
  this->fCurrentIter = 0;
  this->SetConvergence(0, 0.);
  this->SetNewtonConvergence(0, 0.);
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
  for(int i = 1; i < this->fAllSolutions.NElements(); i++){
    this->GetSolution(i).Zero();
  }//i
  
  this->SetAllMaterialsDeltaT();
  for(this->fCurrentIter = 0; this->fCurrentIter < this->fNIter; this->fCurrentIter++){
  
    
  
  }//iterations
  

}//method

template<class T>
void TPZTransientAnalysis::SetLastState(){
  TPZCompMesh * mesh = this->Mesh();
  const int nmat = mesh->MaterialVec().NElements();  
  for(int i = 0; i < nmat; i++){
    TPZMaterial * mat = mesh->MaterialVec()[i];
    TPZTransientMaterial< T > * trans = dynamic_cast<TPZTransientMaterial< T > *>(mat);
    if (trans){
      trans->SetTimeStep( this->TimeStep() );
    }
  }
}
  
void TPZTransientAnalysis::SetCurrentState(){

}
  
void TPZTransientAnalysis::SetAllMaterialsDeltaT(){
  TPZCompMesh * mesh = this->Mesh();
  const int nmat = mesh->MaterialVec().NElements();  
  for(int i = 0; i < nmat; i++){
    TPZMaterial * mat = mesh->MaterialVec()[i];
    TPZTransientMaterial< TPZMatPoisson3d > * trans = dynamic_cast<TPZTransientMaterial< TPZMatPoisson3d > *>(mat);
    if (trans){
      trans->SetTimeStep( this->TimeStep() );
    }
  }
}
 