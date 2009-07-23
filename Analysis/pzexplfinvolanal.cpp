//$Id: pzexplfinvolanal.cpp,v 1.1 2009-07-23 20:37:56 fortiago Exp $

#include "pzexplfinvolanal.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"
#include "checkconv.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"

using namespace std;


TPZExplFinVolAnal::TPZExplFinVolAnal(TPZCompMesh *mesh, std::ostream &out):TPZAnalysis(mesh,out){
  this->fTimeStep = -1.;
  this->fNMaxIter = -1.;
  this->Set(-1.,0.,0.);
  this->SetInitialSolutionAsZero();
  this->SetSaveFrequency(0,0);
}

TPZExplFinVolAnal::~TPZExplFinVolAnal(){
///nothing to be done here
}

void TPZExplFinVolAnal::SetInitialSolution(TPZFMatrix & InitialSol){
  const int nrows = this->Mesh()->Solution().Rows();
  const int ncols = this->Mesh()->Solution().Cols();
  if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
    PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
  }
  else{
    this->fSolution = InitialSol;
  }
  TPZAnalysis::LoadSolution();
}

void TPZExplFinVolAnal::SetInitialSolutionAsZero(){
  TPZFMatrix & MeshSol = this->Mesh()->Solution();
  this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
  this->fSolution.Zero();
}


void TPZExplFinVolAnal::Run(std::ostream &out, bool linesearch){

  ofstream File("Pocos.txt");
  File << "t(mes)\tPi(Pa)\tPp(Pa)\tQwiSC(m3/day)\tQoiSC(m3/day)\tQwpSC(m3/day)\tQopSC(m3/day)\t\tQwiFundo(m3/day)\tQoiFundo(m3/day)\tQwpFundo(m3/day)\tQopFundo(m3/day)\n";

  this->PostProcess(this->fDXResolution);

  TPZFMatrix NextSol;

  REAL steadynorm;
  for(int iter = 0; iter < this->fNMaxIter; iter++){

      this->UpdateSolution(fSolution,fRhs,NextSol);

      ///checking steady state
      steadynorm = 0.;
      for(int i = 0; i < fSolution.Rows(); i++){
        double val = (NextSol(i,0)-fSolution(i,0));
        steadynorm += val*val;
      }
      steadynorm = sqrt(steadynorm);

      fSolution = NextSol;
      TPZAnalysis::LoadSolution();

   if (this->fSaveFrequency){
    if (!(iter % fSaveFrequency)){
      this->fSimulationTime = (iter+1.)*fTimeStep;
      this->PostProcess(this->fDXResolution);
    }
   }

   std::cout << "*********** Steady state error at iteration " << iter << ", " << iter*fTimeStep << " = " << steadynorm << "\n\n";
   if (!fForceAllSteps){
    if (steadynorm < this->fSteadyTol){
      std::cout << "Steady state solution achieved\n\n";
      break;
    }
   }
   std::cout.flush();

  }///time step iterations

}///method

void TPZExplFinVolAnal::PostProcess(int resolution, int dimension){
    REAL T = this->fSimulationTime;
    this->fTime = T;
    TPZAnalysis::PostProcess(resolution, dimension);
}//method


void TPZExplFinVolAnal::Set(REAL timestep, int niter, REAL eps, bool ForceAllSteps){
  this->fTimeStep = timestep;
  this->fNMaxIter = niter;
  this->fSteadyTol = eps;
  this->fForceAllSteps = ForceAllSteps;
}

void TPZExplFinVolAnal::SetSaveFrequency(int SaveFrequency, int resolution){
  this->fSaveFrequency = SaveFrequency;
  this->fDXResolution = resolution;
}

REAL TPZExplFinVolAnal::TimeStep(){
  return this->fTimeStep;
}

void TPZExplFinVolAnal::AssembleFluxes(TPZFMatrix & rhs, std::set<int> *MaterialIds){
  const int nelem = this->Mesh()->NElements();
  TPZElementMatrix ef(this->Mesh(), TPZElementMatrix::EF);

  for(int iel = 0; iel < nelem; iel++){
    TPZCompEl *el = this->Mesh()->ElementVec()[iel];
    if(!el) continue;

    if(MaterialIds){
      TPZAutoPointer<TPZMaterial> mat = el->Material();
      if (!mat) continue;
      int matid = mat->Id();
      if (MaterialIds->find(matid) == MaterialIds->end()) continue;
    }///if

    TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(el);
    if(!face) continue;

    el->CalcResidual(ef);
    ef.ComputeDestinationIndices();
    rhs.AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);

  }///for iel = TPZInterfaceElement

  for(int iel=0; iel < nelem; iel++){
    TPZCompEl *el = this->Mesh()->ElementVec()[iel];
    if(!el) continue;

    if(MaterialIds){
      TPZAutoPointer<TPZMaterial> mat = el->Material();
      if (!mat) continue;
      int matid = mat->Id();
      if (MaterialIds->find(matid) == MaterialIds->end()) continue;
    }///if

    TPZInterpolationSpace* sp = dynamic_cast<TPZInterpolationSpace*>(el);
    if(!sp) continue;

    if(el->Reference()->Dimension() != 3) continue;

    const REAL volume = el->Reference()->Volume();
    sp->InitializeElementMatrix(ef);
    ef.ComputeDestinationIndices();
    const int n = ef.fDestinationIndex.NElements();
    for(int i = 0; i < n; i++){
      int pos = ef.fDestinationIndex[i];
      rhs(pos,0) *= this->fTimeStep/volume;
    }

  }///for iel = TPZInterpolationSpace

}///void

void TPZExplFinVolAnal::UpdateSolution(TPZFMatrix &LastSol, TPZFMatrix & rhs, TPZFMatrix &NextSol){
  int sz = fCompMesh->NEquations();
  rhs.Redim(sz,1);

  ///Euler explicit:
  this->AssembleFluxes(rhs);
  NextSol = LastSol;
  NextSol.operator -=(rhs);

  ///RK2
  TPZFMatrix uEtoile(sz,1,0.);
//   fazer
  
  ///RK3
//   fazer
  

}///void
