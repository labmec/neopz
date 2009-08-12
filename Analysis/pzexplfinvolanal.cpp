//$Id: pzexplfinvolanal.cpp,v 1.3 2009-08-12 21:04:57 fortiago Exp $

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
  this->fSimulationTime = 0.;
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


void TPZExplFinVolAnal::Run(std::ostream &out){

  fSimulationTime = 0.;

  this->PostProcess(this->fDXResolution);

  TPZFMatrix LastSol, NextSol;
  LastSol = fSolution;

  REAL steadynorm;

  for(int iter = 0; iter < this->fNMaxIter; iter++){

      this->UpdateSolution(LastSol,NextSol);

      if(iter == (fNMaxIter -1)){
        ofstream solfile("solucaoOlivier.txt");
        solfile.precision(12);
        for(int is = 0; is < NextSol.Rows(); is++){
          solfile << NextSol(is,0) << "\t";
          if((is+1)%5 == 0) solfile << "\n";
        }
        solfile << "\n";
      }

      ///checking steady state
      steadynorm = 0.;
      for(int i = 0; i < LastSol.Rows(); i++){
        double val = (NextSol(i,0)-LastSol(i,0));
        steadynorm += val*val;
      }
      steadynorm = sqrt(steadynorm);

      LastSol = NextSol;
      fSolution = NextSol;
      TPZAnalysis::LoadSolution();

   if (this->fSaveFrequency){
    if ((!(iter % fSaveFrequency)) || (iter == (this->fNMaxIter-1))){
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

void TPZExplFinVolAnal::AssembleFluxes(const TPZFMatrix & Solution, std::set<int> *MaterialIds){
  this->fSolution = Solution;
  this->LoadSolution();

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
    this->fRhs.AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);

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
      this->fRhs(pos,0) *= this->fTimeStep/volume;
    }

  }///for iel = TPZInterpolationSpace

}///void

void TPZExplFinVolAnal::UpdateSolution(TPZFMatrix &LastSol, TPZFMatrix &NextSol){

  const int order = 1;
  if(order == 1){
    ///Euler explicit:
    int sz = fCompMesh->NEquations();
    fRhs.Redim(sz,1);
    this->AssembleFluxes(LastSol);
    NextSol = LastSol;
    NextSol -= fRhs;
  }
  if(order == 2){
    ///RK2
    int sz = fCompMesh->NEquations();
    this->fRhs.Redim(sz,1);
    TPZFMatrix uEtoile(sz,1,0.);
    ///uEtoile = un - 0.5*rhs(un)
    this->fRhs.Zero();
    this->AssembleFluxes(LastSol);
    uEtoile = LastSol;
    uEtoile.ZAXPY(-0.5, this->fRhs );

    ///un+1 = un - rhs(uEtoile)
    this->fRhs.Zero();
    this->AssembleFluxes(uEtoile);
    NextSol = LastSol;
    NextSol -= this->fRhs;
  }
  if(order == 3){
    ///RK3: un+1 = un + dT (1/6 k1 +2/3 k2 +1/6 k3)
    int sz = fCompMesh->NEquations();
    this->fRhs.Redim(sz,1);
    NextSol = LastSol;
      ///un+1 = un

    this->AssembleFluxes(LastSol);///k1*dT = -rhs
    NextSol.ZAXPY(-1./6.,this->fRhs);
      ///un+1 = un +1/6 k1*dt

        ///keeping values
        TPZFMatrix uEtoile(sz,1,0.);
        uEtoile = LastSol;
        uEtoile += this->fRhs;

      ///fSolution = un - 0.5*rhs(un) = un+0.5 *K1*dt
    TPZFMatrix uStar;
    uStar = LastSol;
    uStar.ZAXPY(-0.5,this->fRhs);
    this->fRhs.Zero();
    this->AssembleFluxes(uStar);///k2*dT = -rhs
    NextSol.ZAXPY(-2./3.,this->fRhs);

      ///un+1 = un +1/6 k1*dt+2/3 k2*dt
    uEtoile.ZAXPY(-2.,this->fRhs);
    this->fRhs.Zero();
    this->AssembleFluxes(uEtoile);///k3*dT = -rhs
    NextSol.ZAXPY(-1./6.,this->fRhs);
  }

}///void
