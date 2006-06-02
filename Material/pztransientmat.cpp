// -*- c++ -*-

//$Id: pztransientmat.cpp,v 1.1 2006-06-02 17:03:59 tiago Exp $

#include "pztransientmat.h"
#include "pzpoisson3d.h"

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(int nummat, int dim, REAL TimeStep):TBASEMAT(nummat, dim){
  this->SetTimeStep(TimeStep);
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::~TPZTransientMaterial(){
  //NOTHING TO BE DONE
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef){
  if (this->fStep == ECurrent){
    TBASEMAT::Contribute(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);
    this->ContributeSolutionRhs(sol, phi, weight, ef);
    this->ContributeTangent(phi, weight, ek);
    return;
  }
  
  if (this->fStep == ELast){
    this->ContributeSolutionRhs(sol, phi, weight, ef);
    return;
  }
  
  PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
  
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                            TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
  if (this->fStep == ECurrent){
    TBASEMAT::ContributeBC(x,sol,weight,axes,phi,ek,ef,bc);
    return;
  }
  
  if (this->fStep == ELast){
    return;
  }
  
  PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeSolutionRhs(TPZVec<REAL> &sol, TPZFMatrix &phi, REAL weight, TPZFMatrix &ef){
  REAL Mult = +1.; //Last solution is added to residual
  if (this->fStep == ECurrent) Mult = -1.; //Current solution is subtracted from residual
  const int phr = phi.Rows();
  const int nstate = this->NStateVariables();
  const REAL DeltaT = this->TimeStep();
  int i, k;
  for(i = 0; i < phr; i++) {
    for(k = 0; k < nstate; k++){
      ef(i*nstate+k, 0) += Mult * weight * sol[k] * phi(i,0) / DeltaT;
    }//k
  }//i
}//method

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeTangent(TPZFMatrix &phi, REAL weight, TPZFMatrix &ek){
  const int phr = phi.Rows();
  const int nstate = this->NStateVariables();
  const REAL DeltaT = this->TimeStep();
  int i, j, k;
  for(i = 0; i < phr; i++) {
    for(j = 0; j < phr; j++){
      for(k = 0; k < nstate; k++){
        ek(i*nstate+k, j*nstate+k) += weight * phi(i,0) * phi(j,0) / DeltaT;
      }//k
    }//j
  }//i
}//method


template class TPZTransientMaterial< TPZMatPoisson3d >;
