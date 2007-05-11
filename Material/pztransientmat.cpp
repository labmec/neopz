// -*- c++ -*-

//$Id: pztransientmat.cpp,v 1.5 2007-05-11 19:15:18 joao Exp $
 
#include "pztransientmat.h"

template<class TBASEMAT>
int TPZTransientMaterial< TBASEMAT >::gTemporalIntegrator = EImplicit;

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::SetExplicit(){
  this->gTemporalIntegrator = EExplicit;
}
template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::SetImplicit(){
  this->gTemporalIntegrator = EImplicit;
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(int nummat, int dim, REAL TimeStep):TBASEMAT(nummat, dim){
  this->SetTimeStep(TimeStep);
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::TPZTransientMaterial(const TPZTransientMaterial &cp):TBASEMAT(cp){
  this->fStep = cp.fStep;
  this->fTimeStep = cp.fTimeStep;
}

template<class TBASEMAT>
TPZTransientMaterial< TBASEMAT >::~TPZTransientMaterial(){
  //NOTHING TO BE DONE
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::Contribute(TPZMaterialData &data,
                                                  REAL weight,
                                                  TPZFMatrix &ek,
                                                  TPZFMatrix &ef){
 if (this->gTemporalIntegrator == EImplicit){
  if (this->fStep == ECurrent){
    TBASEMAT::Contribute(data,weight,ek,ef);
    this->ContributeSolutionRhs(data.sol, data.phi, weight, ef);
    this->ContributeTangent(data.phi, weight, ek);
    return;
  }
  
  if (this->fStep == ELast){
    this->ContributeSolutionRhs(data.sol, data.phi, weight, ef);
    return;
  }
 }//EImplicit
 
 if (this->gTemporalIntegrator == EExplicit){
  if (this->fStep == EMassMatrix){
    this->ContributeTangent(data.phi, weight, ek);
    return;
  }
  if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
    TBASEMAT::Contribute(data,weight,ek,ef);
    return;
  }
 }//EExplicit
  
 PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
  
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeBC(TPZMaterialData &data,
                                                    REAL weight,
                                                    TPZFMatrix &ek,
                                                    TPZFMatrix &ef,
                                                    TPZBndCond &bc){
 if (this->gTemporalIntegrator == EImplicit){
  if (this->fStep == ECurrent){
    TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
    return;
  }
  
  if (this->fStep == ELast){
    return;
  }
 }//EImplicit

 if (this->gTemporalIntegrator == EExplicit){
  if (this->fStep == EMassMatrix){
    return;
  }
  if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
    TBASEMAT::ContributeBC(data,weight,ek,ef,bc);
    return;
  }
 }//EExplicit
  
  PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeInterface(TPZMaterialData &data,
                                                           REAL weight,
                                                           TPZFMatrix &ek,
                                                           TPZFMatrix &ef){
 if (this->gTemporalIntegrator == EImplicit){
  if (this->fStep == ECurrent){
    TBASEMAT::ContributeInterface(data, weight, ek, ef);
    return;
  }
  
  if (this->fStep == ELast){
    return;
  }
 }//EImplicit
 
 if (this->gTemporalIntegrator == EExplicit){
  if (this->fStep == EMassMatrix){
    return;
  }
  if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
    TBASEMAT::ContributeInterface(data, weight, ek, ef);
    return;
  }
 }//EExplicit
  
  PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at LINE " << __LINE__ << std::endl;
}

template<class TBASEMAT>
void TPZTransientMaterial< TBASEMAT >::ContributeBCInterface(TPZMaterialData &data,
                                                             REAL weight, 
                                                             TPZFMatrix &ek,
                                                             TPZFMatrix &ef,
                                                             TPZBndCond &bc){
 if (this->gTemporalIntegrator == EImplicit){
  if (this->fStep == ECurrent){
    TBASEMAT::ContributeBCInterface(data, weight,ek, ef, bc);
    return;
  }
  
  if (this->fStep == ELast){
    return;
  }
 }//EImplicit
 
 if (this->gTemporalIntegrator == EExplicit){
  if (this->fStep == EMassMatrix){
    return;
  }
  if (this->fStep == EFluxOnly){ //Calcula ef = F-ku
    TBASEMAT::ContributeBCInterface(data, weight,  ek, ef, bc);
    return;
  }
 }//EExplicit
  
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

#include "pzpoisson3d.h"
template class TPZTransientMaterial< TPZMatPoisson3d >;

#include "pznonlinearpoisson3d.h"
template class TPZTransientMaterial< TPZNonLinearPoisson3d >;

#include "pzburger.h"
template class TPZTransientMaterial< TPZBurger >;
