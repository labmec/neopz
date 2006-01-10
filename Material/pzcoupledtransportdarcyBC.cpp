
#include "pzcoupledtransportdarcyBC.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"

TPZCoupledTransportDarcyBC::TPZCoupledTransportDarcyBC(TPZCoupledTransportDarcy *material, int id) : TPZBndCond(){
  this->SetId(id);
  this->fMaterial = material;
  this->fMaterials[0] = NULL;
  this->fMaterials[1] = NULL;
}

TPZCoupledTransportDarcyBC::~TPZCoupledTransportDarcyBC(){}

void TPZCoupledTransportDarcyBC::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,
                                            TPZVec<REAL> &sol,TPZFMatrix &dsol,
                                            REAL weight,TPZFMatrix &axes,
                                            TPZFMatrix &phi,TPZFMatrix &dphi,
                                            TPZFMatrix &ek,TPZFMatrix &ef){
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (bc) bc->Contribute(x, jacinv, sol, dsol, weight, axes, phi, dphi, ek, ef);
}

void TPZCoupledTransportDarcyBC::ContributeInterface(TPZVec<REAL> &x,
                                                     TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                                     TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                                     REAL weight,TPZVec<REAL> &normal,
                                                     TPZFMatrix &phiL,TPZFMatrix &phiR,
                                                     TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                     TPZFMatrix &ek,TPZFMatrix &ef) {
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (bc) bc->ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, ek, ef);
}

void TPZCoupledTransportDarcyBC::ContributeInterface(TPZVec<REAL> &x,
                                                     TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                                     TPZFMatrix &dsolL,TPZFMatrix &dsolR,
                                                     REAL weight,TPZVec<REAL> &normal,
                                                     TPZFMatrix &phiL,TPZFMatrix &phiR,
                                                     TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                     TPZFMatrix &ef) {
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (bc) bc->ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, ef);
}

void TPZCoupledTransportDarcyBC::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
			 TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
			 TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
			 TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (bc) bc->ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, ek, ef, LeftPOrder, RightPOrder, faceSize);
}


