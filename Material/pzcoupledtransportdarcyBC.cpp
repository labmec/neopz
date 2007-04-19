
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
  if (!bc) return;
  this->UpdateConvectionDir(dsol);
  TPZMaterialData data;
  data.x = x;
  data.jacinv = jacinv;
  data.sol = sol;
  data.dsol = dsol;
  data.axes = axes;
  data.phi = phi;
  data.dphix = dphi;
  bc->Contribute(data, weight, ek, ef);
}

void TPZCoupledTransportDarcyBC::ContributeInterface(TPZVec<REAL> &x,
                                                     TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                                     TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                                     REAL weight,TPZVec<REAL> &normal,
                                                     TPZFMatrix &phiL,TPZFMatrix &phiR,
                                                     TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                     TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                                     TPZFMatrix &ek,TPZFMatrix &ef) {
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (!bc) return;
  this->UpdateConvectionDirInterface(dsolL, dsolR, phiL, phiR);
  TPZMaterialData data;
  data.x = x;
  data.soll = solL;
  data.solr = solR;
  data.dsoll = dsolL;
  data.dsolr = dsolR;
  data.axesleft = axesleft;
  data.axesright = axesright;
  data.phil = phiL;
  data.phir = phiR;
  data.dphixl = dphiL;
  data.dphixr = dphiR;
  data.normal = normal;
  bc->ContributeInterface(data, weight, ek, ef);
}

void TPZCoupledTransportDarcyBC::ContributeInterface(TPZVec<REAL> &x,
                                                     TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                                     TPZFMatrix &dsolL,TPZFMatrix &dsolR,
                                                     REAL weight,TPZVec<REAL> &normal,
                                                     TPZFMatrix &phiL,TPZFMatrix &phiR,
                                                     TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                     TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                                     TPZFMatrix &ef) {
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (!bc) return;
  this->UpdateConvectionDirInterface(dsolL, dsolR, phiL, phiR);
  TPZMaterialData data;
  data.x = x;
  data.soll = solL;
  data.solr = solR;
  data.dsoll = dsolL;
  data.dsolr = dsolR;
  data.axesleft = axesleft;
  data.axesright = axesright;
  data.phil = phiL;
  data.phir = phiR;
  data.dphixl = dphiL;
  data.dphixr = dphiR;
  data.normal = normal;
  bc->ContributeInterface(data, weight, ef);
}

void TPZCoupledTransportDarcyBC::
  ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                      TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                      TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                      TPZFMatrix &axesleft, TPZFMatrix &axesright,
                      TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
  TPZBndCond * bc = this->GetCurrentMaterial();
  if (!bc) return;
  this->UpdateConvectionDirInterface(dsolL, dsolR, phiL, phiR);
  TPZMaterialData data;
  data.x = x;
  data.soll = solL;
  data.solr = solR;
  data.dsoll = dsolL;
  data.dsolr = dsolR;
  data.axesleft = axesleft;
  data.axesright = axesright;
  data.phil = phiL;
  data.phir = phiR;
  data.dphixl = dphiL;
  data.dphixr = dphiR;
  data.normal = normal;
  data.leftp = LeftPOrder;
  data.rightp = RightPOrder;
  data.HSize = faceSize;
  bc->ContributeInterface(data, weight, ek, ef);
}

void TPZCoupledTransportDarcyBC::UpdateConvectionDir(TPZFMatrix &dsol){
  TPZCoupledTransportDarcy * mat = dynamic_cast<TPZCoupledTransportDarcy*>(this->Material().operator ->());
  if (!mat){
    PZError << __PRETTY_FUNCTION__ << " FATAL ERROR" << std::endl;
    exit(-1);
  }
  mat->UpdateConvectionDir(dsol);
}

void TPZCoupledTransportDarcyBC::UpdateConvectionDirInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                                              TPZFMatrix &phiL, TPZFMatrix &phiR){
  TPZCoupledTransportDarcy * mat = dynamic_cast<TPZCoupledTransportDarcy*>(this->Material().operator ->());
  if (!mat){
    PZError << __PRETTY_FUNCTION__ << " FATAL ERROR" << std::endl;
    exit(-1);
  }
  if (phiL.Rows()) mat->UpdateConvectionDir(dsolL);
  if (phiR.Rows()) mat->UpdateConvectionDir(dsolR);
}
