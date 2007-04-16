// -*- c++ -*-
// $Id: pzdiscgal.cpp,v 1.6 2007-04-16 20:27:36 tiago Exp $

#include "pzdiscgal.h"
#include "pzmaterialdata.h"

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin() : TPZMaterial(){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

TPZDiscontinuousGalerkin::~TPZDiscontinuousGalerkin(){}

char *TPZDiscontinuousGalerkin::Name() { return "TPZDiscontinuousGalerkin"; }

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
  if (data.fNeedsHSize){
    this->ContributeInterface(data.x, data.soll, data.solr, data.dsoll, data.dsolr,
                              weight, data.normal, data.phil, data.phir, data.dphixl, data.dphixr,
                              data.axesleft, data.axesright, ek, ef,
                              data.leftp, data.rightp, data.HSize);
  }
  else{
    this->ContributeInterface(data.x, data.soll, data.solr, data.dsoll, data.dsolr,
                              weight, data.normal, data.phil, data.phir, data.dphixl, data.dphixr,
                              data.axesleft, data.axesright, ek, ef);
  }
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){
  TPZFMatrix fakeek(ef.Rows(), ef.Rows(), 0.);
  this->ContributeInterface(data, weight, fakeek, ef);
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
  if (data.fNeedsHSize){
    this->ContributeBCInterface(data.x, data.soll, data.dsoll,
                                weight, data.normal, data.phil, data.dphixl,
                                data.axesleft, ek, ef, bc, data.leftp, data.HSize);
  }
  else{
    this->ContributeBCInterface(data.x, data.soll, data.dsoll,
                                weight, data.normal, data.phil, data.dphixl,
                                data.axesleft, ek, ef, bc);
  }
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef,TPZBndCond &bc){
  TPZFMatrix fakeek(ef.Rows(), ef.Rows(), 0.);
  this->ContributeBCInterface(data, weight, fakeek, ef, bc);
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                                   TPZFMatrix &ef){
  TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
  ContributeInterface(x, solL, solR, dsolL, dsolR,
                      weight, normal, phiL, phiR,
                      dphiL, dphiR,axesleft,axesright, ek, ef);
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR, TPZFMatrix &axesleft,
                                                   TPZFMatrix &axesright,
                                                   TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
  this->ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, axesleft, axesright, ek, ef);
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			                             TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                                     TPZFMatrix axesleft,
                                                     TPZFMatrix &ef,TPZBndCond &bc){
  TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
  ContributeBCInterface(x, solL,  dsolL, weight, normal, phiL, dphiL, ek, ef, bc);
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                                     TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                                     TPZFMatrix &axesleft,
                                                     TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize){

  this->ContributeBCInterface(x, solL, dsolL, weight, normal, phiL, dphiL, ek, ef, bc);
}

int TPZDiscontinuousGalerkin::IsInterfaceConservative(){
  return 0;
}

void TPZDiscontinuousGalerkin::InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu, TPZVec<REAL> &leftNormalDeriv,
                                              TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                                              TPZVec<REAL> &values){
  PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
}

void TPZDiscontinuousGalerkin::BCInterfaceJumps(TPZVec<REAL> &leftu,
                                                TPZBndCond &bc,
                                                TPZVec<REAL> &values){
  PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
}



