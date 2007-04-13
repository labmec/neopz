// -*- c++ -*-
// $Id: pzdiscgal.cpp,v 1.5 2007-04-13 15:17:54 tiago Exp $

#include "pzdiscgal.h"

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin() : TPZMaterial(){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

TPZDiscontinuousGalerkin::~TPZDiscontinuousGalerkin(){}

char *TPZDiscontinuousGalerkin::Name() { return "TPZDiscontinuousGalerkin"; }

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

void TPZDiscontinuousGalerkin::ContributeInterface(TPZVec<REAL> &x,
                                                   TPZVec<REAL> &sol, TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                                   TPZFMatrix &dsol, TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                                   REAL weight,TPZVec<REAL> &normal,
                                                   TPZFMatrix &phiL, TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                   TPZFMatrix &axes, TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                                   TPZFMatrix &ek,TPZFMatrix &ef){
  ContributeInterface(x,solL,solR,dsolL,dsolR,weight,normal,phiL,phiR,dphiL,dphiR,axesleft,axesright,ek,ef);
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZVec<REAL> &x,
                                                   TPZVec<REAL> &sol, TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                                   TPZFMatrix &dsol, TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                                   REAL weight,TPZVec<REAL> &normal,
                                                   TPZFMatrix &phiL, TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                   TPZFMatrix &axes, TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                                   TPZFMatrix &ef){
  TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
  ContributeInterface(x,solL,solR,dsolL,dsolR,weight,normal,phiL,phiR,dphiL,dphiR,axesleft,axesright,ek,ef);
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



