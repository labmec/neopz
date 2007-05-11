// -*- c++ -*-
// $Id: pzdiscgal.cpp,v 1.7 2007-05-11 19:15:17 joao Exp $

#include "pzdiscgal.h"
#include "pzmaterialdata.h"

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin() : TPZMaterial(){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

TPZDiscontinuousGalerkin::~TPZDiscontinuousGalerkin(){}

char *TPZDiscontinuousGalerkin::Name() { return "TPZDiscontinuousGalerkin"; }

void TPZDiscontinuousGalerkin::FillDataRequirementsInterface(TPZMaterialData &data){
  data.SetAllRequirements(true);
  data.fNeedsSol = false;
}

void TPZDiscontinuousGalerkin::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){
  TPZFMatrix fakeek(ef.Rows(), ef.Rows(), 0.);
  this->ContributeInterface(data, weight, fakeek, ef);
}

void TPZDiscontinuousGalerkin::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef,TPZBndCond &bc){
  TPZFMatrix fakeek(ef.Rows(), ef.Rows(), 0.);
  this->ContributeBCInterface(data, weight, fakeek, ef, bc);
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



