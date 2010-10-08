// -*- c++ -*-
// $Id: pzdiscgal.cpp,v 1.12 2010-10-08 19:18:15 fortiago Exp $

#include "pzdiscgal.h"
#include "pzmaterialdata.h"
#include "pzmaterialid.h"

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin() : TPZMaterial(){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat){}

TPZDiscontinuousGalerkin::TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

TPZDiscontinuousGalerkin::~TPZDiscontinuousGalerkin(){}

std::string TPZDiscontinuousGalerkin::Name() { return "TPZDiscontinuousGalerkin"; }

void TPZDiscontinuousGalerkin::FillDataRequirementsInterface(TPZMaterialData &data){
  data.SetAllRequirements(true);
  data.fNeedsSol = false;
  if(fLinearContext == false){
    data.fNeedsNeighborSol = true;
  }
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

void TPZDiscontinuousGalerkin::InterfaceJump(TPZVec<REAL> &x, 
                                              TPZVec<REAL> &leftu,
                                              TPZVec<REAL> &rightu,
                                              TPZVec<REAL> &jump){
  const int n = leftu.NElements();
  jump.Resize(n);
  for(int i = 0; i < n; i++){
    jump[i] = leftu[i] - rightu[i];
  }
}

void TPZDiscontinuousGalerkin::BCInterfaceJump(TPZVec<REAL> &x, 
                                               TPZVec<REAL> &leftu,
                                               TPZBndCond &bc,
                                               TPZVec<REAL> & jump){
  PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
  DebugStop();
}

int TPZDiscontinuousGalerkin::ClassId() const{
  return TPZDISCONTINUOUSGALERKIN;
}

void TPZDiscontinuousGalerkin::Write(TPZStream &buf, int withclassid){
  TPZMaterial::Write(buf, withclassid);
}

void TPZDiscontinuousGalerkin::Read(TPZStream &buf, void *context){
  TPZMaterial::Read(buf, context);
}
