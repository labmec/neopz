// -*- c++ -*- 

//$Id: pzpoisson3dreferred.cpp,v 1.5 2007-05-11 19:15:18 joao Exp $

#include "pzpoisson3dreferred.h"

using namespace std;

TPZMatPoisson3dReferred::TPZMatPoisson3dReferred(int nummat, int dim):TPZMatPoisson3d(nummat,dim){ 
    this->falpha = -1.; 
}
  
TPZMatPoisson3dReferred::~TPZMatPoisson3dReferred(){ }

void TPZMatPoisson3dReferred::SetConvectionTerm(TPZFMatrix &dsol, TPZFMatrix &axes){
  const int dim = this->Dimension();
  TPZManVector<REAL,3> V(dim);
  const int pos = 1; // first column stores current solution. second column stores previous solution which is requested here. Then pos = 1.
  for(int i = 0; i < dim; i++){
    V[i] = -1. * dsol(i,pos);
  }//for
  
  for(int i = 0; i < 3; i++) this->fConvDir[i] = 0.;
  
  switch(dim) {
    case 1:
      this->fConvDir[0] = axes(0,0) * V[0];
      this->fConvDir[1] = axes(0,1) * V[0];
      this->fConvDir[2] = axes(0,2) * V[0];
    break;
    case 2:
      this->fConvDir[0] = axes(0,0) * V[0] + axes(1,0) * V[1];
      this->fConvDir[1] = axes(0,1) * V[0] + axes(1,1) * V[1];
      this->fConvDir[2] = axes(0,2) * V[0] + axes(1,2) * V[1];
    break;
    case 3:
      this->fConvDir[0] = axes(0,0) * V[0] + axes(1,0) * V[1] + axes(2,0) * V[2];
      this->fConvDir[1] = axes(0,1) * V[0] + axes(1,1) * V[1] + axes(2,1) * V[2];
      this->fConvDir[2] = axes(0,2) * V[0] + axes(1,2) * V[1] + axes(2,2) * V[2];
    break;
    default:
      PZError << "TPZMatPoisson3dReferred::SetConvectionTerm - Dimension error " << fDim << endl;
  }
  
  this->fC = this->falpha;  
}

void TPZMatPoisson3dReferred::SetConvectionTermInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR){
  const int dim = this->Dimension();
  const int pos = 1; // first column stores current solution. second column stores previous solution which is requested here. Then pos = 1.
  for(int i = 0; i < dim; i++){
    this->fConvDir[i] = -1. * ( 0.5 * ( dsolL(i,pos) + dsolR(i,pos) ) );
  }//for
  this->fC = this->falpha;  
}

void TPZMatPoisson3dReferred::Contribute(TPZMaterialData &data,
                                         REAL weight,
                                         TPZFMatrix &ek, 
                                         TPZFMatrix &ef){

  this->SetConvectionTerm(data.dsol, data.axes);
  TPZMatPoisson3d::Contribute(data, weight, ek, ef);
}

void TPZMatPoisson3dReferred::ContributeBC(TPZMaterialData &data,
                                           REAL weight,
                                           TPZFMatrix &ek,
                                           TPZFMatrix &ef,
                                           TPZBndCond &bc){
  TPZMatPoisson3d::ContributeBC(data, weight, ek, ef, bc);
}

void TPZMatPoisson3dReferred::ContributeInterface(TPZMaterialData &data,
                                                  REAL weight,
                                                  TPZFMatrix &ek,
                                                  TPZFMatrix &ef){
TPZFMatrix dsolL=data.dsoll;
TPZFMatrix dsolR=data.dsolr;
  this->SetConvectionTerm(dsolL, dsolR);
  TPZMatPoisson3d::ContributeInterface(data, weight, ek, ef);
}
  
void TPZMatPoisson3dReferred::ContributeBCInterface(TPZMaterialData &data,
                                                    REAL weight,
                                                    TPZFMatrix &ek,
                                                    TPZFMatrix &ef,
                                                    TPZBndCond &bc){
  this->SetConvectionTerm(data.dsoll, data.dsoll);
  TPZMatPoisson3d::ContributeBCInterface(data, weight,  ek, ef, bc);
}
