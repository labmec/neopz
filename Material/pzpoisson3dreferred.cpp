// -*- c++ -*-

//$Id: pzpoisson3dreferred.cpp,v 1.2 2006-07-06 15:57:40 tiago Exp $

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

void TPZMatPoisson3dReferred::SetConvectionTerm(TPZFMatrix &dsolL, TPZFMatrix &axesL, TPZFMatrix &dsolR, TPZFMatrix &axesR){
//tiago
  const int dim = this->Dimension();
  const int pos = 1; // first column stores current solution. second column stores previous solution which is requested here. Then pos = 1.
  for(int i = 0; i < dim; i++){
    this->fConvDir[i] = -1. * ( 0.5 * ( dsolL(i,pos) + dsolR(i,pos) ) );
  }//for
  this->fC = this->falpha;  
}

void TPZMatPoisson3dReferred::Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                                         TPZVec<REAL> &sol, TPZFMatrix &dsol,
                                         REAL weight,TPZFMatrix &axes,
                                         TPZFMatrix &phi, TPZFMatrix &dphi,
                                         TPZFMatrix &ek, TPZFMatrix &ef){
  this->SetConvectionTerm(dsol, axes);
  TPZMatPoisson3d::Contribute(x, jacinv, sol, dsol, weight, axes, phi, dphi, ek, ef);                                          
}

void TPZMatPoisson3dReferred::ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                                           TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
  TPZMatPoisson3d::ContributeBC(x, sol, weight, axes, phi, ek, ef, bc);
}

void TPZMatPoisson3dReferred::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                                  TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                                  TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                                  TPZFMatrix &ek,TPZFMatrix &ef){
  PZError << __PRETTY_FUNCTION__ << " - this method is not working properly. Expect huge problems or better than that abort your simullation.\n";
//   this->SetConvectionTerm(dsolL, dsolR, axesL, axesR);
  TPZMatPoisson3d::ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, ek, ef);
}
  
void TPZMatPoisson3dReferred::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                                    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
  PZError << __PRETTY_FUNCTION__ << " - this method is not working properly. Expect huge problems or better than that abort your simullation.\n";
//   this->SetConvectionTerm(dsolL, axesL);
  TPZMatPoisson3d::ContributeBCInterface(x, solL, dsolL, weight, normal, phiL, dphiL, ek, ef, bc);
}
