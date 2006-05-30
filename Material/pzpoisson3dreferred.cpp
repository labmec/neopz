// -*- c++ -*-

//$Id: pzpoisson3dreferred.cpp,v 1.1 2006-05-30 17:45:00 tiago Exp $

#include "pzpoisson3dreferred.h"

using namespace std;

TPZMatPoisson3dReferred::TPZMatPoisson3dReferred(int nummat, int dim):TPZMatPoisson3d(nummat,dim){ 
    this->falpha = -1.; 
}
  
TPZMatPoisson3dReferred::~TPZMatPoisson3dReferred(){ }

void TPZMatPoisson3dReferred::SetConvectionTerm(TPZFMatrix &dsol){
  const int dim = this->Dimension();
  for(int i = 0; i < dim; i++){
    this->fConvDir[i] = dsol(i,0);
  }//for
  this->fC = this->falpha;  
}

void TPZMatPoisson3dReferred::SetConvectionTerm(TPZFMatrix &dsolL, TPZFMatrix &dsolR){
  const int dim = this->Dimension();
  const int pos = 1; // first column stores current solution. second column stores previous solution which is requested here. Then pos = 1.
  for(int i = 0; i < dim; i++){
    this->fConvDir[i] = 0.5 * ( dsolL(i,pos) + dsolR(i,pos) );
  }//for
  this->fC = this->falpha;  
}

void TPZMatPoisson3dReferred::Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                                         TPZVec<REAL> &sol, TPZFMatrix &dsol,
                                         REAL weight,TPZFMatrix &axes,
                                         TPZFMatrix &phi, TPZFMatrix &dphi,
                                         TPZFMatrix &ek, TPZFMatrix &ef){
  this->SetConvectionTerm(dsol);
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
  this->SetConvectionTerm(dsolL, dsolR);
  TPZMatPoisson3d::ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, ek, ef);
}
  
void TPZMatPoisson3dReferred::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                                    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
  this->SetConvectionTerm(dsolL);
  TPZMatPoisson3d::ContributeBCInterface(x, solL, dsolL, weight, normal, phiL, dphiL, ek, ef, bc);
}
