// -*- c++ -*-

//$Id: pznonlinearpoisson3d.h,v 1.1 2006-07-31 13:46:25 tiago Exp $

#ifndef MATNLPOISSON3DH
#define MATNLPOISSON3DH

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

class TPZNonLinearPoisson3d : public TPZMatPoisson3dReferred {

 public:

  TPZNonLinearPoisson3d(int nummat, int dim);
  
  virtual ~TPZNonLinearPoisson3d();
  
  bool IsReferred(){ return this->fIsReferred;}
  
  void SetReferred(bool Is){ this->fIsReferred = Is; }

  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
               
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                            TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
                            
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
                                     
  protected:
    bool fIsReferred;

};

#endif
