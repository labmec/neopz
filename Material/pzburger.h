// -*- c++ -*-

//$Id: pzburger.h,v 1.3 2007-01-27 14:49:27 phil Exp $

#ifndef BURGERH
#define BURGERH

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

class TPZBurger : public TPZMatPoisson3dReferred {

 public:

  enum EStabilizationScheme{ESUPG = 1, EGRADIENT = 2};
  
  static int gStabilizationScheme;

  TPZBurger(int nummat, int dim);
  
  TPZBurger(const TPZBurger &cp);
  
  virtual ~TPZBurger();
  
//   virtual int HasForcingFunction() {
//     return true;
//   }
   
  bool IsReferred(){ return this->fIsReferred;}
  
  void SetReferred(bool Is){ this->fIsReferred = Is; }
  
  REAL fSolRef;

  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef){
    if (TPZBurger::gStabilizationScheme == ESUPG){
      this->ContributeSUPG(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);
    }
    if (TPZBurger::gStabilizationScheme == EGRADIENT){
      this->ContributeGradStab(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);
    }
  }//Contribute
  
  void ContributeGradStab(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  void ContributeSUPG(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
               
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                            TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
                            
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                   TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL, 
                                     TPZFMatrix &axesleft,
                                     TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
                                     
  protected:
    bool fIsReferred;

};

#endif
