// -*- c++ -*-

//$Id: pzburger.h,v 1.4 2007-05-11 19:15:17 joao Exp $

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

  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef){
    if (TPZBurger::gStabilizationScheme == ESUPG){
      this->ContributeSUPG(data.x,data.jacinv,data.sol,data.dsol,weight,data.axes,data.phi,data.dphix,ek,ef);
    }
    if (TPZBurger::gStabilizationScheme == EGRADIENT){
      this->ContributeGradStab(data.x,data.jacinv,data.sol,data.dsol,weight,data.axes,data.phi,data.dphix,ek,ef);
    }
  }//Contribute
  
  void ContributeGradStab(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  void ContributeSUPG(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
               
  virtual void ContributeBC(TPZMaterialData &data,
                              REAL weight,
                              TPZFMatrix &ek,
                              TPZFMatrix &ef,
                              TPZBndCond &bc);

  virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZMaterialData &data,
                                       REAL weight,
                                       TPZFMatrix &ek,
                                       TPZFMatrix &ef,
                                       TPZBndCond &bc);

  protected:
    bool fIsReferred;

};

#endif
