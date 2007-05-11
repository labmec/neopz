// -*- c++ -*-

//$Id: pznonlinearpoisson3d.h,v 1.4 2007-05-11 19:15:18 joao Exp $

#ifndef MATNLPOISSON3DH
#define MATNLPOISSON3DH

#include "pzpoisson3d.h"
#include "pzpoisson3dreferred.h"

class TPZNonLinearPoisson3d : public TPZMatPoisson3dReferred {

 public:

  TPZNonLinearPoisson3d(int nummat, int dim);
  
  TPZNonLinearPoisson3d(const TPZNonLinearPoisson3d &cp);
  
  virtual ~TPZNonLinearPoisson3d();
  
  bool IsReferred(){ return this->fIsReferred;}
  
  void SetReferred(bool Is){ this->fIsReferred = Is; }

  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
               
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
