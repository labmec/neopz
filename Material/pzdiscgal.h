// -*- c++ -*-
// $Id: pzdiscgal.h,v 1.1 2003-11-25 17:08:00 phil Exp $
#ifndef TPZDISCGALHPP
#define TPZDISCGALHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"


class TPZDiscontinuousGalerkin  : public TPZMaterial {

  
  public :

  TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat) {}

  /**copy constructor*/
  TPZDiscontinuousGalerkin(TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

  virtual  ~TPZDiscontinuousGalerkin(){};
  
  char *Name() { return "TPZDiscontinuousGalerkin"; }
  
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef) = 0;
  
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) = 0;

  //  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
  //			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ef);
};
  
#endif
