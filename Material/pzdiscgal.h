// -*- c++ -*-
// $Id: pzdiscgal.h,v 1.2 2004-02-05 16:12:28 tiago Exp $
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

  /**
   * Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces. 
   * Interfaces may be conservative or not conservative. It is important to agglomeration techniques
   * when using multigrid pre-conditioner. Conservative interfaces into agglomerate elements do not
   * need to be computed. However non-conservative interfaces must be computed in all multigrid levels.
   * Default is non-conservative, because of the computation of a conservative interface into an agglomerate
   * does not ruin the solution.
   * @since Feb 05, 2004
   */
  virtual int IsInterfaceConservative(){ return 0; }
};
  
#endif
