// -*- c++ -*-

//$Id: pzpoisson3dreferred.h,v 1.1 2006-05-30 17:45:00 tiago Exp $

#ifndef MATPOISSON3DREFERREDH
#define MATPOISSON3DREFERREDH

#include <iostream>
#include "pzpoisson3d.h"
#include "pzfmatrix.h"

/**
 * This class implements a version of TPZMatPoisson3d where the convection term is given at each integration point
 * from a previous calculation.
 * The convection term fC * fConvDir = fAlpha * grad(sol) where grad(sol) is the gradient of the previous solution.
 */
class TPZMatPoisson3dReferred : public TPZMatPoisson3d {

protected:

  REAL falpha;
  
  void SetConvectionTerm(TPZFMatrix &dsol);
  void TPZMatPoisson3dReferred::SetConvectionTerm(TPZFMatrix &dsolL, TPZFMatrix &dsolR);

public:
 
  TPZMatPoisson3dReferred(int nummat, int dim);

  virtual ~TPZMatPoisson3dReferred();

  TPZMatPoisson3dReferred(const TPZMatPoisson3dReferred &copy) : TPZMatPoisson3d(copy){
    this->falpha = copy.falpha;
  }
  
  virtual TPZMaterial *NewMaterial(){
    return new TPZMatPoisson3dReferred(*this);
  }
  
  void SetAlpha(REAL alpha){
    this->falpha = alpha;
  }
  
  REAL GetAlpha(){
    return this->falpha;
  }
    
  
  virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                          TPZVec<REAL> &sol, TPZFMatrix &dsol,
                          REAL weight,TPZFMatrix &axes,
                          TPZFMatrix &phi, TPZFMatrix &dphi,
                          TPZFMatrix &ek, TPZFMatrix &ef);

  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
			    
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
    std::cout << "Please, implement me - " << __PRETTY_FUNCTION__ << std::endl;        
  }
				   
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
				     TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize){
    std::cout << "Please, implement me - " << __PRETTY_FUNCTION__ << std::endl;
  }    
             
};

#endif
