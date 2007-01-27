// -*- c++ -*-

//$Id: pzpoisson3dreferred.h,v 1.4 2007-01-27 14:49:27 phil Exp $

#ifndef MATPOISSON3DREFERREDH
#define MATPOISSON3DREFERREDH

#include <iostream>
#include "pzpoisson3d.h"
#include "pzfmatrix.h"

/**
 * This class implements a version of TPZMatPoisson3d where the convection term is given at each integration point
 * from a previous calculation.
 * The convection term fC * fConvDir = - fAlpha * grad(sol) where grad(sol) is the gradient of the previous solution.
 */
class TPZMatPoisson3dReferred : public TPZMatPoisson3d {

protected:

  REAL falpha;
  
  /** SetConvectionTerm
   */
  void SetConvectionTerm(TPZFMatrix &dsol, TPZFMatrix &axes);
  
  /** SeConvectionTerm for ContributeInterface methods
   * It expect dsolL and dsolR to be dSol/dX, i.e. the derivatives 
   * with respect to the global coordinates.
   */
  void SetConvectionTermInterface(TPZFMatrix &dsolL, TPZFMatrix &dsolR);

public:
 
  TPZMatPoisson3dReferred(int nummat, int dim);

  virtual ~TPZMatPoisson3dReferred();
  
//  virtual int HasForcingFunction() {return true;}

  TPZMatPoisson3dReferred(const TPZMatPoisson3dReferred &copy) : TPZMatPoisson3d(copy){
    this->falpha = copy.falpha;
  }
  
  virtual TPZAutoPointer<TPZMaterial> NewMaterial(){
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
       TPZFMatrix &axesleft, TPZFMatrix &axesright,
       TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			    TPZFMatrix &phiL,TPZFMatrix &dphiL, 
       TPZFMatrix &axesleft,
       TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
			    
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
       TPZFMatrix &axesleft, TPZFMatrix &axesright,
       TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){
    std::cout << "Please, implement me - " << __PRETTY_FUNCTION__ << std::endl;        
  }
				   
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
				     TPZFMatrix &phiL,TPZFMatrix &dphiL, 
         TPZFMatrix &axesleft,
         TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize){
    std::cout << "Please, implement me - " << __PRETTY_FUNCTION__ << std::endl;
  }    
             
};

#endif
