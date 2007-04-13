// -*- c++ -*-
// $Id: pzdiscgal.h,v 1.11 2007-04-13 14:29:07 cesar Exp $
#ifndef TPZDISCGALHPP
#define TPZDISCGALHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"


/// This class defines the interface which material objects need to implement for discontinuous Galerkin formulations
class TPZDiscontinuousGalerkin  : public TPZMaterial {


  public :

  TPZDiscontinuousGalerkin() : TPZMaterial(){}

  TPZDiscontinuousGalerkin(int nummat) : TPZMaterial(nummat) {}

  /**copy constructor*/
  TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy) : TPZMaterial(copy) {}

  virtual  ~TPZDiscontinuousGalerkin();

  char *Name() { return "TPZDiscontinuousGalerkin"; }

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
       TPZFMatrix &axesleft, TPZFMatrix &axesright,
       TPZFMatrix &ek,TPZFMatrix &ef) = 0;

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
       TPZFMatrix &axesleft, TPZFMatrix &axesright,
				   TPZFMatrix &ef)
				   {
				      TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
				      ContributeInterface(x, solL, solR, dsolL, dsolR,
				                          weight, normal, phiL, phiR,
							  dphiL, dphiR,axesleft,axesright, ek, ef);
				   }

virtual void ContributeInterface(TPZVec<REAL> &x,
                                 TPZVec<REAL> &sol, TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                 TPZFMatrix &dsol, TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                 REAL weight,TPZVec<REAL> &normal,
                                 TPZFMatrix &phiL, TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                 TPZFMatrix &axes, TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                 TPZFMatrix &ek,TPZFMatrix &ef)
{
  ContributeInterface(x,solL,solR,dsolL,dsolR,weight,normal,phiL,phiR,dphiL,dphiR,axesleft,axesright,ek,ef);
}

virtual void ContributeInterface(TPZVec<REAL> &x,
                                 TPZVec<REAL> &sol, TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                 TPZFMatrix &dsol, TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                 REAL weight,TPZVec<REAL> &normal,
                                 TPZFMatrix &phiL, TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                 TPZFMatrix &axes, TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                 TPZFMatrix &ef)
{
  TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
  ContributeInterface(x,solL,solR,dsolL,dsolR,weight,normal,phiL,phiR,dphiL,dphiR,axesleft,axesright,ek,ef);
}

virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR, TPZFMatrix &axesleft,
       TPZFMatrix &axesright,
				   TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){

   this->ContributeInterface(x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, axesleft, axesright, ek, ef);

}

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			    TPZFMatrix &phiL,TPZFMatrix &dphiL,
       TPZFMatrix &axesleft,
       TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) = 0;

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			    TPZFMatrix &phiL,TPZFMatrix &dphiL,
       TPZFMatrix axesleft,
       TPZFMatrix &ef,TPZBndCond &bc)
			    {
                                TPZFMatrix ek(ef.Rows(),ef.Rows(),0.);
				ContributeBCInterface(x, solL,  dsolL, weight, normal,
			                              phiL, dphiL, ek, ef, bc);
			    }

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
				     TPZFMatrix &phiL,TPZFMatrix &dphiL,
         TPZFMatrix &axesleft,
         TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize){

     this->ContributeBCInterface(x, solL, dsolL, weight, normal, phiL, dphiL, ek, ef, bc);

  }

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

  /** Compute interface jumps
   * values[1] = (solleft - solright)^2
   * values[2] = (dsolleft - dsolright)^2
   * values[0] = values[1] + values[2]
   * @since Feb 14, 2006
   */
  virtual void InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu, TPZVec<REAL> &leftNormalDeriv,
                         TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                         TPZVec<REAL> &values){
    PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
  }

  /** Compute interface jumps from element to Dirichlet boundary condition
   * values[1] = (solleft - solright)^2
   * values[2] = (dsolleft - dsolright)^2
   * values[0] = values[1] + values[2]
   * @since Mar 08, 2006
   */
  virtual void BCInterfaceJumps(TPZVec<REAL> &leftu,
                                TPZBndCond &bc,
                                TPZVec<REAL> &values){
    PZError << __PRETTY_FUNCTION__ << " - method not implemented in derived class" << std::endl;
  }

};

#endif
