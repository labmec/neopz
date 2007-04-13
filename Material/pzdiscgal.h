// -*- c++ -*-
// $Id: pzdiscgal.h,v 1.13 2007-04-13 15:17:54 tiago Exp $
#ifndef TPZDISCGALHPP
#define TPZDISCGALHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"


/// This class defines the interface which material objects need to implement for discontinuous Galerkin formulations
class TPZDiscontinuousGalerkin  : public TPZMaterial {


  public :

  TPZDiscontinuousGalerkin();

  TPZDiscontinuousGalerkin(int nummat);

  /**copy constructor*/
  TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy);

  virtual ~TPZDiscontinuousGalerkin();

  char *Name();

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                   TPZFMatrix &ek,TPZFMatrix &ef) = 0;

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                   TPZFMatrix &ef);

virtual void ContributeInterface(TPZVec<REAL> &x,
                                 TPZVec<REAL> &sol, TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                 TPZFMatrix &dsol, TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                 REAL weight,TPZVec<REAL> &normal,
                                 TPZFMatrix &phiL, TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                 TPZFMatrix &axes, TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                 TPZFMatrix &ek,TPZFMatrix &ef);

virtual void ContributeInterface(TPZVec<REAL> &x,
                                 TPZVec<REAL> &sol, TPZVec<REAL> &solL,TPZVec<REAL> &solR,
                                 TPZFMatrix &dsol, TPZFMatrix &dsolL, TPZFMatrix &dsolR,
                                 REAL weight,TPZVec<REAL> &normal,
                                 TPZFMatrix &phiL, TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                 TPZFMatrix &axes, TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                 TPZFMatrix &ef);

virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                 TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                 TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR, TPZFMatrix &axesleft,
                                 TPZFMatrix &axesright,
                                 TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize);

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                     TPZFMatrix &axesleft,
                                     TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) = 0;

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                     TPZFMatrix axesleft,
                                     TPZFMatrix &ef,TPZBndCond &bc);

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                     TPZFMatrix &axesleft,
                                     TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize);

  /**
   * Dicontinuous galerkin materials implement contribution of discontinuous elements and interfaces.
   * Interfaces may be conservative or not conservative. It is important to agglomeration techniques
   * when using multigrid pre-conditioner. Conservative interfaces into agglomerate elements do not
   * need to be computed. However non-conservative interfaces must be computed in all multigrid levels.
   * Default is non-conservative, because of the computation of a conservative interface into an agglomerate
   * does not ruin the solution.
   * @since Feb 05, 2004
   */
  virtual int IsInterfaceConservative();

  /** Compute interface jumps
   * values[1] = (solleft - solright)^2
   * values[2] = (dsolleft - dsolright)^2
   * values[0] = values[1] + values[2]
   * @since Feb 14, 2006
   */
  virtual void InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu, TPZVec<REAL> &leftNormalDeriv,
                         TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                         TPZVec<REAL> &values);

  /** Compute interface jumps from element to Dirichlet boundary condition
   * values[1] = (solleft - solright)^2
   * values[2] = (dsolleft - dsolright)^2
   * values[0] = values[1] + values[2]
   * @since Mar 08, 2006
   */
  virtual void BCInterfaceJumps(TPZVec<REAL> &leftu,
                                TPZBndCond &bc,
                                TPZVec<REAL> &values);

};

#endif
