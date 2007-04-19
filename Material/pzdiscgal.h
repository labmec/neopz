// -*- c++ -*-
// $Id: pzdiscgal.h,v 1.15 2007-04-19 11:43:12 tiago Exp $
#ifndef TPZDISCGALHPP
#define TPZDISCGALHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"

class TPZMaterialData;


/// This class defines the interface which material objects need to implement for discontinuous Galerkin formulations
class TPZDiscontinuousGalerkin  : public TPZMaterial {


  public :

  TPZDiscontinuousGalerkin();

  TPZDiscontinuousGalerkin(int nummat);

  /**copy constructor*/
  TPZDiscontinuousGalerkin(const TPZDiscontinuousGalerkin &copy);

  virtual ~TPZDiscontinuousGalerkin();

  char *Name();

  /**
   * It computes a contribution to stiffness matrix and load vector at one integration point 
   * @param data [in]
   * @param weight [in]
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @since April 16, 2007
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /**
   * It computes a contribution to residual vector at one integration point 
   * @param data [in]
   * @param weight [in]
   * @param ef [out] is the load vector
   * @since April 16, 2007
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef);

  /**
   * It computes a contribution to stiffness matrix and load vector at one BC integration point 
   * @param data [in]
   * @param weight [in]
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition object
   * @since April 16, 2007
   */
  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  /**
   * It computes a contribution to residual vector at one BC integration point 
   * @param data [in]
   * @param weight [in]
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition object
   * @since April 16, 2007
   */
  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef,TPZBndCond &bc);

///deprecated method
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                   TPZFMatrix &ek,TPZFMatrix &ef){
        PZError << "\nThis method ( " << __PRETTY_FUNCTION__ 
                << " ) was pure virtual before TPZMaterialData\n";
      }

///deprecated method
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                   TPZFMatrix &axesleft, TPZFMatrix &axesright,
                                   TPZFMatrix &ef);

///deprecated method
virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                 TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                 TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR, TPZFMatrix &axesleft,
                                 TPZFMatrix &axesright,
                                 TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize);

///deprecated method
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                     TPZFMatrix &axesleft,
                                     TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
        PZError << "\nThis method ( " << __PRETTY_FUNCTION__ 
                << " ) was pure virtual before TPZMaterialData\n";
      }

///deprecated method
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                     TPZFMatrix &phiL,TPZFMatrix &dphiL,
                                     TPZFMatrix axesleft,
                                     TPZFMatrix &ef,TPZBndCond &bc);

///deprecated method
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
