// -*- c++ -*-
// $Id: pzdiscgal.h,v 1.17 2007-06-08 00:09:31 cesar Exp $
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

  /** Fill material data parameter with necessary requirements for the
    * ContributeInterface method. Here, in base class, all requirements are considered
    * as necessary. Each derived class may optimize performance by selecting
    * only the necessary data.
    * @since April 10, 2007
    */
  virtual void FillDataRequirementsInterface(TPZMaterialData &data);

  /**
   * It computes a contribution to stiffness matrix and load vector at one integration point 
   * @param data [in]
   * @param weight [in]
   * @param ek [out] is the stiffness matrix
   * @param ef [out] is the load vector
   * @since April 16, 2007
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef) = 0;

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
  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) = 0;

  /**
   * It computes a contribution to residual vector at one BC integration point 
   * @param data [in]
   * @param weight [in]
   * @param ef [out] is the load vector
   * @param bc [in] is the boundary condition object
   * @since April 16, 2007
   */
  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef,TPZBndCond &bc);

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


  virtual int NStateVariables() = 0;


virtual void ContributeInterfaceErrors(TPZMaterialData &data,
                                        REAL weight,
                                        TPZVec<REAL> &nkL,
                                        TPZVec<REAL> &nkR,
                                        int &errorid){
  PZError << "Method not implemented\n";
}

virtual void ContributeInterfaceBCErrors(TPZMaterialData &data,
                                           REAL weight,
                                           TPZVec<REAL> &nk,
                                           TPZBndCond &bc,
                                           int &errorid){
  PZError << "Method not implemented\n";
}


};

#endif
