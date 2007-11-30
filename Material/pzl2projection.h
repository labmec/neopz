//$Id: pzl2projection.h,v 1.7 2007-11-30 11:36:43 phil Exp $

#ifndef PZL2PROJECTION_H
#define PZL2PROJECTION_H

#include "pzmaterial.h"
#include "pzdiscgal.h"

/**
 * Implements an L2 projection to constant solution values.
 * @since April 23, 2007
*/
class TPZL2Projection : public TPZDiscontinuousGalerkin{

private:

  /** Problem dimension */
  int fDim;

 /** Number of state variables */
  int fNStateVars;

 /** Constant solution vector */
  TPZVec<REAL> fSol;

 /** Argument defining this material is a referred material */
  bool fIsReferred;

public:

  /** Class constructor 
   * @param id material id
   * @param dim problem dimension
   * @param nstate number of state variables
   * @param sol constant solution vector
   */
  TPZL2Projection(int id, int dim, int nstate, TPZVec<REAL> &sol);

  /** Class destructor
   */
  ~TPZL2Projection();

  /** Copy constructor
   */
  TPZL2Projection(const TPZL2Projection &cp);

  /** Solution indices of post-processing
    */
  enum ESolutionVars { ENone = 0, ESolution = 1 };

  /** Contribute method
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /** To satisfy base class interface.
   */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc){
    //NOTHING TO BE DONE HERE
  }

  /** To satisfy base class interface.
   */
  virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
    //NOTHING TO BE DONE HERE
  }

  /** To satisfy base class interface.
   */
  virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
    //NOTHING TO BE DONE HERE
  }

  /** Returns problem dimension 
   */
  virtual int Dimension(){ return this->fDim; }

  /** Returns number of state variables
   */
  virtual int NStateVariables(){ return this->fNStateVars; }
      /**
   * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
   * @param data[in] stores all input data
   * @param weight[in] is the weight of the integration rule
   * @param ek[out] is the stiffness matrix
   * @param ef[out] is the load vector
   * @param bc[in] is the boundary condition material
   * @since April 16, 2007
       */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc)
  {
  }


  /** Define if material is referred or not */
  void SetIsReferred(bool val);

  /** To create another material of the same type */
  virtual TPZAutoPointer<TPZMaterial> NewMaterial();

  /** It returns the variable index associated with the name */
  virtual int VariableIndex(char *name);

  /** It returns the number of variables associated with the variable
   *  indexed by var.  
   * var is obtained by calling VariableIndex
   */
  virtual int NSolutionVariables(int var);

  /** It returns the solution associated with the var index based on
   * the finite element approximation
   */
  virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                        TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);

};

#endif
