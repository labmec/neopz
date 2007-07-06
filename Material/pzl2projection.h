//$Id: pzl2projection.h,v 1.6 2007-07-06 18:27:42 tiago Exp $

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
