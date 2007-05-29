//$Id: pzl2projection.h,v 1.4 2007-05-29 21:28:40 tiago Exp $

#ifndef PZL2PROJECTION_H
#define PZL2PROJECTION_H

#include <pzmaterial.h>

/**
 * Implements an L2 projection to constant solution values.
 * @since April 23, 2007
*/
class TPZL2Projection : public TPZMaterial{

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

  /** Contribute method
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /** To satisfy base class interface.
   */
  virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc){
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

  /**To create another material of the same type*/
  virtual TPZAutoPointer<TPZMaterial> NewMaterial();

};

#endif
