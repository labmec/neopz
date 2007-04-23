//$Id: pzl2projection.h,v 1.1 2007-04-23 19:02:56 tiago Exp $

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

  /** Contribute method
   */
  virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  /** Returns problem dimension 
   */
  virtual int Dimension(){ return this->fDim; }

  /** Returns number of state variables
   */
  virtual int NStateVariables(){ return this->fNStateVars; }

};

#endif
