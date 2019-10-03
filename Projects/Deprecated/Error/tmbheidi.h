#ifndef TMBHEIDI_H
#define TMBHEIDI_H

#include "pzerror_ind.h"
#include "TPZCompElDisc.h"
#include "pzvec.h"
#include "pzstack.h"

class TPZCompMesh;

/// Error indicator for finite volume solutions
class TMBHeidi : public TPZErrorIndicator{

 public:

  /**
   * Simple constructor
   */
  TMBHeidi (int nstate, TPZVec<int> &dimstate, TPZVec<int> &statetoanlyse, TPZVec<REAL> &sol);

  /**
   * Default destructor
   */
 virtual  ~TMBHeidi ();

  /**
   * Returns the indexes of the elements indicated to be refined
   */
  virtual void MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side,int sidedim = -1, int sidestate = -1);


 protected:

  void FindMaxMin(TPZVec<REAL> &max, TPZVec<REAL> &min, TPZVec<REAL> &delta);



};
#endif
