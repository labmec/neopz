/***************************************************************************
                          tmbxubinshih.h  -  description
                             -------------------
    begin                : Thu Oct 23 2003
    copyright            : (C) 2003 by cesar
    email                : cesar@labmec.fec.unicamp.br
 ***************************************************************************/
#ifndef TMBXUBINSHIH_H
#define TMBXUBINSHIH_H

#include "pzerror_ind.h"

/**Implements the Xubin-Shih Error indicator
  *@author Edimar Cesar Rylo
  *@since october 2003
  */

class TMBXubinShih : public TPZErrorIndicator  {

public: 
  /**
   * Creates a TMBXubinShih object
   * @param cmesh: a PZ computational mesh
   * @param nstate: number of state variables in solution vector
   * @param dimstate: dimension of each state variable
   * @param statetoanalyse: indexes of the state variables to\
                            evaluates the indicator
   * @param sol: solution vector
   * @param source: source error indicator = 1; location error = 0
   */
  TMBXubinShih(TPZCompMesh *cmesh, int nstate, TPZVec<int> &dimstate,
                            TPZVec<int> &statetoanalyse,TPZVec<REAL> &sol,
                            int source = 0);

  /**
   * @see TPZErrorIndicator class documentation
   */
  virtual void MarkedElListH(TPZVec<int> &eleindex, TPZVec<int> &side, int sidedim = -1, int sidestate = -1);

  /**
   * Destructor
   */
  ~TMBXubinShih();


protected: // Protected methods
  /**
   * Returns the angle between two given vectors
   * @param u first vector
   * @param v second vector
   */
  REAL EvaluateCosAngle(TPZVec<REAL> u, TPZVec<REAL> v);

  /**
   * Evaluates the Xubin-Shih location error indicator over the mesh
   */
  void ElLocationIndicator(TPZVec<int> &statetoanalyse, TPZVec<REAL> &results);
  /**
   * Evaluates the Xubin-Shih location error indicator over an element
   */
  void ElLocationIndicator(TPZCompEl *cel, TPZVec<int> &statetoanalyse, TPZVec<REAL> &results);

  /**
   * Evaluates the source error indicator over the mesh
   */
  void SourceIndicator(TPZVec<int> &statetoanalyse, TPZVec<REAL> &results);

  /**
   * Evaluates the source error indicator over an element
   */
  void ElSourceIndicator(TPZCompEl *cel, TPZVec<int> &statetoanalyse, TPZVec<REAL> &results);

  /**
   * Evaluates the vector L2 norm
   */
  REAL L2Norm(TPZVec<REAL> &vec);
 
  REAL cos45;

  /**
   * Set the source error indicator if defined. Elsewhere defines the
   * location error indicator
   */
  int fSource;

};

#endif
