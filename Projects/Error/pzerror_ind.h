#ifndef TPZERROR_IND_H
#define TPZERROR_IND_H

#include "pzcmesh.h"
#include "pzvec.h"
#include "pzstack.h"


//class TPZCompMesh;


enum MErrorIndicator {
  EHeidi,     // 0
  EXubin_Shih // 1
};


class TPZErrorIndicator {

 public:

  /**
   * Simple constructor
   * @param nstate number of state variables in solution vector
   * @param dimstate dimension of each state variable
   * @param statetoanalyse index of the state variables to analyse
   * @param sol vector containing solution that will be analysed
   * @param cmesh compuatational PZ mesh. Used when topological information is required
   */
  TPZErrorIndicator (int nstate, TPZVec<int>& dimstate, TPZVec<int>& statetoanalyse,  TPZVec<REAL> &sol, TPZCompMesh *cmesh );

  /**
   * Default destructor
   */
  virtual ~TPZErrorIndicator ();

  /**
   * Defines the computational mesh 
   * @param cmseh computational PZ mesh 
   */
  void SetMesh(TPZCompMesh *cmesh);

  /**
   * Defines the solution vector and the number of state variables
   * @param sol solution vector
   * @param nstate number of state variables
   * @param dimstate dimension of each state variable
   */
  void SetSolution(TPZVec<REAL> &sol, int nstate, TPZVec<int> &dimstate, TPZVec<int> &statetoanalyse);

  /**
   * Defines the maximum error for each state variable
   */
  void SetError (TPZVec<REAL> &maxerror, TPZVec<REAL> &minerror, TPZVec<int> &erantype);

  /**
   * Returns the indexes of the elements indicated to be refined
   * and their orders of refinement
   */
//  virtual void MarkedElListHP(TPZStack<TPZGeoEl*> &geoelvec, TPZStack<int> porders)=0;

  /**
   * Returns the indexes of the elements indicated to be refined
   */
  virtual void MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side,int sidedim, int sidestate) = 0;

 
 protected:

  /**
   * Computational mesh to analyse error
   */
  TPZCompMesh *fMesh;

  /**
   * Solution
   */
  TPZVec<REAL> fSolution;

  /**
   * Number of state variables in solution vector
   */
  int fNState;

  /**
   * Indexes of the state variables to anlyse
   */
  TPZVec<int> fState;

  /**
   * Dimension of the elements per state variable
   */
  TPZVec<int> fDim;

  /**
   * Maximum admissible error per state variable
   * If the error obtained in one element is greater than this value
   * the element will be marked to refine
   */
  TPZVec<REAL> fMaxError;

  /**
   * Minimum error -->> Used to aglomerate the lower error elements
   */
  TPZVec<REAL> fMinError;

  /**
   * Error analysis type per state variable
   */
  TPZVec<int> fErAnType;

  /**
   * Number of elements to analyse
   */
  int fNElements;

  /**
   * Contains the number of data per element. Note that each element has a \
   number of state variables and each state variable has a dimension.
  */
  int fNDataEl;
  
  /**
   * Return int the vector perm the ascending  index ordering of the vector sol
   */
  void Sort(TPZFMatrix &error, TPZFMatrix &perm);

  /**
   * Position index in the solution vector for a given dimension of a given state\
   * variable of the specified element
   * @param elem element
   * @param state state variable index
   * @param dim dimension
   */
  int Index(int elem, int state, int dim);


  /** Returns the side to refine.
   * @param cel element to analyse
   * @param sidedim dimension of the sides which will be analysed
   * @param sidestate statevariable thats define the side to refine
   */
  int GetRefSide(TPZCompEl *cel, int sidedim, int sidestate, TPZMatrix *mat);

 private:

};
#endif
