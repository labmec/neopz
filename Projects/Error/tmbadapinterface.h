#ifndef TMBADAPTINTERFACE_H
#define TMBADAPTINTERFACE_H

#include "pzcmesh.h"
#include "pzvec.h"
#include "pzstack.h"

/// Interface for refining meshes based on a simple error estimator
class TMBAdaptInterface {

 public:

  /**
   * Simple constructor
   @param cmesh Computational PZ mesh
   @param nstate number of state variables in the solution vector
   @param dimstate dimension of each state variable
   @param statetoanalyse indexes of the state variables to analyse
   @param sol solution vector
   */
  TMBAdaptInterface (	TPZCompMesh *cmesh, int nstate, TPZVec<int> &dimstate,
  										TPZVec<int> &statetoanalyse, TPZVec<REAL> &sol);

  /**
   * Defines the state varibles to drive the adaptive process
   * @param state Ids of the state variables that will be used in error analysis -\
   * the ids must have the format of the enum of CFDK state variables
   */
	//  void SetStateVariables(TPZVec<int> &state);
  
  /**
   * Defines the maximum error per state variable
   */
  void SetMaxMinError (TPZVec<REAL> &maxervec, TPZVec<REAL> &minervec);

  /**
   * Defines the minimum dimension whose elements could have.\
   * The elements shortest than this value are not divided
   */
  void SetMinElSize (REAL size) {fMinElSize = size;}

  /**
   * Defines the maximum refinement level of each element
   */
   void SetMaxRefLevel (int level) {fMaxLevel = level;}

  /**
   * Default destructor
   */
  virtual ~TMBAdaptInterface ();

  /**
   * Return the adatpted mesh
   */
   TPZCompMesh *GetAdaptedMesh(TPZVec<int> &erind, TPZVec<int> &erantype, bool level_check,
                                int reftype = 0, int side_ref_dim = -1, int side_ref_state = -1 );
  /** Returns the side to refine.
@param cel element to analyse
@param sidedim dimension of the sides which will be analysed */
  int GetRefSide(TPZCompEl *cel, int sidedim);

 protected:
  /**
   * Computational mesh to analyse error
   */
  TPZCompMesh *fMesh;
  /**
   * Solution
   */
  TPZVec<REAL> *fSolution;
  /**
   * State variables in solution vector
   */
  TPZVec<int> fDimState;
  /**
   * StateVariables to analyse
   */
  TPZVec<int> fStateAnalyse;   
  /**
   * Error indicator per state variable
   */
  TPZVec<int> fErInd;
  /**
   * Type of anlaysis 	0 - absolute error value comparison
  	 										1 - relative error value comparison
	 */
	TPZVec<int> fErAnType;
  /**
   * Maximum admissible error per state variable
   * If the error obtained in one element is greater than this value
   * the element will be marked to refine
   */
  TPZVec<REAL> fMaxError;
  /**
   * Minimum admissible error per state variable
   * If the error obtained in one element is lower than this value
   * the element will be marked to aglomerate
   */
  TPZVec<REAL> fMinError;
  /**
   * Minimum size of an element. If the element is shortest than fMinSize
   * the element is not refined even the error indicator select it.
   */
  REAL fMinElSize;
  /**
   * Maximum level of refinement for all elements
   */
	int fMaxLevel;
	
  /**
   * Store the status from the error indicator for each element
   * The status should be:
   * - 0 : nothing to do.
   * - 1 : marked to refine.
   * - (-1) : marked to aglomerate.
   */
  TPZVec<int> fMark;

  /**
   * Type of Refinement Pattern
   * - 0 : Uniform
   * - 1 : Side Based (Rib)
   */
   //int fRefType;


 protected:
  /**
   * Transport the CFDK solution to PZ mesh
   */
  //void SetPZMeshSolution(TPZVec<REAL> &sol);
  /**
   * Find the volumes on PZ mesh and fill the stack with their indexes
   */
  //void FillVolIndex(TPZStack<int> &volindex);
   /**
   * Returns the indexes of the elements indicated to be refined and their orders of refinement.
   * This method is not to be used in CFD Embraer Project.
   */
  //virtual void MarkedElListHP(TPZStack<TPZGeoEl*> &geoelvec, TPZStack<int> porders);

  /**
   * Returns the indexes of the elements indicated to be refined
   */
  void MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side, int sidedim = -1, int sidestate = -1);

 private:
	/**
	 * Refine the specified element if the verifications of level, size and neighbors level are done.
	 */
 	void Refine (TPZCompEl *el, bool leve_check);

  /**
   * Used to make the correct cast and change the refinement pattern
   */
  bool ChangeRefPattern(TPZGeoEl * gel, int side);

  /**
   * Aglomerate the elements if the verifications of subelements refinement and neighbors are done.
   */
  void Aglomerate (TPZCompEl *el);

  int MaxLevel(TPZCompMesh *mesh);

  /** Implements side based refinement for the specified element. */
  void SideRefine(TPZCompEl *cel,int side);
};
#endif
