/**
 * @file
 * @brief Contains the declaration of TMBAdaptInterface class
 */

#ifndef TMBADAPTINTERFACE_H
#define TMBADAPTINTERFACE_H

#include "pzcmesh.h"
#include "pzvec.h"
#include "pzstack.h"

/**
 * @brief Interface for refining meshes based on a simple error estimator
 */
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
	TMBAdaptInterface(TPZCompMesh *cmesh, int nstate, TPZVec<int> &dimstate,
					  TPZVec<int> &statetoanalyse, TPZVec<REAL> &sol);
	
	/** @brief Defines the maximum error per state variable */
	void SetMaxMinError (TPZVec<REAL> &maxervec, TPZVec<REAL> &minervec);
	
	/**
	 * @brief Defines the minimum dimension whose elements could have. \n
	 * The elements shortest than this value are not divided
	 */
	void SetMinElSize (REAL size) {fMinElSize = size;}
	
	/** @brief Defines the maximum refinement level of each element */
	void SetMaxRefLevel (int level) {fMaxLevel = level;}
	
	/** @brief Default destructor */
	virtual ~TMBAdaptInterface ();
	
	/** @brief Return the adatpted mesh */
	TPZCompMesh *GetAdaptedMesh(TPZVec<int> &erind, TPZVec<int> &erantype, bool level_check,
                                int reftype = 0, int side_ref_dim = -1, int side_ref_state = -1 );
	
	/**
	 * @brief Returns the side to refine.
	 * @param cel element to analyse
	 * @param sidedim dimension of the sides which will be analysed 
	 */
	int GetRefSide(TPZCompEl *cel, int sidedim);
	
protected:
	/** @brief Computational mesh to analyse error */
	TPZCompMesh *fMesh;
	/** @brief Solution vector */
	TPZVec<REAL> *fSolution;
	/** @brief State variables in solution vector */
	TPZVec<int> fDimState;
	/** @brief StateVariables to analyse */
	TPZVec<int> fStateAnalyse;   
	/** @brief Error indicator per state variable */
	TPZVec<int> fErInd;
	/** @brief Type of anlaysis 	
	 * \li 0 - absolute error value comparison
	 * \li 1 - relative error value comparison
	 */
	TPZVec<int> fErAnType;
	/**
	 * @brief Maximum admissible error per state variable. \n
	 * If the error obtained in one element is greater than this value
	 * the element will be marked to refine
	 */
	TPZVec<REAL> fMaxError;
	/**
	 * @brief Minimum admissible error per state variable. \n
	 * If the error obtained in one element is lower than this value
	 * the element will be marked to aglomerate
	 */
	TPZVec<REAL> fMinError;
	/**
	 * @brief Minimum size of an element. If the element is shortest than fMinSize
	 * the element is not refined even the error indicator select it.
	 */
	REAL fMinElSize;
	/** @brief Maximum level of refinement for all elements */
	int fMaxLevel;
	
	/**
	 * @brief Store the status from the error indicator for each element. \n
	 * The status should be:
	 * \li 0 : nothing to do.
	 * \li 1 : marked to refine.
	 * \li (-1) : marked to aglomerate.
	 */
	TPZVec<int> fMark;
	
protected:
	/** @brief Returns the indexes of the elements indicated to be refined */
	void MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side, int sidedim = -1, int sidestate = -1);
	
private:
	/** @brief Refine the specified element if the verifications of level, size and neighbors level are done. */
 	void Refine (TPZCompEl *el, bool leve_check);
	
	/** @brief Used to make the correct cast and change the refinement pattern */
	bool ChangeRefPattern(TPZGeoEl * gel, int side);
	
	/** @brief Aglomerate the elements if the verifications of subelements refinement and neighbors are done. */
	void Aglomerate (TPZCompEl *el);
	
	int MaxLevel(TPZCompMesh *mesh);
	
	/** @brief Implements side based refinement for the specified element. */
	void SideRefine(TPZCompEl *cel,int side);
};

#endif
