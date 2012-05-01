/**
 * @file
 * @brief Contains the declaration of the TPZErrorIndicator class
 */

#ifndef TPZERROR_IND_H
#define TPZERROR_IND_H

#include "pzcmesh.h"
#include "pzvec.h"
#include "pzstack.h"


enum MErrorIndicator {
	EHeidi,     // 0
	EXubin_Shih // 1
};


class TPZErrorIndicator {
	
public:
	
	/**
	 * @brief Simple constructor
	 * @param nstate number of state variables in solution vector
	 * @param dimstate dimension of each state variable
	 * @param statetoanalyse index of the state variables to analyse
	 * @param sol vector containing solution that will be analysed
	 * @param cmesh compuatational PZ mesh. Used when topological information is required
	 */
	TPZErrorIndicator (int nstate, TPZVec<int>& dimstate, TPZVec<int>& statetoanalyse,  TPZVec<REAL> &sol, TPZCompMesh *cmesh );
	
	/** @brief Default destructor */
	virtual ~TPZErrorIndicator ();
	
	/**
	 * @brief Defines the computational mesh 
	 * @param cmesh computational PZ mesh 
	 */
	void SetMesh(TPZCompMesh *cmesh);
	
	/**
	 * Defines the solution vector and the number of state variables
	 * @param sol solution vector
	 * @param nstate number of state variables
	 * @param dimstate dimension of each state variable
	 * @param statetoanalyse [out]
	 */
	void SetSolution(TPZVec<REAL> &sol, int nstate, TPZVec<int> &dimstate, TPZVec<int> &statetoanalyse);
	
	/** @brief Defines the maximum error for each state variable */
	void SetError (TPZVec<REAL> &maxerror, TPZVec<REAL> &minerror, TPZVec<int> &erantype);
	
	/** @brief Returns the indexes of the elements indicated to be refined */
	virtual void MarkedElListH(TPZVec<int> &elindex, TPZVec<int> &side,int sidedim, int sidestate) = 0;
	
	
protected:
	
	/** @brief Computational mesh to analyse error */
	TPZCompMesh *fMesh;
	
	/** @brief Solution */
	TPZVec<REAL> fSolution;
	
	/** @brief Number of state variables in solution vector */
	int fNState;
	
	/** @brief Indexes of the state variables to anlyse */
	TPZVec<int> fState;
	
	/** @brief Dimension of the elements per state variable */
	TPZVec<int> fDim;
	
	/**
	 * @brief Maximum admissible error per state variable
	 * @note If the error obtained in one element is greater than this value
	 * the element will be marked to refine
	 */
	TPZVec<REAL> fMaxError;
	
	/** @brief Minimum error -->> Used to aglomerate the lower error elements */
	TPZVec<REAL> fMinError;
	
	/** @brief Error analysis type per state variable */
	TPZVec<int> fErAnType;
	
	/** @brief Number of elements to analyse */
	int fNElements;
	
	/**
	 * @brief Contains the number of data per element. Note that each element has a 
	 * number of state variables and each state variable has a dimension.
	 */
	int fNDataEl;
	
	/** @brief Return int the vector perm the ascending  index ordering of the vector sol */
	void Sort(TPZFMatrix &error, TPZFMatrix &perm);
	
	/**
	 * @brief Position index in the solution vector for a given dimension of a given state\
	 * variable of the specified element
	 * @param elem element
	 * @param state state variable index
	 * @param dim dimension
	 */
	int Index(int elem, int state, int dim);
	
	
	/**
	 * @brief Returns the side to refine.
	 * @param cel element to analyse
	 * @param sidedim dimension of the sides which will be analysed
	 * @param sidestate statevariable thats define the side to refine
	 * @param mat matrix with result
	 */
	int GetRefSide(TPZCompEl *cel, int sidedim, int sidestate, TPZMatrix *mat);
	
};
#endif
