// TPZSubMeshFrontalAnalysis.h: interface for the TPZSubMeshFrontalAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZSUBMESHFRONTANALYSIS_H
#define TPZSUBMESHFRONTANALYSIS_H


#include "pzanalysis.h"
#include "pzmatred.h"
#include "TPZFrontMatrix.h"
class TPZSubCompMesh;
class TPZFront;

#include "pzfmatrix.h"	// Added by ClassView

/**
 * @brief Analysis for substructuring. Use a frontal matrix.
 * @ingroup analysis
 */
class TPZSubMeshFrontalAnalysis : public TPZAnalysis  
{
private:
	TPZFMatrix fReferenceSolution;
	
	TPZSubCompMesh *fMesh;
	
	TPZFront *fFront;
	
public:
	virtual void LoadSolution(TPZFMatrix &sol);
	/**
	 * @brief Constructor: create an object analysis from one mesh
	 **/
	TPZSubMeshFrontalAnalysis(TPZSubCompMesh *mesh);
	
	/**
	 * @brief Destructor
	 **/
	virtual ~TPZSubMeshFrontalAnalysis();
	
	
	/**
	 * @brief Run: assemble the stiffness matrix
	 **/
	void Run(std::ostream &out);
	
	/**
	 * @brief CondensedSolution: returns the condensed stiffness
	 *matrix - ek - and the condensed solution vector - ef
	 */
	void CondensedSolution(TPZFMatrix &ek, TPZFMatrix &ef);
	
	/**
	 * @brief Sets the front matrix
	 */
	void SetFront(TPZFront &front) { fFront = &front;}
	
};

#endif 
