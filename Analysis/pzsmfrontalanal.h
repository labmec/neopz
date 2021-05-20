/**
 * @file
 * @brief Contains TPZSubMeshFrontalAnalysis class which implements the analysis for substructuring.
 */

#ifndef TPZSUBMESHFRONTANALYSIS_H
#define TPZSUBMESHFRONTANALYSIS_H


#include "TPZLinearAnalysis.h"
#include "pzmatred.h"
#include "TPZFrontMatrix.h"
class TPZSubCompMesh;
template<class TVar>
class TPZFront;
template<class TVar>
class TPZFMatrix;
#include "TPZSolutionMatrix.h"

/**
 * @brief Analysis for substructuring. Use a frontal matrix. \ref analysis "Analysis"
 * @ingroup analysis
 */
class TPZSubMeshFrontalAnalysis : public TPZLinearAnalysis  
{
private:
	/** @brief Solution vector */
	TPZSolutionMatrix fReferenceSolution;
	/** @brief The computational sub mesh */
	TPZSubCompMesh *fMesh;
	
	/** @brief The decomposition process and frontal matrix */
	TPZFront<STATE> *fFront;
protected:
	template<class TVar>
	void LoadSolutionInternal(TPZFMatrix<TVar> &mySol,
							  const TPZFMatrix<TVar> &myRhs,
							  const TPZFMatrix<TVar> &myRefSol,
							  const TPZFMatrix<TVar> &sol);
public:
	
	virtual void LoadSolution(const TPZFMatrix<STATE> &sol) override;
	
	/** @brief Constructor: create an object analysis from one mesh */
	TPZSubMeshFrontalAnalysis(TPZSubCompMesh *mesh);
	
	/** @brief Destructor */
	virtual ~TPZSubMeshFrontalAnalysis();
	
	/** @brief Run: assemble the stiffness matrix */
	void Run(std::ostream &out) override;
	
	/** @brief CondensedSolution: returns the condensed stiffness matrix - ek - and the condensed solution vector - ef */
	template<class TVar>
	void CondensedSolution(TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef);
	
	/** @brief Sets the front matrix */
	void SetFront(TPZFront<STATE> &front) { fFront = &front;}
	
};

extern template
void TPZSubMeshFrontalAnalysis::CondensedSolution<STATE>(TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
//TODOCOMPLEX
// extern template
// void TPZSubMeshFrontalAnalysis::CondensedSolution<CSTATE>(TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef);
#endif 
