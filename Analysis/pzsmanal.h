/**
 * @file
 * @brief Contains TPZSubMeshAnalysis class which implements the analysis procedure to computational sub mesh.
 */

#ifndef TPZSUBMESHANALYSIS_H
#define TPZSUBMESHANALYSIS_H


#include "TPZLinearAnalysis.h"
#include "pzmatred.h"
class TPZSubCompMesh;

#include "pzfmatrix.h"

/** 
 * @brief Analysis procedure to computational sub mesh. \ref analysis "Analysis"
 * @ingroup analysis
 */
class TPZSubMeshAnalysis : public TPZLinearAnalysis  
{
private:
	/** @brief Solution vector */
	TPZFMatrix<STATE> fReferenceSolution;
	
	/** @brief Stiffness matrix to sub mesh */
	TPZAutoPointer<TPZMatrix<STATE> > fReducableStiff;
	/** @brief The computational sub mesh */
	TPZSubCompMesh *fMesh;
	
public:
	virtual void LoadSolution(const TPZFMatrix<STATE> &sol) override;
	/** @brief Constructor: create an object analysis from one mesh */
	TPZSubMeshAnalysis(TPZSubCompMesh *mesh = 0);
	
	/** @brief Destructor */
	virtual ~TPZSubMeshAnalysis();
	
	TPZAutoPointer<TPZMatrix<STATE> > Matrix()
	{
		return fReducableStiff;
	}
    
    /** @brief Set the computational mesh of the analysis. */
    virtual void SetCompMesh(TPZCompMesh * mesh, bool mustOptimizeBandwidth) override;
    
	
	/** @brief Run: assemble the stiffness matrix */
	void Run(std::ostream &out) override;
	
	/** @brief CondensedSolution: returns the condensed stiffness matrix - ek - and the condensed solution vector - ef */
	void CondensedSolution(TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/** @brief Assemble the global stiffness matrix and put it into the reducable stiffness matrix */
	virtual void Assemble() override;

        int ClassId() const override;

    /** @brief compute the reduced right hand side using the current stiffness. Abort if there is no stiffness computed */
    void ReducedRightHandSide(TPZFMatrix<STATE> &rhs);
private:
	template<class TVar>
	void AssembleInternal();
};

#endif 
