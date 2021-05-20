/**
 * @file
 */

#ifndef POROANALYSIS_H
#define POROANALYSIS_H

#include "ConsLaw/TPZConsLaw.h"
#include "pzelastoplasticanalysis.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include <iostream>

class TPZPoroElastoPlasticAnalysis : public TPZElastoPlasticAnalysis {
	
public:
	
	TPZPoroElastoPlasticAnalysis(TPZCompMesh *mesh,std::ostream &out);
	
	TPZPoroElastoPlasticAnalysis();
	
	virtual ~TPZPoroElastoPlasticAnalysis();
	
	virtual void Run(std::ostream &out,REAL tol, int numiter,
					 TPZPostProcAnalysis * ppAnalysis, int res);
	/** @brief Calls the appropriate sequence of methods to build a solution or a time stepping sequence */
	virtual void Run(std::ostream &out = std::cout) {
        TPZLinearAnalysis::Run(out);
    }
	
	/**
	 * @brief Informs the pertinent poroelastoplastic materials the current timestep
	 * @param deltaT [in] timestep
	 */
	virtual void SetDeltaT(const REAL deltaT);
	
	/**
	 * @brief Assembles and solves a L2 matrix to generate a constant pore pressure
	 * initial solution state
	 * Should be called only after assigning a solve method
	 */
	void SetInitialPorePressure(REAL initPress);
	
protected:
	
	/**
	 * @brief Reimplemented in order to manage the different kinds of varibles.
	 * Displacements are incremental and pressure is replaceable.
	 */
	virtual REAL AcceptSolution(const int ResetOutputDisplacements = 0);
	
	/** @brief Searches in the computational mesh for all the porous material instantiations */
	int FindPorousMaterials();
	
	/**
	 * @brief Sets the solution vector to be the one
	 * representing the Advanced State 
	 * or the advanced implicit solution.
	 */
	void SetAdvancedState();
	
	/**
	 * @brief Sets the solution vector to be the one
	 * representing the Current State 
	 * or the last explicit solution.
	 */
	void SetLastState();
	
	/**
	 * @brief Informs the Analysis class the time at which
	 * the current solution in the computational
	 * mesh belongs, so that the materials can
	 * choose whether to contribute implicitly
	 * or explicitly.
	 */
	void SetContributionTime(TPZContributeTime time);
	
protected:
	
	/** @brief Stores a list of Ids of all porous materials in the computational mehs */
	TPZStack<int> fPorousMaterialIds;
	
	/**
	 * @brief Stores the last state RHS assembled state, sparing computational time by avoiding
	 * evaluating the same contributions in each assemble() operation.
	 */
	TPZFMatrix<REAL> fRhsLast;
	
};

#endif

