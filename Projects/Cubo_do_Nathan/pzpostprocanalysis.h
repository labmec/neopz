//$Id: pzpostprocanalysis.h,v 1.3 2009-05-16 04:25:34 erick Exp $
#ifndef PZPOSTPROCANALYSIS_H
#define PZPOSTPROCANALYSIS_H

#include "pzanalysis.h"
#include "pzcompel.h"
#include "TPZGeoElement.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <iostream>
#include <string>


/**
 * The Post Processing Analysis makes the interface among the computational
 * analysis and itself, being also in charge of getting the current solution
 * from the main analysis.
 */

class TPZPostProcAnalysis : public TPZAnalysis {

public:

TPZPostProcAnalysis(TPZAnalysis * pRef);

TPZPostProcAnalysis();

virtual ~TPZPostProcAnalysis();
	
/**
 *	Assemble() blank implementation in order to avoid its usage. In such an Analysis
 * class the Assemble() method is useless.
 */
virtual  void Assemble();

/**
 *	Solve() blank implementation in order to avoid its usage. In such an Analysis
 * class the Solve() method is useless.
 */

virtual void Solve();

/** 
 * TransferSolution is in charge of transferring the solution from the base analysis/mesh
 * to the current post processing mesh, solving for the dof to provide a LSM extrapolation.
 */
void TransferSolution();

/**
 * Informs the material IDs and the variable names that are to be postprocessed, if
 * matching. Should be called right after the class instantiation
 */
void SetPostProcessVariables(TPZVec<int> & matIds, TPZVec<std::string> &varNames);
		
static void SetAllCreateFunctionsPostProc();
		
protected:

	TPZAnalysis * fpMainAnalysis;

/**
 * TPZCompElPostProc<TCOMPEL> creation function setup
 */

public:
	
static TPZCompEl * CreatePostProcDisc(  TPZGeoEl *gel, TPZCompMesh &mesh, int &index);

};

#endif

