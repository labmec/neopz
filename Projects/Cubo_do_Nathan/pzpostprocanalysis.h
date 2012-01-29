//$Id: pzpostprocanalysis.h,v 1.5 2010-10-20 18:41:37 diogo Exp $
#ifndef PZPOSTPROCANALYSIS_H
#define PZPOSTPROCANALYSIS_H

#include "pzanalysis.h"
#include "pzcompel.h"
#include "TPZGeoElement.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzstrmatrix.h"


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
		
static void SetAllCreateFunctionsPostProc(TPZCompMesh *cmesh);
//static void SetAllCreateFunctionsContinuous(TPZCompMesh *cmesh);
void AutoBuildDisc();
    
protected:
	
	TPZAnalysis * fpMainAnalysis;

/**
 * TPZCompElPostProc<TCOMPEL> creation function setup
 */

public:

	
static TPZCompEl * CreatePostProcDisc(  TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
	
static TPZCompEl * CreatePointEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreateLinearEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreateQuadEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreateTriangleEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreateCubeEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreatePyramEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreateTetraEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);
static TPZCompEl * CreatePrismEl( TPZGeoEl *gel, TPZCompMesh &mesh, int &index);

};

#endif

