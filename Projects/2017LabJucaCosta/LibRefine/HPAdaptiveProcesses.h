//
//  SimilarUniformRefinements.hpp
//  PZ
//
//  Created by labmec on 08/10/17.
//
//

#ifndef HPADAPTIVEPROCESSES_HPP
#define HPADAPTIVEPROCESSES_HPP

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"

#include "pzgeoelbc.h"

#include "pzlog.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZFrontStructMatrix.h"
#include "TPZParFrontStructMatrix.h"

#include "TPZParSkylineStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzfstrmatrix.h"

#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzpoisson3d.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZRefPatternTools.h"

#include <time.h>
#include <stdio.h>
#include <math.h>

#include <fstream>
#include <cmath>


/* 1. Functions to calculate Errors from Approximate solution and Referenced solution
 */
 /**
 * Get Global L2 Error for solution and the L2 error for each element.
 * Return the maxime L2 error by elements. Also return in MinErrorByElement argument the minime L2 error for all elements of the mesh.
 */
bool ProcessingErrorUAndDUKnowingExactSol(TPZAnalysis &an, TPZVec<REAL> &ErrorVecByIteration, int nref, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU);

bool ProcessingErrorUKnowingExactSol(TPZAnalysis &an, TPZVec<REAL> &ErrorVecByIteration, int nref, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU);


/* 2. Functions to apply strategies based on defined tables based on measure of solution error and/or gradient error
 */
 // HP adaptive for strategies in specific tables
bool ApplyingHPAdaptiveStrategyBasedOnU_I(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol,int MaxPOrder,int MaxHLevel,std::ostream &out=std::cout);
bool ApplyingHPAdaptiveStrategyBasedOnU_II(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);
bool ApplyingHPAdaptiveStrategyBasedOnU_III(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_III(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_IV(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_V(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);
bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_VI(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);

bool ApplyingHPAdaptiveStrategyBasedOnUAndDU_XI(TPZCompMesh *cmesh, TPZVec<STATE> &ErrorU, TPZVec<STATE> &ErrorDU, TPZVec<REAL> &Tol, int MaxPOrder, int MaxHLevel,std::ostream &out=std::cout);

void ApplyHPRefinement(TPZCompMesh *cmesh, TPZVec<int64_t> &PRef, int MaxPOrder, TPZVec<int64_t> &HRef, int MaxHLevel);

/* 3. Functions contructing computational meshes.
 */

void PrintNRefinementsByType(int64_t nels,int64_t newnels,int64_t hrefcounter,int64_t prefcounter,std::ostream &out = std::cout);



#endif
