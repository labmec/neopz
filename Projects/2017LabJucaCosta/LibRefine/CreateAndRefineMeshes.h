//
//  SimilarUniformRefinements.hpp
//  PZ
//
//  Created by labmec on 08/10/17.
//
//

#ifndef SimilarUniformRefinements_hpp
#define SimilarUniformRefinements_hpp

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

//#include "TPZCreateHDivMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;

/**  Global variables
std::ofstream out("OutPoissonArcTan.txt",std::ios::app);             // To store output of the console
// ABOUT H P ADAPTIVE
int MaxPOrder = 7;     // Maximum order for p refinement allowed
int MaxHLevel = 4;      // Maximum level for h refinement allowed
int MaxHUsed = 0;
int MaxPUsed = 0;
// Poisson problem
STATE ValueK = 1;
STATE F = sqrt(ValueK);
int ModelDimension;
// Circunference with high gradient - data
TPZManVector<REAL,3> CCircle(3,0.5);
REAL RCircle = 0.25;

*/

/* 1. Functions contructing geometrical meshes.
 Projects:
 Poisson3D_Shock
 */

TPZGeoMesh *CreateGeomMesh(int typeel,int mat,int bc0,int bc1=0,int bc2=0);
TPZGeoMesh *CreateGeomMesh(std::string &nome);

TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,int typeel,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel,int id_bc0,int id_bc1=0,int id_bc2=0);

int MaxLevelReached(TPZCompMesh *cmesh);

void PrintNRefinementsByType(int nref, int64_t nels,int64_t newnels,TPZVec<int64_t> &counter,std::ostream &out = std::cout);


/* 2. Functions to uniform refinement of the geometric meshes.
 Projects:
 Poisson3D_Shock
 */

/** Fucntions to apply refinement. */
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh);
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, TPZVec<int> *MatIdsVec=NULL);

void RegularizeMesh(TPZGeoMesh *gmesh,int dimension);


/* 3. Functions contructing computational meshes.
 Projects:
 Poisson3D_Shock
 */





#endif /* SimilarUniformRefinements_hpp */
