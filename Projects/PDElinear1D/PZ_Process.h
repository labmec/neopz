/*
 * @file
 * @brief Contains declaration of functions showing the process to generate geometric and computational meshes and the assemble and solve linear system associated.
 */

#ifndef PZPROCESSESDECLARATIONHPP
#define PZPROCESSESDECLARATIONHPP

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzmat1dlin.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"

#include <iostream>
#include <string>
#include <math.h>

/* GENERIC FUNCTIONS TO APPLY FEM USING PZ */

// Creates a Geometrical Mesh. It is related with one-dimensional domain - interval (into space)
TPZGeoMesh *GeomMesh(int h,TPZVec<int> &matId,TPZVec<int> &bc, TPZVec<REAL> &xL, TPZVec<REAL> &xR);
// Creates a Computational Mesh. It is related as one-dimensional linear diff. eq.
TPZCompMesh *CompMesh(TPZGeoMesh *gmesh, int p, TPZMaterial *material, TPZVec<int> &bc, TPZVec<int> &bcType);
// Assemble and Solve the generated linear System
void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh);

// Refines geometric mesh level by level from given mesh
void UniformRefinement(int h,TPZGeoMesh *gmesh);


/* AUXILIAR FUNCTIONS TO HELP US VISUALIZATION */

// Output as VTK (Visualization Tool Kit) format
void OutputVTK(std::string &outVTK, TPZCompMesh *cmesh,TPZAnalysis &an);
// Output as Mathematica format
void OutputMathematica(std::ofstream &outMath,int var,int pointsByElement,TPZCompMesh *cmesh);

#endif