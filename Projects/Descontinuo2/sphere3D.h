#ifndef SPHERE3DHPP
#define SPHERE3DHPP

#include "pzreal.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pzgmesh.h"
#include "pzflowcmesh.h"
#include "pzeulerconslaw.h"


void SpherePoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms, int nSubdiv);
TPZGeoMesh * CreateSphereGeoMesh(TPZVec< TPZVec< REAL > > & nodes,
                                 TPZVec< TPZVec< int > > & elms,
                                 MElementType ElType, int matId,
                                 TPZVec<TPZGeoEl *> & gEls,
                                 int nSubdiv);
// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh *
SphereCompMesh(REAL CFL, REAL delta,
               int degree, int nSubdiv,
               TPZArtDiffType DiffType,
               TPZTimeDiscr Diff_TD,
               TPZTimeDiscr ConvVol_TD,
               TPZTimeDiscr ConvFace_TD);


#endif