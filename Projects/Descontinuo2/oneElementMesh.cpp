
#include "oneElementMesh.h"

#include "pzeuleranalysis.h"
#include "pzconslaw.h"
#include "pzmaterial.h"
#include "pzeulerconslaw.h"
#include "pzartdiff.h"
#include "pzreal.h"
#include "pzvec.h"
#include "pzflowcmesh.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include <iostream>
#include <fstream>
#include "TPZGeoElement.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzgeoquad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzbndcond.h"
#include "pzblock.h"
// creates an one-quadrilateral element mesh

using namespace std;

 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
void OneElMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
{
   REAL x1 = 0.,
        x2 = 5.,
	y1 = 0.,
	y2 = 3.1;

   pt.Resize(4);
   TPZVec<REAL> coord(3);

   coord[0] = x1;
   coord[1] = y1;
   coord[2] = 0.;
   pt[0] = coord;

   coord[0] = x2;
   coord[1] = y1;
   coord[2] = 0.;
   pt[1] = coord;

   coord[0] = x2;
   coord[1] = y2;
   coord[2] = 0.;
   pt[2] = coord;

   coord[0] = x1;
   coord[1] = y2;
   coord[2] = 0.;
   pt[3] = coord;

// quadrilateral data

   TPZVec< int > nodes(4);

   elms.Resize(1);

   nodes[0] = 0;
   nodes[1] = 1;
   nodes[2] = 2;
   nodes[3] = 3;
   elms[0] = nodes;

}

TPZGeoMesh * CreateOneElGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls,
			   int nSubdiv)
{
//   TPZGeoMesh * gmesh = new TPZGeoMesh;

   gEls.Resize(elms.NElements());
   gmesh->NodeVec().Resize(nodes.NElements());
   int i;
   for(i = 0; i < nodes.NElements(); i++)
   {
      gmesh->NodeVec()[i].Initialize(nodes[i],*gmesh);
   }

   for( i = 0; i < elms.NElements(); i++)
   {
      gEls[i] = gmesh->CreateGeoElement(ElType, elms[i], matId, i);
   }

// Constructing neighborhood

   gmesh->BuildConnectivity();

   if(nSubdiv == 1)
   {
// Dividing elements to create a mesh of 4 elems.

     TPZVec< TPZGeoEl * > firstDivision;
     gEls[0]->Divide(firstDivision);
   }

   if(nSubdiv > 1)PZError << "CreateOneElGeoMesh unsupported number of subdivisions";

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh *
   OneElCompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
                 int degree, int nSubdiv,
		 TPZArtDiffType DiffType,
		 TPZTimeDiscr Diff_TD,
		 TPZTimeDiscr ConvVol_TD,
		 TPZTimeDiscr ConvFace_TD)
{
   TPZCompEl::SetgOrder(degree);
   REAL gamma = 1.4;

// Configuring the PZ to generate discontinuous elements
//    TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>
//                 ::SetCreateFunction(TPZCompElDisc::CreateDisc);
//
//    TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>
//                 ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   const int dim = 2;
//   int interfdim = dim -1;
//   TPZCompElDisc::gInterfaceDimension = interfdim;


// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   OneElMeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateOneElGeoMesh(cmesh->Reference(), nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

   //TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);

// Creating the materials
   TPZEulerConsLaw * matp = new TPZEulerConsLaw(1/*nummat*/,
                                            0/*timeStep*/,
					    gamma /*gamma*/,
					    dim /* dim*/,
					    DiffType);
// Setting initial solution
   matp->SetForcingFunction(NULL);
   // Setting the time discretization method
   matp->SetTimeDiscr(Diff_TD,
                     ConvVol_TD,
		     ConvFace_TD);
   //mat->SetDelta(0.1); // Not necessary, since the artDiff
   // object computes the delta when it equals null.

   matp->SetCFL(CFL);
/*
   REAL us = sqrt(5.5 * 5.5 + 3.3 * 3.3);
   REAL press = 2.;
   REAL cspeed = sqrt(1.4*press/1.7);
   REAL lambdaMax = us + cspeed;
*/
   matp->SetDelta(delta);
   TPZAutoPointer<TPZMaterial> mat(matp);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZAutoPointer<TPZMaterial>  bc;
   TPZFMatrix<REAL> val1(4,4), val2(4,1);
   REAL ro = 1.7,
	//u = 5.5,
	//v = 0,//3.3,
	p = 2.;
	//vel2 = u*u + v*v;

// copiado do Cedric
   //CC Todas as arestas: DIRICHLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = ro;
   val2(1,0) = /*ro */ .5;
   val2(2,0) = /*ro */ 0.;
   val2(3,0) = p;///(gamma-1.0) + 0.5 * ro * vel2;
/*
   val2(0,0) = ro;
   val2(1,0) = ro * u;
   val2(2,0) = ro * v;
   val2(3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;*/
   TPZGeoElBC((TPZGeoEl *)gElem[0],4,-1);
   TPZGeoElBC((TPZGeoEl *)gElem[0],5,-1);
   TPZGeoElBC((TPZGeoEl *)gElem[0],6,-1);
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-1);
   bc = mat->CreateBC(mat,-1,9,val1,val2);
   cmesh->InsertMaterialObject(bc);

   cmesh->AutoBuild();
   //cmesh->AdjustBoundaryElements();

// printing meshes

   ofstream geoOut("geomesh.out");
   gmesh->Print(geoOut);
   geoOut.close();

   ofstream compOut("compmesh.out");
   cmesh->Print(compOut);
   compOut.close();

// generating initial guess for the mesh solution
   TPZFMatrix<REAL> Solution = cmesh->Solution();
   Solution.Zero();

   int nVars = Solution.Rows();
   int k;
   for(k = 0; k < nVars; k++)Solution(k,0) = 0;//.05;
   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      REAL ro = 1.7,
	   u = .5,
	   v = 0.01,//0.0821588,//3.5,
	   p = 2.,
	   c = sqrt(1.4 * p / ro);
      REAL vel2 = (u*u + v*v) * c*c;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u * c;
      Solution(blockOffset+2,0) = ro * v * c;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
 /*     ro = 1.4;
      u = 3.3;
      v = 5.5;
      p = 4.;
      vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = (p/(gamma-1.0) + 0.5 * ro * vel2);
*/
/*      int nVars = Solution.Rows();
      int k;
      for(k = 4; k < nVars; k++)Solution(k,0) = .05;*/
   }

   cmesh->LoadSolution(Solution);

   return cmesh;
}
