
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
#include "pzgeoquad.h"
#include "pzrefquad.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzblock.h"
// creates an one-quadrilateral element mesh

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

TPZGeoMesh * CreateOneElGeoMesh(TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls)
{
   TPZGeoMesh * gmesh = new TPZGeoMesh;

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

   TPZVec< TPZGeoEl * > firstDivision;
   gEls[0]->Divide(firstDivision);

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * OneElCompMesh()
{
   TPZCompElDisc::gDegree = 1;
   REAL gamma = 1.4;

// Configuring the PZ to generate discontinuous elements
   TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>
                ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>
                ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   int interfdim = 1;
   TPZCompElDisc::gInterfaceDimension = interfdim;


// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   OneElMeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateOneElGeoMesh(nodes, elms, EQuadrilateral, 1, gElem);

   TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);

// Creating the materials
   TPZEulerConsLaw2 * mat = new TPZEulerConsLaw2(1/*nummat*/,
                                            0/*timeStep*/,
					    gamma /*gamma*/,
					    2 /* dim*/,
					    LeastSquares_AD /*pzartdiff.h*/);
// Setting initial solution
   mat->SetForcingFunction(NULL);
   // Setting the time discretization method
   mat->SetTimeDiscr(Implicit_TD/*Diff*/,
                     Implicit_TD/*ConvVol*/,
		     Implicit_TD/*ConvFace*/);
   //mat->SetDelta(0.1); // Not necessary, since the artDiff
   // object computes the delta when it equals null.

   mat->SetCFL(1000.);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZBndCond * bc;
   TPZFMatrix val1(4,4), val2(4,1);
   REAL ro = 1.7,
	u = 5.5,
	v = 3.3,
	p = 2.,
	vel2 = u*u + v*v;

// copiado do Cedric
   //CC Todas as arestas: DIRICHLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = ro;
   val2(1,0) = ro * u;
   val2(2,0) = ro * v;
   val2(3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
   TPZGeoElBC((TPZGeoEl *)gElem[0],4,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[0],5,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[0],6,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-1,*gmesh);
   bc = mat->CreateBC(-1,3,val1,val2);
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
   TPZFMatrix Solution = cmesh->Solution();
   Solution.Zero();

   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j);

      REAL ro = 1.7,
	   u = 5.5,
	   v = 3.3,
	   p = 2.,
	   vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
      ro = 1.4;
      u = 3.3;
      v = 5.5;
      p = 4.;
      vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = (p/(gamma-1.0) + 0.5 * ro * vel2);
   }

   cmesh->LoadSolution(Solution);

   return cmesh;
}
