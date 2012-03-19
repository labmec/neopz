#include "SimpleShock.h"
#include "pzeuleranalysis.h"
#include "pzconslaw.h"
#include "pzmaterial.h"
#include "pzbndcond.h"
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

using namespace std;

 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
const int nEl = 1;

// Creates a mesh for the simple shock problem

void SSMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
{
   REAL y1=0,
	y2=.5;

   pt.Resize(2*nEl + 2);
   TPZVec<REAL> coord(3);

   int i;
   int nPts = nEl*2 + 2;
   for(i = 0; i < nPts/2; i++)
   {
      coord[0] = ((double)i)/((double)nEl);
      coord[1] = y1;
      coord[2] = 0.;
      pt[i] = coord;

      coord[0] = ((double)i)/((double)nEl);
      coord[1] = y2;
      coord[2] = 0.;
      pt[nPts/2 + i] = coord;
   }

// quadrilateral data

   TPZVec< int > nodes(4);

   elms.Resize(nEl);

   for(i = 0; i < nEl; i++)
   {
      nodes[0] = i;
      nodes[1] = i+1;
      nodes[2] = nPts/2 + i + 1;
      nodes[3] = nPts/2 + i;
      elms[i] = nodes;
   }
}

TPZGeoMesh * CreateSSGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls,
			   int nSubdiv)
{
   //TPZGeoMesh * gmesh = new TPZGeoMesh;

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
      TPZVec< TPZGeoEl * > firstDivision;
      for(i = 0; i < gEls.NElements();i++)gEls[i]->Divide(firstDivision);
   }


   if(nSubdiv > 1)PZError << "CreateSSGeoMesh unsupported number of subdivisions";

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * SSCompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
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

   int dim = 2;
//   int interfdim = dim -1;
//   TPZCompElDisc::gInterfaceDimension = interfdim;


// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   SSMeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateSSGeoMesh(cmesh->Reference(),nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

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

/*
   REAL us = sqrt(2.6 * 2.6 + .51 * .51);
   REAL press = 1.52819;
   REAL cspeed = sqrt(1.4*press/1.7);
   REAL lambdaMax = us + cspeed;
*/
   matp->SetCFL(CFL);
   matp->SetDelta(delta);

   TPZAutoPointer<TPZMaterial> mat(matp);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZBndCond *bc;
   TPZFMatrix<REAL> val1(4,4), val2(4,1);
   REAL rhol  = 1.,
        rhoul = 1.,
	rhovl = 0.,
	rhoel = 0.9463,
	rhor  = 2.66709340161093,
	rhour = 1.,
	rhovr = 0.,
	rhoer = 2.19642;


   //CC ARESTA INFERIOR E SUPERIOR : PAREDE
   val1.Zero();
   val2.Zero();
   int i;
   for(i = 0; i < nEl; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],4,-1);
      TPZGeoElBC((TPZGeoEl *)gElem[i],6,-1);
   }
   bc = mat->CreateBC(mat,-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA DIREITA : DIRICHLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = rhol;
   val2(1,0) = rhoul;
   val2(2,0) = rhovl;
   val2(3,0) = rhoel;
//   TPZGeoElBC((TPZGeoEl *)gElem[0],5,-2);
   TPZGeoElBC((TPZGeoEl *)gElem[nEl-1],5,-2);
   bc = mat->CreateBC(mat,-2,3,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA ESQUERDA : DIRICHLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = rhol;
   val2(1,0) = rhoul;
   val2(2,0) = rhovl;
   val2(3,0) = rhoel;
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-3);
//   TPZGeoElBC((TPZGeoEl *)gElem[1],7,-3);
   bc = mat->CreateBC(mat,-3,3,val1,val2);
   cmesh->InsertMaterialObject(bc);

/*
   //CC ARESTA ESQUERDA : INFLOW
   val1.Zero();
   val2.Zero();
   ro = 1.0;
   u = 2.9;
   v = 0.0;
   p = 0.714286;
   val2(0,0) = ro;
   val2(1,0) = ro * u;
   val2(2,0) = ro * v;
   vel2 = u*u+v*v;
   val2(3,0) = p/(gamma-1.0) +  0.5 * ro * vel2;
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-4);
   TPZGeoElBC((TPZGeoEl *)gElem[4],7,-4);
   bc = mat->CreateBC(-4,3,val1,val2);
   cmesh->InsertMaterialObject(bc);
*/
   cmesh->AutoBuild();
//   cmesh->AdjustBoundaryElements();

// printing meshes

   ofstream geoOut("geomesh.out");
   gmesh->Print(geoOut);
   geoOut.close();

   ofstream compOut("compmesh.out");
   cmesh->Print(compOut);
   compOut.close();

// generating initial guess for the mesh solution
   TPZFMatrix<REAL> Solution = cmesh->Solution();

   int nVars = Solution.Rows();
   for(int k = 0; k < nVars; k++)Solution(k)=-.01;

   Solution.Zero();
   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      Solution(blockOffset  ,0) = (rhor+rhol);
      Solution(blockOffset+1,0) = (rhour+rhoul);
      Solution(blockOffset+2,0) = (rhovr+rhovl);
      Solution(blockOffset+3,0) = (rhoer+rhoel);
/*
      if(j < NSolutionBlocks/2)
      {
         Solution(blockOffset  ,0) = rhol;
         Solution(blockOffset+1,0) = rhoul;
         Solution(blockOffset+2,0) = rhovl;
         Solution(blockOffset+3,0) = rhoel;
      }else{
         Solution(blockOffset  ,0) = rhor;
         Solution(blockOffset+1,0) = rhour;
         Solution(blockOffset+2,0) = rhovr;
         Solution(blockOffset+3,0) = rhoer;
      }

*/
/*
      Solution(blockOffset  ,0) = 2.01136;
      Solution(blockOffset+1,0) = 1.11521;
      Solution(blockOffset+2,0) = 0;
      Solution(blockOffset+3,0) = 1.65903;

//      for(int k = 4; k < 12; k++)Solution(blockOffset+k)=-1;*/

   }
   cmesh->LoadSolution(Solution);

/*
   {
      //Solution.Zero();

      for(i = 0; i < nEl/2; i ++)
      {
         int blockOffset = cmesh->Block().Position(i) + lastShapeFun;
         Solution(blockOffset  ,0) = rhol;
         Solution(blockOffset+1,0) = rhoul;
         Solution(blockOffset+2,0) = rhovl;
         Solution(blockOffset+3,0) = rhoel;
      }

      for(i; i < nEl; i++)
      {
         int blockOffset = cmesh->Block().Position(i) + lastShapeFun;
         Solution(blockOffset  ,0) = rhor;
         Solution(blockOffset+1,0) = rhour;
         Solution(blockOffset+2,0) = rhovr;
         Solution(blockOffset+3,0) = rhoer;
      }

   cmesh->LoadSolution(Solution);
   }
*/
   return cmesh;
}
