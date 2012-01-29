
#include "reflectedShock.h"

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

// Creates a mesh for the reflected shock problem
 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
void RSMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
{
   REAL p1=0.7,
        p2=0.902026,
	p3=1.80405,
	p4=2.96598,
	p5=3.2,
	p6=4.12791,
	p7=0.22;

   pt.Resize(15);
   TPZVec<REAL> coord(3);

   coord[0] = 0.;
   coord[1] = 0.;
   coord[2] = 0.;
   pt[0] = coord;

   coord[0] = p1+p1/4;
   coord[1] = 0.;
   coord[2] = 0.;
   pt[1] = coord;

   coord[0] = p3;
   coord[1] = 0.;
   coord[2] = 0.;
   pt[2] = coord;

   coord[0] = p5-p7/1.5;
   coord[1] = 0.;
   coord[2] = 0.;
   pt[3] = coord;

   coord[0] = p6;
   coord[1] = 0.;
   coord[2] = 0.;
   pt[4] = coord;

   coord[0] = 0.;
   coord[1] = p7;
   coord[2] = 0.;
   pt[5] = coord;

   coord[0] = p1;
   coord[1] = p7;
   coord[2] = 0.;
   pt[6] = coord;

   coord[0] = p2;
   coord[1] = .5;
   coord[2] = 0.;
   pt[7] = coord;

   coord[0] = p3;
   coord[1] = .75;
   coord[2] = 0.;
   pt[8] = coord;

   coord[0] = p4;
   coord[1] = .5;
   coord[2] = 0.;
   pt[9] = coord;

   coord[0] = p5;
   coord[1] = p7;
   coord[2] = 0.;
   pt[10] = coord;

   coord[0] = p6;
   coord[1] = p7;
   coord[2] = 0.;
   pt[11] = coord;

   coord[0] = 0.;
   coord[1] = 1.;
   coord[2] = 0.;
   pt[12] = coord;

   coord[0] = p3;
   coord[1] = 1.;
   coord[2] = 0.;
   pt[13] = coord;

   coord[0] = p6;
   coord[1] = 1.;
   coord[2] = 0.;
   pt[14] = coord;

// quadrilateral data

   TPZVec< int > nodes(4);

   elms.Resize(9);

   nodes[0] = 0;
   nodes[1] = 1;
   nodes[2] = 6;
   nodes[3] = 5;
   elms[0] = nodes;

   nodes[0] = 1;
   nodes[1] = 2;
   nodes[2] = 7;
   nodes[3] = 6;
   elms[1] = nodes;

   nodes[0] = 2;
   nodes[1] = 3;
   nodes[2] = 10;
   nodes[3] = 9;
   elms[2] = nodes;

   nodes[0] = 3;
   nodes[1] = 4;
   nodes[2] = 11;
   nodes[3] = 10;
   elms[3] = nodes;

   nodes[0] = 5;
   nodes[1] = 6;
   nodes[2] = 7;
   nodes[3] = 12;
   elms[4] = nodes;

   nodes[0] = 7;
   nodes[1] = 2;
   nodes[2] = 9;
   nodes[3] = 8;
   elms[5] = nodes;

   nodes[0] = 10;
   nodes[1] = 11;
   nodes[2] = 14;
   nodes[3] = 9;
   elms[6] = nodes;

   nodes[0] = 7;
   nodes[1] = 8;
   nodes[2] = 13;
   nodes[3] = 12;
   elms[7] = nodes;

   nodes[0] = 8;
   nodes[1] = 9;
   nodes[2] = 14;
   nodes[3] = 13;
   elms[8] = nodes;
}

TPZGeoMesh * CreateRSGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
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

   if(nSubdiv > 1)PZError << "CreateRSGeoMesh unsupported number of subdivisions";

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * RSCompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
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
   RSMeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateRSGeoMesh(cmesh->Reference(), nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

//   TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);

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
   REAL us = sqrt(2.6 * 2.6 + .51 * .51);
   REAL press = 1.52819;
   REAL cspeed = sqrt(1.4*press/1.7);
   REAL lambdaMax = us + cspeed;

   cout << .22/(2);
*/
   matp->SetDelta(delta);

   TPZAutoPointer<TPZMaterial> mat(matp);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZAutoPointer<TPZMaterial>  bc;
   TPZFMatrix val1(4,4), val2(4,1);
   REAL ro = 1.7,
	u = 2.61934,
	v = -0.50632,
	p = 1.52819,
	vel2 = u*u + v*v;

// copiado do Cedric
   //CC ARESTA INFERIOR : PAREDE
   val1.Zero();
   val2.Zero();
   TPZGeoElBC((TPZGeoEl *)gElem[0],4,-1);
   TPZGeoElBC((TPZGeoEl *)gElem[1],4,-1);
   TPZGeoElBC((TPZGeoEl *)gElem[2],4,-1);
   TPZGeoElBC((TPZGeoEl *)gElem[3],4,-1);
   bc = mat->CreateBC(mat,-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA DIREITA : OUTFLOW
   val1.Zero();
   val2.Zero();
   TPZGeoElBC((TPZGeoEl *)gElem[3],5,-2);
   TPZGeoElBC((TPZGeoEl *)gElem[6],5,-2);
   bc = mat->CreateBC(mat,-2,4,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA SUPERIOR : DIRICHLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = ro;
   val2(1,0) = ro * u;
   val2(2,0) = ro * v;
   val2(3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
   TPZGeoElBC((TPZGeoEl *)gElem[7],6,-3);
   TPZGeoElBC((TPZGeoEl *)gElem[8],6,-3);
   bc = mat->CreateBC(mat,-3,3,val1,val2);
   cmesh->InsertMaterialObject(bc);


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
   bc = mat->CreateBC(mat,-4,3,val1,val2);
   cmesh->InsertMaterialObject(bc);

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
   TPZFMatrix Solution = cmesh->Solution();

   int nVars = Solution.Rows();
   for(int k = 0; k < nVars; k++)Solution(k)=.1;


   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      REAL ro = 1.7,
           u = 2.9,
           v = 0,
           p = 2.714286,
           vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
/*
      for(int k = 4; k < 12; k++)Solution(blockOffset+k)=-1;*/

   }
   cmesh->LoadSolution(Solution);
/*
   // Exact Solution (Lucia Catabriga, Alvaro Coutinho)
   {
    //  Solution.Zero();
      REAL ro = 1.0,
           u = 2.9,
           v = 0,
           p = 0.714286,
           vel2 = u*u + v*v;

      int blockOffset = cmesh->Block().Position(0) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      blockOffset = cmesh->Block().Position(1) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      blockOffset = cmesh->Block().Position(4) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      ro = 1.7;
      u = 2.61934;
      v = -0.50632;
      p = 1.52819;
      vel2 = u*u + v*v;

      blockOffset = cmesh->Block().Position(5) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      blockOffset = cmesh->Block().Position(7) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      blockOffset = cmesh->Block().Position(8) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      ro = 2.68728;
      u = 2.40140;
      v = 0;
      p = 2.93407;
      vel2 = u*u + v*v;
            blockOffset = cmesh->Block().Position(2) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      blockOffset = cmesh->Block().Position(3) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;

      blockOffset = cmesh->Block().Position(6) + lastShapeFun;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;


   cmesh->LoadSolution(Solution);
   }
//*/
   return cmesh;
}
