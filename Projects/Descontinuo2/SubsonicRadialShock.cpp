#include "SubsonicRadialShock.h"
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

// creates an one-quadrilateral element mesh
 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
const int nn = 2;
const int mm = 4;


void SRSPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
{
   REAL ri = 1.8,
        ro = 5.,
	deltar = ro - ri,
	alpha = M_PI / (2. * nn)/8.,
	alphai,
	r;

   pt.Resize((mm+1)*(nn+1));
   TPZVec<REAL> coord(3);

   // defining points
   int i, j, index;
   for(j = 0; j < mm+1; j++)
      for(i = 0; i < nn+1; i++)
         {
            index = i + j * (nn + 1);
            r = ri + deltar/mm*j;
	    alphai = alpha * i;
	    coord[0] = r * cos(alphai);
	    coord[1] = r * sin(alphai);
	    coord[2] = 0.;
	    pt[index] = coord;
	 }

   // definig elements
   // quadrilateral data
   TPZVec< int > nodes(4);
   elms.Resize(mm * nn);

   for(j = 0; j < mm; j++)
      for(i = 0; i < nn; i++)
         {
            index = i + j * nn;
	    nodes[0] = i + j * (nn+1);
	    nodes[1] = i + (j + 1) * (nn+1);
	    nodes[2] = nodes[1] + 1;
	    nodes[3] = nodes[0] + 1;
	    elms[index] = nodes;
	 }

}

TPZGeoMesh * CreateSRSGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls,
			   int nSubdiv)
{
   //TPZGeoMesh * gmesh = new TPZGeoMesh;

   gEls.Resize(elms.NElements());
   gmesh->NodeVec().Resize(nodes.NElements());
   int i;   for(i = 0; i < nodes.NElements(); i++)
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
   SRSCompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
                 int degree, int nSubdiv,
		 TPZArtDiffType DiffType,
		 TPZTimeDiscr Diff_TD,
		 TPZTimeDiscr ConvVol_TD,
		 TPZTimeDiscr ConvFace_TD)
{
   TPZCompEl::SetgOrder(degree);
   REAL gamma = 1.4;
   int i;

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
   SRSPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateSRSGeoMesh(cmesh->Reference(), nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

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

   TPZBndCond *bc;
   TPZFMatrix val1(4,4), val2(4,1);

   //aresta superior: Parede
   val1.Zero();
   val2.Zero();
   for( i = 0; i < mm; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[nn * (i + 1) - 1],6,-1);
   }
   //aresta inferior: Parede
   for( i = 0; i < mm; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[nn * i],4,-1);
   }

   bc = mat->CreateBC(mat,-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //raio interno: INLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = 1. /*rho*/;
   val2(1,0) = .7919128627912566 /*Mach*/;
   val2(3,0) = 5.; //rhoE
   for( i = 0; i < nn; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],7,-2);
   }

   bc = mat->CreateBC(mat,-2,7,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //raio externo: OUTLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = 1.3158190229242444;//0;//.139006621;
   val2(1,0) = .2050900094073414 /*Mach*/;
   val2(3,0) = 6.319340316405589;//.1;//0;//5.;
   for( i = 0; i < nn; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i + (mm - 1) * nn],5,-3);
   }

   bc = mat->CreateBC(mat,-3,8,val1,val2);
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

      REAL ro = 1.0,
	   u = 1e-8,
	   v = 1e-8,
	   p = 5.,
	   vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
   }

   cmesh->LoadSolution(Solution);

   return cmesh;
}
