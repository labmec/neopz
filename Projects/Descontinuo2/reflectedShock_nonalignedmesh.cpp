
#include "reflectedShock_nonalignedmesh.h"

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
#include "pzintel.h"

using namespace std;

// Creates a mesh for the reflected shock problem
 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
void RSNAMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
{
   REAL x1 = 0.,
        x2 = 4.12791,
	y1 = 0.,
	y2 = 1.;

   pt.Resize(6);
   TPZVec<REAL> coord(3);

   coord[0] = x1;
   coord[1] = y1;
   coord[2] = 0.;
   pt[0] = coord;

   coord[0] = x2/2.;
   coord[1] = y1;
   coord[2] = 0.;
   pt[1] = coord;

   coord[0] = x2;
   coord[1] = y1;
   coord[2] = 0.;
   pt[2] = coord;

   coord[0] = x1;
   coord[1] = y2;
   coord[2] = 0.;
   pt[3] = coord;

   coord[0] = x2/2;
   coord[1] = y2;
   coord[2] = 0.;
   pt[4] = coord;

   coord[0] = x2;
   coord[1] = y2;
   coord[2] = 0.;
   pt[5] = coord;

// quadrilateral data

   TPZVec< int > nodes(4);

   elms.Resize(2);

   nodes[0] = 0;
   nodes[1] = 1;
   nodes[2] = 4;
   nodes[3] = 3;
   elms[0] = nodes;

   nodes[0] = 1;
   nodes[1] = 2;
   nodes[2] = 5;
   nodes[3] = 4;
   elms[1] = nodes;

}

TPZGeoMesh * CreateRSNAGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
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

   if(nSubdiv > 2)PZError << "CreateRSNAGeoMesh unsupported number of subdivisions";

   if(nSubdiv > 0)
   {
      int j;
      TPZVec< TPZGeoEl * > firstDivision;
      for(i = 0; i < gEls.NElements(); i++)
      {
         gEls[i]->Divide(firstDivision);
         if(nSubdiv > 1)
         {
            TPZVec< TPZGeoEl * > secondDivision;
            for(j = 0; j < firstDivision.NElements();j++)firstDivision[j]->Divide(secondDivision);
	 }
      }
   }


   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * RSNACompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
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
   RSNAMeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateRSNAGeoMesh(cmesh->Reference(), nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

   //TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);
   cmesh->SetDimModel(2);

// Creating the materials
   TPZEulerConsLaw * matp = new TPZEulerConsLaw(1/*nummat*/,
                                            0/*timeStep*/,
					    gamma /*gamma*/,
					    2 /* dim*/,
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
*/
//   cout << .22/(2/**lambdaMax*/);

   matp->SetDelta(delta);

   TPZAutoPointer<TPZMaterial> mat(matp);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZMaterial *bc;
   TPZFMatrix<REAL> val1(4,4), val2(4,1);
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
   bc = mat->CreateBC(mat,-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA DIREITA : OUTFLOW
   val1.Zero();
   val2.Zero();
   TPZGeoElBC((TPZGeoEl *)gElem[1],5,-2);
   bc = mat->CreateBC(mat,-2,4,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA SUPERIOR : DIRICHLET
   val1.Zero();
   val2.Zero();
   val2(0,0) = ro;
   val2(1,0) = ro * u;
   val2(2,0) = ro * v;
   val2(3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
   TPZGeoElBC((TPZGeoEl *)gElem[0],6,-3);
   TPZGeoElBC((TPZGeoEl *)gElem[1],6,-3);
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
   bc = mat->CreateBC(mat,-4,3,val1,val2);
   cmesh->InsertMaterialObject(bc);

   cmesh->AutoBuild();
   cmesh->AdjustBoundaryElements();

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
   for(int k = 0; k < nVars; k++)Solution(k)=.1;

   int nel = cmesh->NElements();
   int iel;
   for(iel=0; iel<nel; iel++)
   {
     TPZCompEl *cel = cmesh->ElementVec()[iel];
     TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
     TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
     if(intel)
     {
       int nnodes = intel->Reference()->NNodes();
       int in;
       for(in = 0; in<nnodes; in++)
       {
         int blnum = intel->SideConnect(0,in)->SequenceNumber();
         int blockOffset = cmesh->Block().Position(blnum);

         REAL ro = 1.7,
         u = 2.9,
         v = 0,
         p = 2.714286,
         vel2 = u*u + v*v;
         Solution(blockOffset  ,0) = ro;
         Solution(blockOffset+1,0) = ro * u;
         Solution(blockOffset+2,0) = ro * v;
         Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
       }
     } else if(disc)
     {
       int conind = disc->ConnectIndex();
       // skip boundary elements
       if(conind < 0) continue;
       int blnum = cmesh->ConnectVec()[conind].SequenceNumber();
       int blpos = cmesh->Block().Position(blnum);
       int blsize = cmesh->Block().Size(blnum);
       int blockOffset = blpos+blsize-(dim+2);
       REAL ro = 1.7,
       u = 2.9,
       v = 0,
       p = 2.714286,
       vel2 = u*u + v*v;
       Solution(blockOffset  ,0) = ro;
       Solution(blockOffset+1,0) = ro * u;
       Solution(blockOffset+2,0) = ro * v;
       Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
     }
   }

/*   int j, NSolutionBlocks;
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

   }*/
   cmesh->LoadSolution(Solution);
   return cmesh;
}
