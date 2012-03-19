
#include "ShockTube2d.h"

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
const int nSTEl = 2 * 1;

// Creates a mesh for the simple shock problem

void STMeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
{
   cout << "\nAlpha [Deg]?\n";
   REAL alpha;
   cin >> alpha;
   alpha*= 3.14159 / 180.;

   REAL cosA = cos(alpha);
   REAL sinA = sin(alpha);

   REAL y1=0,
	y2=.5;

   pt.Resize(2*nSTEl + 2);
   TPZVec<REAL> coord(3);

   int i;
   int nPts = nSTEl*2 + 2;

   for(i = 0; i < nPts/2; i++)
   {
      coord[0] = ((double)i)/((double)nSTEl) * cosA - y1 * sinA;
      coord[1] = y1 * cosA + ((double)i)/((double)nSTEl) * sinA;
      coord[2] = 0.;
      pt[i] = coord;

      coord[0] = ((double)i)/((double)nSTEl) * cosA - y2 * sinA;
      coord[1] = y2 * cosA + ((double)i)/((double)nSTEl) * sinA;
      coord[2] = 0.;
      pt[nPts/2 + i] = coord;
   }

// quadrilateral data

   TPZVec< int > nodes(4);

   elms.Resize(nSTEl);

   for(i = 0; i < nSTEl; i++)
   {
      nodes[0] = i;
      nodes[1] = i+1;
      nodes[2] = nPts/2 + i + 1;
      nodes[3] = nPts/2 + i;
      elms[i] = nodes;
   }
}

TPZGeoMesh * CreateSTGeoMesh(TPZGeoMesh *gmesh, TPZVec< TPZVec< REAL > > & nodes,
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
// Dividing elements to create a mesh of 4 elems.

      TPZVec< TPZGeoEl * > firstDivision;
      for(i = 0; i < gEls.NElements();i++)gEls[i]->Divide(firstDivision);
   }

   if(nSubdiv > 1)PZError << "CreateOneElGeoMesh unsupported number of subdivisions";

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh * STCompMesh(TPZFlowCompMesh *cmesh, REAL CFL, REAL delta,
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
   STMeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateSTGeoMesh(cmesh->Reference(), nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

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
   matp->SetDelta(delta);

   TPZAutoPointer<TPZMaterial> mat(matp);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZBndCond *bc;

   REAL rhol  = 1.,
        ul = 1.e-8,
	vl = 1.e-8,
	pl = 1.;
   REAL rhoel = pl/(gamma - 1.) + rhol * (ul * ul + vl * vl);

   REAL rhor  = 3.,
        ur = 1.e-8,
	vr = 1.e-8,
	pr = 3.;
   REAL rhoer = pr/(gamma - 1.) + rhor * (ur * ur + vr * vr);

   //CC Todas as arestas : PAREDE
   TPZFMatrix<REAL> val1, val2;
   val1.Zero();
   val2.Zero();
   int i;
   for(i = 0; i < nSTEl; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],4,-1);
      TPZGeoElBC((TPZGeoEl *)gElem[i],6,-1);
   }
   // aresta direita
   TPZGeoElBC((TPZGeoEl *)gElem[nSTEl-1],5,-1);
   // aresta esquerda
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-1);

   bc = mat->CreateBC(mat,-1,5,val1,val2);
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
   TPZFMatrix<REAL> Solution = cmesh->Solution();

   int nVars = Solution.Rows();
   for(int k = 0; k < nVars; k++)Solution(k)=-.1;

   Solution.Zero();
   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks / 2; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      Solution(blockOffset  ,0) = rhol;
      Solution(blockOffset+1,0) = ul * rhol;
      Solution(blockOffset+2,0) = vl * rhol;
      Solution(blockOffset+3,0) = rhoel;

      blockOffset = cmesh->Block().Position(j + NSolutionBlocks / 2) + lastShapeFun;

      Solution(blockOffset  ,0) = rhor;
      Solution(blockOffset+1,0) = ur * rhor;
      Solution(blockOffset+2,0) = vr * rhor;
      Solution(blockOffset+3,0) = rhoer;

   }
   cmesh->LoadSolution(Solution);

   return cmesh;
}
