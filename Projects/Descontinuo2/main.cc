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

void error(char *) {}

void InitialGuess(TPZVec<REAL> &x,TPZVec<REAL> &result){
// makes nothing
   result.Resize(4 /*(2 dimensões + rho + energy)*/);
   result.Fill(0.);
}

void MeshPoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms)
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

TPZGeoMesh * CreateGeoMesh(TPZVec< TPZVec< REAL > > & nodes,
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

   TPZCompElDisc:: gInterfaceDimension = 1;

   return gmesh;
}

int main()
{
   TPZCompElDisc::gDegree = 2;
   REAL gamma = 1.4;

// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   MeshPoints(nodes, elms);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateGeoMesh(nodes, elms, EQuadrilateral, 1, gElem);


// Por que o Cedric faz isto?
   TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>
                ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>
                ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   int interfdim = 1;
   TPZCompElDisc::gInterfaceDimension = interfdim;



// Constructing neighborhood

   gmesh->BuildConnectivity();

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
   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZBndCond * bc;
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
   TPZGeoElBC((TPZGeoEl *)gElem[0],4,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[1],4,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[2],4,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[3],4,-1,*gmesh);
   bc = mat->CreateBC(-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA DIREITA : OUTFLOW
   val1.Zero();
   val2.Zero();
   TPZGeoElBC((TPZGeoEl *)gElem[3],5,-2,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[6],5,-2,*gmesh);
   bc = mat->CreateBC(-2,4,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //CC ARESTA SUPERIOR : DIRICHLET
   val1.Zero();
   val2.Zero();
   ro = 1.7;
   u = 2.61934;
   v = -0.50632;
   p = 1.52819;
   vel2 = u*u + v*v;
   val2(0,0) = ro;
   val2(1,0) = ro * u;
   val2(2,0) = ro * v;
   val2(3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
   TPZGeoElBC((TPZGeoEl *)gElem[7],6,-3,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[8],6,-3,*gmesh);
   bc = mat->CreateBC(-3,3,val1,val2);
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
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-4,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[4],7,-4,*gmesh);
   bc = mat->CreateBC(-4,3,val1,val2);
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

   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j);

      REAL ro = 1.7,
           u = 2.61934,
           v = -0.50632,
           p = 1.52819,
           vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = ro;
      Solution(blockOffset+1,0) = ro * u;
      Solution(blockOffset+2,0) = ro * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * ro * vel2;
   }

   cmesh->LoadSolution(Solution);

// Creating the analysis object

   ofstream anFile("analysis.out");
   TPZEulerAnalysis An(cmesh, anFile);

   // Creating the structural matrix
   TPZBandStructMatrix StrMatrix(cmesh);
   An.SetStructuralMatrix(StrMatrix);

   // Creating the solver for the linearized systems
   TPZStepSolver Solver;
   Solver.SetDirect(ELU);// ECholesky -> simétrica e positiva definida
   An.SetSolver(Solver);

   An.Run(anFile);

   return 0;
}
