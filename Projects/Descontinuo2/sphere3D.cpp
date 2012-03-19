#include "sphere3D.h"
#include "ratio.h"
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
#include "pzshapecube.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzrefquad.h"
#include "TPZRefCube.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "pzblock.h"

using namespace std;

// creates a quarter of a tridimentional sphere
 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
void SpherePoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms, int nSubdiv)
{

   int index;
   int nn = nSubdiv + 1;
   int nLayerPts = nSubdiv * 2;

   REAL ri = 1., ro = 25.;

   // the shape of the quarter of the sphere is made
   // based on the shape below:
   //          _ ____
   //        /  5    2 |
   //      *     \   | |
   //     /      \ A | |
   //   *         \  | |
   //  |          \  | |
   // 3      C    4--1 |
   //  \          |  | |
   //   *        |   | |
   //     \      | B | |
   //      *    |    | |
   //        \_ 7____6 |
   //                  |
   //                [\|/]
   //                 /
   // the points 1, 2 and 3 represents the projections of vectors
   // [0,0,-1], [0,1,0], [0,0,1]
   // Point 4 is the projection of alpha*[1]+beta*[3] onto the sphere
   // Point 5 is the projection of alpha'*[2]+beta'*[3]
   // The final renumbering scheme will differ totally from these
   // construction parameters.
   // Points 6 and 7 are obtained by summetry.
   // The regions A and B are meshed in one line of elements.
   // Regioc C is squared out.

   TPZVec<REAL> v1(3,0.),
                v2(3,0.),
		v3(3,0.),
		v4(3,0.),
		v5(3,0.),
		v6(3,0.),
		v7(3,0.);

   TPZVec<REAL> coord(3), ref1(3), ref2(3);
   REAL alpha = 1. / ((double) nn) * .7;
   REAL beta = 1. - alpha;
   REAL alpha2 = 1. / ((double) nn) * 1.2;
   REAL beta2 = 1. - alpha2;
   REAL alphan, betan;
   int i, j, k;

   v1[2] = 1.;

   v2[1] = 1./sqrt(2.);
   v2[2] = 1./sqrt(2.);

   v3[0] = -1.;

   REAL size4 = 0., size5 = 0.;
   // computing v4 and v5 normalized
   for(i = 0; i < 3; i++)
   {
      v4[i] = v3[i] * alpha  + v1[i] * beta;
      size4 += v4[i] * v4[i];
      v5[i] = v3[i] * alpha2 + v2[i] * beta2;
      size5 += v5[i] * v5[i];
   }
   size5 = sqrt(size5);
   size4 = sqrt(size4);
   for(i = 0; i < 3; i++)
   {
      v4[i] /= size4;
      v5[i] /= size5;
   }

   v6 = v2;
   v7 = v5;


   v7[1] *= -1.;
   v6[1] *= -1.;


   // notable measures
   int nPtsPerLayer = 2 * nn + 1;
   int nPtsPerSphereLayer = nPtsPerLayer * nn + nn - 1;
   int nPts = nPtsPerSphereLayer * nLayerPts;

   // index of reference points
   int refIndex[8];
   refIndex[3] = 0;
   refIndex[5] = nn - 1;
   refIndex[2] = nn;
   refIndex[7] = nPtsPerLayer * (nn-1);
   refIndex[4] = refIndex[7] + nn - 1;
   refIndex[1] = refIndex[4] + 1;
   refIndex[6] = nPtsPerLayer * nn;

   pt.Resize(nPts);


   // initializing ref points
   pt[refIndex[1]] = v1;
   pt[refIndex[2]] = v2;
   pt[refIndex[3]] = v3;
   pt[refIndex[4]] = v4;
   pt[refIndex[5]] = v5;
   pt[refIndex[6]] = v6;
   pt[refIndex[7]] = v7;

   // creating A and B-region points
   for(i = 0; i < nn; i++)
   {
      alphan = xpg(1., i, nn-1);
      betan = 1.-alphan;
      // A -left region Points
      index = i * nPtsPerLayer + nn;
      coord[0] = pt[refIndex[1]][0] * alphan + pt[refIndex[2]][0] * betan;
      coord[1] = pt[refIndex[1]][1] * alphan + pt[refIndex[2]][1] * betan;
      coord[2] = pt[refIndex[1]][2] * alphan + pt[refIndex[2]][2] * betan;
      pt[index] = coord;
   }
   for(i = 0; i < nn-1; i++)
   {
      alphan = xpg(1., i, nn-1);
      betan = 1.-alphan;
      // B -left region Points
      index = i + nPtsPerLayer * nn;
      coord[0] = pt[refIndex[1]][0] * alphan + pt[refIndex[6]][0] * betan;
      coord[1] = pt[refIndex[1]][1] * alphan + pt[refIndex[6]][1] * betan;
      coord[2] = pt[refIndex[1]][2] * alphan + pt[refIndex[6]][2] * betan;
      pt[index] = coord;
   }

   // creating C-region points for both sides (with boundaries)
   for(i = 0; i < nn; i++) // loop in the direction 3-7; 5-4
   {
      int baseindex1 = i * nPtsPerLayer;
      int baseindex2 = baseindex1 + nn - 1;
      alphan = xpg(1., i, nn-1);
      betan = 1. - alphan;

      // computing points at border 3-7
      coord[0] = pt[refIndex[7]][0] * alphan + pt[refIndex[3]][0] * betan;
      coord[1] = pt[refIndex[7]][1] * alphan + pt[refIndex[3]][1] * betan;
      coord[2] = pt[refIndex[7]][2] * alphan + pt[refIndex[3]][2] * betan;
      pt[baseindex1] = coord;
      // right elements
      coord[0] *= -1;
      pt[baseindex1 + nPtsPerLayer - 1] = coord;

      // computing points at border 5-4
      coord[0] = pt[refIndex[4]][0] * alphan + pt[refIndex[5]][0] * betan;
      coord[1] = pt[refIndex[4]][1] * alphan + pt[refIndex[5]][1] * betan;
      coord[2] = pt[refIndex[4]][2] * alphan + pt[refIndex[5]][2] * betan;
      pt[baseindex2] = coord;
      // right elements
      coord[0] *= -1;
      pt[baseindex1 + nn + 1] = coord;

      for(j = 1; j < nn - 1; j++) // loop in the direction 3-5; 7-4
      {
         // left point
         index = baseindex1 + j;
	 alphan = xpg(1., j, nn-1);
	 betan = 1. - alphan;

	 coord[0] = pt[baseindex2][0] * alphan + pt[baseindex1][0] * betan;
	 coord[1] = pt[baseindex2][1] * alphan + pt[baseindex1][1] * betan;
	 coord[2] = pt[baseindex2][2] * alphan + pt[baseindex1][2] * betan;
         pt[baseindex1 + j] = coord;

	 // right point
	 coord[0] *= -1;
	 pt[baseindex1 + nPtsPerLayer - 1 - j] = coord;
      }
   }

   // moving point with ref number = 4 - less deformed elements
   pt[refIndex[4]][0] *= .7;
   pt[refIndex[4]+2][0] *= .7;

   //normalizing vectors, so that they fit in a sphere
   for(i = 0; i < nPtsPerSphereLayer; i++)
   {
      REAL size = pt[i][0] * pt[i][0] +
                  pt[i][1] * pt[i][1] +
		  pt[i][2] * pt[i][2];
      size = sqrt(size);
      pt[i][0] /= size;
      pt[i][1] /= size;
      pt[i][2] /= size;
   }
#ifdef CONIC_EXTERNAL
// conic external points
   //creating internal / external points
   for(j = 0; j < nPtsPerSphereLayer; j++)
   {
      coord[0] = pt[j][0];
      coord[1] = pt[j][1];
      coord[2] = pt[j][2];

      // adjusting points to cones:
      REAL a,b,c, delta, sol1, sol2;
      a = - coord[0] * coord[0] +
            coord[1] * coord[1] +
	    coord[2] * coord[2];
      b = fabs(coord[0]) * 2.;
      c = -1;
      delta = b * b - 4. * a * c;
      delta = sqrt(delta);

      if(fabs(a) < 1e-4)
      {
         sol1 = 1. / b;
      }else{
         sol1 =( -b + delta) / 2. / a;
	 sol2 =( -b - delta) / 2. / a;

	 if(sol1 * sol2 < 0.)
	 {
            sol1 = max(sol1, sol2);
	 }else{
            sol1 = min(sol1, sol2);
	 }
      }

      coord[0] *= ro * sol1;
      coord[1] *= ro * sol1;
      coord[2] *= ro * sol1;
      pt[(nLayerPts - 1) * nPtsPerSphereLayer + j] = coord;

      coord[0] = pt[j][0] * ri;
      coord[1] = pt[j][1] * ri;
      coord[2] = pt[j][2] * ri;
      pt[j] = coord;

      // other points
      for(i = 1; i < nLayerPts-1; i++)
      {
         alpha = xpg(pow((ro-ri),1./((double) nLayerPts-1)), i, nLayerPts-1);
	 beta = 1. - alpha;
         coord[0] = pt[(nLayerPts - 1) * nPtsPerSphereLayer + j][0] * alpha +
	            pt[j][0] * beta;
         coord[1] = pt[(nLayerPts - 1) * nPtsPerSphereLayer + j][1] * alpha +
	            pt[j][1] * beta;
         coord[2] = pt[(nLayerPts - 1) * nPtsPerSphereLayer + j][2] * alpha +
	            pt[j][2] * beta;
         pt[i * nPtsPerSphereLayer + j] = coord;
      }
   }

#else
   //creating all points
   for(i = nLayerPts-1; i >= 0; i--)
   {
      alpha = ri + xpg(pow(ro,(REAL)(1./(nLayerPts-1))), i, nLayerPts-1) * (ro - ri);
      for(j = 0; j < nPtsPerSphereLayer; j++)
      {
         coord[0] = pt[j][0] * alpha;
	 coord[1] = pt[j][1] * alpha;
	 coord[2] = pt[j][2] * alpha;

	 pt[i * nPtsPerSphereLayer + j] = coord;
      }
   }
#endif
// cube data

   int nElsPerLayer = nn * 2;
   int nElsPerSphereLayer = nSubdiv * nElsPerLayer + nSubdiv * 2;
   int nLayerEls = nLayerPts - 1;
   int nEls = nLayerEls * nElsPerSphereLayer;

   TPZVec< int > nodes(8);
   elms.Resize(nEls);

   for(i = 0; i < nLayerEls ; i++) // layer loop
   {
      for(j = 0; j < nSubdiv; j++)
      {
         index = i * nElsPerSphereLayer + j * nn * 2;
         for(k = 0; k < nElsPerLayer; k++)
	 {
	    nodes[0] = i * nPtsPerSphereLayer + j * nPtsPerLayer + k;
	    nodes[1] = nodes[0] + 1;
	    nodes[2] = nodes[0] + nPtsPerLayer + 1;
	    nodes[3] = nodes[2] - 1;
	    nodes[4] = nodes[0] + nPtsPerSphereLayer;
	    nodes[5] = nodes[1] + nPtsPerSphereLayer;
	    nodes[6] = nodes[2] + nPtsPerSphereLayer;
	    nodes[7] = nodes[3] + nPtsPerSphereLayer;
	    elms[index + k] = nodes;
	 }
      }
      index = i * nElsPerSphereLayer + nSubdiv * nElsPerLayer;
      for(j = 0; j < nSubdiv; j++)
      {
         // left elements
	 nodes[0] = i * nPtsPerSphereLayer + nSubdiv * nPtsPerLayer + j;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = (i+1) * nPtsPerSphereLayer - nSubdiv + j + 1;
	 nodes[3] = nodes[2] - 1;
	 if(j == nSubdiv - 1)
	 {
            nodes[2] = nodes[1] + 1;
	 }
	 nodes[4] = nodes[0] + nPtsPerSphereLayer;
	 nodes[5] = nodes[1] + nPtsPerSphereLayer;
	 nodes[6] = nodes[2] + nPtsPerSphereLayer;
	 nodes[7] = nodes[3] + nPtsPerSphereLayer;
	 elms[index + j] = nodes;

	 // right elements
	 nodes[0] = i * nPtsPerSphereLayer + nn * nPtsPerLayer - j - 2;
	 nodes[1] = nodes[0] + 1;
	 int tmp  = nodes[2];
	 nodes[2] = nodes[3];
	 nodes[3] = tmp;
	 if(j == nSubdiv - 1)
	 {
            nodes[3] = nodes[0] - 1;
	 }
	 nodes[4] = nodes[0] + nPtsPerSphereLayer;
	 nodes[5] = nodes[1] + nPtsPerSphereLayer;
	 nodes[6] = nodes[2] + nPtsPerSphereLayer;
	 nodes[7] = nodes[3] + nPtsPerSphereLayer;
	 elms[index + nSubdiv * 2 - j - 1] = nodes;
      }

   }
}

TPZGeoMesh * CreateSphereGeoMesh(TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls,
			   int nSubdiv)
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

   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh *
   SphereCompMesh(REAL CFL, REAL delta,
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
//
//    TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>
//                 ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   const int dim = 3;
//   int interfdim = dim -1;
//   TPZCompElDisc::gInterfaceDimension = interfdim;


// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   SpherePoints(nodes, elms, nSubdiv);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateSphereGeoMesh(nodes, elms, /*EQuadrilateral*/ECube, 1, gElem, nSubdiv);

   TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);

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
   TPZFMatrix<REAL> val1(5,5), val2(5,1);

   int nn = nSubdiv + 1;
   int nElsPerLayer = nn * 2;
   int nElsPerSphereLayer = nSubdiv * nElsPerLayer + nSubdiv * 2;
   int nLayerEls = nSubdiv * 2 - 1;
   int i, j;

   //Sphere faces: Wall
   val1.Zero();
   val2.Zero();
   for( i = 0; i < nElsPerSphereLayer; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],20,-1);
   }

   bc = mat->CreateBC(mat,-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   // Symmetry faces: slipwall
   val1.Zero();
   val2.Zero();
   // upper faces
   for( i = 0; i < nLayerEls; i++)
      for( j = 0; j < nElsPerLayer; j++)
      {
         TPZGeoElBC((TPZGeoEl *)gElem[i*nElsPerSphereLayer + j],21,-2);
      }
   // lower faces
   for( i = 0; i < nLayerEls; i++)
      for( j = 0; j < nn; j++)
      {
         // left elms.
         TPZGeoElBC((TPZGeoEl *)gElem[i*nElsPerSphereLayer + j * nElsPerLayer],24,-2);
	 // right elms.
	 if(j < nn-1)
	 {
            TPZGeoElBC((TPZGeoEl *)gElem[i*nElsPerSphereLayer + (j+1) * nElsPerLayer - 1],22,-2);
	 }else{
            TPZGeoElBC((TPZGeoEl *)gElem[(i+1)*nElsPerSphereLayer - 1],22,-2);
	 }
      }
   // slipwall
   bc = mat->CreateBC(mat,-2,12,val1,val2);
   cmesh->InsertMaterialObject(bc);

   //external sphere faces: inflow/outflow
   REAL Mach;
   cout << "\nMach number\n";
   cin >> Mach;

   val2(0,0) = 1.;// rho
   val2(1,0) = Mach;// Machx
   val2(2,0) = 0.;// Machy
   val2(3,0) = 0.;// Machz
   val2(4,0) = 1.;// pressure

   TPZVec<REAL> param(3,0.), X(3,0.);
   for( i = 0; i < nElsPerSphereLayer; i++)
   {

      TPZGeoEl * pEl = gElem[i +  (nLayerEls - 1) * nElsPerSphereLayer];
      pEl -> X(param, X);

      if(X[0] > 0.)
      {
      // outflow
         TPZGeoElBC(pEl,25,-3);
      }else{
      // inflow
         TPZGeoElBC(pEl,25,-4);
      }
   }
   // outflow
   bc = mat->CreateBC(mat,-3,9,val1,val2);
   cmesh->InsertMaterialObject(bc);
   // inflow
   bc = mat->CreateBC(mat,-4,9,val1,val2);
   cmesh->InsertMaterialObject(bc);

   cmesh->AutoBuild();

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
   for(i = 0; i < nVars; i++)Solution(i,0) = 0;
   int NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      REAL rho = 1.0,
           p = 1.,
	   u = sqrt(1.4 * p / rho) * Mach,
	   v = 0.,
	   w = 0.,
	   vel2 = u*u + v*v + w*w;
      Solution(blockOffset  ,0) = rho;
      Solution(blockOffset+1,0) = rho * u;
      Solution(blockOffset+2,0) = rho * v;
      Solution(blockOffset+3,0) = rho * w;
      Solution(blockOffset+4,0) = p/(gamma-1.0) + 0.5 * rho * vel2;
   }

   cmesh->LoadSolution(Solution);

   return cmesh;
}
