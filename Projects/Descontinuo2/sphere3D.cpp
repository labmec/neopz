#include "ratio.h"
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
// creates a quarter of a tridimentional sphere

void SpherePoints(TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms, int nSubdiv)
{

   int index;
   int nn = nSubdiv + 1;
   int nLayerPts = nn;

   REAL ri = 1., ro = 5.;

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
   //                 /\
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
   REAL alpha = 1. / ((double) nn + 1);
   REAL beta = 1. - alpha;
   REAL alpha2 = 1. / ((double) nn);
   REAL beta2 = 1. - alpha2;
   REAL alphan, betan;
   int i, j, k;

   v1[2] = 1.;
   v2[1] = 1.;
   v3[0] = -1.;

   for(i = 0; i < 3; i++)
   {
      v4[i] = v3[i] * alpha  + v1[i] * beta;
      v5[i] = v3[i] * alpha2 + v2[i] * beta2;
      v6[i] = v2[i];
      v7[i] = v5[i];
   }

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

   //normalizing vectors, so that they fit in a sphere
   for(i = 0; i < nPtsPerLayer; i++)
   {
      REAL size = pt[i][0] * pt[i][0] +
                  pt[i][1] * pt[i][1] +
		  pt[i][2] * pt[i][2];
      size = sqrt(size);
      pt[i][0] /= size;
      pt[i][1] /= size;
      pt[i][2] /= size;
   }

   //creating all points
   for(i = nLayerPts-1; i >= 0; i--)
   {
      alpha = ri + xpg(4., i, nLayerPts-1) * (ri - ro);
      for(j = 0; j < nPtsPerLayer; j++)
      {
         pt[i * nLayerPts + j][0] = pt[i][0] * alpha;
	 pt[i * nLayerPts + j][1] = pt[i][1] * alpha;
	 pt[i * nLayerPts + j][2] = pt[i][2] * alpha;
      }
   }

// quadrilateral data

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
      index = i * nElsPerSphereLayer + nSubDiv * nElsPerLayer;
      for(j = 0; j < nSubdiv; j++)
      {
         // left elements
	 nodes[0] = i * nPtsPerLayer + nSubdiv * nPtsPerLayer + j;
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
	 nodes[0] = i * nPtsPerSphereLayer + nn * nPtsPerLayer - j - 1;
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
	 elms[index + nSubdiv * 2 - j] = nodes;
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
      gEls[i] = gmesh->CreateGeoElement(ElType, elms[i], matId, i, nSubdiv);
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
   TPZCompElDisc::gDegree = degree;
   REAL gamma = 1.4;

// Configuring the PZ to generate discontinuous elements
   TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>
                ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>
                ::SetCreateFunction(TPZCompElDisc::CreateDisc);

   const int dim = 2;
//   int interfdim = dim -1;
//   TPZCompElDisc::gInterfaceDimension = interfdim;


// Retrieving the point coordinates and element references
   TPZVec< TPZVec< REAL > > nodes;
   TPZVec< TPZVec< int  > > elms;
   TPZVec< TPZGeoEl *> gElem;
   SpherePoints(nodes, elms, nSubdiv);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateSphereGeoMesh(nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

   TPZFlowCompMesh * cmesh = new TPZFlowCompMesh(gmesh);

// Creating the materials
   TPZEulerConsLaw2 * mat = new TPZEulerConsLaw2(1/*nummat*/,
                                            0/*timeStep*/,
					    gamma /*gamma*/,
					    dim /* dim*/,
					    DiffType);
// Setting initial solution
   mat->SetForcingFunction(NULL);
   // Setting the time discretization method
   mat->SetTimeDiscr(Diff_TD,
                     ConvVol_TD,
		     ConvFace_TD);
   //mat->SetDelta(0.1); // Not necessary, since the artDiff
   // object computes the delta when it equals null.

   mat->SetCFL(CFL);
/*
   REAL us = sqrt(5.5 * 5.5 + 3.3 * 3.3);
   REAL press = 2.;
   REAL cspeed = sqrt(1.4*press/1.7);
   REAL lambdaMax = us + cspeed;
*/
   mat->SetDelta(delta);

   cmesh -> InsertMaterialObject(mat);

// Boundary conditions

   TPZBndCond * bc;
   TPZFMatrix val1(4,4), val2(4,1);
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
   TPZGeoElBC((TPZGeoEl *)gElem[0],4,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[0],5,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[0],6,-1,*gmesh);
   TPZGeoElBC((TPZGeoEl *)gElem[0],7,-1,*gmesh);
   bc = mat->CreateBC(-1,9,val1,val2);
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
