
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

// This file generates a mesh for the NACA airfoils
double scale = 4.;
double entrance = 4. * scale,
             exitlength = 4. * scale,
	     cord = 5.,
	     height = 10. * scale,
	     q = 1.5,
	     qn = 1.5;
/*const int m = 6,
          n = 2,
	  l = 2;*/

double PI = 3.14159265359;

int l, m, n, k;

const int digits = 0012;

double Spg(double a0, double q, int N)
{
   return a0*(pow(q, N+1.)-1.)/(q-1.);
}

double xpg(double q, int n, int N)
{
   if(q==1.)return ((double)n) / ((double)N);
   return Spg(1., q, n-1)/Spg(1., q, N-1);
}

double P(int digits)
{
   int aux = digits/100;
   aux -= ((int)(aux/10))*10;
   return (double)aux/10.;
}

double M(int digits, double cord)
{
   int aux = digits/1000;
   return (double)aux/100.*cord;
}

double TT(int digits, double cord)
{
   int aux = digits - ((int)(digits/100))*100;
   return (double)aux/100.*cord;
}

// Mean line for the wing
double yc(double x, double c, double P, double M)
{
   if(x/c<P)
   {
      return M/P/P*(2.*P*x/c-x*x/c/c);
   }else
   {
      return M/(1.-P)/(1.-P)*(1.-2.*P+2.*P*x/c-x*x/c/c);
   }
}

double dyc(double x, double c, double P, double M)
{
   if(x/c<P)
   {
      return 2.*M/P/P*(P-x/c)/c;
   }else
   {
      return 2.*M/(1.-P)/(1.-P)*(P-x/c)/c;
   }
}

// thickness
double yt(double x, double TT, double c)
{
   double aux = x/c;
   const double a0 = 1.4845,
                a1 = -.6300,
		a2 = -1.7580,
		a3 = 1.4215,
		a4 = -.5075;

   return TT * (a0*sqrt(aux)+
                 a1*aux+
		 a2*aux*aux+
		 a3*aux*aux*aux+
		 a4*aux*aux*aux*aux);
}

// superior profile
double xu(double x, int digits, double c)
{
   return x-yt(x, TT(digits,c), c)*sin(atan(dyc(x, c, P(digits), M(digits,c))));
}

double yu(double x, int digits, double c)
{
   double MM = M(digits,c);
   double PP = P(digits);
   double TTT = TT(digits,c);
   return yc(x, c, PP, MM) + yt(x, TTT, c)*cos(atan(dyc(x, c, PP, MM)));
}

// inferior profile
double xl(double x, int digits, double c)
{
   return x+yt(x, TT(digits,c), c)*sin(atan(dyc(x, c, P(digits), M(digits,c))));
}

double yl(double x, int digits, double c)
{
   double MM = M(digits,c);
   double PP = P(digits);
   double TTT = TT(digits,c);
   return yc(x, c, PP, MM) - yt(x, TTT, c)*cos(atan(dyc(x, c, PP, MM)));
}

// with attack angle
// superior profile
double xua(double x, int digits, double c, double angle)
{
   return (xu(x, digits, c)-c/2.)*cos(angle) + yu(x, digits, c) * sin(angle) + c/2.;
}

double yua(double x, int digits, double c, double angle)
{
   return yu(x, digits, c)*cos(angle) - (xu(x, digits, c)-c/2.) * sin(angle);
}

// inferior profile
double xla(double x, int digits, double c, double angle)
{
   return (xl(x, digits, c)-c/2.)*cos(angle) + yl(x, digits, c) * sin(angle) + c/2.;
}

double yla(double x, int digits, double c, double angle)
{
   return yl(x, digits, c)*cos(angle) - (xl(x, digits, c)-c/2.) * sin(angle);
}

void NACAPoints(int FourDigits, TPZVec< TPZVec<REAL> > & pt, TPZVec< TPZVec< int> > &elms, int nSubdiv)
{
/*
 cout << "\nNumber of Points along NACA\n";
 cin >> m;
 */

 m = nSubdiv;
/*
 cout << "\nNumber of Points at the input BC\n";
 cin >> l;
 cout << "\nNumber of Elements between NACA and BC\n";
 cin >> n;
 cout << "\nqm\n";
 cin >> q;
 cout << "\nqn\n";
 cin >> qn;

n = m / 4;
l = m / 2;
q  = pow(4., 1./(double)m);
qn = pow(2., 1./(double)n);


 cout << "P" << P(1224) << endl;
 cout << "M" << M(1224,7.) << endl;
 cout << "TT" << TT(1224,7.) << endl;

 cout << "xu" << xu(6.05, 1224, 7) << endl;
 cout << "yu" << yu(6.05, 1224, 7) << endl;
 cout << "xl" << xl(6.05, 1224, 7) << endl;
 cout << "yl" << yl(6.05, 1224, 7) << endl;

*/
 n = 5 * m / 9 * (int) sqrt(scale);
 l = m / 3;
 k = n / 2 + 1;

 q  = pow(5.6 * scale, 1./(double)m);
 qn = pow(2.6 * scale, 1./(double)n);


   int index, indexPt, indexPt2;
   double angle;

   TPZVec<REAL> coord(3), coordBC(3), coordNACA(3);
   elms.Resize(m*n*2);
   pt.Resize(2 * m + n * (2*m+1) + /*exit elements*/ (2*n+1)*k);

   cout << "\nNumber of Points: "<< pt.NElements();

   cout << "\nAirfoil angle [degrees]\n";

   cin >> angle;

   angle *= PI/180.;

   // defining points on the surface of the airfoil
   // i: airfoil index;
   // j: radial index
   int i, j;
   for(i = 1; i < m; i++)
      {
         double x = xpg(q, i, m);
	 coord[0] = xua(x * cord, FourDigits, cord, angle) + entrance;
	 coord[1] = yua(x * cord, FourDigits, cord, angle) + height/2.;
	 coord[2] = 0.;
	 pt[i] = coord;

	 coord[0] = xla(x * cord, FourDigits, cord, angle) + entrance;
	 coord[1] = yla(x * cord, FourDigits, cord, angle) + height/2.;
	 coord[2] = 0.;
	 pt[2*m -i] = coord;
      }
   coord[0] = xua(0., FourDigits, cord, angle) + entrance;///*xu(0., FourDigits, cord) +*/ entrance + ;
   coord[1] = yua(0. * cord, FourDigits, cord, angle) + height/2.;///*yu(0., FourDigits, cord) +*/ height/2.;
   coord[2] = 0.;
   pt[0] = coord;

   coord[0] = xla(1. * cord, FourDigits, cord, angle) + entrance;//cord + /*xu(1., FourDigits, cord) +*/ entrance;
   coord[1] = yla(1. * cord, FourDigits, cord, angle) + height/2.;///*yu(1., FourDigits, cord) +*/ height/2.;
   coord[2] = 0.;
   pt[m] = coord;

   //defining points at the boundary
   // index of the leftmost centered point
   index = 2*m + (n-1) * ((2*m)+1);
   coord[0] = 0.;
   coord[1] = height/2.;
   coord[2] = 0.;
   pt[index] = coord;

   for(i = 1; i < l; i++)
   {
      coord[0] = 0.;
      coord[1] =  xpg(1., i, l) * height/2. + height/2.;
      coord[2] = 0.;
      pt[index + i] = coord;

      coord[0] = 0.;
      coord[1] = -xpg(1., i, l) * height/2. + height/2.;
      coord[2] = 0.;
      pt[index + (2*m+1) - i] = coord;
   }

   for(i = l; i <= m; i++)
   {
      coord[0] = xpg(1., i-l, m-l) * (entrance + cord);
      coord[1] = height;
      coord[2] = 0.;
      pt[index + i] = coord;

      coord[0] = xpg(1., i-l, m-l) * (entrance + cord);
      coord[1] = 0.;
      coord[2] = 0.;
      pt[index + (2*m+1) - i] = coord;
   }

   // defining intermediate points
   for(j = 1; j < n; j++)
   {
      //creating division rule
      double ratio = xpg(qn, j, n);

      // resolving entrance centered point
      index = n * ((2*m)+1) - 1; // index of BC point
      indexPt = j * ((2*m)+1) - 1; // index of leftmost layer point
      coordNACA = pt[0];
      coordBC   = pt[index];
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt] = coord;

      // resolving exit points
      coordNACA = pt[m];

      indexPt2 = indexPt + m;
      coordBC   = pt[index + m];// upper point
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt2] = coord;

      indexPt2 = indexPt + m + 1;
      coordBC   = pt[index + m + 1];// bottom point
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt2] = coord;

      // resolving other points
      for(i = 1; i < m; i++)
      {
      // upper points
      coordNACA = pt[i];
      indexPt2 = indexPt + i;
      coordBC   = pt[index + i];// BC upper point
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt2] = coord;

      // bottom points
      coordNACA = pt[2*m-i];
      indexPt2 = indexPt+2*m+1-i;
      coordBC   = pt[index+2*m+1-i];// BC bottom point
      coord[0] = ratio * coordBC[0] + (1.-ratio) * coordNACA[0];
      coord[1] = ratio * coordBC[1] + (1.-ratio) * coordNACA[1];
      coord[2] = ratio * coordBC[2] + (1.-ratio) * coordNACA[2];
      pt[indexPt2] = coord;

      }

   }

   // Creating exit points
   indexPt = 2 * m + n * (2*m+1); // first exit point index
   // resolving ponts closer to NACA
   for(i = 1; i <= k; i++)
   {
      // center point
      index = m; // index of existent point
      indexPt2 = indexPt + n + (2*n+1)*(i-1) ;
      coord[0] = entrance + cord + exitlength * xpg(qn, i, k);
      coord[1] = pt[index][1];
      coord[2] = 0.;
      pt[indexPt2] = coord;
      for(j = 1; j <= n; j++)
      {
         //upper points
         index = j * (2 * m + 1) + m - 1; // index of existent point
         indexPt2 = indexPt + n - j + (2*n+1)*(i-1);
         coord[0] = entrance + cord + exitlength * xpg(qn, i, k);
         coord[1] = pt[index][1];
         coord[2] = 0.;
         pt[indexPt2] = coord;

         // bottom points
         index = j * (2 * m + 1) + m; // index of existent point
         indexPt2 = indexPt + n + j + (2*n+1)*(i-1);
         coord[0] = entrance + cord + exitlength * xpg(qn, i, k);
         coord[1] = pt[index][1];
         coord[2] = 0.;
         pt[indexPt2] = coord;
      }
   }

   // definig elements
   // quadrilateral data
   TPZVec< int > nodes(4);
   elms.Resize(m * n * 2 + 2 * n * k);

   cout << "\nNumber of Elements: "<< elms.NElements() << endl;

   // the first row (closest to NACA is special)
   //upper elements
   nodes[0] = 0;
   nodes[1] = 1;
   nodes[2] = 1+2*m;
   nodes[3] = 2*m;
   elms[0] = nodes;
   // bottom elements
   nodes[0] = 2*m-1;
   nodes[1] = 0;
   nodes[2] = 2*m;
   nodes[3] = 4*m;
   elms[2*m-1] = nodes;
   for(i = 1; i < m; i++)
   {
   //upper elements
      index = i;
      nodes[0] = index;
      nodes[1] = nodes[0]+1;
      nodes[2] = nodes[1]+2*m;
      nodes[3] = nodes[2]-1;
      elms[index] = nodes;
   // bottom elements
      index = 2*m-i-1;
      nodes[0] = index;
      nodes[1] = nodes[0] + 1;
      nodes[2] = nodes[1] + 2 * m + 1;
      nodes[3] = nodes[2] - 1;
      elms[index] = nodes;
   }

   // elements at rows farther from NACA profile
   for(j = 1; j < n; j++)
   {
         // leftmost upper elements
         index = j * 2*m;
	 nodes[0] = index+ j - 1;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + 2*m + 1;
	 nodes[3] = nodes[2] - 1;
	 elms[index] = nodes;

	 // leftmost lower elements
         index = (j+1) * 2*m - 1;
	 nodes[0] = index + j;
	 nodes[1] = nodes[0] - 2*m;//j * 2*m;//index + 2;
	 nodes[2] = nodes[1] + 2 * m + 1;
	 nodes[3] = nodes[0] + 2 * m + 1;
	 elms[index] = nodes;
      for(i = 1; i < m; i++)
      {
         // upper elements
         index = i + j * 2*m;
	 nodes[0] = index + j - 1;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + 2*m + 1;
	 nodes[3] = nodes[2] - 1;
	 elms[index] = nodes;

	 // lower elements
         index = (j+1) * 2*m - i - 1;
	 nodes[0] = index + j;
	 nodes[1] = nodes[0] + 1;
	 nodes[2] = nodes[1] + 2 * m + 1;
	 nodes[3] = nodes[2] - 1;
	 elms[index] = nodes;
      }
   }

   // exit elements
   indexPt = 2 * m + n * (2*m+1); // first exit point index
   // elements closer to NACA
   index = 2*m * n + n - 1;
   nodes[0] = m;
   nodes[1] = indexPt + n;
   nodes[2] = nodes[1] - 1;
   nodes[3] = nodes[0] + 2*m;
   elms[index] = nodes;

   index = 2*m * n + n;
   nodes[0] = 3 * m + 1;
   nodes[1] = indexPt + n + 1;
   nodes[2] = nodes[1] - 1;
   nodes[3] = m;
   elms[index] = nodes;
   for(i = 1; i < n; i++)
   {
      // upper elements
      index = 2 * m * n + n - i - 1;
      nodes[0] = 3*m + (i-1) * (2*m+1);
      nodes[1] = indexPt + n - i;
      nodes[2] = nodes[1] - 1;
      nodes[3] = nodes[0] + 2 * m + 1;
      elms[index] = nodes;

      // bottom elements
      index = 2 * m * n + n + i;
      nodes[0] = 3*m + i * (2*m+1) + 1;
      //nodes[0] = 3*m + (i-1) * (2*m+1) + 1;
      nodes[1] = indexPt + n + i + 1;
      nodes[2] = nodes[1] - 1;
      nodes[3] = nodes[0] - (2*m + 1);
      elms[index] = nodes;
   }
   // other elements
   for(i = 0; i < 2 * n; i++)
      for(j = 1; j < k; j++)
      {
         index = 2*m*n + j*2*n + i;
	 nodes[0] = indexPt + (2*n+1) * (j-1) + i + 1;
	 nodes[1] = nodes[0] + (2*n+1);
	 nodes[2] = nodes[1] - 1;
	 nodes[3] = nodes[0] - 1;
	 elms[index] = nodes;
      }
}

TPZGeoMesh * CreateNACAGeoMesh(TPZVec< TPZVec< REAL > > & nodes,
                           TPZVec< TPZVec< int > > & elms,
			   MElementType ElType, int matId,
			   TPZVec<TPZGeoEl *> & gEls,
			   int nSubdiv)
{
   TPZGeoMesh * gmesh = new TPZGeoMesh;

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
/*
   if(nSubdiv == 1)
   {
// Dividing elements to create a mesh of 4 elems.

     TPZVec< TPZGeoEl * > firstDivision;
     gEls[0]->Divide(firstDivision);
   }

   if(nSubdiv > 1)PZError << "CreateOneElGeoMesh unsupported number of subdivisions";
*/
   return gmesh;
}


// Creating all the geometric and computational meshes
// for the reflected shock problem.

TPZFlowCompMesh *
   NACACompMesh(REAL CFL, REAL delta,
                 int degree, int nSubdiv,
		 TPZArtDiffType DiffType,
		 TPZTimeDiscr Diff_TD,
		 TPZTimeDiscr ConvVol_TD,
		 TPZTimeDiscr ConvFace_TD)
{
   TPZCompElDisc::gDegree = degree;
   REAL gamma = 1.4;
   int i;
   double Mach;

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
   NACAPoints(digits, nodes, elms, nSubdiv);

// Creating the geometric mesh
   TPZGeoMesh * gmesh = CreateNACAGeoMesh(nodes, elms, EQuadrilateral, 1, gElem, nSubdiv);

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

   //aresta interna NACA: Wall
   val1.Zero();
   val2.Zero();
   for( i = 0; i < 2*m; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[i],4,-1,*gmesh);
   }
   bc = mat->CreateBC(-1,5,val1,val2);
   cmesh->InsertMaterialObject(bc);

   cout << "\nMach number\n";
   cin >> Mach;

   // leftmost bc face: Inlet
   val1.Zero();
   val2.Zero();
   val2(0,0) = 1.;// rho
   val2(1,0) = Mach;// Mach
   val2(3,0) = 2.;// pressure
   for( i = 0; i < l; i++)
   {
      TPZGeoElBC((TPZGeoEl *)gElem[(n-1)*2*m+i],6,-2,*gmesh);
      TPZGeoElBC((TPZGeoEl *)gElem[n*2*m-i-1]  ,6,-2,*gmesh);
   }
   //bc = mat->CreateBC(-2,7,val1,val2);
   bc = mat->CreateBC(-2,10,val1,val2); // inflow/outflow
   cmesh->InsertMaterialObject(bc);

   // upper and lower extern NACA BC faces
   // Wall
   for( i = (n-1)*2*m + l; i < n*2*m - l; i++)
   {
      //TPZGeoElBC((TPZGeoEl *)gElem[i],6,-3,*gmesh);
      TPZGeoElBC((TPZGeoEl *)gElem[i],6,-3,*gmesh);
   }
   // exit upper and bottom faces: Wall
   for(i = 0; i < k; i++)
   {
      /*TPZGeoElBC((TPZGeoEl *)gElem[2*m*n+i*n*2],6,-3,*gmesh);
      TPZGeoElBC((TPZGeoEl *)gElem[2*m*n+(i+1)*n*2-1],4,-3,*gmesh);*/
      TPZGeoElBC((TPZGeoEl *)gElem[2*m*n+i*n*2],6,-3,*gmesh);
      TPZGeoElBC((TPZGeoEl *)gElem[2*m*n+(i+1)*n*2-1],4,-3,*gmesh);
   }
   // rightmost exit face: Exit
   for(i = 0; i < 2*n; i++)
   {
      //TPZGeoElBC((TPZGeoEl *)gElem[2*m*n+(k-1)*n*2+i],5,-4,*gmesh);
      TPZGeoElBC((TPZGeoEl *)gElem[2*m*n+(k-1)*n*2+i],5,-3,*gmesh);
   }

   bc = mat->CreateBC(-3,9,val1,val2); // inflow/outflow
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
   for(i = 0; i < nVars; i++)Solution(i,0) = 0;//.05;
   int j, NSolutionBlocks;
   //TPZBlock * pBlock = cmesh->Block();
   NSolutionBlocks = cmesh->Block().NBlocks();
   int nShape = Solution.Rows() / NSolutionBlocks / (dim + 2);
   int lastShapeFun = (nShape - 1)*(dim+2);
   for(j = 0; j < NSolutionBlocks; j++)
   {
      int blockOffset = cmesh->Block().Position(j) + lastShapeFun;

      REAL rho = 1.0,
           p = 2.,
	   u = sqrt(1.4 * p / rho) * Mach,
	   v = 0.,
	   vel2 = u*u + v*v;
      Solution(blockOffset  ,0) = rho;
      Solution(blockOffset+1,0) = rho * u;
      Solution(blockOffset+2,0) = rho * v;
      Solution(blockOffset+3,0) = p/(gamma-1.0) + 0.5 * rho * vel2;
   }

   cmesh->LoadSolution(Solution);

   return cmesh;
}
