#include <iostream>
#include "pzfmatrix.h"
#include "pzelasmat.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzlog.h"
#include "tpzgeoblend.h"
#include "tpzellipse3d.h"
#include "tpzarc3d.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeopoint.h"
#include "pzgeopyramid.h"
#include "pzgeoquad.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzgmesh.h"
#include "tpzmathtools.h"

#include "TPZPlaneFracture.h"

//#include "tpzchangeel.h"
//#include "tpzquadraticline.h"
//#include "tpzquadratictrig.h"
//#include "tpzquadraticquad.h"
//#include "tpzquadratictetra.h"
//#include "tpzquadraticpyramid.h"
//#include "tpzquadraticprism.h"
//#include "tpzquadraticcube.h"

#include "pzgengrid.h"

using namespace std;

const int fractPlaneId = 1;
const int matPoint = -1;

//** returns an xy plane geoMesh based in quadrilaterals */
TPZGeoMesh * GeneratePlaneMesh(int nrows, int ncols, double deltaX, double deltaY);

//** just for visualize given dots in vtk */
TPZGeoMesh * GeneratePlaneMesh(int nrows, int ncols, TPZVec<REAL> &fractureDots, double deltaX, double deltaY);

void FillFractureDotsExampleEllipse(TPZVec<REAL> &fractureDots);

//----------------------------------------------------------------------------------------------------------------------------------

/*
 int mainSimple(int argc, char * const argv[])
 {	
 gRefDBase.InitializeRefPatterns();
 
 int nrows = 6;
 int ncols = 11;
 TPZGeoMesh * gmesh = GenerateGMesh(nrows, ncols);	
 gmesh->Print();
 
 int QTDdivisions = 8;
 TPZPlaneFracture plfrac(gmesh, 0, 10, 55, QTDdivisions);
 int nnodes = 9;
 TPZVec<REAL> fractureDots(3*nnodes);
 for(int n = 0; n < 3*nnodes; n++)
 {
 fractureDots[n] = 0.;
 }
 int node;
 
 node = 0;
 fractureDots[3*node] = 0.;		fractureDots[3*node+1] = 8.;
 
 node = 1;
 fractureDots[3*node] = 4;	fractureDots[3*node+1] = 8.;
 
 node = 2;
 fractureDots[3*node] = 6.;	fractureDots[3*node+1] = 8.;
 
 node = 3;
 fractureDots[3*node] = 6.;	fractureDots[3*node+1] = 6.;
 
 node = 4;
 fractureDots[3*node] = 10.;	fractureDots[3*node+1] = 4.;
 
 node = 5;
 fractureDots[3*node] = 8.;	fractureDots[3*node+1] = 2.;
 
 node = 6;
 fractureDots[3*node] = 6.;	fractureDots[3*node+1] = 4.;
 
 node = 7;
 fractureDots[3*node] = 2.;	fractureDots[3*node+1] = 4.;
 
 node = 8;
 fractureDots[3*node] = 0.;	fractureDots[3*node+1] = 4.;
 
 //
 TPZGeoMesh * gmesh4VTK = GenerateGMesh(nrows,ncols,fractureDots);
 std::ofstream outDots("FractureDots.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh4VTK, outDots, false);
 //
 
 TPZGeoMesh * fractureMesh = plfrac.GetFractureMesh(fractureDots);
 
 std::ofstream out("Fracture.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(fractureMesh, out, true);
 
 return 0;
 }
 */

//Crack Ellipse

int main(int argc, char * const argv[])
{	
    gRefDBase.InitializeRefPatterns();
    //gRefDBase.InitializeUniformRefPattern(ETriangle);
    
    int nrows = 51;
    int ncols = 101;
    double deltaX = 2.5;
    double deltaY = 4.;
    TPZGeoMesh * gmesh = GeneratePlaneMesh(nrows,ncols, deltaX,deltaY);
        
    int QTDdivisions = 4;
    TPZPlaneFracture plfrac(gmesh, 0, 1, ncols, QTDdivisions);
    
    TPZVec<REAL> fractureDots;
    FillFractureDotsExampleEllipse(fractureDots);
    
    TPZGeoMesh * gmesh4VTK = GeneratePlaneMesh(nrows,ncols,fractureDots, deltaX, deltaY);
    std::ofstream outDots("FractureDots.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh4VTK, outDots, false);
    
    TPZGeoMesh * fractureMesh = plfrac.GetFractureMesh(fractureDots);
    
    std::ofstream out("Fracture.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fractureMesh, out, true);
    
    return 0;
}

//Crack Crazy
/*
 int mainCrazy(int argc, char * const argv[])
 {	
 gRefDBase.InitializeRefPatterns();
 
 int nrows = 51;
 int ncols = 101;
 TPZGeoMesh * gmesh = GenerateGMesh(nrows,ncols);
 
 int QTDdivisions = 15;
 TPZPlaneFracture plfrac(gmesh, 0, 1, ncols, QTDdivisions);
 int nnodes = 62;
 TPZVec<REAL> fractureDots(3*nnodes);
 for(int n = 0; n < 3*nnodes; n++)
 {
 fractureDots[n] = 0.;
 }
 int node;
 
 node = 0;
 double desloc = 0.;
 
 fractureDots[3*node] = 3.; fractureDots[3*node+1] = 80.;
 
 node = 1;
 
 fractureDots[3*node] = 9.02368; fractureDots[3*node+1] = 82.1856;
 
 node = 2;
 
 fractureDots[3*node] = 25.6151; fractureDots[3*node+1] = 82.9122;
 
 node = 3;
 
 fractureDots[3*node] = 48.4505; fractureDots[3*node+1] = 82.5833;
 
 node = 4;
 
 fractureDots[3*node] = 74.1186; fractureDots[3*node+1] = 81.5259;
 
 node = 5;
 
 fractureDots[3*node] = 100. + desloc; fractureDots[3*node+1] = 80.;//<----------------------
 
 node = 6;
 
 fractureDots[3*node] = 124.157; fractureDots[3*node+1] = 78.2079;
 
 node = 7;
 
 fractureDots[3*node] = 145.23; fractureDots[3*node+1] = 76.3027;
 
 node = 8;
 
 fractureDots[3*node] = 162.349; fractureDots[3*node+1] = 74.3954;
 
 node = 9;
 
 fractureDots[3*node] = 175.043; fractureDots[3*node+1] = 72.5625;
 
 node = 10;
 
 fractureDots[3*node] = 183.172; fractureDots[3*node+1] = 70.8516;
 
 node = 11;
 
 fractureDots[3*node] = 186.85; fractureDots[3*node+1] = 69.2873;
 
 node = 12;
 
 fractureDots[3*node] = 186.391; fractureDots[3*node+1] = 67.8763;
 
 node = 13;
 
 fractureDots[3*node] = 182.253; fractureDots[3*node+1] = 66.6118;
 
 node = 14;
 
 fractureDots[3*node] = 174.99; fractureDots[3*node+1] = 65.4771;
 
 node = 15;
 
 fractureDots[3*node] = 165.212; fractureDots[3*node+1] = 64.4495;
 
 node = 16;
 
 fractureDots[3*node] = 153.552; fractureDots[3*node+1] = 63.5025;
 
 node = 17;
 
 fractureDots[3*node] = 140.633; fractureDots[3*node+1] = 62.609;
 
 node = 18;
 
 fractureDots[3*node] = 127.05; fractureDots[3*node+1] = 61.7426;
 
 node = 19;
 
 fractureDots[3*node] = 113.346; fractureDots[3*node+1] = 60.8796;
 
 node = 20;
 
 fractureDots[3*node] = 100. + desloc; fractureDots[3*node+1] = 60.;//<----------------------
 
 node = 21;
 
 fractureDots[3*node] = 87.418; fractureDots[3*node+1] = 59.0884;
 
 node = 22;
 
 fractureDots[3*node] = 75.9253; fractureDots[3*node+1] = 58.1348;
 
 node = 23;
 
 fractureDots[3*node] = 65.7644; fractureDots[3*node+1] = 57.1342;
 
 node = 24;
 
 fractureDots[3*node] = 57.0956; fractureDots[3*node+1] = 56.0873;
 
 node = 25;
 
 fractureDots[3*node] = 50.; fractureDots[3*node+1] = 55.;
 
 node = 26;
 
 fractureDots[3*node] = 44.4857; fractureDots[3*node+1] = 53.8827;
 
 node = 27;
 
 fractureDots[3*node] = 40.495; fractureDots[3*node+1] = 52.75;
 
 node = 28;
 
 fractureDots[3*node] = 37.9143; fractureDots[3*node+1] = 51.6198;
 
 node = 29;
 
 fractureDots[3*node] = 36.5847; fractureDots[3*node+1] = 50.5124;
 
 node = 30;
 
 fractureDots[3*node] = 36.3137; fractureDots[3*node+1] = 49.4496;
 
 node = 31;
 
 fractureDots[3*node] = 36.888; fractureDots[3*node+1] = 48.4535;
 
 node = 32;
 
 fractureDots[3*node] = 38.0857; fractureDots[3*node+1] = 47.5455;
 
 node = 33;
 
 fractureDots[3*node] = 39.6892; fractureDots[3*node+1] = 46.7455;
 
 node = 34;
 
 fractureDots[3*node] = 41.4972; fractureDots[3*node+1] = 46.0702;
 
 node = 35;
 
 fractureDots[3*node] = 43.3363; fractureDots[3*node+1] = 45.5328;
 
 node = 36;
 
 fractureDots[3*node] = 45.0707; fractureDots[3*node+1] = 45.1417;
 
 node = 37;
 
 fractureDots[3*node] = 46.6112; fractureDots[3*node+1] = 44.8996;
 
 node = 38;
 
 fractureDots[3*node] = 47.9218; fractureDots[3*node+1] = 44.8032;
 
 node = 39;
 
 fractureDots[3*node] = 49.0243; fractureDots[3*node+1] = 44.8424;
 
 node = 40;
 
 fractureDots[3*node] = 50.; fractureDots[3*node+1] = 45.;
 
 node = 41;
 
 fractureDots[3*node] = 50.9889; fractureDots[3*node+1] = 45.2519;
 
 node = 42;
 
 fractureDots[3*node] = 52.1852; fractureDots[3*node+1] = 45.5668;
 
 node = 43;
 
 fractureDots[3*node] = 53.8292; fractureDots[3*node+1] = 45.9072;
 
 node = 44;
 
 fractureDots[3*node] = 56.1951; fractureDots[3*node+1] = 46.2297;
 
 node = 45;
 
 fractureDots[3*node] = 59.5747; fractureDots[3*node+1] = 46.4863;
 
 node = 46;
 
 fractureDots[3*node] = 64.2563; fractureDots[3*node+1] = 46.626;
 
 node = 47;
 
 fractureDots[3*node] = 70.4977; fractureDots[3*node+1] = 46.5964;
 
 node = 48;
 
 fractureDots[3*node] = 78.4951; fractureDots[3*node+1] = 46.3461;
 
 node = 49;
 
 fractureDots[3*node] = 88.3449; fractureDots[3*node+1] = 45.8275;
 
 node = 50;
 
 fractureDots[3*node] = 100.; fractureDots[3*node+1] = 45.;
 
 node = 51;
 
 fractureDots[3*node] = 113.22; fractureDots[3*node+1] = 43.8338;
 
 node = 52;
 
 fractureDots[3*node] = 127.511; fractureDots[3*node+1] = 42.3139;
 
 node = 53;
 
 fractureDots[3*node] = 142.068; fractureDots[3*node+1] = 40.4454;
 
 node = 54;
 
 fractureDots[3*node] = 155.692; fractureDots[3*node+1] = 38.2585;
 
 node = 55;
 
 fractureDots[3*node] = 166.721; fractureDots[3*node+1] = 35.8147;
 
 node = 56;
 
 fractureDots[3*node] = 172.932; fractureDots[3*node+1] = 33.2139;
 
 node = 57;
 
 fractureDots[3*node] = 171.452; fractureDots[3*node+1] = 30.601;
 
 node = 58;
 
 fractureDots[3*node] = 158.644; fractureDots[3*node+1] = 28.175;
 
 node = 59;
 
 fractureDots[3*node] = 129.998; fractureDots[3*node+1] = 26.1969;
 
 node = 60;
 
 fractureDots[3*node] = 80.; fractureDots[3*node+1] = 25.;
 
 node = 61;
 
 fractureDots[3*node] = 2.; fractureDots[3*node+1] = 25.;
 
 //
 TPZGeoMesh * gmesh4VTK = GenerateGMesh(nrows,ncols,fractureDots);
 std::ofstream outDots("FractureDots.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh4VTK, outDots, false);
 //
 
 TPZGeoMesh * fractureMesh = plfrac.GetFractureMesh(fractureDots);
 
 std::ofstream out("Fracture.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(fractureMesh, out, true);
 
 return 0;
 }
 */

//Crack Crazy Crossed
/*
 int mainCrazyCrossed(int argc, char * const argv[])
 {	
 gRefDBase.InitializeRefPatterns();
 
 int nrows = 51;
 int ncols = 101;
 TPZGeoMesh * gmesh = GenerateGMesh(nrows,ncols);
 
 int QTDdivisions = 5;
 TPZPlaneFracture plfrac(gmesh, 0, 1, ncols, QTDdivisions);
 int nnodes = 62;
 TPZVec<REAL> fractureDots(3*nnodes);
 for(int n = 0; n < 3*nnodes; n++)
 {
 fractureDots[n] = 0.;
 }
 int node;
 
 node = 0;
 double desloc = 0.;
 
 fractureDots[3*node] = 3.; fractureDots[3*node+1] = 80.;
 
 node = 1;
 
 fractureDots[3*node] = 9.02368; fractureDots[3*node+1] = 82.1856;
 
 node = 2;
 
 fractureDots[3*node] = 25.6151; fractureDots[3*node+1] = 82.9122;
 
 node = 3;
 
 fractureDots[3*node] = 48.4505; fractureDots[3*node+1] = 82.5833;
 
 node = 4;
 
 fractureDots[3*node] = 74.1186; fractureDots[3*node+1] = 81.5259;
 
 node = 5;
 
 fractureDots[3*node] = 100. + desloc; fractureDots[3*node+1] = 80.;//<----------------------
 
 node = 6;
 
 fractureDots[3*node] = 124.157; fractureDots[3*node+1] = 78.2079;
 
 node = 7;
 
 fractureDots[3*node] = 145.23; fractureDots[3*node+1] = 76.3027;
 
 node = 8;
 
 fractureDots[3*node] = 162.349; fractureDots[3*node+1] = 74.3954;
 
 node = 9;
 
 fractureDots[3*node] = 175.043; fractureDots[3*node+1] = 72.5625;
 
 node = 10;
 
 fractureDots[3*node] = 183.172; fractureDots[3*node+1] = 70.8516;
 
 node = 11;
 
 fractureDots[3*node] = 186.85; fractureDots[3*node+1] = 69.2873;
 
 node = 12;
 
 fractureDots[3*node] = 186.391; fractureDots[3*node+1] = 67.8763;
 
 node = 13;
 
 fractureDots[3*node] = 182.253; fractureDots[3*node+1] = 66.6118;
 
 node = 14;
 
 fractureDots[3*node] = 174.99; fractureDots[3*node+1] = 65.4771;
 
 node = 15;
 
 fractureDots[3*node] = 165.212; fractureDots[3*node+1] = 64.4495;
 
 node = 16;
 
 fractureDots[3*node] = 153.552; fractureDots[3*node+1] = 63.5025;
 
 node = 17;
 
 fractureDots[3*node] = 140.633; fractureDots[3*node+1] = 42.609;
 
 node = 18;
 
 fractureDots[3*node] = 127.05; fractureDots[3*node+1] = 41.7426;
 
 node = 19;
 
 fractureDots[3*node] = 113.346; fractureDots[3*node+1] = 40.8796;
 
 node = 20;
 
 fractureDots[3*node] = 100. + desloc; fractureDots[3*node+1] = 60.;//<----------------------
 
 node = 21;
 
 fractureDots[3*node] = 87.418; fractureDots[3*node+1] = 59.0884;
 
 node = 22;
 
 fractureDots[3*node] = 75.9253; fractureDots[3*node+1] = 58.1348;
 
 node = 23;
 
 fractureDots[3*node] = 65.7644; fractureDots[3*node+1] = 57.1342;
 
 node = 24;
 
 fractureDots[3*node] = 57.0956; fractureDots[3*node+1] = 56.0873;
 
 node = 25;
 
 fractureDots[3*node] = 50.; fractureDots[3*node+1] = 55.;
 
 node = 26;
 
 fractureDots[3*node] = 44.4857; fractureDots[3*node+1] = 53.8827;
 
 node = 27;
 
 fractureDots[3*node] = 40.495; fractureDots[3*node+1] = 52.75;
 
 node = 28;
 
 fractureDots[3*node] = 37.9143; fractureDots[3*node+1] = 51.6198;
 
 node = 29;
 
 fractureDots[3*node] = 36.5847; fractureDots[3*node+1] = 50.5124;
 
 node = 30;
 
 fractureDots[3*node] = 36.3137; fractureDots[3*node+1] = 49.4496;
 
 node = 31;
 
 fractureDots[3*node] = 36.888; fractureDots[3*node+1] = 48.4535;
 
 node = 32;
 
 fractureDots[3*node] = 38.0857; fractureDots[3*node+1] = 47.5455;
 
 node = 33;
 
 fractureDots[3*node] = 39.6892; fractureDots[3*node+1] = 46.7455;
 
 node = 34;
 
 fractureDots[3*node] = 41.4972; fractureDots[3*node+1] = 46.0702;
 
 node = 35;
 
 fractureDots[3*node] = 43.3363; fractureDots[3*node+1] = 45.5328;
 
 node = 36;
 
 fractureDots[3*node] = 45.0707; fractureDots[3*node+1] = 45.1417;
 
 node = 37;
 
 fractureDots[3*node] = 46.6112; fractureDots[3*node+1] = 44.8996;
 
 node = 38;
 
 fractureDots[3*node] = 47.9218; fractureDots[3*node+1] = 44.8032;
 
 node = 39;
 
 fractureDots[3*node] = 49.0243; fractureDots[3*node+1] = 44.8424;
 
 node = 40;
 
 fractureDots[3*node] = 50.; fractureDots[3*node+1] = 45.;
 
 node = 41;
 
 fractureDots[3*node] = 50.9889; fractureDots[3*node+1] = 45.2519;
 
 node = 42;
 
 fractureDots[3*node] = 52.1852; fractureDots[3*node+1] = 45.5668;
 
 node = 43;
 
 fractureDots[3*node] = 53.8292; fractureDots[3*node+1] = 45.9072;
 
 node = 44;
 
 fractureDots[3*node] = 56.1951; fractureDots[3*node+1] = 46.2297;
 
 node = 45;
 
 fractureDots[3*node] = 59.5747; fractureDots[3*node+1] = 46.4863;
 
 node = 46;
 
 fractureDots[3*node] = 64.2563; fractureDots[3*node+1] = 46.626;
 
 node = 47;
 
 fractureDots[3*node] = 70.4977; fractureDots[3*node+1] = 46.5964;
 
 node = 48;
 
 fractureDots[3*node] = 78.4951; fractureDots[3*node+1] = 46.3461;
 
 node = 49;
 
 fractureDots[3*node] = 88.3449; fractureDots[3*node+1] = 45.8275;
 
 node = 50;
 
 fractureDots[3*node] = 100.; fractureDots[3*node+1] = 45.;
 
 node = 51;
 
 fractureDots[3*node] = 113.22; fractureDots[3*node+1] = 63.8338;
 
 node = 52;
 
 fractureDots[3*node] = 127.511; fractureDots[3*node+1] = 62.3139;
 
 node = 53;
 
 fractureDots[3*node] = 142.068; fractureDots[3*node+1] = 60.4454;
 
 node = 54;
 
 fractureDots[3*node] = 155.692; fractureDots[3*node+1] = 38.2585;
 
 node = 55;
 
 fractureDots[3*node] = 166.721; fractureDots[3*node+1] = 35.8147;
 
 node = 56;
 
 fractureDots[3*node] = 172.932; fractureDots[3*node+1] = 33.2139;
 
 node = 57;
 
 fractureDots[3*node] = 171.452; fractureDots[3*node+1] = 30.601;
 
 node = 58;
 
 fractureDots[3*node] = 158.644; fractureDots[3*node+1] = 28.175;
 
 node = 59;
 
 fractureDots[3*node] = 129.998; fractureDots[3*node+1] = 26.1969;
 
 node = 60;
 
 fractureDots[3*node] = 80.; fractureDots[3*node+1] = 25.;
 
 node = 61;
 
 fractureDots[3*node] = 2.; fractureDots[3*node+1] = 25.;
 
 //
 TPZGeoMesh * gmesh4VTK = GenerateGMesh(nrows,ncols,fractureDots);
 std::ofstream outDots("FractureDots.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(gmesh4VTK, outDots, false);
 //
 
 TPZGeoMesh * fractureMesh = plfrac.GetFractureMesh(fractureDots);
 
 std::ofstream out("Fracture.vtk");
 TPZVTKGeoMesh::PrintGMeshVTK(fractureMesh, out, true);
 
 return 0;
 }
 //----------------------------------------------------------------------------------------------------------------------------------
 */


TPZGeoMesh * GeneratePlaneMesh(int nrows, int ncols, double deltaX, double deltaY)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxElementId((nrows-1)*(ncols-1) - 1);
	gmesh->SetMaxNodeId(nrows*ncols - 1);
	
	const int Qnodes = nrows*ncols;
	TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
	
	int nodeId = 0;
	
	//setting nodes coords
	for(int r = 0; r < nrows; r++)
	{
		for(int c = 0; c < ncols; c++)
		{
			NodeCoord[nodeId][0] = c*deltaX;
			NodeCoord[nodeId][1] = r*deltaY;
			nodeId++;
		}
	}
	
	//initializing gmesh->NodeVec()
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		gmesh->NodeVec()[n] = Node[n]; 
	}
	
	//inserting quadrilaterals
	int elId = 0;
	TPZVec <int> Topol(4);
	for(int r = 0; r < (nrows-1); r++)
	{
		for(int c = 0; c < (ncols-1); c++)
		{
			Topol[0] = ncols*r+c; Topol[1] = ncols*r+c+1; Topol[2] = ncols*(r+1)+c+1; Topol[3] = ncols*(r+1)+c;
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,fractPlaneId,*gmesh);
			elId++;
		}
	}
	
	gmesh->BuildConnectivity();
	
	return gmesh;
}

//** just for visualize given dots in vtk */
TPZGeoMesh * GeneratePlaneMesh(int nrows, int ncols, TPZVec<REAL> &fractureDots, double deltaX, double deltaY)
{
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->SetMaxElementId((nrows-1)*(ncols-1) - 1);
	gmesh->SetMaxNodeId(nrows*ncols - 1);
	
	const int Qnodes = nrows*ncols + (fractureDots.NElements()/3);
	TPZVec < TPZVec <REAL> > NodeCoord(Qnodes);
	for(int i = 0; i < Qnodes; i++) NodeCoord[i].Resize(3,0.);
	
	int nodeId = 0;
	
	//setting nodes coords
	for(int r = 0; r < nrows; r++)
	{
		for(int c = 0; c < ncols; c++)
		{
			NodeCoord[nodeId][0] = c*deltaX;
			NodeCoord[nodeId][1] = r*deltaY;
			nodeId++;
		}
	}
	for(int r = 0; r < (fractureDots.NElements()/3); r++)
	{
		NodeCoord[nodeId][0] = fractureDots[3*r];
		NodeCoord[nodeId][1] =fractureDots[3*r+1];
		nodeId++;
	}	
	
	//initializing gmesh->NodeVec()
	gmesh->NodeVec().Resize(Qnodes);
	TPZVec <TPZGeoNode> Node(Qnodes);
	for(int n = 0; n < Qnodes; n++)
	{
		Node[n].SetNodeId(n);
		Node[n].SetCoord(NodeCoord[n]);
		gmesh->NodeVec()[n] = Node[n]; 
	}
	
	//inserting quadrilaterals
	int elId = 0;
	TPZVec <int> Topol(4);
	for(int r = 0; r < (nrows-1); r++)
	{
		for(int c = 0; c < (ncols-1); c++)
		{
			Topol[0] = ncols*r+c; Topol[1] = ncols*r+c+1; Topol[2] = ncols*(r+1)+c+1; Topol[3] = ncols*(r+1)+c;
			new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elId,Topol,fractPlaneId,*gmesh);
			elId++;
		}
	}
	Topol.Resize(1);
	for(int n = nrows*ncols; n < Qnodes; n++)
	{
		Topol[0] = n;
		new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (elId,Topol, matPoint,*gmesh);
	}
	
	gmesh->BuildConnectivity();
	
	return gmesh;
}

void FillFractureDotsExampleEllipse(TPZVec<REAL> &fractureDots)
{
    int nnodes = 80;
    
    fractureDots.Resize(3*nnodes, 0.);
    int node;
    
    node = 0;
    
    fractureDots[3*node] = 0.5; fractureDots[3*node+1] = 95.;
    
    node = 1;
    
    fractureDots[3*node] = 5.; fractureDots[3*node+1] = 94.9852;
    
    node = 2;
    
    fractureDots[3*node] = 10.; fractureDots[3*node+1] = 94.9408;
    
    node = 3;
    
    fractureDots[3*node] = 15.; fractureDots[3*node+1] = 94.8667;
    
    node = 4;
    
    fractureDots[3*node] = 20.; fractureDots[3*node+1] = 94.7627;
    
    node = 5;
    
    fractureDots[3*node] = 25.; fractureDots[3*node+1] = 94.6286;
    
    node = 6;
    
    fractureDots[3*node] = 30.; fractureDots[3*node+1] = 94.4643;
    
    node = 7;
    
    fractureDots[3*node] = 35.; fractureDots[3*node+1] = 94.2692;
    
    node = 8;
    
    fractureDots[3*node] = 40.; fractureDots[3*node+1] = 94.0431;
    
    node = 9;
    
    fractureDots[3*node] = 45.; fractureDots[3*node+1] = 93.7854;
    
    node = 10;
    
    fractureDots[3*node] = 50.; fractureDots[3*node+1] = 93.4956;
    
    node = 11;
    
    fractureDots[3*node] = 55.; fractureDots[3*node+1] = 93.173;
    
    node = 12;
    
    fractureDots[3*node] = 60.; fractureDots[3*node+1] = 92.8169;
    
    node = 13;
    
    fractureDots[3*node] = 65.; fractureDots[3*node+1] = 92.4264;
    
    node = 14;
    
    fractureDots[3*node] = 70.; fractureDots[3*node+1] = 92.0006;
    
    node = 15;
    
    fractureDots[3*node] = 75.; fractureDots[3*node+1] = 91.5385;
    
    node = 16;
    
    fractureDots[3*node] = 80.; fractureDots[3*node+1] = 91.0387;
    
    node = 17;
    
    fractureDots[3*node] = 85.; fractureDots[3*node+1] = 90.4998;
    
    node = 18;
    
    fractureDots[3*node] = 90.; fractureDots[3*node+1] = 89.9204;
    
    node = 19;
    
    fractureDots[3*node] = 95.; fractureDots[3*node+1] = 89.2986;
    
    node = 20;
    
    fractureDots[3*node] = 100.; fractureDots[3*node+1] = 88.6323;
    
    node = 21;
    
    fractureDots[3*node] = 105.; fractureDots[3*node+1] = 87.9193;
    
    node = 22;
    
    fractureDots[3*node] = 110.; fractureDots[3*node+1] = 87.1567;
    
    node = 23;
    
    fractureDots[3*node] = 115.; fractureDots[3*node+1] = 86.3416;
    
    node = 24;
    
    fractureDots[3*node] = 120.; fractureDots[3*node+1] = 85.4702;
    
    node = 25;
    
    fractureDots[3*node] = 125.; fractureDots[3*node+1] = 84.5384;
    
    node = 26;
    
    fractureDots[3*node] = 130.; fractureDots[3*node+1] = 83.541;
    
    node = 27;
    
    fractureDots[3*node] = 135.; fractureDots[3*node+1] = 82.4721;
    
    node = 28;
    
    fractureDots[3*node] = 140.; fractureDots[3*node+1] = 81.3243;
    
    node = 29;
    
    fractureDots[3*node] = 145.; fractureDots[3*node+1] = 80.0886;
    
    node = 30;
    
    fractureDots[3*node] = 150.; fractureDots[3*node+1] = 78.7537;
    
    node = 31;
    
    fractureDots[3*node] = 155.; fractureDots[3*node+1] = 77.305;
    
    node = 32;
    
    fractureDots[3*node] = 160.; fractureDots[3*node+1] = 75.7233;
    
    node = 33;
    
    fractureDots[3*node] = 165.; fractureDots[3*node+1] = 73.9822;
    
    node = 34;
    
    fractureDots[3*node] = 170.; fractureDots[3*node+1] = 72.0442;
    
    node = 35;
    
    fractureDots[3*node] = 175.; fractureDots[3*node+1] = 69.8515;
    
    node = 36;
    
    fractureDots[3*node] = 180.; fractureDots[3*node+1] = 67.3077;
    
    node = 37;
    
    fractureDots[3*node] = 185.; fractureDots[3*node+1] = 64.2256;
    
    node = 38;
    
    fractureDots[3*node] = 190.; fractureDots[3*node+1] = 60.125;
    
    node = 39;
    
    fractureDots[3*node] = 195.; fractureDots[3*node+1] = 50.;
    
    node = 40;
    
    fractureDots[3*node] = 195.; fractureDots[3*node+1] = 50.;
    
    node = 41;
    
    fractureDots[3*node] = 190.; fractureDots[3*node+1] = 39.875;
    
    node = 42;
    
    fractureDots[3*node] = 185.; fractureDots[3*node+1] = 35.7744;
    
    node = 43;
    
    fractureDots[3*node] = 180.; fractureDots[3*node+1] = 32.6923;
    
    node = 44;
    
    fractureDots[3*node] = 175.; fractureDots[3*node+1] = 30.1485;
    
    node = 45;
    
    fractureDots[3*node] = 170.; fractureDots[3*node+1] = 27.9558;
    
    node = 46;
    
    fractureDots[3*node] = 165.; fractureDots[3*node+1] = 26.0178;
    
    node = 47;
    
    fractureDots[3*node] = 160.; fractureDots[3*node+1] = 24.2767;
    
    node = 48;
    
    fractureDots[3*node] = 155.; fractureDots[3*node+1] = 22.695;
    
    node = 49;
    
    fractureDots[3*node] = 150.; fractureDots[3*node+1] = 21.2463;
    
    node = 50;
    
    fractureDots[3*node] = 145.; fractureDots[3*node+1] = 19.9114;
    
    node = 51;
    
    fractureDots[3*node] = 140.; fractureDots[3*node+1] = 18.6757;
    
    node = 52;
    
    fractureDots[3*node] = 135.; fractureDots[3*node+1] = 17.5279;
    
    node = 53;
    
    fractureDots[3*node] = 130.; fractureDots[3*node+1] = 16.459;
    
    node = 54;
    
    fractureDots[3*node] = 125.; fractureDots[3*node+1] = 15.4616;
    
    node = 55;
    
    fractureDots[3*node] = 120.; fractureDots[3*node+1] = 14.5298;
    
    node = 56;
    
    fractureDots[3*node] = 115.; fractureDots[3*node+1] = 13.6584;
    
    node = 57;
    
    fractureDots[3*node] = 110.; fractureDots[3*node+1] = 12.8433;
    
    node = 58;
    
    fractureDots[3*node] = 105.; fractureDots[3*node+1] = 12.0807;
    
    node = 59;
    
    fractureDots[3*node] = 100.; fractureDots[3*node+1] = 11.3677;
    
    node = 60;
    
    fractureDots[3*node] = 95.; fractureDots[3*node+1] = 10.7014;
    
    node = 61;
    
    fractureDots[3*node] = 90.; fractureDots[3*node+1] = 10.0796;
    
    node = 62;
    
    fractureDots[3*node] = 85.; fractureDots[3*node+1] = 9.50016;
    
    node = 63;
    
    fractureDots[3*node] = 80.; fractureDots[3*node+1] = 8.96134;
    
    node = 64;
    
    fractureDots[3*node] = 75.; fractureDots[3*node+1] = 8.46154;
    
    node = 65;
    
    fractureDots[3*node] = 70.; fractureDots[3*node+1] = 7.99937;
    
    node = 66;
    
    fractureDots[3*node] = 65.; fractureDots[3*node+1] = 7.57359;
    
    node = 67;
    
    fractureDots[3*node] = 60.; fractureDots[3*node+1] = 7.18313;
    
    node = 68;
    
    fractureDots[3*node] = 55.; fractureDots[3*node+1] = 6.82703;
    
    node = 69;
    
    fractureDots[3*node] = 50.; fractureDots[3*node+1] = 6.50444;
    
    node = 70;
    
    fractureDots[3*node] = 45.; fractureDots[3*node+1] = 6.21462;
    
    node = 71;
    
    fractureDots[3*node] = 40.; fractureDots[3*node+1] = 5.95692;
    
    node = 72;
    
    fractureDots[3*node] = 35.; fractureDots[3*node+1] = 5.73079;
    
    node = 73;
    
    fractureDots[3*node] = 30.; fractureDots[3*node+1] = 5.53573;
    
    node = 74;
    
    fractureDots[3*node] = 25.; fractureDots[3*node+1] = 5.37135;
    
    node = 75;
    
    fractureDots[3*node] = 20.; fractureDots[3*node+1] = 5.23731;
    
    node = 76;
    
    fractureDots[3*node] = 15.; fractureDots[3*node+1] = 5.13333;
    
    node = 77;
    
    fractureDots[3*node] = 10.; fractureDots[3*node+1] = 5.05921;
    
    node = 78;
    
    fractureDots[3*node] = 5.; fractureDots[3*node+1] = 5.0148;
    
    node = 79;
    
    fractureDots[3*node] = 0.5; fractureDots[3*node+1] = 5.;
}
