#include <iostream>

#include "pzgmesh.h"
#include "pzgeopoint.h"
#include "TPZRefPatternDataBase.h"
#include "tpzgeoelrefpattern.h"
#include "TPZPlaneFracture.h"

//to delete
#include "tpzchangeel.h"
#include "pzmaterial.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzstrmatrix.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeotetrahedra.h"
#include "TPZVTKGeoMesh.h"
#include "pzelasmat.h"
#include "PlaneFracture/TPZJIntegral.h"

#include "TPZTimer.h"

using namespace std;

const int matPoint = -3;

//** just for visualize given dots in vtk */

void FillFractureDotsExampleEllipse(TPZVec< std::pair<REAL,REAL> > &fractureDots);
void FillFractureDotsExampleCrazy(TPZVec< std::pair<REAL,REAL> > &fractureDots);
void FillFractureDotsCircle(REAL center, REAL radius, TPZVec< std::pair<REAL,REAL> > &fractureDots);

TPZGeoMesh * PlaneMesh(REAL lf, REAL ldom, REAL hdom, REAL lmax);
TPZCompMesh * PlaneMesh(TPZGeoMesh * gmesh, REAL pressureInsideCrack);

//----------------------------------------------------------------------------------------------------------------------------------

/** Exemplo da utilizacao da quadratura adaptativa (integral adaptativa) */


#include "adapt.h"


int mainCRAZY(int argc, char * const argv[])
{	
    std::cout << "\e";
    TPZTimer readRef("ReadingRefPatterns");
    readRef.start();
    
    //#define writeAgain
    #ifdef writeAgain
        gRefDBase.InitializeRefPatterns();
    #else
        std::ifstream inRefP("RefPatternsUsed.txt");
        gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
    #endif
    
    readRef.stop();
    std::cout << "DeltaT leitura refpatterns = " << readRef.seconds() << " s" << std::endl;
    
    REAL lengthX = 5.;
    REAL lengthY = 7.;
    REAL Lmax = 10.;
    
    REAL lw = 202.;
    REAL bulletDepthTVDIni = 20.;
    REAL bulletDepthTVDFin = 80.;
    TPZVec<TPZLayerProperties> layerVec(5);
    //stretch #0
    
    layerVec[0] = TPZLayerProperties(1.E5, 0.25, 1.E5, 1.E1, 0., 11.);
    layerVec[1] = TPZLayerProperties(1.E5, 0.25, 1.E5, 1.E1, 0., 11.);
    layerVec[2] = TPZLayerProperties(1.E5, 0.25, 1.E5, 1.E1, 0., 11.);
    layerVec[3] = TPZLayerProperties(1.E5, 0.25, 1.E5, 1.E1, 0., 11.);
    layerVec[4] = TPZLayerProperties(1.E5, 0.25, 1.E5, 1.E1, 0., 11.);
    
//    posTVD_stress[0][0.]  = 2.;
//    posTVD_stress[0][11.]  = 3.;
//    //stretch #1
//    posTVD_stress[1][11.]  = 4.;
//    posTVD_stress[1][25.]  = 7.;
//    //stretch #2
//    posTVD_stress[2][25.]  = 3.;
//    posTVD_stress[2][37.]  = 6.;
//    //stretch #3
//    posTVD_stress[3][37.]  = 6.;
//    posTVD_stress[3][63.]  = 8.;
//    //stretch #4
//    posTVD_stress[4][63.] = 8.;
//    posTVD_stress[4][210.] = 10.;
    TPZPlaneFracture plfrac(layerVec, bulletDepthTVDIni, bulletDepthTVDFin, lengthX, lengthY, Lmax);
    
    TPZVec< std::pair<REAL,REAL> > fractureDots(0);
    FillFractureDotsExampleEllipse(fractureDots);
    
    TPZTimer clockIni2("PartyBegins2");
    clockIni2.start();    
    
    REAL pressureInsideCrack = 5.;
    std::string vtkFile = "fracturePconstant.vtk";
    plfrac.RunThisFractureGeometry(fractureDots, pressureInsideCrack, vtkFile, true);
    
    clockIni2.stop();
    std::cout << "DeltaT get fracture cmesh = " << clockIni2.seconds() << " s" << std::endl;
  
    std::ofstream outRefP("RefPatternsUsed.txt");
    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}

int mainCircles(int argc, char * const argv[])
{
    std::cout << "\e";
    TPZTimer readRef("ReadingRefPatterns");
    readRef.start();
    
    //#define writeAgain
#ifdef writeAgain
    gRefDBase.InitializeRefPatterns();
#else
    std::ifstream inRefP("RefPatternsUsed.txt");
    gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
#endif
    
    readRef.stop();
    std::cout << "DeltaT leitura refpatterns = " << readRef.seconds() << " s" << std::endl;
    
    REAL lengthX = 50.;
    REAL lengthY = 1.;
    REAL Lmax = 1.;
    
    REAL bulletDepthIni =  0.;
    REAL bulletDepthFin = 100.;
    
    TPZVec<TPZLayerProperties> layerVec(1);
    layerVec[0] = TPZLayerProperties(1.E5,0.25,1.E5,1.E1,0.,100.);
    TPZPlaneFracture plfrac(layerVec, bulletDepthIni, bulletDepthFin, lengthX, lengthY, Lmax);
    
//    REAL Rini = 10.;
//    REAL Rfin = 40;
//    int nRadius = 15;

//    for(int r = 0; r < nRadius; r++)
//    {
//        std::stringstream nm;
//        nm << "circle" << r << ".vtk";
//        
//        TPZVec< std::pair<REAL,REAL> > fractureDots;
//        FillFractureDotsCircle(50., (Rfin-Rini)/(nRadius-1)*r+Rini, fractureDots);
//        TPZGeoMesh * gmesh = plfrac.GetFractureGeoMesh(fractureDots);
//        
//        std::ofstream outC(nm.str().c_str());
//        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outC, true);
//    }
    
    std::ofstream outRefP("RefPatternsUsed.txt");
    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}


int main/*3D*/(int argc, char * const argv[])
{    
    std::cout << "\e";
    TPZTimer readRef("ReadingRefPatterns");
    readRef.start();
    
    //#define writeAgain
    #ifdef writeAgain
    gRefDBase.InitializeRefPatterns();
    #else
    std::ifstream inRefP("RefPatternsUsed.txt");
    gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
    #endif
    
    readRef.stop();
    std::cout << "DeltaT leitura refpatterns = " << readRef.seconds() << " s" << std::endl;
    
    REAL lengthX = 250.;
    REAL lengthY = 50.;
    REAL Lmax = 10.;
    
    REAL bulletDepthIni =  0.;
    REAL bulletDepthFin = 100.;
    
    TPZVec<TPZLayerProperties> layerVec(1);
    layerVec[0] = TPZLayerProperties(1.E5,0.25,1.E5,1.E1,0.,100.);
    TPZPlaneFracture plfrac(layerVec, bulletDepthIni, bulletDepthFin, lengthX, lengthY, Lmax);
    
    TPZVec< std::pair<REAL,REAL> > fractureDots;
    FillFractureDotsExampleCrazy(fractureDots);
    
    
    TPZTimer clockIni2("PartyBegins2");
    clockIni2.start();    
    
    REAL pressureInsideCrack = 10.;
    std::string vtkFile = "fracturePconstant0.vtk";
    plfrac.RunThisFractureGeometry(fractureDots, pressureInsideCrack, vtkFile, true);
    
    clockIni2.stop();
    std::cout << "DeltaT get fracture cmesh = " << clockIni2.seconds() << " s" << std::endl;
    
//    std::ofstream outRefP("RefPatternsUsed.txt");
//    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}


void FillFractureDotsExampleEllipse(TPZVec<std::pair<REAL,REAL> > &fractureDots)
{
    int nnodes = 80;
    
    fractureDots.Resize(nnodes);
    int node;
    REAL shiftZ = -100.;
    
    node = 0;
    
    fractureDots[node] = std::make_pair(0.5,shiftZ + 95.);
    
    node = 1;
    
    fractureDots[node] = std::make_pair(5.,shiftZ + 94.9852);
    
    node = 2;
    
    fractureDots[node] = std::make_pair(10.,shiftZ + 94.9408);
    
    node = 3;
    
    fractureDots[node] = std::make_pair(15.,shiftZ + 94.8667);
    
    node = 4;
    
    fractureDots[node] = std::make_pair(20.,shiftZ + 94.7627);
    
    node = 5;
    
    fractureDots[node] = std::make_pair(25.,shiftZ + 94.6286);
    
    node = 6;
    
    fractureDots[node] = std::make_pair(30.,shiftZ + 94.4643);
    
    node = 7;
    
    fractureDots[node] = std::make_pair(35.,shiftZ + 94.2692);
    
    node = 8;
    
    fractureDots[node] = std::make_pair(40.,shiftZ + 94.0431);
    
    node = 9;
    
    fractureDots[node] = std::make_pair(45.,shiftZ + 93.7854);
    
    node = 10;
    
    fractureDots[node] = std::make_pair(50.,shiftZ + 93.4956);
    
    node = 11;
    
    fractureDots[node] = std::make_pair(55.,shiftZ + 93.173);
    
    node = 12;
    
    fractureDots[node] = std::make_pair(60.,shiftZ + 92.8169);
    
    node = 13;
    
    fractureDots[node] = std::make_pair(65.,shiftZ + 92.4264);
    
    node = 14;
    
    fractureDots[node] = std::make_pair(70.,shiftZ + 92.0006);
    
    node = 15;
    
    fractureDots[node] = std::make_pair(75.,shiftZ + 91.5385);
    
    node = 16;
    
    fractureDots[node] = std::make_pair(80.,shiftZ + 91.0387);
    
    node = 17;
    
    fractureDots[node] = std::make_pair(85.,shiftZ + 90.4998);
    
    node = 18;
    
    fractureDots[node] = std::make_pair(90.,shiftZ + 89.9204);
    
    node = 19;
    
    fractureDots[node] = std::make_pair(95.,shiftZ + 89.2986);
    
    node = 20;
    
    fractureDots[node] = std::make_pair(100.,shiftZ + 88.6323);
    
    node = 21;
    
    fractureDots[node] = std::make_pair(105.,shiftZ + 87.9193);
    
    node = 22;
    
    fractureDots[node] = std::make_pair(110.,shiftZ + 87.1567);
    
    node = 23;
    
    fractureDots[node] = std::make_pair(115.,shiftZ + 86.3416);
    
    node = 24;
    
    fractureDots[node] = std::make_pair(120.,shiftZ + 85.4702);
    
    node = 25;
    
    fractureDots[node] = std::make_pair(125.,shiftZ + 84.5384);
    
    node = 26;
    
    fractureDots[node] = std::make_pair(130.,shiftZ + 83.541);
    
    node = 27;
    
    fractureDots[node] = std::make_pair(135.,shiftZ + 82.4721);
    
    node = 28;
    
    fractureDots[node] = std::make_pair(140.,shiftZ + 81.3243);
    
    node = 29;
    
    fractureDots[node] = std::make_pair(145.,shiftZ + 80.0886);
    
    node = 30;
    
    fractureDots[node] = std::make_pair(150.,shiftZ + 78.7537);
    
    node = 31;
    
    fractureDots[node] = std::make_pair(155.,shiftZ + 77.305);
    
    node = 32;
    
    fractureDots[node] = std::make_pair(160.,shiftZ + 75.7233);
    
    node = 33;
    
    fractureDots[node] = std::make_pair(165.,shiftZ + 73.9822);
    
    node = 34;
    
    fractureDots[node] = std::make_pair(170.,shiftZ + 72.0442);
    
    node = 35;
    
    fractureDots[node] = std::make_pair(175.,shiftZ + 69.8515);
    
    node = 36;
    
    fractureDots[node] = std::make_pair(180.,shiftZ + 67.3077);
    
    node = 37;
    
    fractureDots[node] = std::make_pair(185.,shiftZ + 64.2256);
    
    node = 38;
    
    fractureDots[node] = std::make_pair(190.,shiftZ + 60.125);
    
    node = 39;
    
    fractureDots[node] = std::make_pair(195.,shiftZ + 50.);
    
    node = 40;
    
    fractureDots[node] = std::make_pair(195.,shiftZ + 42.);
    
    node = 41;
    
    fractureDots[node] = std::make_pair(190.,shiftZ + 39.875);
    
    node = 42;
    
    fractureDots[node] = std::make_pair(185.,shiftZ + 35.7744);
    
    node = 43;
    
    fractureDots[node] = std::make_pair(180.,shiftZ + 32.6923);
    
    node = 44;
    
    fractureDots[node] = std::make_pair(175.,shiftZ + 30.1485);
    
    node = 45;
    
    fractureDots[node] = std::make_pair(170.,shiftZ + 27.9558);
    
    node = 46;
    
    fractureDots[node] = std::make_pair(165.,shiftZ + 26.0178);
    
    node = 47;
    
    fractureDots[node] = std::make_pair(160.,shiftZ + 24.2767);
    
    node = 48;
    
    fractureDots[node] = std::make_pair(155.,shiftZ + 22.695);
    
    node = 49;
    
    fractureDots[node] = std::make_pair(150.,shiftZ + 21.2463);
    
    node = 50;
    
    fractureDots[node] = std::make_pair(145.,shiftZ + 19.9114);
    
    node = 51;
    
    fractureDots[node] = std::make_pair(140.,shiftZ + 18.6757);
    
    node = 52;
    
    fractureDots[node] = std::make_pair(135.,shiftZ + 17.5279);
    
    node = 53;
    
    fractureDots[node] = std::make_pair(130.,shiftZ + 16.459);
    
    node = 54;
    
    fractureDots[node] = std::make_pair(125.,shiftZ + 15.4616);
    
    node = 55;
    
    fractureDots[node] = std::make_pair(120.,shiftZ + 14.5298);
    
    node = 56;
    
    fractureDots[node] = std::make_pair(115.,shiftZ + 13.6584);
    
    node = 57;
    
    fractureDots[node] = std::make_pair(110.,shiftZ + 12.8433);
    
    node = 58;
    
    fractureDots[node] = std::make_pair(105.,shiftZ + 12.0807);
    
    node = 59;
    
    fractureDots[node] = std::make_pair(100.,shiftZ + 11.3677);
    
    node = 60;
    
    fractureDots[node] = std::make_pair(95.,shiftZ + 10.7014);
    
    node = 61;
    
    fractureDots[node] = std::make_pair(90.,shiftZ + 10.0796);
    
    node = 62;
    
    fractureDots[node] = std::make_pair(85.,shiftZ + 9.50016);
    
    node = 63;
    
    fractureDots[node] = std::make_pair(80.,shiftZ + 8.96134);
    
    node = 64;
    
    fractureDots[node] = std::make_pair(75.,shiftZ + 8.46154);
    
    node = 65;
    
    fractureDots[node] = std::make_pair(70.,shiftZ + 7.99937);
    
    node = 66;
    
    fractureDots[node] = std::make_pair(65.,shiftZ + 7.57359);
    
    node = 67;
    
    fractureDots[node] = std::make_pair(56.,shiftZ + 7.18313);
    
    node = 68;
    
    fractureDots[node] = std::make_pair(55.,shiftZ - 7.);//+ 6.82703);
    
    node = 69;
    
    fractureDots[node] = std::make_pair(50.,shiftZ + 6.50444);
    
    node = 70;
    
    fractureDots[node] = std::make_pair(45.,shiftZ + 6.21462);
    
    node = 71;
    
    fractureDots[node] = std::make_pair(40.,shiftZ + 5.95692);
    
    node = 72;
    
    fractureDots[node] = std::make_pair(31.,shiftZ + 5.73079);
    
    node = 73;
    
    fractureDots[node] = std::make_pair(29.,shiftZ - 6.);//5.53573);
    
    node = 74;
    
    fractureDots[node] = std::make_pair(25.,shiftZ + 5.37135);
    
    node = 75;
    
    fractureDots[node] = std::make_pair(20.,shiftZ + 5.23731);
    
    node = 76;
    
    fractureDots[node] = std::make_pair(15.,shiftZ + 5.13333);
    
    node = 77;
    
    fractureDots[node] = std::make_pair(10.,shiftZ + 5.05921);
    
    node = 78;
    
    fractureDots[node] = std::make_pair(5.,shiftZ + 5.0148);
    
    node = 79;
    
    fractureDots[node] = std::make_pair(0.5,shiftZ + 5.);
}


void FillFractureDotsExampleCrazy(TPZVec<std::pair<REAL,REAL> > &fractureDots)
{
    int nnodes = 62;
    
    fractureDots.Resize(nnodes);
    int node;
    
    REAL shiftZ = -100.;
    
    ///
    node = 0;
    REAL desloc = 0.;
    
    fractureDots[node] = std::make_pair(5.,shiftZ + 80.);
    
    node = 1;
    
    fractureDots[node] = std::make_pair(9.02368,shiftZ + 82.1856);
    
    node = 2;
    
    fractureDots[node] = std::make_pair(25.6151,shiftZ + 82.9122);
    
    node = 3;
    
    fractureDots[node] = std::make_pair(48.4505,shiftZ + 82.5833);
    
    node = 4;
    
    fractureDots[node] = std::make_pair(74.1186,shiftZ + 81.5259);
    
    node = 5;
    
    fractureDots[node] = std::make_pair(100. + desloc,shiftZ + 80.);//<----------------------
    
    node = 6;
    
    fractureDots[node] = std::make_pair(124.157,shiftZ + 78.2079);
    
    node = 7;
    
    fractureDots[node] = std::make_pair(145.23,shiftZ + 76.3027);
    
    node = 8;
    
    fractureDots[node] = std::make_pair(162.349,shiftZ + 74.3954);
    
    node = 9;
    
    fractureDots[node] = std::make_pair(175.043,shiftZ + 72.5625);
    
    node = 10;
    
    fractureDots[node] = std::make_pair(183.172,shiftZ + 70.8516);
    
    node = 11;
    
    fractureDots[node] = std::make_pair(186.85,shiftZ + 69.2873);
    
    node = 12;
    
    fractureDots[node] = std::make_pair(186.391,shiftZ + 67.8763);
    
    node = 13;
    
    fractureDots[node] = std::make_pair(182.253,shiftZ + 66.6118);
    
    node = 14;
    
    fractureDots[node] = std::make_pair(174.99,shiftZ + 65.4771);
    
    node = 15;
    
    fractureDots[node] = std::make_pair(165.212,shiftZ + 64.4495);
    
    node = 16;
    
    fractureDots[node] = std::make_pair(153.552,shiftZ + 63.5025);
    
    node = 17;
    
    fractureDots[node] = std::make_pair(140.633,shiftZ + 62.609);
    
    node = 18;
    
    fractureDots[node] = std::make_pair(127.05,shiftZ + 61.7426);
    
    node = 19;
    
    fractureDots[node] = std::make_pair(113.346,shiftZ + 60.8796);
    
    node = 20;
    
    fractureDots[node] = std::make_pair(100. + desloc,shiftZ + 60.);//<----------------------
    
    node = 21;
    
    fractureDots[node] = std::make_pair(87.418,shiftZ + 59.0884);
    
    node = 22;
    
    fractureDots[node] = std::make_pair(75.9253,shiftZ + 58.1348);
    
    node = 23;
    
    fractureDots[node] = std::make_pair(65.7644,shiftZ + 57.1342);
    
    node = 24;
    
    fractureDots[node] = std::make_pair(57.0956,shiftZ + 56.0873);
    
    node = 25;
    
    fractureDots[node] = std::make_pair(50.,shiftZ + 55.);
    
    node = 26;
    
    fractureDots[node] = std::make_pair(44.4857,shiftZ + 53.8827);
    
    node = 27;
    
    fractureDots[node] = std::make_pair(40.495,shiftZ + 52.75);
    
    node = 28;
    
    fractureDots[node] = std::make_pair(37.9143,shiftZ + 51.6198);
    
    node = 29;
    
    fractureDots[node] = std::make_pair(36.5847,shiftZ + 50.5124);
    
    node = 30;
    
    fractureDots[node] = std::make_pair(36.3137,shiftZ + 49.4496);
    
    node = 31;
    
    fractureDots[node] = std::make_pair(36.888,shiftZ + 48.4535);
    
    node = 32;
    
    fractureDots[node] = std::make_pair(38.0857,shiftZ + 47.5455);
    
    node = 33;
    
    fractureDots[node] = std::make_pair(39.6892,shiftZ + 46.7455);
    
    node = 34;
    
    fractureDots[node] = std::make_pair(41.4972,shiftZ + 46.0702);
    
    node = 35;
    
    fractureDots[node] = std::make_pair(43.3363,shiftZ + 45.5328);
    
    node = 36;
    
    fractureDots[node] = std::make_pair(45.0707,shiftZ + 45.1417);
    
    node = 37;
    
    fractureDots[node] = std::make_pair(46.6112,shiftZ + 44.8996);
    
    node = 38;
    
    fractureDots[node] = std::make_pair(47.9218,shiftZ + 44.8032);
    
    node = 39;
    
    fractureDots[node] = std::make_pair(49.0243,shiftZ + 44.8424);
    
    node = 40;
    
    fractureDots[node] = std::make_pair(50.,shiftZ + 45.);
    
    node = 41;
    
    fractureDots[node] = std::make_pair(50.9889,shiftZ + 45.2519);
    
    node = 42;
    
    fractureDots[node] = std::make_pair(52.1852,shiftZ + 45.5668);
    
    node = 43;
    
    fractureDots[node] = std::make_pair(53.8292,shiftZ + 45.9072);
    
    node = 44;
    
    fractureDots[node] = std::make_pair(56.1951,shiftZ + 46.2297);
    
    node = 45;
    
    fractureDots[node] = std::make_pair(59.5747,shiftZ + 46.4863);
    
    node = 46;
    
    fractureDots[node] = std::make_pair(64.2563,shiftZ + 46.626);
    
    node = 47;
    
    fractureDots[node] = std::make_pair(70.4977,shiftZ + 46.5964);
    
    node = 48;
    
    fractureDots[node] = std::make_pair(78.4951,shiftZ + 46.3461);
    
    node = 49;
    
    fractureDots[node] = std::make_pair(88.3449,shiftZ + 45.8275);
    
    node = 50;
    
    fractureDots[node] = std::make_pair(100.,shiftZ + 45.);
    
    node = 51;
    
    fractureDots[node] = std::make_pair(113.22,shiftZ + 43.8338);
    
    node = 52;
    
    fractureDots[node] = std::make_pair(127.511,shiftZ + 42.3139);
    
    node = 53;
    
    fractureDots[node] = std::make_pair(142.068,shiftZ + 40.4454);
    
    node = 54;
    
    fractureDots[node] = std::make_pair(155.692,shiftZ + 38.2585);
    
    node = 55;
    
    fractureDots[node] = std::make_pair(166.721,shiftZ + 35.8147);
    
    node = 56;
    
    fractureDots[node] = std::make_pair(172.932,shiftZ + 33.2139);
    
    node = 57;
    
    fractureDots[node] = std::make_pair(171.452,shiftZ + 30.601);
    
    node = 58;
    
    fractureDots[node] = std::make_pair(158.644,shiftZ + 28.175);
    
    node = 59;
    
    fractureDots[node] = std::make_pair(129.998,shiftZ + 26.1969);
    
    node = 60;
    
    fractureDots[node] = std::make_pair(80.,shiftZ + 25.);
    
    node = 61;
    
    fractureDots[node] = std::make_pair(22.,shiftZ + 25.);
}



void FillFractureDotsCircle(REAL center, REAL radius, TPZVec< std::pair<REAL,REAL> > &fractureDots)
{
    REAL Lmax = 2.;
    int nsteps = M_PI * radius / Lmax;
    if(nsteps < 10) nsteps = 10;
    REAL ang = M_PI / nsteps;
    
    int nnodes = nsteps + 1;
    fractureDots.Resize(nnodes);
    
    for(int node = 0; node < nnodes; node++)
    {
        REAL vx = 0.1 + radius*sin(node*ang);
        REAL vz = radius*cos(node*ang) - center;
        fractureDots[node] = std::make_pair(vx , vz);
    }
    std::cout.flush();
}




#define usingRefUnif
//#define usingRefdir
//#define usingQPoints

//#define writeAgain

int main2D(int argc, char * const argv[])
{
#ifdef usingRefdir
    #ifdef writeAgain
        gRefDBase.InitializeRefPatterns();
    #else
        std::ifstream inRefP("RefPatternsUsed.txt");
        gRefDBase.ReadRefPatternDBase("RefPatternsUsed.txt");
    #endif
#endif
    
    REAL lf = 1.;
    REAL ldom = 2.;
    REAL hdom = 2.;
    REAL lmax = 0.25;
    TPZGeoMesh * gmesh = PlaneMesh(lf, ldom, hdom, lmax);
    
    std::ofstream pppt("ppt.txt");
    for(int ppp = 1; ppp <= 100; ppp++)
    {
        REAL pressure = 250.*ppp;
        TPZCompMesh * cmesh = PlaneMesh(gmesh,pressure);
        
        ////Analysis
        TPZAnalysis an(cmesh);
        
        TPZSkylineStructMatrix skylin(cmesh); //caso simetrico
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        
        an.SetStructuralMatrix(skylin);
        an.SetSolver(step);
        an.Run();
        
        ////Post Processing
        TPZManVector<std::string,10> scalnames(0), vecnames(0);

        vecnames.Resize(1);
        vecnames[0] = "displacement";
        
        scalnames.Resize(3);
        scalnames[0] = "SigmaX";
        scalnames[1] = "SigmaY";
        scalnames[2] = "TauXY";
        
        int div = 0;
        std::stringstream postp;
        postp << "Elastic2dP" << (int)(pressure) << ".vtk";
        an.DefineGraphMesh(2,scalnames,vecnames,postp.str());
        an.PostProcess(div,2);
        
        ////2D J-Integral
//        TPZVec<REAL> Origin(3,0.);
//        for(int el = 0; el < gmesh->NElements(); el++)
//        {
//            if(gmesh->ElementVec()[el]->MaterialId() == __1DcrackTipMat)
//            {
//                Origin[0] = gmesh->ElementVec()[el]->Node(0).Coord(0);
//                break;
//            }
//        }
//        
//        TPZVec<REAL> normalDirection(3,0.);
//        normalDirection[2] = 1.;
//        REAL radius = 0.5;  //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//        Path2D * j2Dpath = new Path2D(cmesh, Origin, normalDirection, radius, pressure);
//        
//        JIntegral2D Jpath;
//        Jpath.PushBackPath2D(j2Dpath);
//        TPZVec<REAL> j2Dintegral = Jpath.IntegratePath2D(0);
//        
//        pppt << j2Dintegral[0] << "\n";
    }
    return 0;
}



TPZGeoMesh * PlaneMesh(REAL lf, REAL ldom, REAL hdom, REAL lmax)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    int ndivfrac = int(lf/lmax + 0.5);
    int ndivoutfrac = int((ldom - lf)/lmax + 0.5);
    int ndivh = int(hdom/lmax + 0.5);
    
    int ncols = ndivfrac + ndivoutfrac + 1;
    int nrows = ndivh + 1;
    int nnodes = nrows*ncols;
    
    gmesh->NodeVec().Resize(nnodes);
    
    REAL deltadivfrac = lf/ndivfrac;
    REAL deltadivoutfrac = (ldom-lf)/ndivoutfrac;
    REAL deltandivh = hdom/ndivh;
    
    int nid = 0;
    for(int r = 0; r < nrows; r++)
    {
        for(int c = 0; c < ncols; c++)
        {
            REAL x, y;
            if(c <= ndivfrac)
            {
                x = c*deltadivfrac;
            }
            else
            {
                x = ndivfrac*deltadivfrac + (c-ndivfrac)*deltadivoutfrac;
            }
            y = r*deltandivh;
            
            TPZVec<REAL> coord(3,0.);
            coord[0] = x;
            coord[1] = y;
            gmesh->NodeVec()[r*ncols + c].SetCoord(coord);
            gmesh->NodeVec()[r*ncols + c].SetNodeId(nid);
            nid++;
        }
    }
    
    TPZGeoEl * gel = NULL;
    TPZVec<long> topol(4);
    long indx = 0;
    for(int r = 0; r < nrows-1; r++)
    {
        for(int c = 0; c < ncols-1; c++)
        {
            topol[0] = r*(ncols) + c;
            topol[1] = r*(ncols) + c + 1;
            topol[2] = r*(ncols) + c + 1 + ncols;
            topol[3] = r*(ncols) + c + ncols;
            
            gel = gmesh->CreateGeoElement(EQuadrilateral, topol, globMaterialIdStruct.__3DrockMat_linear, indx);
            gel->SetId(indx);
            indx++;
        }
    }
    
    gmesh->BuildConnectivity();
    
    int nelem = gmesh->NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        
        //south BC
        TPZGeoElSide sideS(gel,4);
        TPZGeoElSide neighS(sideS.Neighbour());
        if(sideS == neighS)
        {
            if(el < ndivfrac)
            {
                gel->CreateBCGeoEl(4, globMaterialIdStruct.__2DfractureMat_inside);
            }
            else
            {
                gel->CreateBCGeoEl(4, globMaterialIdStruct.__2DfractureMat_outside);
            }
        }
        
        //east BC
        TPZGeoElSide sideE(gel,5);
        TPZGeoElSide neighE(sideE.Neighbour());
        if(sideE == neighE)
        {
            gel->CreateBCGeoEl(5, globMaterialIdStruct.__2DrightMat);
        }
        
        //north BC
        TPZGeoElSide sideN(gel,6);
        TPZGeoElSide neighN(sideN.Neighbour());
        if(sideN == neighN)
        {
            gel->CreateBCGeoEl(6, globMaterialIdStruct.__2DfarfieldMat);
        }
        
        //west BC
        TPZGeoElSide sideW(gel,7);
        TPZGeoElSide neighW(sideW.Neighbour());
        if(sideW == neighW)
        {
            gel->CreateBCGeoEl(7, globMaterialIdStruct.__2DleftMat);
        }
    }
    
    topol.Resize(1);
    topol[0] = ndivfrac;
    gel = gmesh->CreateGeoElement(EPoint, topol, globMaterialIdStruct.__1DcrackTipMat, indx);
    
    gmesh->BuildConnectivity();
    
#ifdef usingQPoints
    TPZGeoElSide pt(gel,0);
    TPZGeoElSide ptneigh(pt.Neighbour());
    while(pt != ptneigh)
    {
        if(ptneigh.Element()->HasSubElement() == false)
        {
            int neighSide = ptneigh.Side();
            TPZGeoEl * ptneighEl = TPZChangeEl::ChangeToQuarterPoint(gmesh, ptneigh.Element()->Id(), neighSide);
            ptneigh = ptneighEl->Neighbour(neighSide);
        }
        else
        {
            ptneigh = ptneigh.Neighbour();
        }
    }
#endif
    
#ifdef usingRefUnif
    int nrefUnif = 3;
    for(int ref = 0; ref < nrefUnif; ref++)
    {
        nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
            if(gmesh->ElementVec()[el]->MaterialId() == globMaterialIdStruct.__2DfractureMat_inside)
            {
                TPZVec<TPZGeoEl*> sons;
                gmesh->ElementVec()[el]->Divide(sons);
                continue;
            }
            for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
            {
                TPZGeoElSide gelside(gmesh->ElementVec()[el],s);
                TPZGeoElSide neighside(gelside.Neighbour());
                bool refinedAlready = false;
                while(neighside != gelside)
                {
                    if(neighside.Element()->MaterialId() == globMaterialIdStruct.__2DfractureMat_inside)
                    {
                        TPZVec<TPZGeoEl*> sons;
                        gmesh->ElementVec()[el]->Divide(sons);
                        refinedAlready = true;
                        break;
                    }
                    neighside = neighside.Neighbour();
                }
                if(refinedAlready == true)
                {
                    break;
                }
            }
        }
    }
#endif
    
#ifdef usingRefdir
    std::set<int> matDir;
    //matDir.insert(__2DfractureMat_inside);
    matDir.insert(__1DcrackTipMat);
    int nrefDir = 1;
    for(int ref = 0; ref < nrefDir; ref++)
    {
        nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            if(!gmesh->ElementVec()[el]) continue;
            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
            TPZRefPatternTools::RefineDirectional(gmesh->ElementVec()[el], matDir);
        }
    }
#endif

    
    return gmesh;
}


TPZCompMesh * PlaneMesh(TPZGeoMesh * gmesh, REAL pressureInsideCrack)
{
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDimModel(2);
    cmesh->SetDefaultOrder(2);
    cmesh->SetAllCreateFunctionsContinuous();
    
    STATE young = 0.3e5;
    STATE poisson = 0.25;
    
    int planeStrain = 0;
    //int planeStress = 1;
    int planeWhat = planeStrain;
    
    TPZMaterial * materialLin = new TPZElasticityMaterial(globMaterialIdStruct.__3DrockMat_linear, young, poisson, 0., 0., planeWhat);
    cmesh->InsertMaterialObject(materialLin);
    
    ////BCs
    TPZFMatrix<STATE> k(3,3,0.), f(3,1,0.);
    int newmann = 1, directionalNullDirich = 3;
    
    {
        f(0,0) = 1.;
        TPZBndCond * nullXleft = materialLin->CreateBC(materialLin, globMaterialIdStruct.__2DleftMat, directionalNullDirich, k, f);
        cmesh->InsertMaterialObject(nullXleft);
        TPZBndCond * nullXright = materialLin->CreateBC(materialLin, globMaterialIdStruct.__2DrightMat, directionalNullDirich, k, f);
        cmesh->InsertMaterialObject(nullXright);
        
        f.Zero();
        f(1,0) = 1.;
        TPZBndCond * nullYoutside = materialLin->CreateBC(materialLin, globMaterialIdStruct.__2DfractureMat_outside, directionalNullDirich, k, f);
        cmesh->InsertMaterialObject(nullYoutside);
        
        k.Zero();
//        f(1,0) = traction;
//        TPZBndCond * newmanFarfield = materialLin->CreateBC(materialLin, __2DfarfieldMat, newmann, k, f);
//        cmesh->InsertMaterialObject(newmanFarfield);
        f(1,0) = 1.;
        TPZBndCond * newmanFarfield = materialLin->CreateBC(materialLin, globMaterialIdStruct.__2DfarfieldMat, directionalNullDirich, k, f);
        cmesh->InsertMaterialObject(newmanFarfield);
        
        f(1,0) = pressureInsideCrack;
        TPZBndCond * newmanInside = materialLin->CreateBC(materialLin, globMaterialIdStruct.__2DfractureMat_inside, newmann, k, f);
        cmesh->InsertMaterialObject(newmanInside);
    }
    
    cmesh->AutoBuild();
    
    return cmesh;
}
