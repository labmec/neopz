#include <iostream>

#include "pzgmesh.h"
#include "pzgeopoint.h"
#include "TPZRefPatternDataBase.h"
#include "tpzgeoelrefpattern.h"
#include "TPZPlaneFractureMesh.h"

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
    int nstripes = 3;
    
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
    TPZPlaneFractureMesh plfrac(layerVec, bulletDepthTVDIni, bulletDepthTVDFin, lengthX, lengthY, Lmax, nstripes);
    
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
    int nstripes = 3;
    
    TPZVec<TPZLayerProperties> layerVec(1);
    layerVec[0] = TPZLayerProperties(1.E5,0.25,1.E5,1.E1,0.,100.);
    TPZPlaneFractureMesh plfrac(layerVec, bulletDepthIni, bulletDepthFin, lengthX, lengthY, Lmax, nstripes);
    
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
    int nstripes = 5;
    
    TPZVec<TPZLayerProperties> layerVec(3);
    layerVec[0] = TPZLayerProperties(1.E5,0.25,0.,0.,0.,33.333);
    layerVec[1] = TPZLayerProperties(1.E5,0.25,0.,0.,33.333,66.666);
    layerVec[2] = TPZLayerProperties(1.E5,0.25,0.,0.,66.666,100.0);
    TPZPlaneFractureMesh plfrac(layerVec, bulletDepthIni, bulletDepthFin, lengthX, lengthY, Lmax, nstripes);
    
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
