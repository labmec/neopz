#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;

const int matPoint = -3;

//** just for visualize given dots in vtk */

void FillFractureDotsExampleEllipse(TPZVec< std::pair<REAL,REAL> > &fractureDots);
//---------------------------------------------------------------------------------------------------------------------------------


int main(int argc, char * const argv[])
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
    
    //Transient data
    REAL Ttot = 50.;//3600.; /** em segundos */
    REAL maxDeltaT = 10.;//600.; /** em segundos */
    int nTimes = 1; /** quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes) */
    globTimeControl.SetTimeControl(Ttot, maxDeltaT, nTimes);
    
    //Geometry data
    REAL lengthX = 50.;
    REAL lengthY = 8.;
    REAL Lmax = 2.;
    
    REAL bulletTVDIni = 40.;
    REAL bulletTVDFin = 60.;
    int nstripes = 1;

    //Material data
    TPZVec<TPZLayerProperties> layerVec(3);
    
    REAL Young = 1.E5;
    REAL Poisson = 0.25;
    REAL SigMax  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin  = -100.;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL TVDi0 = 0.;
    REAL TVDf0 = 30.;
    REAL TVDi1 = TVDi0;
    REAL TVDf1 = 70.;
    REAL TVDi2 = TVDi1;
    REAL TVDf2 = 100.;
    
    REAL KIc = 300.;
    
    REAL Cl = 1.E-4;
    REAL Pe = 100.;//Sempre positivo
    REAL gradPref = 100.;
    REAL vsp = 1.E-8;
    
    layerVec[0] = TPZLayerProperties(Young, Poisson, SigMax, SigMin, SigConf, TVDi0, TVDf0, KIc, Cl, Pe, gradPref, vsp);
    layerVec[1] = TPZLayerProperties(Young, Poisson, SigMax, SigMin, SigConf, TVDi1, TVDf1, KIc, Cl, Pe, gradPref, vsp);
    layerVec[2] = TPZLayerProperties(Young, Poisson, SigMax, SigMin, SigConf, TVDi2, TVDf2, KIc, Cl, Pe, gradPref, vsp);

    //Fluid injection data
    REAL QinjWell = -2.;//m3/s
    REAL visc = 0.001E-6;
    
    //J-Integral data
    REAL Jradius = 2.0;

    //Simulation p-order data
    int porder = 1;
    
    TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax, nstripes,
                                                                 QinjWell, visc,
                                                                 Jradius,
                                                                 porder);

    plfrac->Run();
    
//    std::ofstream outRefP("RefPatternsUsed.txt");
//    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}


void FillFractureDotsExampleEllipse(TPZVec<std::pair<REAL,REAL> > &fractureDots)
{
    int nnodes = 80;
    
    fractureDots.Resize(nnodes);
    int node;
    REAL shiftZ = -110.;
    
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
    
    fractureDots[node] = std::make_pair(55.,shiftZ + 6.82703);
    
    node = 69;
    
    fractureDots[node] = std::make_pair(50.,shiftZ + 6.50444);
    
    node = 70;
    
    fractureDots[node] = std::make_pair(45.,shiftZ + 6.21462);
    
    node = 71;
    
    fractureDots[node] = std::make_pair(40.,shiftZ + 5.95692);
    
    node = 72;
    
    fractureDots[node] = std::make_pair(31.,shiftZ + 5.73079);
    
    node = 73;
    
    fractureDots[node] = std::make_pair(29.,shiftZ + 5.53573);
    
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

