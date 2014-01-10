#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;



#define PDoreTest
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
    
#ifdef PDoreTest
    
    //Transient data
    REAL Ttot = 25.159243081728423 * 60.; /** em segundos */
    REAL maxDeltaT = Ttot/30.; /** em segundos */
    int nTimes = 1; /** quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes) */
    globTimeControl.SetTimeControl(Ttot, maxDeltaT, nTimes);
    
    //Geometry data
    REAL lengthX = 50.;
    REAL lengthY = 8.;
    REAL Lmax = 2.;
    
    REAL bulletTVDIni = 2100.;
    REAL bulletTVDFin = 2120.;
    int nstripes = 1;
    
    //Material data
    TPZVec<TPZLayerProperties> layerVec(3);
    
    REAL Young0 = 34473785900.;
    REAL Poisson0 = 0.28;
    REAL SigMax0  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin0  = -31351839.9;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf0 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young1 = 20684271600.;
    REAL Poisson1 = 0.25;
    REAL SigMax1  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin1  = -30212136.5;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf1 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young2 = 34473785900.;
    REAL Poisson2 = 0.29;
    REAL SigMax2  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin2  = -32845244.3;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf2 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL TVDi0 = 2080.;
    REAL TVDf0 = 2100.;
    REAL TVDi1 = TVDi0;
    REAL TVDf1 = 2120.;
    REAL TVDi2 = TVDi1;
    REAL TVDf2 = 2140.;
    
    REAL KIc0 = 5.83553E6;
    REAL KIc1 = 2.26115E6;
    REAL KIc2 = 7.04121E6;
    
    REAL maxDisplacement = 5.;
    
    REAL Cl = 0.00019674755398733676; //<<<<<<<<<< ??????
    REAL Pe = 3447378.64/5.*3.;//Sempre positivo
    REAL gradPref = 2.*3447378.64;//n*500 psi
    REAL vsp = 0.;
    
    layerVec[0] = TPZLayerProperties(Young0, Poisson0, SigMax0, SigMin0, SigConf0, TVDi0, TVDf0, 1.6*KIc0, Cl, Pe, gradPref, vsp);
    layerVec[1] = TPZLayerProperties(Young1, Poisson1, SigMax1, SigMin1, SigConf1, TVDi1, TVDf1, 1.4*KIc1, Cl, Pe, gradPref, vsp);
    layerVec[2] = TPZLayerProperties(Young2, Poisson2, SigMax2, SigMin2, SigConf2, TVDi2, TVDf2, 1.4*KIc2, Cl, Pe, gradPref, vsp);
    
    //Fluid injection data
    REAL QinjWell = /*2.*(-0.0067920697287366);// <-- Q tirando leakoff || --> Qtot = */2.*(-0.052995764976);//m3/s
    REAL visc = 200.02E-3;//N.s/m2
    
    //J-Integral data
    REAL Jradius = 1.0;
    
#else
    
    //Transient data
    REAL Ttot = 50.; /** em segundos */
    REAL maxDeltaT = 10.; /** em segundos */
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
    REAL SigMax  =    0.;   //<<<<<<<============= PRE-STRESS XX
    REAL SigMin  = -100.;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf =    0.;   //<<<<<<<============= PRE-STRESS ZZ
    
    REAL TVDi0 = 0.;
    REAL TVDf0 = 40.;
    REAL TVDi1 = TVDi0;
    REAL TVDf1 = 60.;
    REAL TVDi2 = TVDi1;
    REAL TVDf2 = 100.;
    
    REAL KIc = 190.;
    
    REAL maxDisplacement = 5.;
    
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
    REAL Jradius = 1.0;
#endif

    //Simulation p-order data
    int porder = 1;
    
    TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax, nstripes,
                                                                 QinjWell, visc,
                                                                 Jradius,
                                                                 maxDisplacement,
                                                                 porder);

    plfrac->Run();
    
//    std::ofstream outRefP("RefPatternsUsed.txt");
//    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}
