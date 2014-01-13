#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;



#define PDoreTest1
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
    
#ifdef PDoreTest1
    
    //Transient data
    REAL Ttot = 25.159243081728423 * 60.; /** em segundos */
    REAL maxDeltaT = Ttot/30.; /** em segundos */
    int nTimes = 1; /** quantidade de divisao do maxDeltaT para definir minDeltaT (minDeltaT = maxDeltaT/nTimes) */
    globTimeControl.SetTimeControl(Ttot, maxDeltaT, nTimes);
    
    //Geometry data
    REAL lengthX = 70.;
    REAL lengthY = 10.;
    REAL Lmax = 2.0;
    
    REAL bulletTVDIni = 2100.;
    REAL bulletTVDFin = 2120.;
    int nstripes = 1;
    
    //Material data
    TPZVec<TPZLayerProperties> layerVec(3);
    
    REAL Young0 = 4.2747495136E10;
    REAL Poisson0 = 0.25;
    REAL SigMax0  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin0  = -35627279.293;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf0 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young1 = 4.1368543680E10;
    REAL Poisson1 = 0.15;
    REAL SigMax1  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin1  = -34528254.983;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf1 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL Young2 = 4.2747495136E10;
    REAL Poisson2 = 0.25;
    REAL SigMax2  = 0.;      //<<<<<<<============= PRE-STRESS XX
    REAL SigMin2  = -37324079.06;   //<<<<<<<============= PRE-STRESS YY
    REAL SigConf2 = 0.;      //<<<<<<<============= PRE-STRESS ZZ
    
    REAL TVDi0 = 2070.;
    REAL TVDf0 = 2100.;
    REAL TVDi1 = TVDi0;
    REAL TVDf1 = 2120.;
    REAL TVDi2 = TVDi1;
    REAL TVDf2 = 2160.;
    
    REAL KIc0 = 1.09884E6;
    REAL KIc1 = 1.09884E6;
    REAL KIc2 = 1.09884E6;
    
    REAL Cl = 0.00019674755398733676;

    REAL Pe = 10.E6;//Sempre positivo
    REAL gradPref = 10.*3447378.64;//n*500 psi
    REAL vsp = 0.;
    
    layerVec[0] = TPZLayerProperties(Young0, Poisson0, SigMax0, SigMin0, SigConf0, TVDi0, TVDf0, KIc0, Cl, Pe, gradPref, vsp);
    layerVec[1] = TPZLayerProperties(Young1, Poisson1, SigMax1, SigMin1, SigConf1, TVDi1, TVDf1, KIc1, Cl, Pe, gradPref, vsp);
    layerVec[2] = TPZLayerProperties(Young2, Poisson2, SigMax2, SigMin2, SigConf2, TVDi2, TVDf2, KIc2, Cl, Pe, gradPref, vsp);
    
    //Fluid injection data
    REAL QinjWell = /*2.*(-0.0063084833733305975);*/2.*(-0.052995764976);//m3/s
    REAL visc = 200.02E-3;//N.s/m2
    
    //J-Integral data
    REAL Jradius = 1.5;
    
#else
    
#endif

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
