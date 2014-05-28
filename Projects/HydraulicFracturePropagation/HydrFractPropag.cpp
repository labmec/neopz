#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;

//P. Dore
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
    REAL Ttot = 25. * 60.; /** em segundos */
    globTimeControl.SetTimeControl(Ttot);
    
    //Geometry data
    REAL lengthX = 80.;
    REAL lengthY = 60.;
    REAL Lmax = 2.0;
    
    //Material data
    TPZVec<LayerProperties> layerVec(3);
    
    REAL Young0   =  4.2747495136E10 * globStressScale;
    REAL Poisson0 =  0.25;
    REAL SigYY0   = -35.627279293E6  * globStressScale;
    
    REAL Young1   =  4.1368543680E10 * globStressScale;
    REAL Poisson1 =  0.15;
    REAL SigYY1   = -34.528254983E6  * globStressScale;
    
    REAL Young2   =  4.2747495136E10 * globStressScale;
    REAL Poisson2 =  0.25;
    REAL SigYY2   = -37.32407906E6   * globStressScale;
    
    REAL TVDi0 = 2070.;
    REAL TVDf0 = 2100.;
    REAL TVDi1 = TVDf0;
    REAL TVDf1 = 2120.;
    REAL TVDi2 = TVDf1;
    REAL TVDf2 = 2160.;
    
    REAL bulletTVDIni = 2100.;
    REAL bulletTVDFin = 2120.;
    
    //KIc
    REAL KIc0 = 1.09884E6 * globStressScale;
    REAL KIc1 = 1.09884E6 * globStressScale;
    REAL KIc2 = 1.09884E6 * globStressScale;
    
    //Leakoff
    REAL Cl0 = 0.;
    REAL Pe0 = 0.;
    REAL gradPref0 = 1.;
    REAL vsp0 = 0.;
    
    REAL Cl1 = 0.00019674755398733676;
    REAL Pe1 = 34.528254983E6 * globStressScale;
    REAL gradPref1 = 34.528254983E6 * globStressScale;
    REAL vsp1 = 0.;
    
    REAL Cl2 = 0.;
    REAL Pe2 = 0.;
    REAL gradPref2 = 1.;
    REAL vsp2 = 0.;
    
    layerVec[0] = LayerProperties(Young0, Poisson0, SigYY0, TVDi0, TVDf0, KIc0, Cl0, Pe0, gradPref0, vsp0);
    layerVec[1] = LayerProperties(Young1, Poisson1, SigYY1, TVDi1, TVDf1, KIc1, Cl1, Pe1, gradPref1, vsp1);
    layerVec[2] = LayerProperties(Young2, Poisson2, SigYY2, TVDi2, TVDf2, KIc2, Cl2, Pe2, gradPref2, vsp2);
    
    //Fluid injection data
    REAL QinjWell = 1.*(-0.0533333333333);//m3/s
    
    REAL visc = 200.02E-3 * globStressScale;
    
    //J-Integral data
    REAL Jradius = 1.0;

    REAL MaxDispl = 2.0;
    
    bool pressureINdependent = true;//If true, Carter Leakoff Coefficient is pressure independent
    TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax,
                                                                 QinjWell, visc,
                                                                 Jradius,
                                                                 MaxDispl,
                                                                 pressureINdependent);
    
//    globLeakoffStorage.DisableLeakoff();
    plfrac->Run();
    
//    std::ofstream outRefP("RefPatternsUsed.txt");
//    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}

//Abou Sayed
int mainSayed(int argc, char * const argv[])
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
    REAL Ttot = 200. * 60.; /** em segundos */ //ok
    globTimeControl.SetTimeControl(Ttot);
    
    //Geometry data
    REAL lengthX = 80.;//ok
    REAL lengthY = 60.;//ok
    REAL Lmax = 2.0;//ok
    
    //Material data
    TPZVec<LayerProperties> layerVec(3);
    
    REAL Young0   =  4.481592E10 * globStressScale;//ok
    REAL Poisson0 =  0.30;//ok
    REAL SigYY0   = -4.929751255E7  * globStressScale;//ok
    
    REAL Young1   =  5.860543E10 * globStressScale;//ok
    REAL Poisson1 =  0.21;//ok
    REAL SigYY1   = -3.930011490E7  * globStressScale;//ok
    
    REAL Young2   =  3.792116E10 * globStressScale;//ok
    REAL Poisson2 =  0.29;//ok
    REAL SigYY2   = -5.067646395E7   * globStressScale;//ok
    
    REAL TVDi0 = 2740.16099;//ok
    REAL TVDf0 = 2795.02517;//ok
    REAL TVDi1 = TVDf0;
    REAL TVDf1 = 2846.84134;//ok
    REAL TVDi2 = TVDf1;
    REAL TVDf2 = 2941.32965;//ok
    
    REAL bulletTVDIni = 2795.02517;//ok
    REAL bulletTVDFin = 2846.84134;//ok
    
    //KIc
    REAL KIc0 = 2.19769E6 * globStressScale;//ok
    REAL KIc1 = 2.19769E6 * globStressScale;//ok
    REAL KIc2 = 2.19769E6 * globStressScale;//ok
    
    //Leakoff
    REAL Cl0 = 0.000590245;//ok
    REAL Pe0 = 2.482112520E7 * globStressScale;//ok
    REAL gradPref0 = 689475700E6;//ok
    REAL vsp0 = 0.;
    
    REAL Cl1 = 0.000590245;//ok
    REAL Pe1 = 2.482112520E7 * globStressScale;//ok
    REAL gradPref1 = 689475700E6;//ok
    REAL vsp1 = 0.;
    
    REAL Cl2 = 0.000590245;//ok
    REAL Pe2 = 2.482112520E7 * globStressScale;//ok
    REAL gradPref2 = 689475700E6;//ok
    REAL vsp2 = 0.;
    
    layerVec[0] = LayerProperties(Young0, Poisson0, SigYY0, TVDi0, TVDf0, KIc0, Cl0, Pe0, gradPref0, vsp0);
    layerVec[1] = LayerProperties(Young1, Poisson1, SigYY1, TVDi1, TVDf1, KIc1, Cl1, Pe1, gradPref1, vsp1);
    layerVec[2] = LayerProperties(Young2, Poisson2, SigYY2, TVDi2, TVDf2, KIc2, Cl2, Pe2, gradPref2, vsp2);
    
    //Fluid injection data
    REAL QinjWell = 1.*(-0.1324894);//m3/s //ok
    
    REAL visc = 200.0E-3 * globStressScale;//ok
    
    //J-Integral data
    REAL Jradius = 1.0;
    
    REAL MaxDispl = 1.0;
    
    bool pressureINdependent = true;//If true, Carter Leakoff Coefficient is pressure independent
    TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax,
                                                                 QinjWell, visc,
                                                                 Jradius,
                                                                 MaxDispl,
                                                                 pressureINdependent);
    
    //    globLeakoffStorage.DisableLeakoff();
    plfrac->Run();
    
    //    std::ofstream outRefP("RefPatternsUsed.txt");
    //    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}