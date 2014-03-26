#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;


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
    REAL SigYY0   = -3.5627279293E7  * globStressScale;
    
    REAL Young1   =  4.1368543680E10 * globStressScale;
    REAL Poisson1 =  0.15;
    REAL SigYY1   = -3.4528254983E7  * globStressScale;
    
    REAL Young2   =  4.2747495136E10 * globStressScale;
    REAL Poisson2 =  0.25;
    REAL SigYY2   = -3.732407906E7   * globStressScale;
    
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
    REAL Pe1 = 3.4528254983E7 * globStressScale;
    REAL gradPref1 = 3.4528254983E7 * globStressScale;
    REAL vsp1 = 0.;
    
    REAL Cl2 = 0.;
    REAL Pe2 = 0.;
    REAL gradPref2 = 1.;
    REAL vsp2 = 0.;
    
    layerVec[0] = LayerProperties(Young0, Poisson0, SigYY0, TVDi0, TVDf0, KIc0, Cl0, Pe0, gradPref0, vsp0);
    layerVec[1] = LayerProperties(Young1, Poisson1, SigYY1, TVDi1, TVDf1, KIc1, Cl1, Pe1, gradPref1, vsp1);
    layerVec[2] = LayerProperties(Young2, Poisson2, SigYY2,TVDi2, TVDf2, KIc2, Cl2, Pe2, gradPref2, vsp2);
    
    //Fluid injection data
    REAL QinjWell = 1.*(-0.0533333333333);//m3/s (1.* pois os 80 bpm jah eh no poco e nao 1 wing)
    REAL visc = 200.02E-3 * globStressScale;
    
    //J-Integral data
    REAL Jradius = 0.5;

    //Simulation p-order data
    int porder = 1;
    REAL MaxDispl = 2.0;
    
    bool pressureINdependent = true;//If true, Carter Leakoff Coefficient is pressure independent
    TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax,
                                                                 QinjWell, visc,
                                                                 Jradius,
                                                                 porder,
                                                                 MaxDispl,
                                                                 pressureINdependent);
    
//    globLeakoffStorage.DisableLeakoff();
    plfrac->Run();
    
//    std::ofstream outRefP("RefPatternsUsed.txt");
//    gRefDBase.WriteRefPatternDBase(outRefP);
    
    return 0;
}
