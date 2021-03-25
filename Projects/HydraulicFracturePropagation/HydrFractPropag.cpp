#include <iostream>

#include "TPZRefPatternDataBase.h"
#include "TPZPlaneFractureKernel.h"

#include "TPZTimer.h"

using namespace std;


std::string path = PZ_REFPATTERN_DIR;


int main(int argc, char * const argv[])
{
    std::string desiredPath = PZ_REFPATTERN_DIR;
    desiredPath.erase(desiredPath.length()-18,desiredPath.length());
    desiredPath.append("Projects/HydraulicFracturePropagation/PlaneFracture/RefPatternsUsed.txt");
    gRefDBase.ReadRefPatternDBase(desiredPath);
    
    {
        //Transient data
        REAL Ttot = 2. * 60.; /** em segundos */
        globTimeControl.SetTimeControl(Ttot);
        
        //Geometry data
        REAL lengthX = 70.;
        REAL lengthY = 10.;
        REAL Lmax = 2.0;
        
        //Material data
        TPZVec<LayerProperties> layerVec(3);

        ///////////////////////////////////layer 0
        REAL Young0   =  4.3E10 * globStressScale;
        REAL Poisson0 =  0.25;
        REAL SigYY0   = -35.E6 * globStressScale;
        //
        REAL TVDi0 = 0.;
        REAL TVDf0 = 40.;
        //
        REAL KIc0 = 1.1E6 * globStressScale;
        //
        REAL Cl0 = 0.00025;
        REAL Pe0 = 0.;
        REAL gradPref0 = 0.;
        REAL vsp0 = 0.;

        ///////////////////////////////////layer 1
        REAL Young1   =  4.3E10 * globStressScale;
        REAL Poisson1 =  0.25;
        REAL SigYY1   = -35.E6  * globStressScale;
        //
        REAL TVDi1 = TVDf0;
        REAL TVDf1 = 60.;
        //
        REAL KIc1 = 1.1E6 * globStressScale;
        //
        REAL Cl1 = 0.00025;
        REAL Pe1 = 0.;
        REAL gradPref1 = 0.;
        REAL vsp1 = 0.;
        
        ///////////////////////////////////layer 2
        REAL Young2   =  4.3E10 * globStressScale;
        REAL Poisson2 =  0.25;
        REAL SigYY2   = -35.E6 * globStressScale;
        //
        REAL TVDi2 = TVDf1;
        REAL TVDf2 = 100.;
        //
        REAL KIc2 = 1.1E6 * globStressScale;
        //
        REAL Cl2 = 0.00025;
        REAL Pe2 = 0.;
        REAL gradPref2 = 0.;
        REAL vsp2 = 0.;
        
        ///////////////////////////////////canhoneado
        REAL bulletTVDIni = 40.;
        REAL bulletTVDFin = 60.;
        
        
        bool pressureIndependent = true;//If true, Carter Leakoff Coefficient is pressure independent
        
        layerVec[0] = LayerProperties(Young0, Poisson0, SigYY0, TVDi0, TVDf0, KIc0, Cl0, Pe0, gradPref0, vsp0, pressureIndependent);
        layerVec[1] = LayerProperties(Young1, Poisson1, SigYY1, TVDi1, TVDf1, KIc1, Cl1, Pe1, gradPref1, vsp1, pressureIndependent);
        layerVec[2] = LayerProperties(Young2, Poisson2, SigYY2, TVDi2, TVDf2, KIc2, Cl2, Pe2, gradPref2, vsp2, pressureIndependent);
        
        //Fluid injection data
        REAL QinjWell = -0.05;//m3/s
        
        REAL visc = 200.0E-3 * globStressScale;
        
        //J-Integral data
        REAL Jradius = 1.0;
        
        REAL MaxDispl = 2.0;
        
        bool just1stripe = true;
        bool layerStripesToo = false;
        TPZPlaneFractureKernel * plfrac = new TPZPlaneFractureKernel(layerVec, bulletTVDIni, bulletTVDFin, lengthX, lengthY, Lmax,
                                                                     QinjWell, visc,
                                                                     Jradius,
                                                                     MaxDispl,
                                                                     pressureIndependent,
                                                                     just1stripe, layerStripesToo);
        
        //globLeakoffStorage.DisableLeakoff();
        plfrac->Run();
    }
    
    std::ofstream outf(desiredPath);
    gRefDBase.WriteRefPatternDBase(outf);
    return 0;
}
