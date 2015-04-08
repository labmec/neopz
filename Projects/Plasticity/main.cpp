
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#include "TPZPlasticityTest.h"
#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZTensor.h"
#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "tpzgeoelrefpattern.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZTensor.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "GeoMeshClass.h"
#include "pzelastoplastic2D.h"
#include <pzmathyperelastic.h>
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"
#include "clock_timer.h"
#include "TPBrAcidFunc.h"

using namespace pzshape; // needed for TPZShapeCube and related classes



#include "pzlog.h"
//#include "tpztimer.h"
#include "TPZTimer.h"
#include "WellBoreAnalysis.h"
#include "pzbfilestream.h"
#include "TPZProjectEllipse.h"
#include "arglib.h"
#include "run_stats_table.h"

#define MACOS
#ifdef MACOS

#include <iostream>
#include <math.h>
#include <signal.h>
#include <fenv.h>
#include <xmmintrin.h>

#define ENABLE_FPO_EXCEPTIONS _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);


#define DECLARE_FPO_HANDLER_FUNC void InvalidFPOHandler(int signo) {\
switch(signo) {\
case SIGFPE: std::cout << "ERROR : Invalid Arithmetic operation." << std::endl; break;\
}\
exit(signo);\
}

#define ATTACH_FPO_SIGNAL struct sigaction act = {};\
act.sa_handler = InvalidFPOHandler;\
sigaction(SIGFPE, &act, NULL);


DECLARE_FPO_HANDLER_FUNC;
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.main"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse(Logger::getLogger("LogEllipse"));
#endif

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
using namespace tbb;
// If you have issues with: dyld: Library not loaded: libtbb.dylib
// try setting the LD path. Ex:
//   export DYLD_FALLBACK_LIBRARY_PATH=/Users/borin/Desktop/neopz/tbb40_297oss/lib/
#endif


RunStatsTable plast_tot("-tpz_plast_tot", "Raw data table statistics for the main execution.");
clarg::argInt NumberOfThreads("-nt", "Number of threads for WellBoreAnalysis", 1);


void Config1()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    // tottotototo
    TPZWellBoreAnalysis::TConfig::gNumThreads = 20;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config1.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    int Startfrom=1;
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        well.SetInnerLiner(innerradius-0.02);
        
        well.SetLinerAxialStress(-WellPressure);
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        well.GetCurrentConfig()->fHasCompletion = 1;
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->fHasCompletion = 0;
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        
        
        well.PostProcess(0);
        
        TPZBFileStream save;
        save.OpenWrite("Config1-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        std::stringstream out;
        out << "First pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        
        {
            std::ofstream out("gmesh_good.txt");
            well.GetCurrentConfig()->fGMesh.Print(out);
            std::ofstream co("cmesh_good.txt");
            well.GetCurrentConfig()->fCMesh.Print(co);
        }
        
        
        int substepsloc  = 1;
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            int substepsloc = 1;
            well.ExecuteSimulation(substepsloc,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config1-1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying second pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        int substepsloc = 1;
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            int substepsloc = 1;
            well.ExecuteSimulation(substepsloc,out);
            well.PostProcess(0);
        }
        well.AppendExecutionLog(out);
        TPZBFileStream save;
        save.OpenWrite("Config1-2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying third pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        int substepsloc =1;
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(substepsloc,out);
            well.PostProcess(0);
        }
        well.AppendExecutionLog(out);
        TPZBFileStream save;
        save.OpenWrite("Config1-3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying fourth pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        int substepsloc =1;
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(substepsloc,out);
            well.PostProcess(0);
        }
        well.AppendExecutionLog(out);
        TPZBFileStream save;
        save.OpenWrite("Config1-4.bin");
        well.Write(save);
        
        
    }
    well.PrintExecutionLog(std::cout);
}


void Config2()
{
    //EHorizontalWellalongH
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config2.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    const int nsubsteps = 5;
    int Startfrom=0;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EHorizontalWellalongH;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config2-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 1st pass\n";
        
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        std::stringstream out;
        out << "Applying 2nd pass\n";
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 3nd pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 4th pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-4.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==5)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-4.bin");
        well.Read(read);
    }
    
    if (Startfrom <=5)
    {
        
        std::cout << "\n ------- 5 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 5th pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config2-5.bin");
        well.Write(save);
        
        
    }
    
}


void Config3()
{
    //EHorizontalWellalongh
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config3.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    const int nsubsteps = 5;
    int Startfrom=0;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EHorizontalWellalongh;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config3-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        std::stringstream out;
        out << "Applying 1st pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        std::stringstream out;
        out << "Applying 2nd pass\n";
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 3rd pass\n";
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        std::stringstream out;
        out << "Applying 4th pass\n";
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config3-4.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==5)
    {
        TPZBFileStream read;
        read.OpenRead("Config3-4.bin");
        well.Read(read);
    }
    
    if (Startfrom <=5)
    {
        
        std::stringstream out;
        out << "Applying 5th pass\n";
        std::cout << "\n ------- 5 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b,out);
            well.ExecuteSimulation(nsubsteps,out);
            well.PostProcess(0);
        }
        
        TPZBFileStream save;
        save.OpenWrite("Config3-5.bin");
        well.Write(save);
        
        
    }
    
}

void Config4()
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config4.vtk";
    const int nsubsteps = 5;
    well.SetVtkOutPutName(output);
    EPlasticModel Emodel = EElastic;
    if (Emodel == EMohrCoulomb)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 13.;
        REAL Phi = 0.52;
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
    }
    else if (Emodel == ESandler)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
    }
    else if (Emodel == EElastic){ // Mohr-Coulomb with a VERY far way surface
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 1.e8; // Very very big
        REAL Phi = 1.5533430342749532; // 89 degrees
        
#ifdef PlasticPQP
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
        well.GetCurrentConfig()->fModel = EElastic;
#else
        well.SetElasticParameters(poisson, elast);
        well.GetCurrentConfig()->fModel = EElastic;
#endif
    }
    
    
    int Startfrom=0;
    if (Startfrom == 0)
    {
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        std::cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        //well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("wellbore0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::stringstream out;
        out << "Applying 1st pass\n";
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        well.PRefineElementAbove(0.000001, 2,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        
        cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.SetFluidModel(EPenetrating);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        cout << "Penetrating fluid Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        
        well.EvolveBothPressures(4, WellPressure*1.2,reservoirPressure);
        cout << "Higher well pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.EvolveBothPressures(4,WellPressure,reservoirPressure*0.8);
        cout << "Lower reservoir pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        TPZBFileStream save;
        save.OpenWrite("wellbore1.bin");
        well.Write(save);
    }
    
}

void Config5()
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config5.vtk";
    const int nsubsteps = 5;
    well.SetVtkOutPutName(output);
    EPlasticModel Emodel = ESandler;
    if (Emodel == EMohrCoulomb)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 13.;
        REAL Phi = 0.52;
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
    }
    else if (Emodel == ESandler)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
    }
    else if (Emodel == EElastic){ // Mohr-Coulomb with a VERY far way surface
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 1.e8; // Very very big
        REAL Phi = 1.5533430342749532; // 89 degrees
        
#ifdef PlasticPQP
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
        well.GetCurrentConfig()->fModel = EElastic;
#else
        well.SetElasticParameters(poisson, elast);
        well.GetCurrentConfig()->fModel = EElastic;
#endif
    }
    
    
    int Startfrom=0;
    if (Startfrom == 0)
    {
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        std::cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        //well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("wellbore0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::stringstream out;
        out << "Applying 1st pass\n";
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        well.PRefineElementAbove(0.0001, 2,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        
        cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.SetFluidModel(EPenetrating);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        cout << "Penetrating fluid Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        
        //         well.EvolveBothPressures(4, WellPressure*1.2,reservoirPressure);
        //         cout << "Higher well pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        //
        //         well.EvolveBothPressures(4,WellPressure,reservoirPressure*0.8);
        //         cout << "Lower reservoir pressure Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        //
        TPZBFileStream save;
        save.OpenWrite("wellbore1.bin");
        well.Write(save);
    }
    
    //QUEBRA AKI
    {
        TPZBFileStream read;
        read.OpenRead("wellbore1.bin");
        well.Read(read);
    }
    
}


//Reprodução do bug do erick para Philippe
void Config6()
{
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
    STATE biotcoef = 0.6666666;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 70.0; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 100.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config6.vtk";
    const int nsubsteps = 1;
    well.SetVtkOutPutName(output);
    REAL poisson = 0.2;
    REAL elast = 30000.;
    REAL cohesion = 1.e8; // Very very big
    REAL Phi = 1.5533430342749532; // 89 degrees
    
    well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
    well.GetCurrentConfig()->fModel = EElastic;
    
    int Startfrom=0;
    
    if(Startfrom == 0){
        int porder = 1;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        std::cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        well.PostProcess(0);
        
    }
    
    if(Startfrom ==0){
        
        int nsteps = 1;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        TPZWellBoreAnalysis::TConfig::gNumThreads = 0;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("wellbore0.bin");
        well.Write(save);
        
    }
    
}

void Config7()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPBrAcidFunc acidfunc;
    acidfunc.StandardParameters();
    acidfunc.CalculaDerivados();
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    TPZVec<REAL> x(2,0.);
    for (REAL r=innerradius; r <= outerradius; r += (outerradius-innerradius)/15.) {
        x[0] = r;
        TPZManVector<REAL,3> func(1);
        acidfunc.Execute(x, func);
        std::cout << "r = " << r << " Elast " << func[0] << std::endl;
    }
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config7.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    int Startfrom=1;
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        std::stringstream out;
        out << "Applying 0th pass\n";
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config7-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config7-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 1st pass\n";
        
        
        well.GetCurrentConfig()->fAcidParameters.StandardParameters();
        well.GetCurrentConfig()->ActivateAcidification();
        
        well.ExecuteSimulation(1,out);
        
        well.PostProcess(1);
        
        TPZBFileStream save;
        save.OpenWrite("Config7-1.bin");
        well.Write(save);
        
        
    }
}

void Config8()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    //    Sh=-57.2;
    //    SH=-57.2;
    //    SV=-57.2;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 57.2; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config8.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    int Startfrom=0;
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        
        well.SetElasticParameters(poisson, elast);
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        {
            std::stringstream sout;
            well.PrintInitialConfiguration(sout);
            well.AppendExecutionLog(sout);
        }
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config8-0.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config8-0.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        std::stringstream out;
        out << "Applying 1st pass\n";
        
        well.PRefineElementAbove(sqj2_refine, 3,out);
        well.DivideElementsAbove(sqj2_refine,out);
        well.ExecuteSimulation(nsubsteps,out);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Config1-1.bin");
        well.Write(save);
        
        well.AppendExecutionLog(out);
    }
    
    well.PrintExecutionLog(std::cout);
    
}

void Config9()
{
    //EVertical
    //ENonPenetrating
    InitializePZLOG();
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(ETriangle);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    std::cout << std::setprecision(15);
    
    TPZWellBoreAnalysis well;
    
#ifndef USING_TBB
    if (NumberOfThreads.get_value()>=0) {
        TPZWellBoreAnalysis::TConfig::gNumThreads=NumberOfThreads.get_value();
    }
#else
    int number_tbb=NumberOfThreads.get_value();
    if(number_tbb<=0)number_tbb=1;
    task_scheduler_init init(number_tbb);
#endif
    TPZWellBoreAnalysis::TConfig::gNumThreads = 0;
    
    STATE biotcoef = 0.659;
    well.SetBiotCoefficient(biotcoef);
    
    REAL reservoirPressure=57.2;
    well.SetReservoirPressure(reservoirPressure);
    
    
    REAL SH,Sh,SV;
    Sh=-83.5;
    SH=-99.8;
    SV=-85.9;
    TPZManVector<STATE,3> confinementTotal(3,0.);
    confinementTotal[0] = Sh;
    confinementTotal[1] = SH;
    confinementTotal[2] = SV;
    REAL WellPressure = 66.1; //66.6 61.1 57.2
    well.SetConfinementTotalStresses(confinementTotal, WellPressure);
    
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    
    std::string output = "Config-9.vtk";
    well.SetVtkOutPutName(output);
    
    
    REAL sqj2_refine=0.0001;
    
    int Startfrom=0;
    
    const int nsubsteps = 5;
    if (Startfrom == 0)
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        
        bool modelMC =false;
        
        if (modelMC)
        {
            REAL cohesion = 13.;
            REAL Phi = 0.52;
            well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
            
        }
        else
        {
            well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
            
            
        }
        
        int porder = 2;
        int nrad=20;
        int ncircle = 40;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        
        well.SetInnerLiner(innerradius-0.02);
        
        well.SetLinerAxialStress(-WellPressure);
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        well.GetCurrentConfig()->fHasCompletion = 0;
        
        well.GetCurrentConfig()->CreateGeometricMesh();
        
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        
        //Garantir pos processamento da pressao de poros (criacao do biotforce)
        well.GetCurrentConfig()->SetWellPressure();
        
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        
        
        well.PostProcess(0);
        
        TPZBFileStream save;
        save.OpenWrite("Config1-9.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Config1-9.bin");
        well.Read(read);
    }
    
    
    if (Startfrom <=1)
    {
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        std::stringstream out;
        out << "First pass\n";
        
        well.GetCurrentConfig()->fHasCompletion = 1;
        well.GetCurrentConfig()->fPostprocess.SetCompMesh(0);
        
        TPZElasticResponse ER;
        ER.SetUp(210000., 0.27);
        well.GetCurrentConfig()->fCompletionElastic.SetElasticResponse(ER);
        well.GetCurrentConfig()->AddLinerComputationalElements(well.GetCurrentConfig()->fInitialPOrderLiner);
        
        well.GetCurrentConfig()->SetWellPressure();
        
        int substepsloc  = 1;
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        
        TPZBFileStream save;
        save.OpenWrite("Config2-9.bin");
        well.Write(save);
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("Config2-9.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        std::stringstream out;
        int substepsloc  = 1;
        
        well.GetCurrentConfig()->fHasCompletion = 1;
        well.GetCurrentConfig()->fPostprocess.SetCompMesh(0);
        
        TPZElasticResponse ER;
        ER.SetUp(210000., 0.27);
        well.GetCurrentConfig()->fCompletionElastic.SetElasticResponse(ER);
        
        
        WellPressure = 54.; //66.6 61.1 57.2
        well.GetCurrentConfig()->fWellborePressure = WellPressure;
        well.GetCurrentConfig()->SetWellPressure();
        
        
        
        well.ExecuteSimulation(substepsloc,out);
        well.PostProcess(0);
        
        TPZBFileStream save;
        save.OpenWrite("Config3-9.bin");
        well.Write(save);
        
        
        /*
         std::cout << "\n ------- 2 -------- " << std::endl;
         REAL a, b;
         well.ComputeAandB(sqj2_refine, a,b);
         if(a >innerradius )
         {
         well.AddEllipticBreakout(a, b,std::cout);
         int substepsloc = 1;
         well.ExecuteSimulation(substepsloc,std::cout);
         well.PostProcess(0);
         }
         TPZBFileStream save;
         save.OpenWrite("Config3-9.bin");
         well.Write(save);
         }
         if (Startfrom ==3)
         {
         TPZBFileStream read;
         read.OpenRead("Config3-9.bin");
         well.Read(read);
         }
         if (Startfrom <=3)
         {
         well.GetCurrentConfig()->fHasCompletion = 1;
         well.GetCurrentConfig()->fPostprocess.SetCompMesh(0);
         
         TPZElasticResponse ER;
         ER.SetUp(210000., 0.27);
         well.GetCurrentConfig()->fCompletionElastic.SetElasticResponse(ER);
         well.GetCurrentConfig()->AddLinerComputationalElements(well.GetCurrentConfig()->fInitialPOrderLiner);
         
         well.ExecuteSimulation(1, std::cout);
         
         well.PostProcess(1);
         
         TPZBFileStream save;
         save.OpenWrite("Config4-9.bin");
         well.Write(save);
         */
    }
}

int main(int argc, char **argv)
{
    clarg::parse_arguments(argc, argv);
    
    plast_tot.start();
    Config1();
    //    Config2();
    //    Config3();
    //    Config8();
    //    Config6();
    //    Config7();
//    Config9();
    plast_tot.stop();
    
    
    return 0;
}



