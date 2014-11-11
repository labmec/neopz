
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

using namespace pzshape; // needed for TPZShapeCube and related classes



#include "pzlog.h"
//#include "tpztimer.h"
#include "TPZTimer.h"
#include "WellBoreAnalysis.h"
#include "pzbfilestream.h"
#include "TPZProjectEllipse.h"
#include "arglib.h"
#include "run_stats_table.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.main"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerEllipse(Logger::getLogger("LogEllipse"));
#endif


RunStatsTable plast_tot("-tpz_plast_tot", "Raw data table statistics for the main execution.");

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
        
        
        well.GetCurrentConfig()->fWellConfig = EVerticalWell;
        
        
        well.GetCurrentConfig()->CreateMesh();
        
        
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
        
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
            well.PostProcess(0);
        }
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
            well.PostProcess(0);
        }
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
            well.PostProcess(0);
        }
        TPZBFileStream save;
        save.OpenWrite("Config1-4.bin");
        well.Write(save);
        
        
    }
    
    
    
    
 
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
        
        
        well.GetCurrentConfig()->CreateMesh();
        
        
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
        
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        
        
        well.GetCurrentConfig()->CreateMesh();
        
        
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
        
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
        
        std::cout << "\n ------- 5 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        well.PostProcess(0);
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        if(a >innerradius )
        {
            well.AddEllipticBreakout(a, b);
            well.ExecuteSimulation();
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
    
    
    std::string output = "Config3.vtk";
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
      well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
      well.GetCurrentConfig()->fModel = EElastic;
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
        well.GetCurrentConfig()->CreateMesh();
        
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
        std::cout << "\n ------- 1 -------- "<<std::endl;
        
        well.PRefineElementAbove(0.0001, 2);
        well.ExecuteSimulation();
        
        cout << "Average vertical stress " << well.GetCurrentConfig()->AverageVerticalStress() << std::endl;
        
        well.SetFluidModel(EPenetrating);
        well.ExecuteSimulation();
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

int main(int argc, char **argv)
{
    clarg::parse_arguments(argc, argv);
    
    plast_tot.start();
    Config1();
    Config2();
    Config3();
    Config4();
    plast_tot.stop();
    
    
    return 0;
}



