
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


RunStatsTable plast_tot("-tpz_plast_tot", "Raw data table statistics for Plasticity::FindBug");

#ifdef USING_TBB
#include "tbb/task_scheduler_init.h"
#endif
void Config2()
{
#ifdef USING_TBB
    tbb::task_scheduler_init init(12);
#endif
    
    InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();
    TPZWellBoreAnalysis well;
    REAL innerradius = 0.10795;
    REAL outerradius = 3.;
    
    REAL sqj2_refine=1.e-7;
    REAL sqj2_ellips = 1.e-7;
    std::cout << std::setprecision(15);
    TPZManVector<STATE,3> confinement(3,0.);
    REAL SH,Sh,SV;
    SH=-62.1;
    Sh=-45.9;
    SV=-48.2;
    
    confinement[0] = Sh;
    confinement[1] = SH;
    confinement[2] = SV;
    REAL effectivePressure = 19.5; // 19.5 ou 23.4 ou 28.9
    well.SetConfinementEffectiveStresses(confinement, effectivePressure);
    well.SetInnerOuterRadius(innerradius, outerradius);
    
    bool modelMC =false;
    if (modelMC)
    {
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL cohesion = 13.;
        REAL Phi = 0.52;
        well.SetMohrCoulombParameters(poisson, elast, cohesion, Phi, Phi);
        
    }
    else
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
    
    
    int Startfrom=0;
    if (Startfrom == 0)
    {
        int porder = 2;
        int nrad=40;
        int ncircle = 20;
        REAL delx = 0.5*innerradius*M_PI_2/ncircle;
        TPZManVector<int,2> numdiv(2);
        numdiv[0] = nrad;
        numdiv[1] = ncircle;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->CreateMesh();
        
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        //REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 80;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
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
        
        well.PRefineElementAbove(sqj2_refine, 4);
        well.DivideElementsAbove(sqj2_refine);
        //well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_ellips, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        save.OpenWrite("wellbore1.bin");
        well.Write(save);
        
        
    }
    
    
    for(int i=0;i<1;i++)
    {
        std::cout << "\n ------- "<< i+2 <<"-------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 4);
        well.DivideElementsAbove(sqj2_refine);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_ellips, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        stringstream outfile;
        outfile<< "wellbore"<<i+2<<".bin";
        save.OpenWrite(outfile.str());
        well.Write(save);
        
    }
    
    
}

void Config1()
{
    TPZTimer time1,time2;
    time1.start();
    InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();
    
    time2.start();
    TPZWellBoreAnalysis well;
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    
    REAL computedquarter = 7.05761678496926;
    REAL sqj2_refine=0.0001;
    std::cout << std::setprecision(15);
    TPZManVector<STATE,3> confinement(3,0.);
    REAL SH,Sh,SV;
    SH=-62.1;
    Sh=-45.9;
    SV=-48.2;
    
    confinement[0] = SH;
    confinement[1] = Sh;
    confinement[2] = SV;
    REAL effectivePressure = 19.5; // 19.5 ou 23.4 ou 28.9
    well.SetConfinementEffectiveStresses(confinement, effectivePressure);
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
        
        int divisions = 20;
        REAL delx = 0.2*innerradius*M_PI_2/divisions;
        TPZManVector<int,2> numdiv(2,divisions);
        numdiv[1] = 40;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->CreateMesh();
        int porder = 2;
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        
//        well.TestLinearMaterial();
        
        well.GetCurrentConfig()->CreatePostProcessingMesh();
        REAL farfieldwork = well.GetCurrentConfig()->ComputeFarFieldWork();
        //well.LinearConfiguration(1);
        well.PostProcess(0);
        
    }
    if (Startfrom ==0)
    {
        
        int nsteps = 5;
        int numnewton = 90;
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteInitialSimulation(nsteps, numnewton);
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
        
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        save.OpenWrite("wellbore1.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==2)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore1.bin");
        well.Read(read);
    }
    
    if (Startfrom <=2)
    {
        
        std::cout << "\n ------- 2 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        save.OpenWrite("wellbore2.bin");
        well.Write(save);
        
    }
    
    
    
    if (Startfrom ==3)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore2.bin");
        well.Read(read);
    }
    
    if (Startfrom <=3)
    {
        
        std::cout << "\n ------- 3 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        save.OpenWrite("wellbore3.bin");
        well.Write(save);
        
        
    }
    
    
    
    if (Startfrom ==4)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore3.bin");
        well.Read(read);
    }
    
    if (Startfrom <=4)
    {
        
        std::cout << "\n ------- 4 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        save.OpenWrite("wellbore4.bin");
        well.Write(save);
        
        
    }
    
    if (Startfrom ==5)
    {
        TPZBFileStream read;
        read.OpenRead("wellbore4.bin");
        well.Read(read);
    }
    
    if (Startfrom <=5)
    {
        
        std::cout << "\n ------- 5 -------- "<<std::endl;
        well.PRefineElementAbove(sqj2_refine, 3);
        well.DivideElementsAbove(sqj2_refine);
        well.ExecuteSimulation();
        REAL a, b;
        well.ComputeAandB(sqj2_refine, a,b);
        well.AddEllipticBreakout(a, b);
        well.ExecuteSimulation();
        TPZBFileStream save;
        save.OpenWrite("wellbore5.bin");
        well.Write(save);
        
        
    }
    
    
    
 
}



int main(int argc, char **argv)
{
    clarg::parse_arguments(argc, argv);
    
    plast_tot.start();
    Config2();
    plast_tot.stop();
    
    
    return 0;
}



