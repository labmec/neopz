#include "pzlog.h"
#include "TPZMaterial.h"
#include "TPZDarcyAnalysis.h"
#include "TPZFracAnalysis.h"
#include "pzshapelinear.h"

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.frac"));
#endif

int main()
{
    
    TPZMaterial::gBigNumber = 1.e7; // so flow imposition in fracture can have low residuum
    
    // ---------- Parameters ----------
    
    // Reservoir Data
    const REAL phi = 0.1;
    TPZFMatrix<STATE> Kabolute(2,2,0.0);
    Kabolute(0,0) = 20.0*(1.0e6)*4.93466e-14;// 4.93466e-14 m2 -> 50 md
    Kabolute(1,1) = 20.0*(1.0e6)*4.93466e-14;// 4.93466e-14 m2 -> 50 md
    
    const REAL nu = 0.2;
    const REAL E = 1.e4;
    const REAL SigmaConf = 20.; // Sigma min
    
    // Simulation Data
    const REAL theta = 1.; // MUST BE ALWAYS 1
    const REAL InitTime = 0.; // Of course 0
    const REAL timeStep = 1.;
    const REAL Ttot = 1000.;
    const int pOrdQDarcy = 1;
    const int pOrdPDarcy = 1;
    const int pOrdQFrac = 1;
    const int pOrdPFrac = 0;
    const int nel = 100; // Deprecated
    const int npropag = 40; // use for darcy gmesh n elements in x direction
    const int hy = 8; // Number of refinements in y direction
    const REAL Ly = 2000.0; // Reservoir width
    const REAL elsize = 5.; // mm
    const int nthreadsForAssemble = 16; // if 0, program is serial in assemble
    bool isCoupled = false;
    bool plotVTK = true; // it is faster is false
    int hrefvtk = 0;
    std::string PostProcessFileName = "NoCoupledk.vtk";
    if (isCoupled) { PostProcessFileName = "Coupledk.vtk";}
    
    // Fluid Data
    const REAL mu = 1.0e-9; // N/mm2 * s = 1.0e-3 Pa*s -> 1 Centipoise cp ->
    const REAL Density = 1000.0; // 1000.0 kg/m3 -> 1.0 gr/cm3
    
    // Fracture Data
    const REAL Lfrac = 1000.; // mm // Deprecated
    const REAL hf = 40000.; // mm
    const REAL Q = 578.7; // mm3/s/mm - 0.023148 m3/s using 40000 mm
    
    // Leak off data
    const REAL Cl = 10.01; // Carter coefficient
    const REAL Pe = 15.0; // Should come from darcy simulation
    const REAL Pref = 5.0; // pressure where Cl was measured
    const REAL vsp = 0.*0.002; // spurt loss - Better keep 0 for now. I havent debugged it
    
    // Filling TPZFracData structure. Every object has an autopointer to this structure called fData
    TPZAutoPointer<TPZFracData> Data = new TPZFracData;
    Data->SetPostProcessFileName(PostProcessFileName);
    Data->SetK(Kabolute);
    Data->SetTheta(theta);
    Data->SetPorosity(phi);
    Data->SetTime(InitTime);
    Data->SetTotalTime(Ttot);
    Data->SetTimeStep(timeStep);
    Data->SetViscosity(mu);
    Data->SetDensity(Density);
    Data->SetLfrac(Lfrac);
    Data->SetHf(hf);
    Data->SetPoisson(nu);
    Data->SetE(E);
    Data->SetNelFrac(nel);
    Data->SetPorderFlow(pOrdQFrac);
    Data->SetPorderPressure(pOrdPFrac);
    Data->SetQ(Q);
    Data->SetSigmaConf(SigmaConf);
    Data->SetCl(Cl);
    Data->SetPe(Pe);
    Data->SetPref(Pref);
    Data->SetVsp(vsp);
    Data->SetElSize(elsize);
    Data->SetDwDp(); // Has to be called after setting frac parameters
    Data->SetNPropagations(npropag);
    Data->SetHy(hy);
    Data->SetLy(Ly);
    Data->SetIsCoupled(isCoupled);
    Data->SetIsPlotVTK(plotVTK);
    Data->SetHrefPostPro(hrefvtk);
    
    Data->SetNThreadsForAssemble(nthreadsForAssemble);
    Data->SetPorderFlow(pOrdQFrac);
    Data->SetPorderPressure(pOrdPFrac);
    Data->SetPorderDarcyFlow(pOrdQDarcy);
    Data->SetPorderDarcyPressure(pOrdPDarcy);
    
    // Fracture Simulation uncoupled
//    TPZFracAnalysis fracAn(Data);
//    fracAn.Run();
    
    // Reservoir Simulation coupled
      TPZDarcyAnalysis DarcyAn(Data);
      DarcyAn.Run();
    
    return 0;
}
