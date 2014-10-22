#include "pzlog.h"
#include "pzmaterial.h"
#include "TPZDarcyAnalysis.h"
#include "TPZFracAnalysis.h"
#include "pzshapelinear.h"

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.frac"));
#endif

int main()
{
#ifdef LOG4CXX
  std::string dirname = PZSOURCEDIR;
  std::string FileName = dirname;
  FileName = dirname + "/Projects/Frac1DHdiv/";
  FileName += "FracLog.cfg";
  InitializePZLOG(FileName);
#endif
  
  TPZMaterial::gBigNumber = 1.e9;
  
  // ---------- Parametros ----------
  // Reservoir Data
  const REAL phi = 0.1;
  TPZFMatrix<STATE> Kabolute(2,2,0.0);
  Kabolute(0,0) = 4.93466e-14;// 4.93466e-14 m2 -> 50 md
  Kabolute(1,1) = 4.93466e-14;// 4.93466e-14 m2 -> 50 md
  REAL day    = 86400.0;
  REAL year    = 365.0*day;
  
  const REAL nu = 0.2;
  const REAL E = 1.e4;
  const REAL SigmaConf = 20.;
  
  // Simulation Data
  const REAL theta = 1.;
  const REAL InitTime = 0.;
  const REAL timeStep = 1.;
  const REAL Ttot = 100.0;
  const int pOrdQDarcy = 1;
  const int pOrdPDarcy = 1;
  const int pOrdQFrac = 1;
  const int pOrdPFrac = 0;
  const int nel = 100;
  const int npropag = 3;
  const int hy = 2;
  const REAL Ly = 100.;
  const REAL elsize = 50.0;
  std::string PostProcessFileName = "SemiCoupled.vtk";
  const int nthreadsForAssemble = 0; // if 0, program is serial
  
  // Fluid Data
  const REAL mu = 1.0e-8; // 1.0e-3 Pa*s -> 1 Centipoise cp
  const REAL Density = 1000.0; //1000.0 kg/m3 -> 1.0 gr/cm3
  
  // Fracture Data
  const REAL Lfrac = 1000.; // mm
  const REAL hf = 50000.; // mm
  const REAL Q = 100.; // mm3/s/mm
  
  // Leak off data
  const REAL Cl = 0.3; // Carter coefficient
  const REAL Pe = 15.; // Should come from darcy simulation
  const REAL Pref = 10.; // pressure where Cl was measured
  const REAL vsp = 0.*0.002; // spurt loss
  
  // Preenchendo estrutura TPZFracData
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
  Data->SetDwDp();
  Data->SetNPropagations(npropag);
  Data->SetHy(hy);
  Data->SetLy(Ly);
  
  
  Data->SetNThreadsForAssemble(nthreadsForAssemble);
  Data->SetPorderFlow(pOrdQFrac);
  Data->SetPorderPressure(pOrdPFrac);
  Data->SetPorderDarcyFlow(pOrdQDarcy);
  Data->SetPorderDarcyPressure(pOrdPDarcy);
  
  // Fracture Simulation uncoupled
//  TPZFracAnalysis fracAn(Data);
//  fracAn.Run();
  
  // Reservoir Simulation uncoupled
  TPZDarcyAnalysis DarcyAn(Data);
  DarcyAn.Run();
  
  return 0;
}
