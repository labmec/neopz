#include "pzlog.h"
#include "pzmaterial.h"
#include "TPZDarcyAnalysis.h"
#include "TPZFracAnalysis.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif

int main()
{
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  TPZMaterial::gBigNumber = 1.e8;
  
  // ---------- Parametros ----------
  // Reservoir Data
  const REAL phi = 0.1;
  const REAL Density = 1000.0;
  TPZFMatrix<STATE> Kabolute(2,2,0.0);
  Kabolute(0,0) = 1.0e-13;
  Kabolute(1,1) = 1.0e-13;
  const REAL nu = 0.2;
  const REAL E = 1.e4;
  const REAL SigmaConf = 20.;
  
  // Simulation Data
  const REAL theta = 1.;
  const REAL InitTime = 0.;
  const REAL timeStep = 1.;
  const REAL Ttot = 100.;
  const int pOrdQDarcy = 1;
  const int pOrdPDarcy = 1;
  const int pOrdQFrac = 1;
  const int pOrdPFrac = 0;
  const int nel = 100;
  const REAL elsize = 50.;
  std::string PostProcessFileName = "Propag.vtk";
  
  // Fluid Data
  const REAL mu = 1.e-8; // N/mm2 . s
  
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
  
  // Fracture Simulation uncoupled
  Data->SetPorderFlow(pOrdQFrac);
  Data->SetPorderPressure(pOrdPFrac);
  TPZFracAnalysis fracAn(Data);
  fracAn.Run();

  // Reservoir Simulation uncoupled
//  Data->SetPorderFlow(pOrdQDarcy);
//  Data->SetPorderPressure(pOrdQDarcy);
//  TPZDarcyAnalysis DarcyAn(Data);
//  DarcyAn.Run();
  
  return 0;
}
