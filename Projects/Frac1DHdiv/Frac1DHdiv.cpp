#include "pzlog.h"
#include "pzmaterial.h"
#include "TPZFracAnalysis.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif

int main()
{
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  // Parametros
  /*
  const REAL mu = 1./12., InitTime = 0., timeStep = 0.1, Lfrac = 1., Ttot = 1.;
  const REAL hf = 1., nu = 0.2, E = 1.e3;
  const REAL Q = 10., SigmaConf = 90.0;
   */
  
  /*
  const REAL mu = 1.e-8, InitTime = 0., timeStep = 0.1, Lfrac = 1000., Ttot = 1.;
  const REAL hf = 1000., nu = 0.2, E = 1.e3;
  const REAL Q = 0.01e6, SigmaConf = 10.;
   */

  const REAL mu = 1.e-8, InitTime = 0., timeStep = 1., Lfrac = 1000., Ttot = 2.;
  const REAL hf = 50000., nu = 0.2, E = 1.e4;
  const REAL Q = 100., SigmaConf = 0.;
  
  
  const int pOrdQ = 1;
  const int pOrdP = 0;
  const int nel = 100;
  std::string PostProcessFileName = "TransientMathematica.vtk";
  TPZAutoPointer<TPZFracData> Data = new TPZFracData;
  Data->SetPostProcessFileName(PostProcessFileName);
  Data->SetTotalTime(Ttot);
  Data->SetTime(InitTime);
  Data->SetTimeStep(timeStep);
  Data->SetViscosity(mu);
  Data->SetLfrac(Lfrac);
  Data->SetHf(hf);
  Data->SetPoisson(nu);
  Data->SetE(E);
  Data->SetNelFrac(nel);
  Data->SetPorderFlow(pOrdQ);
  Data->SetPorderPressure(pOrdP);
  Data->SetQ(Q);
  Data->SetSigmaConf(SigmaConf);
  
  TPZFracAnalysis fracAn(Data);
  fracAn.Run();
  
  return 0;
}
