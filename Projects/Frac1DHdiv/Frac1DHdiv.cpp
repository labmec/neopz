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

  const REAL phi= 0.1;
  const REAL Density = 1000.0;
  const REAL theta = 1.0;
  
  TPZFMatrix<STATE> Kabolute(2,2,0.0);
  Kabolute(0,0) = 1.0e-13;
  Kabolute(1,1) = 1.0e-13;

  const REAL mu = 1.0e-3;
  const REAL InitTime = 0.;
  const REAL    timeStep = 0.1;
  const REAL    Lfrac = 1000.;
  const REAL Ttot = 1.;

  const REAL hf = 50000.;
  const REAL    nu = 0.2;
  const REAL    E = 1.e4;

  const REAL Q = 100.;
  const REAL    SigmaConf = 0.;



  const int pOrdQDarcy = 1;
  const int pOrdPDarcy = 1;

  const int pOrdQFrac = 1;
  const int pOrdPFrac = 0;
  const int nel = 100;
  std::string PostProcessFileName = "TransientMathematica.vtk";

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
  
//    Data->SetPorderFlow(pOrdQFrac);
//    Data->SetPorderPressure(pOrdPFrac);
//  TPZFracAnalysis fracAn(Data);
//  fracAn.Run();

    Data->SetPorderFlow(pOrdQDarcy);
    Data->SetPorderPressure(pOrdQDarcy);
    TPZDarcyAnalysis DarcyAn(Data);
    DarcyAn.Run();
    
  return 0;
}
