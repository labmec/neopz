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
  const REAL mu = 1./12., time = 0., timeStep = TPZMaterial::gBigNumber, Lfrac = 1.;
  const int pOrdQ = 4;
  const int pOrdP = 3;
  TPZAutoPointer<TPZFracData> Data = new TPZFracData;
  Data->SetTime(time);
  Data->SetTimeStep(timeStep);
  Data->SetViscosity(mu);
  Data->SetLfrac(Lfrac);
  Data->SetPorderFlow(pOrdQ);
  Data->SetPorderPressure(pOrdP);
  
  TPZFracAnalysis fracAn(Data);
  fracAn.Run();
  
  return 0;
}
