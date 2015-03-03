
#include "tpzautopointer.h"
#include "ReservoirData.h"
#include "SimulationData.h"
#include "TPZDarcyAnalysis.h"
#include "pzlog.h"


#include <time.h>


#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.DarcyFlow"));
#endif

int main()
{

    // Simulation Data SI units
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 5;
    bool broyden    = false;
    
    REAL hour       = 3600;
    REAL day        = hour * 24;
    REAL dt         = 1.0*day;
    REAL time       = 100.0*day;
    REAL TolDeltaX  = 1.0*10e-5;
    REAL TolRes     = 1.0*10e-5;

    Dataset->SetDeltaT(dt);
    Dataset->SetTime(time);
    Dataset->SetToleranceDX(TolDeltaX);
    Dataset->SetToleranceRes(TolRes);
    Dataset->SetMaxiterations(maxiter);
    Dataset->SetIsBroyden(broyden);
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer = new ReservoirData;
    
    REAL porosityref    = 0.1;
    REAL densityref     = 900.0;
    REAL viscosityref   = 0.001;
    REAL pressureref    = 1.0*10e6;
    REAL lengthref      = 1.0;
    REAL kref           = 1.0;
    REAL crock          = 1.0*10e-10;
    REAL cfluid         = 1.0*10e-9;
    
    TPZVec<int> MatIds(5);
    MatIds[0]=1;
    MatIds[1]=2;
    MatIds[2]=3;
    MatIds[3]=4;
    MatIds[4]=5;
    
    TPZFMatrix<STATE> Kabsolute(2,2);
    Kabsolute.Zero();
    Kabsolute(0,0) = 1.0;
    Kabsolute(1,1) = 1.0;
    
    Layer->SetLref(lengthref);
    Layer->SetKref(kref);
    Layer->SetEtaref(viscosityref);
    Layer->Rhoref(densityref);
    Layer->SetPhiRef(porosityref);
    Layer->SetPref(pressureref);
    Layer->SetcRock(crock);
    Layer->SetcFluid(cfluid);
    Layer->SetKabsolute(Kabsolute);
    Layer->SetMatIDs(MatIds);
    
    // Creating the analysis
    
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers);
    SandStone.Run();
    
	return 0;
}


