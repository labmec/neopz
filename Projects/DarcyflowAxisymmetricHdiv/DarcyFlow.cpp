#include "pzlog.h"
#include "tpzautopointer.h"
#include "ReservoirData.h"
#include "SimulationData.h"
#include "TPZDarcyAnalysis.h"
#include "pzlog.h"


#include <time.h>

// #ifdef LOG4CXX
// static LoggerPtr logdata(Logger::getLogger("pz.darcy"));
// #endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.DarcyFlow"));
#endif

int main()
{
    
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    // Simulation Data SI units
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 10;
    bool broyden    = false;    // Use this when more than 10000 DOF are required
    bool h1         = false ;
    bool IsDirect   = false;     //  Not Used broyden with Iterative !!!
    bool IsCG       = true;
    int fixedJac    = 1;
    
    int qorder      = 2;
    int porder      = 2;
    int hrefinement = 2;
    
    REAL hour       = 3600;
    REAL day        = hour * 24;
    REAL dt         = 1.0*day;
    REAL time       = 100.0*day;
    REAL TolDeltaX  = 1.0*1e-5;
    REAL TolRes     = 1.0*1e-5;
    
    Dataset->SetIsH1approx(h1);
    Dataset->SetIsDirect(IsDirect);
    Dataset->SetIsCG(IsCG);
    Dataset->Setqorder(qorder);
    Dataset->Setporder(porder);
    Dataset->SetHrefinement(hrefinement);
    Dataset->SetDeltaT(dt);
    Dataset->SetTime(time);
    Dataset->SetToleranceDX(TolDeltaX);
    Dataset->SetToleranceRes(TolRes);
    Dataset->SetMaxiterations(maxiter);
    Dataset->SetFixediterations(fixedJac);
    Dataset->SetIsBroyden(broyden);
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer = new ReservoirData;
    
    bool isGIDGeom      = true;
    REAL porosityref    = 1.0;//0.1;
    REAL densityref     = 1.0;//1000.0;
    REAL viscosityref   = 1.0;//0.001;
    REAL pressureref    = 1.0*10e6;
    REAL lengthref      = 1.0;
    REAL kref           = 1.0;
    REAL crock          = 0.0*1e-10;
    REAL cfluid         = 0.0*1e-8;
    REAL Hres           = 100.0;
    REAL Rres           = 1000.0;
    REAL Top            = -3000.0;
    REAL Rw             = 0.0;
    
    TPZVec<int> MatIds(5);
    MatIds[0]=1;
    MatIds[1]=2;
    MatIds[2]=3;
    MatIds[3]=4;
    MatIds[4]=5;
    
    TPZFMatrix<STATE> Kabsolute(2,2);
    Kabsolute.Zero();
    Kabsolute(0,0) = 1.0;//1.0e-14;
    Kabsolute(1,1) = 1.0;//1.0e-14;
    
    Layer->SetIsGIDGeometry(isGIDGeom);
    Layer->SetLayerTop(Top);
    Layer->SetLayerrw(Rw);
    Layer->SetLayerh(Hres);
    Layer->SetLayerr(Rres);
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
    
    
    std::cout << " Process complete normally." << std::endl;
    return 0;
}


