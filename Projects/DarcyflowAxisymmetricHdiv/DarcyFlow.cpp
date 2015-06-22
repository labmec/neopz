  #include "pzlog.h"
#include "tpzautopointer.h"
#include "ReservoirData.h"

#include "ReducedPVT.h"
#include "PetroPhysicData.h"
#include "SimulationData.h"
#include "TPZDarcyAnalysis.h"
#include "pzlog.h"


#include <time.h>

void LinearTracer();

int main()
{

  LinearTracer();
  
  std::cout << " Process complete normally." << std::endl;
  return 0; 
}

void LinearTracer()
{
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!
    
  
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 0;
    
    // Simulation Data SI units
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 40;
    bool broyden    = false;    // Use this when more than 10000 DOF are required don't used for now!
    bool GR         = false;    // Use Gradient Reconstruction
    bool IsDirect   = true;     // Not Used broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = false;    // Band optimization
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 0;
    
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 10.0*day;

    REAL maxtime    = 500.0*day;
    REAL t0         = 0.0*day;
    REAL TolDeltaX  = 1.0*1e-4;
    REAL TolRes     = 1.0*1e-4;
    
    int  nelemX     =3;
    REAL lengthX    =100.0;
    
    int nelemY      =3;
    REAL lengthY    =100.0;


    
    Dataset->SetGR(GR);
    Dataset->SetIsDirect(IsDirect);
    Dataset->SetIsCG(IsCG);
    Dataset->SetOptband(OptBand);
    Dataset->Setqorder(qorder);
    Dataset->Setporder(porder);
    Dataset->Setsorder(sorder);
    Dataset->SetHrefinement(hrefinement);
    Dataset->SetHPostrefinement(hpostref);
    Dataset->SetDeltaT(dt);
    Dataset->SetMaxTime(maxtime);    
    Dataset->SetTime(t0);
    Dataset->SetToleranceDX(TolDeltaX);
    Dataset->SetToleranceRes(TolRes);
    Dataset->SetMaxiterations(maxiter);
    Dataset->SetFixediterations(fixedJac);
    Dataset->SetIsBroyden(broyden);
    Dataset->SetnElementsx(nelemX);
    Dataset->SetnElementsy(nelemY);
    Dataset->SetLengthElementx(lengthX);
    Dataset->SetLengthElementy(lengthY);
    
    
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel = new PetroPhysicData;
    TPZAutoPointer<ReducedPVT> FluidModel = new ReducedPVT;
    
    // Complete data set
    
    // Reservoir Description
    bool isGIDGeom      = false;
    REAL porosityref    = 0.1;
    REAL pressureref    = 1.0*10e6;
    REAL lengthref      = 1.0;
    REAL kref           = 1.0;
    REAL crock          = 0.0*1e-10;
    REAL Hres           = 100.0;
    REAL Rres           = 1000.0;
    REAL Top            = -3000.0;
    REAL Rw             = 0.0;

    // Reservoir Description linear tracer configuration
    REAL waterdensity       = 1000.0;
    REAL waterviscosity     = 0.001;
    REAL cwater             = 0.0*1e-8;
    REAL oildensity         = 800.0;
    REAL oilviscosity       = 0.001;
    REAL coil               = 0.0*1e-8;
    REAL gasdensity         = 0.0;
    REAL gasviscosity       = 0.0;
    REAL cgas               = 0.0;
    
    
    TPZVec<int> MatIds(5);
    MatIds[0]=1;
    MatIds[1]=2;
    MatIds[2]=3;
    MatIds[3]=4;
    MatIds[4]=5;
    
    TPZFMatrix<STATE> Kabsolute(2,2);
    Kabsolute.Zero();
    Kabsolute(0,0) = 1.0e-14;
    Kabsolute(1,1) = 1.0e-14;
    
    Layer->SetIsGIDGeometry(isGIDGeom);
    Layer->SetLayerTop(Top);
    Layer->SetLayerrw(Rw);
    Layer->SetLayerh(Hres);
    Layer->SetLayerr(Rres);
    Layer->SetLref(lengthref);
    Layer->SetKref(kref);
    Layer->SetPhiRef(porosityref);
    Layer->SetPref(pressureref);
    Layer->SetcRock(crock);
    Layer->SetKabsolute(Kabsolute);
    Layer->SetMatIDs(MatIds);
    
    FluidModel->SetRhoWater(waterdensity);
    FluidModel->SetMuWater(waterviscosity);
    FluidModel->SetcWater(cwater);
    FluidModel->SetRhoOil(oildensity);
    FluidModel->SetMuOil(oilviscosity);
    FluidModel->SetcWater(coil);
    FluidModel->SetRhoGas(gasdensity);
    FluidModel->SetMuGas(gasviscosity);
    FluidModel->SetcGas(cgas);
    
    
    
    
    // Creating the analysis
    
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers,Rocks,FluidModel);
    SandStone.Run();
    
 
  
}
