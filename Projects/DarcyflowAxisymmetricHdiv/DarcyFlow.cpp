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
void NonlinearTracer();
void NonlinearTracerDimensionless();

int main()
{
  
    TPZMaterial::gBigNumber = 1.0e10; // Use this for check of convergence using neumann
//    TPZMaterial::gBigNumber = 1.0e15;
    
//  LinearTracer();
//  NonlinearTracer();
  NonlinearTracerDimensionless();
  
    
  std::cout << " Process complete normally." << std::endl;
  return 0; 
}

void NonlinearTracerDimensionless()
{
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!
    
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    // Simulation Data SI units
    
    // Characteristic Data
    
    REAL Kstr           = 1.0e-15;
    REAL Pstr           = 1.0e7;
    REAL Lstr           = 1000.0;
    REAL gcstr          = 9.81;
    REAL Mustr          = 0.001;
    REAL Rhostr         = ((Pstr)/(Lstr*gcstr));
    REAL Lambdastr      = Rhostr/Mustr;
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 20;
    bool broyden    = false;    // Use this when more than 10000 DOF are required don't used for now!
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = true;    // false means GMRES
    bool OptBand    = false;    // Band optimization
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 4;
    
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    REAL dt         = 1000.0*day*((Kstr*Lambdastr*gcstr)/(Lstr));
    REAL maxtime    = 10000.0*day*((Kstr*Lambdastr*gcstr)/(Lstr));
    REAL t0         = 0.0*day*((Kstr*Lambdastr*gcstr)/(Lstr));
    
    REAL TolDeltaX  = 1.0*1e-5;
    REAL TolRes     = 1.0*1e-5;
    
    int  nelemX     =4;
    REAL lengthX    =250.0/Lstr;
    
    int nelemY      =4;
    REAL lengthY    =250.0/Lstr;
    
    Gravity(0,0)= -0.0;
    Gravity(1,0)= -0.0;

    
    
    
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
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
    Dataset->SetGravity(Gravity);
//    int typeFluxin = 1, typePressurein = 0;
//    int typeFluxout = 3, typePressureout = 2;
    // BCs
    
    TPZVec<REAL> bottombc(4,0.0);
    bottombc[0] = 1;
    bottombc[1] = 0;
    bottombc[2] = 0;
    bottombc[3] = 0;
    
    TPZVec<REAL> rightbc(4,0.0);
    rightbc[0] = 2;
    rightbc[1] = (1.0*1e6)/(Pstr);
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = 0;
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = -10.0;
    leftbc[2] = 1;
    leftbc[3] = 0;
    
    Dataset->SetBottomBC(bottombc);
    Dataset->SetRightBC(rightbc);
    Dataset->SetTopBC(topbc);
    Dataset->SetLeftBC(leftbc);
    
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel = new PetroPhysicData;
    TPZAutoPointer<ReducedPVT> FluidModel = new ReducedPVT;
    
    // Complete data set
    
    // Reservoir Description
    bool isGIDGeom      = false;
    REAL porosityref    = 0.1;
    REAL pressureref    = (1.0*1e6)/(Pstr);
    REAL lengthref      = 1.0;
    REAL kref           = 1.0;
    REAL crock          = (0.0*1e-10)*Pstr;
    REAL Hres           = 100.0/Lstr;
    REAL Rres           = 1000.0/Lstr;
    REAL Top            = -3000.0/Lstr;
    REAL Rw             = 0.0/Lstr;
    
    // Reservoir Description linear tracer configuration
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (1.0*1e-9)*Pstr;
    REAL oildensity         = 1000.0/Rhostr;
    REAL oilviscosity       = 0.001/Mustr;
    REAL coil               = (1.0*1e-9)*Pstr;
    REAL gasdensity         = 0.0/Rhostr;
    REAL gasviscosity       = 0.0/Mustr;
    REAL cgas               = (0.0)*Pstr;
    
    
    TPZVec<int> MatIds(5);
    MatIds[0]=1;
    MatIds[1]=2;
    MatIds[2]=3;
    MatIds[3]=4;
    MatIds[4]=5;
    
    TPZFMatrix<STATE> Kabsolute(2,2);
    Kabsolute.Zero();
    Kabsolute(0,0) = (1.0e-15)/Kstr;
    Kabsolute(1,1) = (1.0e-15)/Kstr;
    
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
    FluidModel->SetcOil(coil);
    FluidModel->SetRhoGas(gasdensity);
    FluidModel->SetMuGas(gasviscosity);
    FluidModel->SetcGas(cgas);
    FluidModel->SetPref(pressureref);
    
    
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

void NonlinearTracer()
{
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!

    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    // Simulation Data SI units
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    TPZFMatrix<REAL> Gravity(2,1);
	 
    int maxiter     = 40;
    bool broyden    = false;    // Use this when more than 10000 DOF are required don't used for now!
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = true;     // Use Static Condensation
    bool IsDirect   = true;     // Not Used broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = false;    // Band optimization
    int fixedJac    = 0;
    
    int qorder      = 2;
    int porder      = 2;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 3;
    
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 1.0*day;
    
    REAL maxtime    = 2.0*day;
    REAL t0         = 0.0*day;
    REAL TolDeltaX  = 1.0*1e-5;
    REAL TolRes     = 1.0*1e-5;
    
    int  nelemX     =4;
    REAL lengthX    =250.0;
    
    int nelemY      =1;
    REAL lengthY    =100.0;
    
    Gravity(0,0)= -0.0;
    Gravity(1,0)= -0.0;

    
    
    
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
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
    Dataset->SetGravity(Gravity);
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    // BCs
    
    TPZVec<REAL> bottombc(4,0.0);
    bottombc[0] = 0;
    bottombc[1] = 0;
    bottombc[2] = 0;
    bottombc[3] = 0;
    
    TPZVec<REAL> rightbc(4,0.0);
    rightbc[0] = 2;
    rightbc[1] = (1.0*1e6);
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 0;
    topbc[1] = 0;
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 0;
    leftbc[1] = (10.0*1e6);
    leftbc[2] = 1;
    leftbc[3] = 0;
    
    Dataset->SetBottomBC(bottombc);
    Dataset->SetRightBC(rightbc);
    Dataset->SetTopBC(topbc);
    Dataset->SetLeftBC(leftbc);
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel = new PetroPhysicData;
    TPZAutoPointer<ReducedPVT> FluidModel = new ReducedPVT;
    
    // Complete data set
    
    // Reservoir Description
    bool isGIDGeom      = false;
    REAL porosityref    = 0.1;
    REAL pressureref    = 1.0*1e6;
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
    REAL cwater             = 1.0*1e-5;
    REAL oildensity         = 1000.0;
    REAL oilviscosity       = 0.001;
    REAL coil               = 1.0*1e-5;
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
    Kabsolute(0,0) = (1.0e-15);
    Kabsolute(1,1) = (1.0e-15);
    
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
    FluidModel->SetcOil(coil);
    FluidModel->SetRhoGas(gasdensity);
    FluidModel->SetMuGas(gasviscosity);
    FluidModel->SetcGas(cgas);
    FluidModel->SetPref(pressureref);
    
    
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

void LinearTracer()
{
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!
    
  
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
//    TPZMaterial::gBigNumber = 1.0e15;
    
    // Simulation Data SI units
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    TPZMaterial::gBigNumber = 1.0e12;
    TPZFMatrix<REAL> Gravity(2,1);
    
    int maxiter     = 40;
    bool broyden    = false;    // Use this when more than 10000 DOF are required don't used for now!
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = true;     // Use Static Condensation
    bool IsDirect   = true;     // Not Used broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = false;    // Band optimization
    bool IsDimen    = true;     // Compute solution with corresponding dimensionless PDE system
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 0;
    
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 0.5*day;

    REAL maxtime    = 10.0*day;
    REAL t0         = 0.0*day;
    REAL TolDeltaX  = 1.0*1e-5;
    REAL TolRes     = 1.0*1e-5;
    
    int  nelemX     =10;
    REAL lengthX    =25.0;
    
    int nelemY      =6;
    REAL lengthY    =25.0;

    Gravity(0,0)= -0.0;
    Gravity(1,0)= -0.0;

    
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
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
    
    
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    // BCs
    

    
    TPZVec<REAL> bottombc(4,0.0);
    bottombc[0] = 1;
    bottombc[1] = 0;
    bottombc[2] = 0;
    bottombc[3] = 0;
    
    TPZVec<REAL> rightbc(4,0.0);
    rightbc[0] = 2;
    rightbc[1] = (1.0*1e6);
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = 0;
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 0;
    leftbc[1] = (5.0*1e6);
    leftbc[2] = 1;
    leftbc[3] = 0;
    
    Dataset->SetBottomBC(bottombc);
    Dataset->SetRightBC(rightbc);
    Dataset->SetTopBC(topbc);
    Dataset->SetLeftBC(leftbc);
    Dataset->SetGravity(Gravity);
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel = new PetroPhysicData;
    TPZAutoPointer<ReducedPVT> FluidModel = new ReducedPVT;
    
    // Complete data set
    
    // Reservoir Description
    bool isGIDGeom      = false;
    REAL porosityref    = 0.1;
    REAL pressureref    = 1.0*1e6;
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
    REAL oildensity         = 500.0;
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
    FluidModel->SetPref(pressureref);
    
    
    
    
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
