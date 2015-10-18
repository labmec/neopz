#include "pzlog.h"
#include "tpzautopointer.h"
#include "ReservoirData.h"

#include "ReducedPVT.h"

#include "WaterPhase.h"
#include "OilPhase.h"
#include "GasPhase.h"

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
    
    TPZMaterial::gBigNumber = 1.0e9; // Use this for check of convergence using neumann
    //    TPZMaterial::gBigNumber = 1.0e15;
    
    NonlinearTracerDimensionless();
    
    
    std::cout << "Program successfully executed." << std::endl;
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
    REAL Tstr           = 355.37;
    REAL Tres           = 355.37;
    REAL Lstr           = 1000.0;
//    REAL gcstr          = 9.81;
    REAL Mustr          = 0.001;
    REAL Rhostr         = 1000.0;
//    REAL Lambdastr      = Rhostr/Mustr;
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 40;
    bool broyden    = false;    // Use this when more than 10000 DOF are required don't used for now!
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;    // Band optimization
    int fixedJac    = 0;
    
    int qorder      = 5;
    int porder      = 5;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 0;
    
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    REAL dt         = 10000000.*day*((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL maxtime    = 10000000.0*day*((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL t0         = 0.0*day*((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    
    REAL TolDeltaX  = 1.0*1e-7;
    REAL TolRes     = 1.0*1e-7;
    
    int  nelemX     =4;
    REAL lengthX    =250.0/Lstr;
    
    int nelemY      =4;
    REAL lengthY    =250.0/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -0.0*((Lstr*Rhostr)/Pstr);
    
    TPZStack<std::string> system;
    system.Push("Oil");
//    system.Push("Water");
//    system.Push("Gas");
    
    Dataset->SetsystemType(system);
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
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> bottombcini(4,0.0);
    bottombcini[0] = 1;
    bottombcini[1] = 0.0;
    bottombcini[2] = 0;
    bottombcini[3] = 1;
    
    TPZVec<REAL> rightbcini(4,0.0);
    rightbcini[0] = 2;
    rightbcini[1] = (1.0*1e6)/(Pstr);
    rightbcini[2] = 0;
    rightbcini[3] = 0;
    
    TPZVec<REAL> topbcini(4,0.0);
    topbcini[0] = 1;
    topbcini[1] = 0.0;
    topbcini[2] = 0;
    topbcini[3] = 1;
    
    TPZVec<REAL> leftbcini(4,0.0);
    leftbcini[0] = 2;
    leftbcini[1] = (1.0*1e6)/(Pstr);
    leftbcini[2] = 0;
    leftbcini[3] = 0;
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Boundary Value Problem
    
    TPZVec<REAL> bottombc(4,0.0);
    bottombc[0] = 1;
    bottombc[1] = 0;
    bottombc[2] = 0;
    bottombc[3] = 0;
    
    TPZVec<REAL> rightbc(4,0.0);
    rightbc[0] = 2;
    rightbc[1] = log(2.0)+0.1;
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = 0;
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = 0.00718196;
    leftbc[2] = 1;
    leftbc[3] = 0;
    
    Dataset->SetBottomBC(bottombcini, bottombc);
    Dataset->SetRightBC(rightbcini, rightbc);
    Dataset->SetTopBC(topbcini, topbc);
    Dataset->SetLeftBC(leftbcini, leftbc);
    
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> Layer         = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel   = new PetroPhysicData;
    TPZAutoPointer<OilPhase> Oil          = new OilPhase;   // alpha
    TPZAutoPointer<WaterPhase> Water        = new WaterPhase;   // beta
    TPZAutoPointer<GasPhase> Gas          = new GasPhase;   // gamma
    
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
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
    REAL oildensity         = 1000.0/Rhostr;
    REAL oilviscosity       = 0.001/Mustr;
    REAL coil               = (1.0*1e-7)*Pstr;
    REAL p_g_ref            = Pstr;
    REAL gasdensity         = Rhostr;
    REAL gasviscosity       = Mustr;
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
    
    
    Water->SetRho(waterdensity);
    Water->SetMu(waterviscosity);
    Water->Setc(cwater);
    Water->SetPRef(p_w_ref);
    Water->SetTRef(Tstr);
    Water->SetTRes(Tres);
    
    Oil->SetRho(oildensity);
    Oil->SetMu(oilviscosity);
    Oil->Setc(coil);
    Oil->SetPRef(p_o_ref);
    Oil->SetTRef(Tstr);
    Oil->SetTRes(Tres);
    
    Gas->SetRho(gasdensity);
    Gas->SetMu(gasviscosity);
    Gas->Setc(cgas);
    Gas->SetPRef(p_g_ref);
    Gas->SetTRef(Tstr);
    Gas->SetTRes(Tres);
    
    TPZVec< TPZAutoPointer<ReducedPVT> > PVTData(3);
    PVTData[2] = Oil.operator->();      // alpha
    PVTData[1] = Water.operator->();    // beta
    PVTData[0] = Gas.operator->();      // gamma
    
//    TPZManVector<REAL> rho(5,0),sv(4,0);
//    sv[1] = 0.167128;
//    Gas->Z(rho, sv);
//    Gas->Bg(rho, sv);
//    Gas->Density(rho, sv);
    
    // Creating the analysis
    
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers,Rocks);
    SandStone.SetFluidData(PVTData);
    SandStone.Run();
    
}