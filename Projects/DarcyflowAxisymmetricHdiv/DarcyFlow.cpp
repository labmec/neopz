#include "pzlog.h"
#include "tpzautopointer.h"
#include "ReservoirData.h"

#include "Phase.h"

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
    
    REAL Kstr           = 1.0e-14;
    REAL Pstr           = 1.0e7;
    REAL Tstr           = 355.37;
    REAL Tres           = 355.37;
    REAL Lstr           = 100.0;
    REAL Mustr          = 0.001;
    REAL Rhostr         = 1000.0;
//    REAL lambdastr      = Rhostr/Mustr;
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 20;
    int nthread     = 8;
    bool broyden    = false;    // Use this when more than 10000 DOF are required don't used for now!
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;    // Band optimization
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 0;
    
    int n_times = 11;
    int n_sub_dt = 50;
    TPZManVector<REAL,10> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 10*day*scale;
    REAL t0         = 0.0*day*scale;

    Reporting_times[0] =  1.0*day*scale;
    Reporting_times[1] = 100.0*day*scale;
    Reporting_times[2] = 200.0*day*scale;
    Reporting_times[3] = 300.0*day*scale;
    Reporting_times[4] = 400.0*day*scale;
    Reporting_times[5] = 500.0*day*scale;
    Reporting_times[6] = 600.0*day*scale;
    Reporting_times[7] = 700.0*day*scale;
    Reporting_times[8] = 800.0*day*scale;
    Reporting_times[9] = 900.0*day*scale;
    Reporting_times[10] = 1000.0*day*scale;
    
    REAL maxtime    = Reporting_times[n_times-1];
    
    REAL TolRes     = 1.0*1e-5;
    REAL TolDeltaX  = 1.0*1e-10;
    
    int  nelemX     =1000;
    REAL dxD        =0.1/Lstr;
    
    int nelemY      =1;
    REAL dyD        =100.0/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -0.0*((Lstr*Rhostr)/Pstr);
    
    REAL angle = 0.0;
    
    REAL S_w_r              = 0.0;
    REAL S_nw_r             = 0.0;
    
    TPZStack<std::string> system;
    system.Push("Oil");
    system.Push("Water");
    //    system.Push("Gas");

    Dataset->SetTimes(Reporting_times);
    Dataset->SetNSubSteps(n_sub_dt);
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
    Dataset->SetNthreads(nthread);
    Dataset->SetIsBroyden(broyden);
    Dataset->SetnElementsx(nelemX);
    Dataset->SetnElementsy(nelemY);
    Dataset->SetLengthElementx(dxD);
    Dataset->SetLengthElementy(dyD);
    Dataset->SetGravity(Gravity);
    Dataset->SetRotationAngle(angle);
    
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
    rightbc[1] = 0.1;
    rightbc[2] = S_w_r;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = 0;
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = -(0.1);
    leftbc[2] = 1.0;
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
    REAL porosityref    = 0.2;
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
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e7)/(Pstr);
    REAL oildensity         = 1000.0/Rhostr;
    REAL oilviscosity       = 0.001/Mustr;
    REAL coil               = (0.0*1.0*1e-9)*Pstr;
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
    Kabsolute(0,0) = (1.0e-14)/Kstr;
    Kabsolute(1,1) = (1.0e-14)/Kstr;
    
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
    Layer->SetS_wett_r(S_w_r);
    Layer->SetS_nwett_r(S_nw_r);
    
    
    Water->SetRho(waterdensity);
    Water->SetMu(waterviscosity);
    Water->Setc(cwater);
    Water->SetPRef(p_w_ref);
    Water->SetTRef(Tstr);
    Water->SetTRes(Tres);
    Water->SetS_wett_r(S_w_r);
    Water->SetS_nwett_r(S_nw_r);
    
    Oil->SetRho(oildensity);
    Oil->SetMu(oilviscosity);
    Oil->Setc(coil);
    Oil->SetPRef(p_o_ref);
    Oil->SetTRef(Tstr);
    Oil->SetTRes(Tres);
    Oil->SetS_wett_r(S_w_r);
    Oil->SetS_nwett_r(S_nw_r);
    
    Gas->SetRho(gasdensity);
    Gas->SetMu(gasviscosity);
    Gas->Setc(cgas);
    Gas->SetPRef(p_g_ref);
    Gas->SetTRef(Tstr);
    Gas->SetTRes(Tres);
    Gas->SetS_wett_r(S_w_r);
    Gas->SetS_nwett_r(S_nw_r);
    
    TPZVec< TPZAutoPointer<Phase> > PVTData(3);
    PVTData[0] = Water.operator->();      // alpha
    PVTData[1] = Oil.operator->();    // beta
    PVTData[2] = Gas.operator->();      // gamma
    
    
    TPZManVector<REAL> krw(5,0.0);
    TPZManVector<REAL> kro(5,0.0);
    TPZManVector<REAL> s_vars(5,0.0);
    
    s_vars[2] = 0.8;
    Water->Kr(krw, s_vars);
    s_vars[2] = 1.0-0.8;
    Oil->Kr(kro, s_vars);
    
    std::cout << "krw = " << krw << std::endl;
    std::cout << "kro = " << kro << std::endl;
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    REAL VD = 100.0;
    REAL cfl = (VD*dt*dyD)/(dxD*dyD*porosityref);
    std::cout << "Starting with CFL = " << cfl << std::endl;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers,Rocks);
    SandStone.SetFluidData(PVTData);
    SandStone.RunAnalysis();
    
    std::cout << "Finished with CFL = " << cfl << std::endl;
    
}