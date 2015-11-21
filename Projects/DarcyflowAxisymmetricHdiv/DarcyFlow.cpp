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
void NonlinearTracer(bool IsDimensionlessQ);
void LinearWithOutReconstruction(bool IsDimensionlessQ);
void LinearWithReconstruction(bool IsDimensionlessQ);

int main()
{
    
    bool IsDimensionlessQ = true;
//    LinearWithOutReconstruction(IsDimensionlessQ);
//    LinearWithReconstruction(IsDimensionlessQ);

    
    NonlinearTracer(IsDimensionlessQ);
    
    std::cout << "Program successfully executed." << std::endl;
    return 0;
}

void LinearWithOutReconstruction(bool IsDimensionlessQ){
    
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!
    
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    // Simulation Data SI units
    
    // Characteristic Data
    REAL Kstr           = 1.0;
    REAL Pstr           = 1.0;
    REAL Tstr           = 355.37;
    REAL Tres           = 355.37;
    REAL Lstr           = 1.0;
    REAL Mustr          = 1.0;
    REAL Rhostr         = 1.0;
    TPZMaterial::gBigNumber = 1.0e15; // Use this for check of convergence using neumann
    
    REAL TolRes     = 1.0*1e-4;
    REAL TolDeltaX  = 1.0*1e-6;
    
    if (IsDimensionlessQ)
    {
        Kstr           = 1.0e-13;
        Pstr           = 1.0e7;
        Tstr           = 355.37;
        Tres           = 355.37;
        Lstr           = 100.0;
        Mustr          = 0.001;
        Rhostr         = 1000.0;
        TPZMaterial::gBigNumber = 1.0e12; // Use this for check of convergence using neumann
        TolRes     = 1.0*1e-5;
        TolDeltaX  = 1.0*1e-7;
    }
    
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 40;
    int nthread     = 8;
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;    // Band optimization
    bool IsAxisy    = false;    // Axisymmetric analysis
    bool IsImpes    = false;    // Impes analysis
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 0;
    
    // Time control parameters
    int n_times  = 10;
    int n_sub_dt = 100;
    int which_dt = n_times;
    TPZManVector<REAL,20> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = (1.0e7/day) * day * scale;
    REAL t0         = 0.0  * day * scale;
    
    for (int it = 0 ; it < n_times; it++) {
        Reporting_times[it] = REAL(it+1)*dt;
    }
    REAL maxtime    = Reporting_times[which_dt-1];
    std::cout << "Reporting times = " << Reporting_times << std::endl;
    std::cout << "Maximum simulation time = " << maxtime <<std::endl;
    
    int  nelemX     =1000;
    REAL dxD        =(100.0/nelemX)/Lstr;
    
    int nelemY      =2;
    REAL dyD        =(20.0/nelemY)/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -0.0*((Lstr*Rhostr)/Pstr);
    bool LinearSegregation = true;
    
    REAL angle = 0.0;
    
    REAL S_w_r              = 0.0;
    REAL S_nw_r             = 0.0;
    
    TPZStack<std::string> system;
    system.Push("Water");
    system.Push("Oil");
    
    
    Dataset->SetTimes(Reporting_times);
    Dataset->SetNSubSteps(n_sub_dt);
    Dataset->SetsystemType(system);
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
    Dataset->SetIsDirect(IsDirect);
    Dataset->SetIsCG(IsCG);
    Dataset->SetOptband(OptBand);
    Dataset->SetAxisymmetricQ(IsAxisy);
    Dataset->SetImpesQ(IsImpes);
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
    Dataset->SetnElementsx(nelemX);
    Dataset->SetnElementsy(nelemY);
    Dataset->SetLengthElementx(dxD);
    Dataset->SetLengthElementy(dyD);
    Dataset->SetGravity(Gravity);
    Dataset->SetLinearSegregationQ(LinearSegregation);
    Dataset->SetRotationAngle(angle);
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> bottombcini(4,0.0);
    bottombcini[0] = 1;
    bottombcini[1] = 0;
    bottombcini[2] = 0;
    bottombcini[3] = 0;
    
    TPZVec<REAL> rightbcini(4,0.0);
    rightbcini[0] = 1;
    rightbcini[1] = 0.0*(1.0*1e6)/(Pstr);
    rightbcini[2] = 0;
    rightbcini[3] = 0;
    
    TPZVec<REAL> topbcini(4,0.0);
    topbcini[0] = 0;
    topbcini[1] = (1.0*1e6)/(Pstr);
    topbcini[2] = 0;
    topbcini[3] = 0;
    
    TPZVec<REAL> leftbcini(4,0.0);
    leftbcini[0] = 1;
    leftbcini[1] = 0.0*(1.0*1e6)/(Pstr);
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
    rightbc[1] = (1.0*1e6)/(Pstr);
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = (0.0*1e6)/(Pstr);
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = -(1.0e-4)*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    leftbc[2] = (1.0 - S_nw_r);
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
    REAL Top            = 0.0/Lstr;
    REAL Rw             = 0.127/Lstr;
    
    // Reservoir Description linear tracer configuration
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
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
    Kabsolute(0,0) = (1.0e-13)/Kstr;
    Kabsolute(1,1) = (1.0e-13)/Kstr;
    
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
    PVTData[0] = Water.operator->();
    PVTData[1] = Oil.operator->();
    PVTData[2] = Gas.operator->();
    
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    // CFL
    
    REAL u = (1.0e-4)*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    REAL cfl = ((dyD*(u/(waterdensity)))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
    std::cout << "cfl = " << cfl << std::endl;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers,Rocks);
    SandStone.SetFluidData(PVTData);
    SandStone.RunAnalysis();
    
    std::cout << std::endl;
    std::cout << "cfl = " << cfl << std::endl;
    
    
}

void LinearWithReconstruction(bool IsDimensionlessQ){
    
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!
    
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    // Simulation Data SI units
    
    // Characteristic Data
    REAL Kstr           = 1.0;
    REAL Pstr           = 1.0;
    REAL Tstr           = 355.37;
    REAL Tres           = 355.37;
    REAL Lstr           = 1.0;
    REAL Mustr          = 1.0;
    REAL Rhostr         = 1.0;
    TPZMaterial::gBigNumber = 1.0e15; // Use this for check of convergence using neumann
    
    REAL TolRes     = 1.0*1e-4;
    REAL TolDeltaX  = 1.0*1e-6;
    
    if (IsDimensionlessQ)
    {
        Kstr           = 1.0e-13;
        Pstr           = 1.0e7;
        Tstr           = 355.37;
        Tres           = 355.37;
        Lstr           = 100.0;
        Mustr          = 0.001;
        Rhostr         = 1000.0;
        TPZMaterial::gBigNumber = 1.0e12; // Use this for check of convergence using neumann
        TolRes     = 1.0*1e-5;
        TolDeltaX  = 1.0*1e-7;
    }
    
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 40;
    int nthread     = 8;
    bool GR         = true;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;    // Band optimization
    bool IsAxisy    = false;    // Axisymmetric analysis
    bool IsImpes    = false;    // Impes analysis
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 1;
    int hrefinement = 0;
    int hpostref    = 0;
    
    // Time control parameters
    int n_times  = 10;
    int n_sub_dt = 100;
    int which_dt = n_times;
    TPZManVector<REAL,20> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = (1.0e7/day) * day * scale;
    REAL t0         = 0.0  * day * scale;
    
    for (int it = 0 ; it < n_times; it++) {
        Reporting_times[it] = REAL(it+1)*dt;
    }
    REAL maxtime    = Reporting_times[which_dt-1];
    std::cout << "Reporting times = " << Reporting_times << std::endl;
    std::cout << "Maximum simulation time = " << maxtime <<std::endl;
    
    int  nelemX     =1000;
    REAL dxD        =(100.0/nelemX)/Lstr;
    
    int nelemY      =2;
    REAL dyD        =(20.0/nelemY)/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -0.0*((Lstr*Rhostr)/Pstr);
    bool LinearSegregation = true;
    
    REAL angle = 0.0;
    
    REAL S_w_r              = 0.0;
    REAL S_nw_r             = 0.0;
    
    TPZStack<std::string> system;
    system.Push("Water");
    system.Push("Oil");
    
    
    Dataset->SetTimes(Reporting_times);
    Dataset->SetNSubSteps(n_sub_dt);
    Dataset->SetsystemType(system);
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
    Dataset->SetIsDirect(IsDirect);
    Dataset->SetIsCG(IsCG);
    Dataset->SetOptband(OptBand);
    Dataset->SetAxisymmetricQ(IsAxisy);
    Dataset->SetImpesQ(IsImpes);
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
    Dataset->SetnElementsx(nelemX);
    Dataset->SetnElementsy(nelemY);
    Dataset->SetLengthElementx(dxD);
    Dataset->SetLengthElementy(dyD);
    Dataset->SetGravity(Gravity);
    Dataset->SetLinearSegregationQ(LinearSegregation);
    Dataset->SetRotationAngle(angle);
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> bottombcini(4,0.0);
    bottombcini[0] = 1;
    bottombcini[1] = 0;
    bottombcini[2] = 0;
    bottombcini[3] = 0;
    
    TPZVec<REAL> rightbcini(4,0.0);
    rightbcini[0] = 1;
    rightbcini[1] = 0.0*(1.0*1e6)/(Pstr);
    rightbcini[2] = 0;
    rightbcini[3] = 0;
    
    TPZVec<REAL> topbcini(4,0.0);
    topbcini[0] = 0;
    topbcini[1] = (1.0*1e6)/(Pstr);
    topbcini[2] = 0;
    topbcini[3] = 0;
    
    TPZVec<REAL> leftbcini(4,0.0);
    leftbcini[0] = 1;
    leftbcini[1] = 0.0*(1.0*1e6)/(Pstr);
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
    rightbc[1] = (1.0*1e6)/(Pstr);
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = (0.0*1e6)/(Pstr);
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = -(1.0e-4)*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    leftbc[2] = (1.0 - S_nw_r);
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
    REAL Top            = 0.0/Lstr;
    REAL Rw             = 0.127/Lstr;
    
    // Reservoir Description linear tracer configuration
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
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
    Kabsolute(0,0) = (1.0e-13)/Kstr;
    Kabsolute(1,1) = (1.0e-13)/Kstr;
    
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
    PVTData[0] = Water.operator->();
    PVTData[1] = Oil.operator->();
    PVTData[2] = Gas.operator->();
    
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    // CFL
    
    REAL u = (1.0e-4)*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    REAL cfl = ((dyD*(u/(waterdensity)))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
    std::cout << "cfl = " << cfl << std::endl;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers,Rocks);
    SandStone.SetFluidData(PVTData);
    SandStone.RunAnalysis();
    
    std::cout << std::endl;
    std::cout << "cfl = " << cfl << std::endl;

    
}

void NonlinearTracer(bool IsDimensionlessQ)
{
    // Important notes:
    // This code consider a homogeneus absolute permeability when Gravitational segregational function is active!
    
    // This code use piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    // Simulation Data SI units
    
    // Characteristic Data
    REAL Kstr           = 1.0;
    REAL Pstr           = 1.0;
    REAL Tstr           = 355.37;
    REAL Tres           = 355.37;
    REAL Lstr           = 1.0;
    REAL Mustr          = 1.0;
    REAL Rhostr         = 1.0;
    TPZMaterial::gBigNumber = 1.0e15; // Use this for check of convergence using neumann
    
    REAL TolRes     = 1.0*1e-4;
    REAL TolDeltaX  = 1.0*1e-6;
    
    if (IsDimensionlessQ)
    {
        Kstr           = 1.0e-13;
        Pstr           = 1.0e7;
        Tstr           = 355.37;
        Tres           = 355.37;
        Lstr           = 100.0;
        Mustr          = 0.001;
        Rhostr         = 1000.0;
        TPZMaterial::gBigNumber = 1.0e12; // Use this for check of convergence using neumann
        TolRes     = 1.0*1e-2;
        TolDeltaX  = 1.0*1e-3;
    }
    
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 40;
    int nthread     = 0;
    bool GR         = true;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;    // Band optimization
    bool IsAxisy    = false;    // Axisymmetric analysis
    bool IsImpes    = false;    // Impes analysis
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 1;
    int hrefinement = 0;
    int hpostref    = 0;
    
    // Time control parameters
    int n_times  = 40;
    int n_sub_dt = 1;
    int which_dt = n_times;
    TPZManVector<REAL,20> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 1.0*(1.0e7/day) * day * scale;
    REAL t0         = 0.0  * day * scale;

    for (int it = 0 ; it < n_times; it++) {
        Reporting_times[it] = REAL(it+1)*dt;
    }
    REAL maxtime    = Reporting_times[which_dt-1];
    std::cout << "Reporting times = " << Reporting_times << std::endl;
    std::cout << "Maximum simulation time = " << maxtime <<std::endl;
    
    int  nelemX     =2;
    REAL dxD        =(100.0/nelemX)/Lstr;
    
    int nelemY      =10;
    REAL dyD        =(100.0/nelemY)/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -10.0*((Lstr*Rhostr)/Pstr);
    bool LinearSegregation = true;
    
    REAL angle = 0.0;
    
    REAL S_w_r              = 0.0;
    REAL S_nw_r             = 0.0;
    
    TPZStack<std::string> system;
    system.Push("Water");
    system.Push("Oil");


    Dataset->SetTimes(Reporting_times);
    Dataset->SetNSubSteps(n_sub_dt);
    Dataset->SetsystemType(system);
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
    Dataset->SetIsDirect(IsDirect);
    Dataset->SetIsCG(IsCG);
    Dataset->SetOptband(OptBand);
    Dataset->SetAxisymmetricQ(IsAxisy);
    Dataset->SetImpesQ(IsImpes);
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
    Dataset->SetnElementsx(nelemX);
    Dataset->SetnElementsy(nelemY);
    Dataset->SetLengthElementx(dxD);
    Dataset->SetLengthElementy(dyD);
    Dataset->SetGravity(Gravity);
    Dataset->SetLinearSegregationQ(LinearSegregation);
    Dataset->SetRotationAngle(angle);
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> bottombcini(4,0.0);
    bottombcini[0] = 1;
    bottombcini[1] = 0;
    bottombcini[2] = 0;
    bottombcini[3] = 0;
    
    TPZVec<REAL> rightbcini(4,0.0);
    rightbcini[0] = 1;
    rightbcini[1] = 0.0*(1.0*1e6)/(Pstr);
    rightbcini[2] = 0;
    rightbcini[3] = 0;
    
    TPZVec<REAL> topbcini(4,0.0);
    topbcini[0] = 0;
    topbcini[1] = (1.0*1e6)/(Pstr);
    topbcini[2] = 0;
    topbcini[3] = 0;
    
    TPZVec<REAL> leftbcini(4,0.0);
    leftbcini[0] = 1;
    leftbcini[1] = 0.0*(1.0*1e6)/(Pstr);
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
    rightbc[0] = 1;
    rightbc[1] = (0.0*1e6)/(Pstr);
    rightbc[2] = 0;
    rightbc[3] = 0;
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 0;
    topbc[1] = (1.0*1e6)/(Pstr);
    topbc[2] = 0;
    topbc[3] = 0;
    
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = -(0.0e-4)*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    leftbc[2] = 0*(1.0 - S_nw_r);
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
    REAL Top            = 0.0/Lstr;
    REAL Rw             = 0.127/Lstr;
    
    // Reservoir Description linear tracer configuration
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
    REAL oildensity         = 800.0/Rhostr;
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
    Kabsolute(0,0) = (1.0e-13)/Kstr;
    Kabsolute(1,1) = (1.0e-13)/Kstr;
    
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
    PVTData[0] = Water.operator->();
    PVTData[1] = Oil.operator->();
    PVTData[2] = Gas.operator->();
    
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = Layer;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    // CFL
    
    REAL u = (1.0e-4)*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    REAL cfl = ((dyD*(u/(waterdensity)))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
    std::cout << "cfl = " << cfl << std::endl;
    
    TPZDarcyAnalysis SandStone(Dataset,Layers,Rocks);
    SandStone.SetFluidData(PVTData);
    SandStone.RunAnalysis();
    
    std::cout << std::endl;
    std::cout << "cfl = " << cfl << std::endl;
}