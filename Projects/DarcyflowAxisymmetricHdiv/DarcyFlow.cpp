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

#include <fstream>
void LinearTracer();
void NonlinearTracer(bool IsDimensionlessQ);
void TwoUnitsModel(bool IsDimensionlessQ);
void DiffusiveModel(bool IsDimensionlessQ);

int main()
{

    bool IsDimensionlessQ = true;
//    TwoUnitsModel(IsDimensionlessQ);
    DiffusiveModel(IsDimensionlessQ);
//    NonlinearTracer(IsDimensionlessQ);
    
    std::cout << "Program successfully executed." << std::endl;
    return 0;
}

void DiffusiveModel(bool IsDimensionlessQ)
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
    TPZMaterial::gBigNumber = 1.0e20; // Use this for check of convergence using neumann
    
    REAL TolRes     = 1.0*1e-2;
    REAL TolDeltaX  = 1.0*1e-5;
    
    if (IsDimensionlessQ)
    {
        Kstr           = 1.0e-13;
        Pstr           = 2.0e7;
        Tstr           = 355.37;
        Tres           = 355.37;
        Lstr           = 500.0;
        Mustr          = 0.001;
        Rhostr         = 1000.0;
        TPZMaterial::gBigNumber = 1.0e8;
        TolRes     = 1.0*1e-4;
        TolDeltaX  = 1.0*1e-6;
    }
    
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 30;
    int nthread     = 6;
    bool GR         = true;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;     // Band optimization
    bool IsAxisy    = false;     // Axisymmetric analysis 1.0/s;
    bool IsTMesh    = false;    // Triangular mesh
    bool IsImpes    = false;    // Impes analysis
    bool IsHydro    = true;    // Hydrostatic bc
    bool IsPGMesh   = true;     // Geometric Progression mesh
    bool IsHetero   = false;     // Heterogeneous k model
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 1;
    int hrefinement = 0;
    int hpostref    = 0;
    
    // Time control parameters
    int n_times  = 10;
    int n_sub_dt = 10;
    int which_dt = n_times;
    TPZManVector<REAL,20> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 1.0*(1000) * day * scale;
    REAL t0         = 0.0  * day * scale;
    
    for (int it = 0 ; it < n_times; it++) {
        Reporting_times[it] = REAL(it+1)*dt;
    }
    REAL maxtime    = Reporting_times[which_dt-1];
    std::cout << "Reporting times = " << Reporting_times << std::endl;
    std::cout << "Maximum simulation time = " << maxtime <<std::endl;
    
    REAL x_l = 500.0;
    REAL y_l = 100.0;
    REAL ratio = 0.95;
    
    int  nelemX     =20;
    if (GR && nelemX == 1 && IsTMesh) {
        nelemX++;
    }
    REAL dxD        =(x_l/nelemX)/Lstr;
    
    int nelemY      =8;
    if (GR && nelemY == 1 && IsTMesh ) {
        nelemY++;
    }
    REAL dyD        =(y_l/nelemY)/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -10.0*((Lstr*Rhostr)/Pstr);
    bool LinearSegregation = false;
    
    REAL angle = 0.0;
    
    REAL S_w_r              = 0.2;
    REAL S_nw_r             = 0.2;
    
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
    Dataset->SetTriangularMesh(IsTMesh);
    Dataset->SetMeshwithPGQ(IsPGMesh);
    Dataset->SetHeterogeneousQ(IsHetero);
    Dataset->SetPGRatio(ratio);
    Dataset->SetHydrostaticBCQ(IsHydro);
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
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> HydraulicUnit1         = new ReservoirData;
    TPZAutoPointer<ReservoirData> HydraulicUnit2         = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel   = new PetroPhysicData;
    TPZAutoPointer<OilPhase> Oil          = new OilPhase;   // alpha
    TPZAutoPointer<WaterPhase> Water        = new WaterPhase;   // beta
    TPZAutoPointer<GasPhase> Gas          = new GasPhase;   // gamma
    
    /////////////////////////////////////////////////
    // Complete data set for HydraulicUnit1
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> top(4,0.0);
    top[0] = 2;
    top[1] = (1.0*1e7)/(Pstr);
    top[2] = 0;
    top[3] = 0;
    
    TPZVec<REAL> nonflux(4,0.0);
    nonflux[0] = 1;
    nonflux[1] = 0;
    nonflux[2] = 0;
    nonflux[3] = 0;
    
    TPZStack< TPZVec<REAL> > initial_bcs;
    initial_bcs.Push(nonflux);
    initial_bcs.Push(nonflux);
    initial_bcs.Push(top);
    initial_bcs.Push(nonflux);    
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    //  typeImpervious = 4;
    
    //  Boundary Value Problem
    
    
    TPZVec<REAL> impervious(4,0.0);
    impervious[0] = 1;
    impervious[1] = 0;
    impervious[2] = 0;
    impervious[3] = 0;
    
    TPZVec<REAL> inlet(4,0.0);
    inlet[0] = 1;
    inlet[1] = -0.01;
    inlet[2] = 0.8;
    inlet[3] = 0;
    
    TPZVec<REAL> outlet(4,0.0);
    outlet[0] = 2;
    outlet[1] = (1.0*1e7)/(Pstr);
    outlet[2] = 0;
    outlet[3] = 0;
    
    TPZStack< TPZVec<REAL> > bcs;
    bcs.Push(impervious);
    bcs.Push(outlet);
    bcs.Push(impervious);
    bcs.Push(inlet);

    
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
    REAL Rw             = 1.0*0.127/Lstr;
    
    // Reservoir Description configuration
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
    REAL oildensity         = 800.0/Rhostr;
    REAL oilviscosity       = 0.001/Mustr;
    REAL coil               = (0.0*1.0*1e-8)*Pstr;
    REAL p_g_ref            = Pstr;
    REAL gasdensity         = Rhostr;
    REAL gasviscosity       = Mustr;
    REAL cgas               = (0.0)*Pstr;
    
    REAL pc_max             = (0.0*1e5)/(Pstr);
    
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
    
    HydraulicUnit1->SetIsGIDGeometry(isGIDGeom);
    HydraulicUnit1->SetBC(initial_bcs, bcs);
    HydraulicUnit1->SetLayerTop(Top);
    HydraulicUnit1->SetLayerrw(Rw);
    HydraulicUnit1->SetLayerh(Hres);
    HydraulicUnit1->SetLayerr(Rres);
    HydraulicUnit1->SetLref(lengthref);
    HydraulicUnit1->SetKref(kref);
    HydraulicUnit1->SetPhiRef(porosityref);
    HydraulicUnit1->SetPref(pressureref);
    HydraulicUnit1->SetcRock(crock);
    HydraulicUnit1->SetKabsolute(Kabsolute);
    HydraulicUnit1->SetMatIDs(MatIds);
    HydraulicUnit1->SetS_wett_r(S_w_r);
    HydraulicUnit1->SetS_nwett_r(S_nw_r);
    
    
    Water->SetRho(waterdensity);
    Water->SetMu(waterviscosity);
    Water->Setc(cwater);
    Water->SetPRef(p_w_ref);
    Water->SetTRef(Tstr);
    Water->SetTRes(Tres);
    Water->SetS_wett_r(S_w_r);
    Water->SetS_nwett_r(S_nw_r);
    Water->SetPc_max(pc_max);
    
    Oil->SetRho(oildensity);
    Oil->SetMu(oilviscosity);
    Oil->Setc(coil);
    Oil->SetPRef(p_o_ref);
    Oil->SetTRef(Tstr);
    Oil->SetTRes(Tres);
    Oil->SetS_wett_r(S_w_r);
    Oil->SetS_nwett_r(S_nw_r);
    Oil->SetPc_max(pc_max);
    
    Gas->SetRho(gasdensity);
    Gas->SetMu(gasviscosity);
    Gas->Setc(cgas);
    Gas->SetPRef(p_g_ref);
    Gas->SetTRef(Tstr);
    Gas->SetTRes(Tres);
    Gas->SetS_wett_r(S_w_r);
    Gas->SetS_nwett_r(S_nw_r);
    Gas->SetPc_max(pc_max);
    
    TPZVec< TPZAutoPointer<Phase> > PVTData(3);
    PVTData[0] = Water.operator->();
    PVTData[1] = Oil.operator->();
    PVTData[2] = Gas.operator->();
    
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = HydraulicUnit1;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
    TPZDarcyAnalysis *SandStone = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
    SandStone->SetFluidData(PVTData);
    SandStone->RunAnalysis();
    delete SandStone;
    
    std::cout << std::endl;
    
}


void TwoUnitsModel(bool IsDimensionlessQ)
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
    TPZMaterial::gBigNumber = 1.0e20; // Use this for check of convergence using neumann
    
    REAL TolRes     = 1.0*1e-2;
    REAL TolDeltaX  = 1.0*1e-5;
    
    if (IsDimensionlessQ)
    {
        Kstr           = 1.0e-13;
        Pstr           = 2.0e7;
        Tstr           = 355.37;
        Tres           = 355.37;
        Lstr           = 500.0;
        Mustr          = 0.001;
        Rhostr         = 1000.0;
        TPZMaterial::gBigNumber = 1.0e8;
        TolRes     = 1.0*1e-4;
        TolDeltaX  = 1.0*1e-6;
    }
    
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 30;
    int nthread     = 0;
    bool GR         = true;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;     // Band optimization
    bool IsAxisy    = false;     // Axisymmetric analysis 1.0/s;
    bool IsTMesh    = false;    // Triangular mesh
    bool IsImpes    = false;    // Impes analysis
    bool IsHydro    = false;    // Hydrostatic bc
    bool IsPGMesh   = true;     // Geometric Progression mesh
    bool IsHetero   = false;     // Heterogeneous k model
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 1;
    int hrefinement = 0;
    int hpostref    = 0;
    
    // Time control parameters
    int n_times  = 10;
    int n_sub_dt = 10;
    int which_dt = n_times;
    TPZManVector<REAL,20> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 1.0*(10) * day * scale;
    REAL t0         = 0.0  * day * scale;
    
    for (int it = 0 ; it < n_times; it++) {
        Reporting_times[it] = REAL(it+1)*dt;
    }
    REAL maxtime    = Reporting_times[which_dt-1];
    std::cout << "Reporting times = " << Reporting_times << std::endl;
    std::cout << "Maximum simulation time = " << maxtime <<std::endl;
    
    REAL x_l = 500.0;
    REAL y_l = 100.0;
    REAL ratio = 0.95;
    
    int  nelemX     =31;
    if (GR && nelemX == 1 && IsTMesh) {
        nelemX++;
    }
    REAL dxD        =(x_l/nelemX)/Lstr;
    
    int nelemY      =15;
    if (GR && nelemY == 1 && IsTMesh ) {
        nelemY++;
    }
    REAL dyD        =(y_l/nelemY)/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -0.0*((Lstr*Rhostr)/Pstr);
    bool LinearSegregation = false;
    
    REAL angle =
    0.0;
    
    REAL S_w_r              = 0.2;
    REAL S_nw_r             = 0.2;
    
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
    Dataset->SetTriangularMesh(IsTMesh);
    Dataset->SetMeshwithPGQ(IsPGMesh);
    Dataset->SetHeterogeneousQ(IsHetero);
    Dataset->SetPGRatio(ratio);
    Dataset->SetHydrostaticBCQ(IsHydro);
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
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> HydraulicUnit1         = new ReservoirData;
    TPZAutoPointer<ReservoirData> HydraulicUnit2         = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel   = new PetroPhysicData;
    TPZAutoPointer<OilPhase> Oil          = new OilPhase;   // alpha
    TPZAutoPointer<WaterPhase> Water        = new WaterPhase;   // beta
    TPZAutoPointer<GasPhase> Gas          = new GasPhase;   // gamma
    
    /////////////////////////////////////////////////
    // Complete data set for HydraulicUnit1
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> top(4,0.0);
    top[0] = 2;
    top[1] = (1.0*1e7)/(Pstr);
    top[2] = 0;
    top[3] = 0;
    
    TPZVec<REAL> nonflux(4,0.0);
    nonflux[0] = 1;
    nonflux[1] = 0;
    nonflux[2] = 0;
    nonflux[3] = 0;
    
    TPZStack< TPZVec<REAL> > initial_bcs;
    initial_bcs.Push(top);
    initial_bcs.Push(nonflux);
    initial_bcs.Push(nonflux);
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    //  typeImpervious = 4;
    
    //  Boundary Value Problem
    
    
    TPZVec<REAL> impervious(4,0.0);
    impervious[0] = 1;
    impervious[1] = 0;
    impervious[2] = 0;
    impervious[3] = 0;
    
    TPZVec<REAL> inlet(4,0.0);
    inlet[0] = 1;
    inlet[1] = -10.0;
    inlet[2] = 0.8;
    inlet[3] = 0;
    
    TPZVec<REAL> outlet(4,0.0);
    outlet[0] = 2;
    outlet[1] = (1.0*1e6)/(Pstr);
    outlet[2] = 0;
    outlet[3] = 0;
    
    TPZStack< TPZVec<REAL> > bcs;
    bcs.Push(impervious);
    bcs.Push(inlet);
    bcs.Push(outlet);
    
    // Reservoir Description
    bool isGIDGeom      = true;
    REAL porosityref    = 0.2;
    REAL pressureref    = (1.0*1e6)/(Pstr);
    REAL lengthref      = 1.0;
    REAL kref           = 1.0;
    REAL crock          = (0.0*1e-10)*Pstr;
    REAL Hres           = 100.0/Lstr;
    REAL Rres           = 1000.0/Lstr;
    REAL Top            = 0.0/Lstr;
    REAL Rw             = 1.0*0.127/Lstr;
    
    // Reservoir Description configuration
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
    REAL oildensity         = 1000.0/Rhostr;
    REAL oilviscosity       = 0.001/Mustr;
    REAL coil               = (0.0*1.0*1e-8)*Pstr;
    REAL p_g_ref            = Pstr;
    REAL gasdensity         = Rhostr;
    REAL gasviscosity       = Mustr;
    REAL cgas               = (0.0)*Pstr;
    
    REAL pc_max             = (0.0*1e5)/(Pstr);
    
    TPZVec<int> MatIds(4);
    MatIds[0]=1;
    MatIds[1]=3;
    MatIds[2]=5;
    MatIds[3]=7;
    
    TPZFMatrix<STATE> Kabsolute(2,2);
    Kabsolute.Zero();
    Kabsolute(0,0) = (1.0e-13)/Kstr;
    Kabsolute(1,1) = (1.0e-13)/Kstr;
    
    HydraulicUnit1->SetIsGIDGeometry(isGIDGeom);
    HydraulicUnit1->SetBC(initial_bcs, bcs);
    HydraulicUnit1->SetLayerTop(Top);
    HydraulicUnit1->SetLayerrw(Rw);
    HydraulicUnit1->SetLayerh(Hres);
    HydraulicUnit1->SetLayerr(Rres);
    HydraulicUnit1->SetLref(lengthref);
    HydraulicUnit1->SetKref(kref);
    HydraulicUnit1->SetPhiRef(porosityref);
    HydraulicUnit1->SetPref(pressureref);
    HydraulicUnit1->SetcRock(crock);
    HydraulicUnit1->SetKabsolute(Kabsolute);
    HydraulicUnit1->SetMatIDs(MatIds);
    HydraulicUnit1->SetS_wett_r(S_w_r);
    HydraulicUnit1->SetS_nwett_r(S_nw_r);
    
    /////////////////////////////////////////////////
    // Complete data set for HydraulicUnit2
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    
    //  Initial Boundary Value Problem
    
    TPZVec<REAL> nonflux2(4,0.0);
    nonflux2[0] = 1;
    nonflux2[1] = 0;
    nonflux2[2] = 0;
    nonflux2[3] = 0;
    
    TPZStack< TPZVec<REAL> > initial_bcs2;
    initial_bcs2.Push(nonflux2);
    initial_bcs2.Push(nonflux2);
    initial_bcs2.Push(nonflux2);
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    //  typeImpervious = 4;
    
    //  Boundary Value Problem
    
    
    TPZVec<REAL> impervious2(4,0.0);
    impervious2[0] = 1;
    impervious2[1] = 0;
    impervious2[2] = 0;
    impervious2[3] = 0;
    
    TPZVec<REAL> inlet2(4,0.0);
    inlet2[0] = 1;
    inlet2[1] = -10.0;
    inlet2[2] = 0.8;
    inlet2[3] = 0;
    
    TPZVec<REAL> outlet2(4,0.0);
    outlet2[0] = 2;
    outlet2[1] = (1.0*1e6)/(Pstr);;
    outlet2[2] = 0;
    outlet2[3] = 0;
    
    TPZStack< TPZVec<REAL> > bcs2;
    bcs2.Push(impervious2);
    bcs2.Push(inlet2);
    bcs2.Push(outlet2);
    
    // Reservoir Description
    isGIDGeom      = true;
    porosityref    = 0.2;
    pressureref    = (1.0*1e6)/(Pstr);
    lengthref      = 1.0;
    kref           = 1.0;
    crock          = (0.0*1e-10)*Pstr;
    Hres           = 100.0/Lstr;
    Rres           = 1000.0/Lstr;
    Top            = 0.0/Lstr;
    Rw             = 1.0*0.127/Lstr;
    
    // Reservoir Description configuration
    p_w_ref            = (1.0*1e6)/(Pstr);
    waterdensity       = 1000.0/Rhostr;
    waterviscosity     = 0.001/Mustr;
    cwater             = (0.0*1.0*1e-10)*Pstr;
    p_o_ref            = (1.0*1e6)/(Pstr);
    oildensity         = 1000.0/Rhostr;
    oilviscosity       = 0.001/Mustr;
    coil               = (0.0*1.0*1e-8)*Pstr;
    p_g_ref            = Pstr;
    gasdensity         = Rhostr;
    gasviscosity       = Mustr;
    cgas               = (0.0)*Pstr;
    
    pc_max             = (0.0*1e5)/(Pstr);
    
    MatIds.Resize(4, 0);
    
    MatIds[0]=2;
    MatIds[1]=4;
    MatIds[2]=6;
    MatIds[3]=8;
    
    Kabsolute(0,0) = (1.0e-13)/Kstr;
    Kabsolute(1,1) = (1.0e-13)/Kstr;
    
    HydraulicUnit2->SetIsGIDGeometry(isGIDGeom);
    HydraulicUnit2->SetBC(initial_bcs2, bcs2);
    HydraulicUnit2->SetLayerTop(Top);
    HydraulicUnit2->SetLayerrw(Rw);
    HydraulicUnit2->SetLayerh(Hres);
    HydraulicUnit2->SetLayerr(Rres);
    HydraulicUnit2->SetLref(lengthref);
    HydraulicUnit2->SetKref(kref);
    HydraulicUnit2->SetPhiRef(porosityref);
    HydraulicUnit2->SetPref(pressureref);
    HydraulicUnit2->SetcRock(crock);
    HydraulicUnit2->SetKabsolute(Kabsolute);
    HydraulicUnit2->SetMatIDs(MatIds);
    HydraulicUnit2->SetS_wett_r(S_w_r);
    HydraulicUnit2->SetS_nwett_r(S_nw_r);
    
    
    Water->SetRho(waterdensity);
    Water->SetMu(waterviscosity);
    Water->Setc(cwater);
    Water->SetPRef(p_w_ref);
    Water->SetTRef(Tstr);
    Water->SetTRes(Tres);
    Water->SetS_wett_r(S_w_r);
    Water->SetS_nwett_r(S_nw_r);
    Water->SetPc_max(pc_max);
    
    Oil->SetRho(oildensity);
    Oil->SetMu(oilviscosity);
    Oil->Setc(coil);
    Oil->SetPRef(p_o_ref);
    Oil->SetTRef(Tstr);
    Oil->SetTRes(Tres);
    Oil->SetS_wett_r(S_w_r);
    Oil->SetS_nwett_r(S_nw_r);
    Oil->SetPc_max(pc_max);
    
    Gas->SetRho(gasdensity);
    Gas->SetMu(gasviscosity);
    Gas->Setc(cgas);
    Gas->SetPRef(p_g_ref);
    Gas->SetTRef(Tstr);
    Gas->SetTRes(Tres);
    Gas->SetS_wett_r(S_w_r);
    Gas->SetS_nwett_r(S_nw_r);
    Gas->SetPc_max(pc_max);
    
    TPZVec< TPZAutoPointer<Phase> > PVTData(3);
    PVTData[0] = Water.operator->();
    PVTData[1] = Oil.operator->();
    PVTData[2] = Gas.operator->();
    
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(2);
    Layers[0] = HydraulicUnit1;
    Layers[1] = HydraulicUnit2;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;

    TPZDarcyAnalysis *SandStone = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
    SandStone->SetFluidData(PVTData);
    SandStone->RunAnalysis();
    delete SandStone;

    std::cout << std::endl;
    
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
    TPZMaterial::gBigNumber = 1.0e20; // Use this for check of convergence using neumann
    
    REAL TolRes     = 1.0*1e-2;
    REAL TolDeltaX  = 1.0*1e-5;
    
    if (IsDimensionlessQ)
    {
        Kstr           = 1.0e-13;
        Pstr           = 2.0e7;
        Tstr           = 355.37;
        Tres           = 355.37;
        Lstr           = 500.0;
        Mustr          = 0.001;
        Rhostr         = 1000.0;
        TPZMaterial::gBigNumber = 1.0e8;
        TolRes     = 1.0*1e-4;
        TolDeltaX  = 1.0*1e-6;
    }
    
    TPZFMatrix<REAL> Gravity(2,1);
    
    TPZAutoPointer<SimulationData> Dataset  = new SimulationData;
    
    int maxiter     = 30;
    int nthread     = 8;
    bool GR         = false;    // Use Gradient Reconstruction
    bool SC         = false;    // Use Static Condensation not working for nonlinear and transient problems
    bool IsDirect   = true;     // No Use broyden with Iterative !!!
    bool IsCG       = false;    // false means GMRES
    bool OptBand    = true;     // Band optimization
    bool IsAxisy    = false;     // Axisymmetric analysis 1.0/s;
    bool IsTMesh    = false;    // Triangular mesh
    bool IsImpes    = false;    // Impes analysis
    bool IsHydro    = false;    // Hydrostatic bc
    bool IsPGMesh   = true;     // Geometric Progression mesh
    bool IsHetero   = false;     // Heterogeneous k model
    int fixedJac    = 0;
    
    int qorder      = 1;
    int porder      = 1;
    int sorder      = 0;
    int hrefinement = 0;
    int hpostref    = 0;
    
    // Time control parameters
    int n_times  = 10;
    int n_sub_dt = 20;
    int which_dt = n_times;
    TPZManVector<REAL,20> Reporting_times(n_times,0.0);
    REAL scale = ((Kstr*Pstr)/(Lstr*Lstr*Mustr));
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    REAL dt         = 1.0*(500) * day * scale;
    REAL t0         = 0.0  * day * scale;
    
    for (int it = 0 ; it < n_times; it++) {
        Reporting_times[it] = REAL(it+1)*dt;
    }
    REAL maxtime    = Reporting_times[which_dt-1];
    std::cout << "Reporting times = " << Reporting_times << std::endl;
    std::cout << "Maximum simulation time = " << maxtime <<std::endl;
    
    REAL x_l = 500.0;
    REAL y_l = 100.0;
    REAL ratio = 0.95;
    
    int  nelemX     =31;
    if (GR && nelemX == 1 && IsTMesh) {
        nelemX++;
    }
    REAL dxD        =(x_l/nelemX)/Lstr;
    
    int nelemY      =15;
    if (GR && nelemY == 1 && IsTMesh ) {
        nelemY++;
    }
    REAL dyD        =(y_l/nelemY)/Lstr;
    
    Gravity(0,0)= -0.0*((Lstr*Rhostr)/Pstr);
    Gravity(1,0)= -0.0*((Lstr*Rhostr)/Pstr);
    bool LinearSegregation = false;
    
    REAL angle =
    0.0;
    
    REAL S_w_r              = 0.2;
    REAL S_nw_r             = 0.2;
    
    TPZStack<std::string> system;
    system.Push("Water");
//    system.Push("Oil");
    
    
    Dataset->SetTimes(Reporting_times);
    Dataset->SetNSubSteps(n_sub_dt);
    Dataset->SetsystemType(system);
    Dataset->SetGR(GR);
    Dataset->SetSC(SC);
    Dataset->SetIsDirect(IsDirect);
    Dataset->SetIsCG(IsCG);
    Dataset->SetOptband(OptBand);
    Dataset->SetAxisymmetricQ(IsAxisy);
    Dataset->SetTriangularMesh(IsTMesh);
    Dataset->SetMeshwithPGQ(IsPGMesh);
    Dataset->SetHeterogeneousQ(IsHetero);
    Dataset->SetPGRatio(ratio);
    Dataset->SetHydrostaticBCQ(IsHydro);
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
    rightbcini[1] = 0.0*(1.0*1e7)/(Pstr);
    rightbcini[2] = 0;
    rightbcini[3] = 0;
    
    TPZVec<REAL> topbcini(4,0.0);
    topbcini[0] = 2;
    topbcini[1] = (1.0*1e7)/(Pstr);
    topbcini[2] = 0;
    topbcini[3] = 0;
    
    TPZVec<REAL> leftbcini(4,0.0);
    leftbcini[0] = 1;
    leftbcini[1] = 0.0*(2.0*1e7)/(Pstr);
    leftbcini[2] = 0;
    leftbcini[3] = 0;
    
    TPZStack< TPZVec<REAL> > initial_bcs;
    initial_bcs.Push(bottombcini);
    initial_bcs.Push(rightbcini);
    initial_bcs.Push(topbcini);
    initial_bcs.Push(leftbcini);
    
    // BCs
    //    int typeFluxin = 1, typePressurein = 0;
    //    int typeFluxout = 3, typePressureout = 2;
    //  typeImpervious = 4;
    
    //  Boundary Value Problem
    
    
    TPZVec<REAL> bottombc(4,0.0);
    bottombc[0] = 1;
    bottombc[1] = (0.0*1e6)/(Pstr);
    bottombc[2] = 0;
    bottombc[3] = 0;
    
    REAL Q = -0.397475/day;//-158.99/day;//158.99/day;
    REAL u = Q/(2.0*M_PI*0.127*y_l);
    REAL rho = 1000.0;
    REAL mu = 0.001;
    REAL k = 1.0e-13;
    REAL muD = mu/mu;
    REAL rhoD = rho/rho;
    REAL kD = k/k;
    
    REAL m = rho * u ;

    REAL mD = m*(Lstr*Mustr/(Kstr*Pstr*Rhostr));
    REAL pr = 0.0;//(2.0*1e7)/(Pstr);
    TPZVec<REAL> rightbc(4,0.0);
    rightbc[0] = 2;
    rightbc[1] = (2.0*1e7)/(Pstr);//0.0;//sin(M_PI*pr);//1.0*log(pr)+0.0*pow(pr,4.0);
    rightbc[2] = 0.0*(1.0 - S_nw_r);
    rightbc[3] = 0;
    
    
    TPZVec<REAL> topbc(4,0.0);
    topbc[0] = 1;
    topbc[1] = (0.0*1e6)/(Pstr);
    topbc[2] = 0;
    topbc[3] = 0;
    
    REAL pl = 0.0;
    TPZVec<REAL> leftbc(4,0.0);
    leftbc[0] = 1;
    leftbc[1] = mD;//sin(M_PI*pl);//;1.0*log(pl)+0.0*pow(pl,4.0);//(1.0*1e7)/(Pstr);
    leftbc[2] = 1.0*(1.0 - S_nw_r);
    leftbc[3] = 0;
    
    TPZStack< TPZVec<REAL> > bcs;
    bcs.Push(bottombc);
    bcs.Push(rightbc);
    bcs.Push(topbc);
    bcs.Push(leftbc);
    
//    Dataset->SetBC(initial_bcs, bcs);
    
    // Reservoir Data SI units
    
    TPZAutoPointer<ReservoirData> HydraulicUnit1         = new ReservoirData;
    TPZAutoPointer<ReservoirData> HydraulicUnit2         = new ReservoirData;
    TPZAutoPointer<PetroPhysicData> RockModel   = new PetroPhysicData;
    TPZAutoPointer<OilPhase> Oil          = new OilPhase;   // alpha
    TPZAutoPointer<WaterPhase> Water        = new WaterPhase;   // beta
    TPZAutoPointer<GasPhase> Gas          = new GasPhase;   // gamma
    
    /////////////////////////////////////////////////
    // Complete data set for HydraulicUnit1
    
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
    REAL Rw             = 1.0*0.127/Lstr;
    
    // Reservoir Description configuration
    REAL p_w_ref            = (1.0*1e6)/(Pstr);
    REAL waterdensity       = 1000.0/Rhostr;
    REAL waterviscosity     = 0.001/Mustr;
    REAL cwater             = (0.0*1.0*1e-10)*Pstr;
    REAL p_o_ref            = (1.0*1e6)/(Pstr);
    REAL oildensity         = 800.0/Rhostr;
    REAL oilviscosity       = 0.001/Mustr;
    REAL coil               = (1.0*1.0*1e-8)*Pstr;
    REAL p_g_ref            = Pstr;
    REAL gasdensity         = Rhostr;
    REAL gasviscosity       = Mustr;
    REAL cgas               = (0.0)*Pstr;
    
    REAL pc_max             = (0.0*1e5)/(Pstr);
    
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
    
    HydraulicUnit1->SetIsGIDGeometry(isGIDGeom);
    HydraulicUnit1->SetLayerTop(Top);
    HydraulicUnit1->SetLayerrw(Rw);
    HydraulicUnit1->SetLayerh(Hres);
    HydraulicUnit1->SetLayerr(Rres);
    HydraulicUnit1->SetLref(lengthref);
    HydraulicUnit1->SetKref(kref);
    HydraulicUnit1->SetPhiRef(porosityref);
    HydraulicUnit1->SetPref(pressureref);
    HydraulicUnit1->SetcRock(crock);
    HydraulicUnit1->SetKabsolute(Kabsolute);
    HydraulicUnit1->SetMatIDs(MatIds);
    HydraulicUnit1->SetS_wett_r(S_w_r);
    HydraulicUnit1->SetS_nwett_r(S_nw_r);
    
    /////////////////////////////////////////////////
    // Complete data set for HydraulicUnit1
    
    // Reservoir Description
    isGIDGeom      = false;
    porosityref    = 0.2;
    pressureref    = (1.0*1e6)/(Pstr);
    lengthref      = 1.0;
    kref           = 1.0;
    crock          = (0.0*1e-10)*Pstr;
    Hres           = 100.0/Lstr;
    Rres           = 1000.0/Lstr;
    Top            = 0.0/Lstr;
    Rw             = 1.0*0.127/Lstr;
    
    // Reservoir Description configuration
    p_w_ref            = (1.0*1e6)/(Pstr);
    waterdensity       = 1000.0/Rhostr;
    waterviscosity     = 0.001/Mustr;
    cwater             = (0.0*1.0*1e-10)*Pstr;
    p_o_ref            = (1.0*1e6)/(Pstr);
    oildensity         = 800.0/Rhostr;
    oilviscosity       = 0.001/Mustr;
    coil               = (1.0*1.0*1e-8)*Pstr;
    p_g_ref            = Pstr;
    gasdensity         = Rhostr;
    gasviscosity       = Mustr;
    cgas               = (0.0)*Pstr;
    
    pc_max             = (0.0*1e5)/(Pstr);
    
    MatIds[0]=1;
    MatIds[1]=2;
    MatIds[2]=3;
    MatIds[3]=4;
    MatIds[4]=5;
    
    Kabsolute(0,0) = (1.0e-13)/Kstr;
    Kabsolute(1,1) = (1.0e-13)/Kstr;
    
    HydraulicUnit2->SetIsGIDGeometry(isGIDGeom);
    HydraulicUnit2->SetLayerTop(Top);
    HydraulicUnit2->SetLayerrw(Rw);
    HydraulicUnit2->SetLayerh(Hres);
    HydraulicUnit2->SetLayerr(Rres);
    HydraulicUnit2->SetLref(lengthref);
    HydraulicUnit2->SetKref(kref);
    HydraulicUnit2->SetPhiRef(porosityref);
    HydraulicUnit2->SetPref(pressureref);
    HydraulicUnit2->SetcRock(crock);
    HydraulicUnit2->SetKabsolute(Kabsolute);
    HydraulicUnit2->SetMatIDs(MatIds);
    HydraulicUnit2->SetS_wett_r(S_w_r);
    HydraulicUnit2->SetS_nwett_r(S_nw_r);
    
    
    Water->SetRho(waterdensity);
    Water->SetMu(waterviscosity);
    Water->Setc(cwater);
    Water->SetPRef(p_w_ref);
    Water->SetTRef(Tstr);
    Water->SetTRes(Tres);
    Water->SetS_wett_r(S_w_r);
    Water->SetS_nwett_r(S_nw_r);
    Water->SetPc_max(pc_max);
    
    Oil->SetRho(oildensity);
    Oil->SetMu(oilviscosity);
    Oil->Setc(coil);
    Oil->SetPRef(p_o_ref);
    Oil->SetTRef(Tstr);
    Oil->SetTRes(Tres);
    Oil->SetS_wett_r(S_w_r);
    Oil->SetS_nwett_r(S_nw_r);
    Oil->SetPc_max(pc_max);
    
    Gas->SetRho(gasdensity);
    Gas->SetMu(gasviscosity);
    Gas->Setc(cgas);
    Gas->SetPRef(p_g_ref);
    Gas->SetTRef(Tstr);
    Gas->SetTRes(Tres);
    Gas->SetS_wett_r(S_w_r);
    Gas->SetS_nwett_r(S_nw_r);
    Gas->SetPc_max(pc_max);
    
    TPZVec< TPZAutoPointer<Phase> > PVTData(3);
    PVTData[0] = Water.operator->();
    PVTData[1] = Oil.operator->();
    PVTData[2] = Gas.operator->();
    
    
    // Creating the analysis
    TPZVec<TPZAutoPointer<ReservoirData> > Layers;
    Layers.Resize(1);
    Layers[0] = HydraulicUnit1;
//    Layers[1] = HydraulicUnit2;
    
    TPZVec<TPZAutoPointer<PetroPhysicData> > Rocks;
    Rocks.Resize(1);
    Rocks[0] = RockModel;
    
//  Computing the approximation rate for each refinement for  the order
    TPZFMatrix<REAL> rates(5,5,0.0);
    TPZManVector<int,5> el_sizesx(5,0);
    TPZManVector<int,5> el_sizesy(5,0);

    el_sizesx[0] = nelemX;
    el_sizesx[1] = nelemX*2;
    el_sizesx[2] = nelemX*4;
    el_sizesx[3] = nelemX*8;
    el_sizesx[4] = nelemX*16;
    
    el_sizesy[0] = nelemY;
    el_sizesy[1] = nelemY*2;
    el_sizesy[2] = nelemY*4;
    el_sizesy[3] = nelemY*8;
    el_sizesy[4] = nelemY*16;
    
//    qorder = 1;
//    porder = 1;
//    sorder = 0;
//    GR = false;
    REAL cfl = 0.0;
    mD = fabs(mD);
    
//    Dataset->SetGR(GR);
//    Dataset->Setqorder(qorder);
//    Dataset->Setporder(porder);
//    Dataset->Setsorder(sorder);
    
    // CFL dimensionless
    cfl = ((dyD*(mD))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
    std::cout << "cfl =" << cfl << std::endl;
    TPZDarcyAnalysis *SandStone = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
    SandStone->SetFluidData(PVTData);
    SandStone->RunAnalysis();
    rates(0,0) = dxD;
    rates(0,1) = SandStone->fL2_norm[0];
    rates(0,2) = SandStone->fl2_norm_flux[0];
    rates(0,3) = SandStone->fL2_norm_s[0];
    rates(0,4) = fabs(cfl);
    delete SandStone;

    std::cout << "cfl =" << cfl << std::endl;    
    return;
    
//    nelemX = el_sizesx[1];
//    dxD        =(x_l/nelemX)/Lstr;
//    nelemY = el_sizesy[1];
//    dyD        =(y_l/nelemY)/Lstr;
//    Dataset->SetnElementsx(nelemX);
//    Dataset->SetLengthElementx(dxD);
//    Dataset->SetnElementsy(nelemY);
//    Dataset->SetLengthElementy(dyD);
//    Dataset->SetNSubSteps(n_sub_dt);
//    Dataset->SetDeltaT(dt);
//    Dataset->SetMaxTime(maxtime);
//    Dataset->SetTime(t0);
//    
//    // CFL dimensionless
//    cfl = ((dyD*(mD))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
//    
//    TPZDarcyAnalysis *SandStone2 = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
//    SandStone2->SetFluidData(PVTData);
//    SandStone2->RunAnalysis();
//    rates(1,0) = dxD;
//    rates(1,1) = SandStone2->fL2_norm[0];
//    rates(1,2) = SandStone2->fl2_norm_flux[0];
//    rates(1,3) = SandStone2->fL2_norm_s[0];
//    rates(1,4) = fabs(cfl);
//    delete SandStone2;
//
//    nelemX = el_sizesx[2];
//    dxD        =(x_l/nelemX)/Lstr;
//    nelemY = el_sizesy[2];
//    dyD        =(y_l/nelemY)/Lstr;
//    Dataset->SetnElementsx(nelemX);
//    Dataset->SetLengthElementx(dxD);
//    Dataset->SetnElementsy(nelemY);
//    Dataset->SetLengthElementy(dyD);
//    Dataset->SetDeltaT(dt);
//    Dataset->SetMaxTime(maxtime);
//    Dataset->SetTime(t0);
//    
//    // CFL dimensionless
//    cfl = ((dyD*(mD))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
//    
//    TPZDarcyAnalysis *SandStone3 = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
//    SandStone3->SetFluidData(PVTData);
//    SandStone3->RunAnalysis();
//    rates(2,0) = dxD;
//    rates(2,1) = SandStone3->fL2_norm[0];
//    rates(2,2) = SandStone3->fl2_norm_flux[0];
//    rates(2,3) = SandStone3->fL2_norm_s[0];
//    rates(2,4) = fabs(cfl);
//    delete SandStone3;
//
//    
//    nelemX = el_sizesx[3];
//    dxD        =(x_l/nelemX)/Lstr;
//    nelemY = el_sizesy[3];
//    dyD        =(y_l/nelemY)/Lstr;
//    Dataset->SetnElementsx(nelemX);
//    Dataset->SetLengthElementx(dxD);
//    Dataset->SetnElementsy(nelemY);
//    Dataset->SetLengthElementy(dyD);
//    Dataset->SetDeltaT(dt);
//    Dataset->SetMaxTime(maxtime);
//    Dataset->SetTime(t0);
//    
//    // CFL dimensionless
//    cfl = ((dyD*(mD))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
//    
//    TPZDarcyAnalysis *SandStone4 = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
//    SandStone4->SetFluidData(PVTData);
//    SandStone4->RunAnalysis();
//    rates(3,0) = dxD;
//    rates(3,1) = SandStone4->fL2_norm[0];
//    rates(3,2) = SandStone4->fl2_norm_flux[0];
//    rates(3,3) = SandStone4->fL2_norm_s[0];
//    rates(3,4) = fabs(cfl);
//    delete SandStone4;
//    
//    nelemX = el_sizesx[4];
//    dxD        =(x_l/nelemX)/Lstr;
//    nelemY = el_sizesy[4];
//    dyD        =(y_l/nelemY)/Lstr;
//    Dataset->SetnElementsx(nelemX);
//    Dataset->SetLengthElementx(dxD);
//    Dataset->SetnElementsy(nelemY);
//    Dataset->SetLengthElementy(dyD);
//    Dataset->SetDeltaT(dt);
//    Dataset->SetMaxTime(maxtime);
//    Dataset->SetTime(t0);
//    
//    // CFL dimensionless
//    cfl = ((dyD*(mD))*(dt/REAL(n_sub_dt)))/(porosityref*dxD*dyD);
//    
//    TPZDarcyAnalysis *SandStone5 = new TPZDarcyAnalysis(Dataset,Layers,Rocks);
//    SandStone5->SetFluidData(PVTData);
//    SandStone5->RunAnalysis();
//    rates(4,0) = dxD;
//    rates(4,1) = SandStone5->fL2_norm[0];
//    rates(4,2) = SandStone5->fl2_norm_flux[0];
//    rates(4,3) = SandStone5->fL2_norm_s[0];
//    rates(4,4) = fabs(cfl);
//    delete SandStone5;
//    
//    rates *= 1000000.0;
//    rates.Print("data = ",std::cout,EMathematicaInput);
    
    std::cout << std::endl;
    
}