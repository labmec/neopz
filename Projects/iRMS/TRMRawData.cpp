//
//  TRMRawData.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//  Implemented by Omar Duran since 7/25/15.
//

#include "TRMRawData.h"


TRMRawData::TRMRawData()
{
    fLw = 0.;
    fHasLiner = true;
    fHasCasing = true;
    
    fReservoirWidth = 0.;
    fReservoirLength = 0.;
    fReservoirHeight = 0.;
    fProdVertPosition = 0.;
    fWellDiam = 0.;
    
    /** @brief Material identifier for interfaces */
    fInterface_mat_Id = 10000;
    
    /** @brief vector that stores all material ids associated with omega domain */
    fOmegaIds.Resize(0);
    
    /** @brief vector that stores all material ids associated with gamma domain */
    fGammaIds.Resize(0);
    
    /** @brief vector that stores all material ids associated with skeleton domain */
    fSkeletonIds.Resize(0);
    
    /** @brief vector that stores pointers to L2 function associated with with gamma domain at intial conditions */
    fIntial_bc_data.Resize(0);
    
    /** @brief vector that stores pointers to L2 function associated with with gamma domain at given conditions */
    fRecurrent_bc_data.Resize(0);
    
    /** @brief Definition of the flow system one - two or three phase */
    fSystemType.Resize(0);
    
    /** @brief Definition of the gravity vector flied */
    fg.Resize(0);
    
    /** @brief ntime steps */
    fn_steps = 0;
    
    /** @brief Store time values to be reported */
    fReportingTimes.Resize(0);
    
    /** @brief Time step */
    fdt = 0.0;
    
    /** @brief Min time step */
    fdt_min = 0.0;
    
    /** @brief Max time step */
    fdt_max = 0.0;
    
    /** @brief Increment dt factor */
    fdt_up = 0.0;
    
    /** @brief Decrement dt factor */
    fdt_down = 0.0;
    
    /** @brief number of corrections steps */
    fn_corrections = 0;
    
    /** @brief residue overal tolerance */
    fepsilon_res = 0.0;
    
    /** @brief correction overal tolerance */
    fepsilon_cor = 0.0;
    
    /** @brief set the use of quasi newton method */
    fIsQuasiNewtonQ = false;
    
    /** @brief phases = {alpha, beta, gamma} */
    fPhases.resize(0);
    
    /** @brief Porperties map */
    fMap = NULL;
    
}

TRMRawData::~TRMRawData()
{
    
}


/**
 * @ingroup Configuration Cases
 * @brief one-phase water flow configuration
 * @since May 08, 2016
 */


/** @brief Define the materials for a primitive mono-phasic example */
void TRMRawData::WaterReservoirBox(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    water->SetRhoModel(0);
    fPhases.Push(water);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[1] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
//    fReportingTimes.Push(std::make_pair(500.0*day,true));
//    fReportingTimes.Push(std::make_pair(450.0*day,false));
//    fReportingTimes.Push(std::make_pair(400.0*day,false));
//    fReportingTimes.Push(std::make_pair(350.0*day,false));
//    fReportingTimes.Push(std::make_pair(300.0*day,false));
//    fReportingTimes.Push(std::make_pair(250.0*day,false));
//    fReportingTimes.Push(std::make_pair(200.0*day,false));
//    fReportingTimes.Push(std::make_pair(150.0*day,false));
//    fReportingTimes.Push(std::make_pair(100.0*day,true));
//    fReportingTimes.Push(std::make_pair(50.0*day,false));
//    fReportingTimes.Push(std::make_pair(40.0*day,false));
//    fReportingTimes.Push(std::make_pair(30.0*day,false));
//    fReportingTimes.Push(std::make_pair(20.0*day,false));
//    fReportingTimes.Push(std::make_pair(10.0*day,false));
    fReportingTimes.Push(std::make_pair(1.0*day,false));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt = 10.0*day;
    fdt_max = 100.0*day;
    fdt_min = 0.1*day;
    fdt_up = 1.5;
    fdt_down = 0.1;
    
    // Numeric controls
    fn_corrections = 20;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.0001;
    fIsQuasiNewtonQ = true;
    
    
    // Rock materials ids
    int Rock = 4;
    fOmegaIds.Push(Rock);
    
    int bc_W = 10;
    int bc_E = 8;
    int bc_S = 7;
    int bc_N = 9;
    int bc_B = 5;
    int bc_T = 6;
    
    if (!Is3DGeometryQ) {
        bc_W = 6;
        bc_E = 8;
        bc_S = 5;
        bc_N = 7;
        bc_B = 100;
        bc_T = 100;
    }
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > W(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > E(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > S(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > N(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_W);
    W[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(E);
    E[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(S);
    S[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(N);
    N[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(N);
    
    fGammaIds.Push(bc_B);
    B[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(B);
    B[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(B);
    
    fGammaIds.Push(bc_T);
    T[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(T);
    T[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(T);

    int bc_lids = 1;
    int bc_Prod = 2;
    int bc_Inj  = 3;

    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);

    
    fGammaIds.Push(bc_lids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Flux));
    fRecurrent_bc_data.Push(WInj);

    
}

void TRMRawData::Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP)
{
    REAL p = 1.0e+7;// 1.0342e+7; // 1500 psi
    P[0] = p;
    return;
}

void TRMRawData::Flux(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    REAL flux_b = -0.0184;
    
    REAL day = 86400;
    REAL flux = flux_b + 0.00*(sin((time/day)/100));
    
    F[0] = flux;
    return;
}

void TRMRawData::FluxAquifer(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    
    REAL flux_aquifer = -0.000000184;
    F[0] = flux_aquifer;
    return;
}

void TRMRawData::Impervious(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    REAL f = 0.0;
    F[0] = f;
    return;
}

// @}


/**
 * @ingroup Configuration Cases
 * @brief one-phase water flow configuration
 * @since May 08, 2016
 */


/** @brief Define the materials for a primitive mono-phasic example */
void TRMRawData::WaterReservoirCircle(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    gas->SetRhoModel(0);
    fPhases.Push(gas);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[1] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    fReportingTimes.Push(std::make_pair(100.0*day,true));
    fReportingTimes.Push(std::make_pair(50.0*day,false));
    fReportingTimes.Push(std::make_pair(40.0*day,false));
    fReportingTimes.Push(std::make_pair(30.0*day,false));
    fReportingTimes.Push(std::make_pair(20.0*day,false));
    fReportingTimes.Push(std::make_pair(10.0*day,false));
    fReportingTimes.Push(std::make_pair(1.0*day,false));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt = 0.1*day;
    fdt_max = 30.0*day;
    fdt_min = 0.1*day;
    fdt_up = 1.5;
    fdt_down = 0.1;
    
    // Numeric controls
    fn_corrections = 50;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);
    
    int bc_Outlet   = 2;
    int bc_Inlet    = 3;
    int bc_Noflux   = 4;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Noflux(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Inlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Outlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_Noflux);
    Noflux[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(Noflux);
    Noflux[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(Noflux);
    
    fGammaIds.Push(bc_Inlet);
    Inlet[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(Inlet);
    Inlet[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Flux));
    fRecurrent_bc_data.Push(Inlet);
    
    fGammaIds.Push(bc_Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fIntial_bc_data.Push(Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fRecurrent_bc_data.Push(Outlet);
    
}


// @}



/** @brief Define the materials for a primitive two-phase flow example and their functions associated */
void TRMRawData::WaterOilReservoirBox(bool Is3DGeometryQ){
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    fSystemType.Push("water");
    water->SetRhoModel(0);
    water->SetRhoModel(0);
    fPhases.Push(water);
    fPhases.Push(water);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[1] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;

    fReportingTimes.Push(std::make_pair(5000.0*day,true));
    fReportingTimes.Push(std::make_pair(4500.0*day,false));
    fReportingTimes.Push(std::make_pair(4000.0*day,false));
    fReportingTimes.Push(std::make_pair(3500.0*day,false));
    fReportingTimes.Push(std::make_pair(3000.0*day,false));
    fReportingTimes.Push(std::make_pair(2500.0*day,false));
    fReportingTimes.Push(std::make_pair(2000.0*day,false));
    fReportingTimes.Push(std::make_pair(1500.0*day,false));
    fReportingTimes.Push(std::make_pair(1000.0*day,true));
    fReportingTimes.Push(std::make_pair(500.0*day,false));
    fReportingTimes.Push(std::make_pair(400.0*day,false));
    fReportingTimes.Push(std::make_pair(300.0*day,false));
    fReportingTimes.Push(std::make_pair(200.0*day,false));
    fReportingTimes.Push(std::make_pair(100.0*day,false));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt = 1.0*day;
    fdt_max = 100.0*day;
    fdt_min = 0.1*day;
    fdt_up = 1.5;
    fdt_down = 0.1;
    
    // Numeric controls
    fn_corrections = 150;
    fepsilon_res = 0.1;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true;
    
    
    // Rock materials ids
    int Rock = 4;
    fOmegaIds.Push(Rock);
    
    int bc_W = 10;
    int bc_E = 8;
    int bc_S = 7;
    int bc_N = 9;
    int bc_B = 5;
    int bc_T = 6;
    
    if (!Is3DGeometryQ) {
        bc_W = 6;
        bc_E = 8;
        bc_S = 100;
        bc_N = 100;
        bc_B = 5;
        bc_T = 7;
    }
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > W(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > E(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > S(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > N(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_W);
    W[0] = std::make_pair(3,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(3,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(E);
    E[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(S);
    S[0] = std::make_pair(3,new TPZDummyFunction<REAL>(Aquifer_2p));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(N);
    N[0] = std::make_pair(3,new TPZDummyFunction<REAL>(Aquifer_2p));
    fRecurrent_bc_data.Push(N);
    
    fGammaIds.Push(bc_B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(B);
    
    fGammaIds.Push(bc_T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(T);
    
    int bc_lids = 1;
    int bc_Prod = 2;
    int bc_Inj  = 3;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);
    
    
    fGammaIds.Push(bc_lids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(3,new TPZDummyFunction<REAL>(FluxInlet_2p));
    fRecurrent_bc_data.Push(WInj);
    
    
}

void TRMRawData::PressureOutlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL p = 1.0e+7;// 1.0342e+7; // 1500 psi
    f[0] = p;
    return;
    
}

void TRMRawData::FluxInlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux_b = -0.0058569, S = 1.0;
    
    REAL day = 86400;
    REAL flux = flux_b + 0.00*(sin((time/day)/100));
    
    f[0] = flux;
    f[1] = S;
    return;
    
}

void TRMRawData::Aquifer_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = -0.000001, S = 1.0;
    
    f[0] = flux;
    f[1] = S;
    return;
    
}

void TRMRawData::Impervious_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = 0.0, S = 0.0;
    f[0] = flux;
    f[1] = S;
    return;
}


// @}


/** @brief Define the materials for a primitive two-phase flow example and their functions associated */
void TRMRawData::WaterOilReservoirVertical(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    fSystemType.Push("water");
    water->SetRhoModel(0);
    water->SetRhoModel(0);
    fPhases.Push(water);
    fPhases.Push(water);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[2] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    //    fn_steps  = 40;
    //    fdt = 1.0*day;
    //    fdt_max = 10.0*day;
    //    fdt_min = 0.5*day;
    //    fdt_up = 1.5;
    //    fdt_down = 0.5;
    
    fn_steps  = 40;
    fdt = 1.0*day;
    fdt_max = 30.0*day;
    fdt_min = 0.5*day;
    fdt_up = 1.5;
    fdt_down = 0.5;
    
    // Numeric controls
    fn_corrections = 2;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = false;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);
    
    int bc_Noflux   = 2;
    int bc_Inlet    = 3;
    int bc_Outlet   = 4;
    int bc_B        = 5;
    int bc_T        = 6;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Noflux(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Inlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Outlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_Noflux);
    Noflux[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(Noflux);
    Noflux[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(Noflux);
    
    fGammaIds.Push(bc_Inlet);
    Inlet[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(Inlet);
    Inlet[0] = std::make_pair(3,new TPZDummyFunction<REAL>(FluxInlet_2p));
    fRecurrent_bc_data.Push(Inlet);
    
    fGammaIds.Push(bc_Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fIntial_bc_data.Push(Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fRecurrent_bc_data.Push(Outlet);
    
    fGammaIds.Push(bc_B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(B);
    
    fGammaIds.Push(bc_T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(T);
    
    
}


// @}


/** @brief Define the materials for a primitive two-phase flow example and their functions associated */
void TRMRawData::WaterOilReservoirCircular(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    fSystemType.Push("oil");
    water->SetRhoModel(1);
    oil->SetRhoModel(1);
    fPhases.Push(water);
    fPhases.Push(oil);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[2] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;

    fReportingTimes.Push(std::make_pair(500.0*day,true));
    fReportingTimes.Push(std::make_pair(450.0*day,false));
    fReportingTimes.Push(std::make_pair(400.0*day,false));
    fReportingTimes.Push(std::make_pair(350.0*day,false));
    fReportingTimes.Push(std::make_pair(300.0*day,false));
    fReportingTimes.Push(std::make_pair(250.0*day,false));
    fReportingTimes.Push(std::make_pair(200.0*day,false));
    fReportingTimes.Push(std::make_pair(150.0*day,false));
    fReportingTimes.Push(std::make_pair(100.0*day,true));
    fReportingTimes.Push(std::make_pair(50.0*day,false));
    fReportingTimes.Push(std::make_pair(40.0*day,false));
    fReportingTimes.Push(std::make_pair(30.0*day,false));
    fReportingTimes.Push(std::make_pair(20.0*day,false));
    fReportingTimes.Push(std::make_pair(10.0*day,false));
    fReportingTimes.Push(std::make_pair(1.0*day,false));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt = 10.0*day;
    fdt_max = 100.0*day;
    fdt_min = 0.1*day;
    fdt_up = 2.0;
    fdt_down = 0.1;
    
    // Numeric controls
    fn_corrections = 50;
    fepsilon_res = 0.1;
    fepsilon_cor = 0.01;
    fIsQuasiNewtonQ = true;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);
    
    int bc_Outlet   = 3;
    int bc_Inlet    = 2;
    int bc_Noflux   = 4;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Noflux(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Inlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Outlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_Noflux);
    Noflux[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(Noflux);
    Noflux[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(Noflux);
    
    fGammaIds.Push(bc_Inlet);
    Inlet[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(Inlet);
    Inlet[0] = std::make_pair(3,new TPZDummyFunction<REAL>(FluxInlet_2p));
    fRecurrent_bc_data.Push(Inlet);
    
    fGammaIds.Push(bc_Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fIntial_bc_data.Push(Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fRecurrent_bc_data.Push(Outlet);
    
}


// @}

/** @brief Define the materials for a primitive three-phase flow example and their functions associated */
void TRMRawData::WaterOilGasReservoirBox(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    fSystemType.Push("oil");
    fSystemType.Push("gas");
    water->SetRhoModel(0);
    oil->SetRhoModel(0);
    gas->SetRhoModel(0);
    fPhases.Push(water);
    fPhases.Push(oil);
    fPhases.Push(gas);
    
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[2] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    fn_steps  = 20;
    fdt = 1.0*day;
    fdt_max = 30.0*day;
    fdt_min = 0.5*day;
    fdt_up = 1.0;
    fdt_down = 1.0;
    
    // Numeric controls
    fn_corrections = 50;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.01;
    fIsQuasiNewtonQ = true;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);
    
    int bc_W = 11;
    int bc_E = 12;
    int bc_S = 13;
    int bc_N = 14;
    int bc_B = 15;
    int bc_T = 16;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > W(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > E(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > S(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > N(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_W);
    W[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(3,new TPZDummyFunction<REAL>(FluxInlet_3p));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fIntial_bc_data.Push(E);
    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fRecurrent_bc_data.Push(N);
    
    if (!Is3DGeometryQ) {
        return;
    }
    
    fGammaIds.Push(bc_B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fRecurrent_bc_data.Push(B);
    
    fGammaIds.Push(bc_T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fRecurrent_bc_data.Push(T);
    
    
}

void TRMRawData::PressureOutlet_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL p = 1.0e+7;// 1.0342e+7; // 1500 psi
    f[0] = p;
    return;
    
}

void TRMRawData::FluxInlet_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = -1.84, S_w = 1.0, S_o = 0.0;
    f[0] = flux;
    f[1] = S_w;
    f[2] = S_o;
    return;
    
}

void TRMRawData::Impervious_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = 0.0, S = 0.0;
    f[0] = flux;
    f[1] = S;
    f[2] = S;
    return;
}


// @}


/** @brief Define the materials for a primitive two-phase flow example and their functions associated */
void TRMRawData::WaterOilGasReservoirCircular(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    fSystemType.Push("water");
    fSystemType.Push("water");
    water->SetRhoModel(0);
    water->SetRhoModel(0);
    water->SetRhoModel(0);
    fPhases.Push(water);
    fPhases.Push(water);
    fPhases.Push(water);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    //fg[2] = -9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    
    fn_steps  = 50;
    fdt = 1.0*day;
    fdt_max = 30.0*day;
    fdt_min = 1.0*day;
    fdt_up = 1.0;
    fdt_down = 1.0;
    
    // Numeric controls
    fn_corrections = 50;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);
    
    int bc_Inlet    = 2;
    int bc_Outlet   = 3;
    int bc_Noflux   = 4;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Noflux(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Inlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > Outlet(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_Noflux);
    Noflux[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(Noflux);
    Noflux[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fRecurrent_bc_data.Push(Noflux);
    
    fGammaIds.Push(bc_Inlet);
    Inlet[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(Inlet);
    Inlet[0] = std::make_pair(3,new TPZDummyFunction<REAL>(FluxInlet_3p));
    fRecurrent_bc_data.Push(Inlet);
    
    fGammaIds.Push(bc_Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fIntial_bc_data.Push(Outlet);
    Outlet[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fRecurrent_bc_data.Push(Outlet);
    
}


// @}
