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
    
    /** @brief Use, level and resolution of MHM process */
    fMHMResolutionQ.first = false;
    fMHMResolutionQ.second.first = 0;
    fMHMResolutionQ.second.second = 0;
    
    /** @brief Use of increased transpor resolution transfers operators */
    fIncreaseTransporResolutionQ.first = false;
    fIncreaseTransporResolutionQ.second = 0;
    
    /** @brief Use of RB method that surrogates */
    fReduceBasisQ.first = false;
    fReduceBasisQ.second.first = false;
    fReduceBasisQ.second.second.Resize(0);
    
    /** @brief Gmsh grid file */
    fGridName = "";
    
    /** @brief Set SPE10 fields file */
    fPermPorFields.first  = "";
    fPermPorFields.second = "";
    
    /** @brief number of blocks i, j and k  */
    fNBlocks.resize(0);
    
    /** @brief size of blocks dx, dy and dz  */
    fBlocks_sizes.resize(0);
    
    /** @brief phases = {alpha, beta, gamma} */
    fPhases.resize(0);
    
    /** @brief Porperties map */
    fMap = NULL;
    
}

TRMRawData::~TRMRawData()
{
    
}

/** @brief Define the materials for a primitive mono-phasic example */
void TRMRawData::SinglePhaseReservoirHMM(bool Is3DGeometryQ){
    
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
    
    int map_model = 1; // constant -> 0, function -> 1, SPE10 interpolation -> 2
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    fGridName = "ch_fem_thiem/reservoir_3D_H.msh";
//    fGridName = "ch_fem_thiem/reservoir_3D_P.msh";
//    fGridName = "ch_fem_thiem/reservoir_3D_T.msh";
    fPermPorFields.first = "ch_fem_thiem/spe_perm.dat";
    fPermPorFields.second = "ch_fem_thiem/spe_phi.dat";
    fNBlocks.Push(60);
    fNBlocks.Push(220);
    fNBlocks.Push(2);
    fBlocks_sizes.Push(1.6666666667);
    fBlocks_sizes.Push(4.5454545455);
    fBlocks_sizes.Push(50.0);
    fMap->SetSpatialFields(fNBlocks, fBlocks_sizes, fPermPorFields);
    fMap->LoadSPE10Map(false);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    //    fReportingTimes.Push(std::make_pair(1000.0*day,true));
    //    fReportingTimes.Push(std::make_pair(900.0*day,true));
    //    fReportingTimes.Push(std::make_pair(800.0*day,true));
    //    fReportingTimes.Push(std::make_pair(700.0*day,true));
    //    fReportingTimes.Push(std::make_pair(600.0*day,true));
    //    fReportingTimes.Push(std::make_pair(400.0*day,true));
    //    fReportingTimes.Push(std::make_pair(300.0*day,true));
    //    fReportingTimes.Push(std::make_pair(200.0*day,true));
    fReportingTimes.Push(std::make_pair(100.0*day,true));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt = 100.0*day;
    fdt_max = 100.0*day;
    fdt_min = 0.1*day;
    fdt_up = 1.5;
    fdt_down = 0.1;
    
    // Numeric controls
    fn_corrections = 10;
    fepsilon_res = 0.1;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true;
    fIsAdataptedQ = true;
    fEnhancedPressureQ = false;
    fMHMResolutionQ.first = true;
    fMHMResolutionQ.second.first = 0; // level
    fMHMResolutionQ.second.second = 2; // fine
    
    
    // Rock materials ids
    int Rock = 5;
    int wellbore_p = 6;
    int wellbore_i = 7;
    fOmegaIds.Push(Rock);
    fOmegaIds.Push(wellbore_p);
    fOmegaIds.Push(wellbore_i);
    
    int bc_W = 13;
    int bc_N = 12;
    int bc_E = 11;
    int bc_S = 10;
    int bc_T = 9;
    int bc_B = 8;
    
    if (!Is3DGeometryQ) {
        bc_W = 9;
        bc_E = 11;
        bc_S = 8;
        bc_N = 10;
        bc_B = 100; // inserted but not being used
        bc_T = 100; // inserted but not being used
    }
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > W(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > E(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > S(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > N(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_W);
    W[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fIntial_bc_data.Push(E);
    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fIntial_bc_data.Push(S);
    S[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fIntial_bc_data.Push(N);
    N[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
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
    
    int bc_p_lids = 1;
    int bc_i_lids = 2;
    int bc_Prod = 3;
    int bc_Inj  = 4;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);
    
    
    fGammaIds.Push(bc_p_lids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_i_lids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureThiem));
    fRecurrent_bc_data.Push(WInj);
    
    
}

void TRMRawData::PressureThiem(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP){
    
    REAL Q = -0.01;
    REAL mu = 0.001;
    REAL h = 10.0;
    REAL Kappa =  1.0e-13;
    REAL p0 = 25.0e6;
    REAL r0 = 50.0;
    REAL r;
    r = sqrt(pt[0]*pt[0] + pt[1]*pt[1]);
    REAL p = p0 - (mu*Q)/(2.0*M_PI*h*Kappa)*log(r/r0);
    
    P[0] = p;
    return;
    
}

void TRMRawData::FluxThiem(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF){
    
    REAL r;
    r = sqrt(pt[0]*pt[0] + pt[1]*pt[1]);
    REAL flux_b = (0.159155/r);//*((coordsX/r)*iHat+(coordsY/r)*jHat)
    F[0] = flux_b;
    return;
}


/** @brief Define the materials for a primitive mono-phasic example */
void TRMRawData::SinglePhaseReservoir(bool Is3DGeometryQ){
    
//    std::pair< int, TPZFunction<REAL> * > bc;
    
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
    fg[1] = -9.81;
    
    int map_model = 0; // constant -> 0, function -> 1, SPE10 interpolation -> 2
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    fGridName = "Meshes/Gmsh/reservoir.msh";
    fPermPorFields.first = "case_2/spe_perm.dat";
    fPermPorFields.second = "case_2/spe_phi.dat";
    fNBlocks.Push(60);
    fNBlocks.Push(220);
    fNBlocks.Push(2);
    fBlocks_sizes.Push(1.6666666667);
    fBlocks_sizes.Push(4.5454545455);
    fBlocks_sizes.Push(50.0);
    fMap->SetSpatialFields(fNBlocks, fBlocks_sizes, fPermPorFields);
    fMap->LoadSPE10Map(false);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
//    fReportingTimes.Push(std::make_pair(1000.0*day,true));
//    fReportingTimes.Push(std::make_pair(900.0*day,true));
//    fReportingTimes.Push(std::make_pair(800.0*day,true));
//    fReportingTimes.Push(std::make_pair(700.0*day,true));
//    fReportingTimes.Push(std::make_pair(600.0*day,true));
//    fReportingTimes.Push(std::make_pair(400.0*day,true));
    fReportingTimes.Push(std::make_pair(500.0*day,true));
    fReportingTimes.Push(std::make_pair(100.0*day,true));
    fReportingTimes.Push(std::make_pair(50.0*day,true));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt       = 50.0*day;
    fdt_max   = 100.0*day;
    fdt_min   = 0.1*day;
    fdt_up    = 1.5;
    fdt_down  = 0.1;
    
    // Numeric controls
    fn_corrections = 20;
    fepsilon_res = 0.5;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true; // Deprecated fixed due to secant method
    fIsAdataptedQ = false;
    fEnhancedPressureQ = false;
    fMHMResolutionQ.first = false;
    fMHMResolutionQ.second.first = 0; // level
    fMHMResolutionQ.second.second = 0; // fine
    
    // RB controls
    fReduceBasisQ.first = false;
    fReduceBasisQ.second.first = false;
    fReduceBasisQ.second.second.Push(1); // x
    fReduceBasisQ.second.second.Push(1); // y
    fReduceBasisQ.second.second.Push(1); // z
    
    // Rock materials ids
    int Rock = 5;
    int wellbore_p = 6;
    int wellbore_i = 7;
    fOmegaIds.Push(Rock);
    fOmegaIds.Push(wellbore_p);
    fOmegaIds.Push(wellbore_i);
    
    int bc_W = 13;
    int bc_N = 12;
    int bc_E = 11;
    int bc_S = 10;
    int bc_T = 9;
    int bc_B = 8;
    
    if (!Is3DGeometryQ) {
        bc_W = 9;
        bc_E = 11;
        bc_S = 8;
        bc_N = 10;
        bc_B = 100; // inserted but not being used
        bc_T = 100; // inserted but not being used
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

    int bc_p_lids = 1;
    int bc_i_lids = 2;
    int bc_Prod = 3;
    int bc_Inj  = 4;

    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);

    
    fGammaIds.Push(bc_p_lids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_i_lids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(2,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure_inj));
    fRecurrent_bc_data.Push(WInj);

    
}

void TRMRawData::Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP)
{
    REAL p = 10.0e+6;// 1.0342e+7; // 1500 psi
    P[0] = p;
    return;
}

void TRMRawData::Pressure_inj(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP)
{
    REAL p = 2.0e+7;// 1.0342e+7; // 1500 psi
    P[0] = p;
    return;
}

void TRMRawData::Flux(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    REAL flux_b = -0.195775;//-0.58569;
    
    REAL day = 86400;
    REAL flux = 0.0;
    REAL c_time = time/day;
    if (c_time > 300.0) {
        flux = flux_b;
    }
    
    F[0] = flux_b;
    return;
}

void TRMRawData::Aquifer(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    
    REAL MPa = 1.0e6;
    
    // Aquifer properties
    REAL pressure_aquifer = 25.0*MPa;
    REAL mu_w = 0.001;
    REAL k = 1.0e-13;
    REAL h = 100.0;
    REAL phi = 0.25;
    REAL ct = 
    
    
    F[0] = pressure_aquifer;
    
    return;
}

void TRMRawData::Impervious(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    REAL f = 0.0;
    F[0] = f;
    return;
}



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




/** @brief Define the materials for a primitive two-phase flow example and their functions associated */
void TRMRawData::TwoPhaseWaterOilReservoir(bool Is3DGeometryQ){
    
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
    fg[1] = -9.81;
//    fg[2] = -9.81;
    
    int map_model = 0; // constant -> 0, function -> 1, SPE10 interpolation -> 2
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    fGridName = "Meshes/Gmsh/reservoir_cad.msh";
    fPermPorFields.first = "case_2/spe_perm.dat";
    fPermPorFields.second = "case_2/spe_phi.dat";
    fNBlocks.Push(60);
    fNBlocks.Push(220);
    fNBlocks.Push(1); // for 2D cases
    fBlocks_sizes.Push(1.6666666667);
    fBlocks_sizes.Push(4.5454545455);
    fBlocks_sizes.Push(50.0);
    fMap->SetSpatialFields(fNBlocks, fBlocks_sizes, fPermPorFields);
    fMap->LoadSPE10Map(false);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    //    fReportingTimes.Push(std::make_pair(1000.0*day,true));
    //    fReportingTimes.Push(std::make_pair(900.0*day,true));
    //    fReportingTimes.Push(std::make_pair(800.0*day,true));
    //    fReportingTimes.Push(std::make_pair(700.0*day,true));
    //    fReportingTimes.Push(std::make_pair(600.0*day,true));
    //    fReportingTimes.Push(std::make_pair(400.0*day,true));
    fReportingTimes.Push(std::make_pair(1000.0*day,false));
    fReportingTimes.Push(std::make_pair(500.0*day,false));
    fReportingTimes.Push(std::make_pair(100.0*day,true));
    fReportingTimes.Push(std::make_pair(0.0*day,true));
    
    fn_steps  = 100;
    fdt       = 100.0*day;
    fdt_max   = 100.0*day;
    fdt_min   = 0.1*day;
    fdt_up    = 1.0;
    fdt_down  = 1.0;
    
    // Numeric controls
    fn_corrections = 20;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true; // Deprecated fixed due to secant method
    fIsAdataptedQ = false;
    fEnhancedPressureQ = false;
    fMHMResolutionQ.first = false;
    fMHMResolutionQ.second.first = 0; // level
    fMHMResolutionQ.second.second = 0; // fine
    fIncreaseTransporResolutionQ.first = true;
    fIncreaseTransporResolutionQ.second = 0;
    
    // RB controls
    fReduceBasisQ.first = true;
    fReduceBasisQ.second.first = false;
    fReduceBasisQ.second.second.Push(2); // x
    fReduceBasisQ.second.second.Push(2); // y
    fReduceBasisQ.second.second.Push(2); // z
    
    // Rock materials ids
    int Rock = 5;
    int wellbore_p = 6;
    int wellbore_i = 7;
    fOmegaIds.Push(Rock);
    fOmegaIds.Push(wellbore_p);
    fOmegaIds.Push(wellbore_i);
    
    int bc_W = 13;
    int bc_N = 12;
    int bc_E = 11;
    int bc_S = 10;
    int bc_T = 9;
    int bc_B = 8;
    
    if (!Is3DGeometryQ) {
        bc_W = 9;
        bc_E = 11;
        bc_S = 8;
        bc_N = 10;
        bc_B = 100; // inserted but not being used
        bc_T = 100; // inserted but not being used
    }
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > W(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > E(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > S(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > N(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > B(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > T(n_data);
    
    fGammaIds.Push(bc_W);
    W[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(E);
    E[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
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
    
    int bc_p_lids = 1;
    int bc_i_lids = 2;
    int bc_Prod = 3;
    int bc_Inj  = 4;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);
    
    
    fGammaIds.Push(bc_p_lids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_2p));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_i_lids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_2p));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(2,new TPZDummyFunction<REAL>(PressureInlet_2p));
    fRecurrent_bc_data.Push(WInj);
    
}

void TRMRawData::PressureOutlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL p = 1.0e+7;// 1.0342e+7; // 1500 psi
    REAL S = 0.0;
    f[0] = p;
    f[1] = S;
    return;
    
}

void TRMRawData::PressureInlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL p = 2.0e+7;// 1.0342e+7; // 3000 psi
    REAL S = 1.0;
    f[0] = p;
    f[1] = S;
    return;
    
}

void TRMRawData::FluxInlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux_b = -0.19523, S = 1.0;
    
    REAL day = 86400;
    REAL flux = flux_b + 0.00*(sin((time/day)/100));
    
    f[0] = flux;
    f[1] = S;
    return;
    
}

void TRMRawData::Aquifer_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = -0.001, S = 1.0;
    
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


/** @brief Define the materials for a primitive two-phase flow example and their functions associated */
void TRMRawData::CaseTracerTransport(bool Is3DGeometryQ){
    
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
    
    int map_model = 0; // constant -> 0, function -> 1, SPE10 interpolation -> 2
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
//    fGridName = "case_2/reservoir_2D_T.msh";
//    fGridName = "case_2/reservoir_2D_Q.msh";
//    fGridName = "case_2/reservoir_3D_T.msh";
    fGridName = "Meshes/Gmsh/reservoir.msh";
//    fGridName = "case_2/reservoir_3D_H.msh";
    fPermPorFields.first = "case_2/spe_perm.dat";
    fPermPorFields.second = "case_2/spe_phi.dat";
    fNBlocks.Push(60);
    fNBlocks.Push(220);
    fNBlocks.Push(2);
    fBlocks_sizes.Push(1.6666666667);
    fBlocks_sizes.Push(4.5454545455);
    fBlocks_sizes.Push(50.0);
    fMap->SetSpatialFields(fNBlocks, fBlocks_sizes, fPermPorFields);
    fMap->LoadSPE10Map(false);
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    fReportingTimes.Push(std::make_pair(1000.0*day,true));
    fReportingTimes.Push(std::make_pair(900.0*day,false));
    fReportingTimes.Push(std::make_pair(800.0*day,false));
    fReportingTimes.Push(std::make_pair(700.0*day,false));
    fReportingTimes.Push(std::make_pair(600.0*day,false));
    fReportingTimes.Push(std::make_pair(500.0*day,false));
    fReportingTimes.Push(std::make_pair(400.0*day,false));
    fReportingTimes.Push(std::make_pair(300.0*day,false));
    fReportingTimes.Push(std::make_pair(200.0*day,false));
    fReportingTimes.Push(std::make_pair(100.0*day,false));
    fReportingTimes.Push(std::make_pair(0.000*day,true));
    
    fn_steps  = 500;
    fdt = 50.0*day;
    fdt_max = 50.0*day;
    fdt_min = 0.01*day;
    fdt_up = 1.0;
    fdt_down = 1.0;
    
    // Numeric controls
    fn_corrections = 50;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.001;
    fIsQuasiNewtonQ = true;
    fMHMResolutionQ.first = false;
    fMHMResolutionQ.second.first = 0;
    fMHMResolutionQ.second.second = 0;
    fIncreaseTransporResolutionQ.first = true;
    fIncreaseTransporResolutionQ.second = 2;
    
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
    W[0] = std::make_pair(0,new TPZDummyFunction<REAL>(CTracer_PressureOutlet_2p));
//    W[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(2,new TPZDummyFunction<REAL>(CTracer_PressureInlet_2p));
//    W[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
//    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(CTracer_PressureOutlet_2p));
    E[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(E);
//    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(CTracer_PressureOutlet_2p));
    E[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(S);
    S[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(N);
    N[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(N);
    
    fGammaIds.Push(bc_B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(B);
    B[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(B);
    
    fGammaIds.Push(bc_T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(T);
    T[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(T);
    
    int bc_lids = 1;
    int bc_Prod = 2;
    int bc_Inj  = 3;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);
    
    fGammaIds.Push(bc_lids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_Impervious_2p));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(CTracer_PressureOutlet_2p));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(CTracer_PressureOutlet_2p));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(4,new TPZDummyFunction<REAL>(CTracer_PressureInlet_2p));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(3,new TPZDummyFunction<REAL>(CTracer_PressureInlet_2p));
    fRecurrent_bc_data.Push(WInj);
    
}

void TRMRawData::CTracer_PressureOutlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL p = 1.0e+6;// 150.0 psi
    f[0] = p;
    return;
    
}

void TRMRawData::CTracer_PressureInlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL p = 1.0e+7;// 1500 psi
    REAL S = 1.0;
    f[0] = p;
    f[1] = S;
    return;
    
}

void TRMRawData::CTracer_Impervious_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = 0.0, S = 0.0;
    f[0] = flux;
    f[1] = S;
    return;
}


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



/** @brief Define the materials for a primitive three-phase flow example and their functions associated */
void TRMRawData::WaterOilGasReservoirBox(bool Is3DGeometryQ){
    
    std::pair< int, TPZFunction<REAL> * > bc;
    
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
    
    int map_model = 0; // constant -> 0, function -> 1, SPE10 interpolation -> 2
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);
    
    
    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    fReportingTimes.Push(std::make_pair(5000.0*day,false));
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
    fdt = 10.0*day;
    fdt_max = 500.0*day;
    fdt_min = 0.1*day;
    fdt_up = 10.0;
    fdt_down = 0.1;
    
    // Numeric controls
    fn_corrections = 20;
    fepsilon_res = 0.01;
    fepsilon_cor = 0.0001;
    fIsQuasiNewtonQ = true;
    fIncreaseTransporResolutionQ.first = true;
    fIncreaseTransporResolutionQ.second = 0;
    
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
    W[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(W);
    W[0] = std::make_pair(3,new TPZDummyFunction<REAL>(Aquifer_3p));
    //    W[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
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
    
    int bc_lids = 1;
    int bc_Prod = 2;
    int bc_Inj  = 3;
    
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WLids(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WPro(n_data);
    TPZVec< std::pair< int, TPZFunction<REAL> * > > WInj(n_data);
    
    
    fGammaIds.Push(bc_lids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fIntial_bc_data.Push(WLids);
    WLids[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious_3p));
    fRecurrent_bc_data.Push(WLids);
    
    fGammaIds.Push(bc_Prod);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fIntial_bc_data.Push(WPro);
    WPro[0] = std::make_pair(0,new TPZDummyFunction<REAL>(PressureOutlet_3p));
    fRecurrent_bc_data.Push(WPro);
    
    fGammaIds.Push(bc_Inj);
    WInj[0] = std::make_pair(4,new TPZDummyFunction<REAL>(Impervious));
    fIntial_bc_data.Push(WInj);
    WInj[0] = std::make_pair(3,new TPZDummyFunction<REAL>(FluxInlet_3p));
    fRecurrent_bc_data.Push(WInj);
    
    
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

void TRMRawData::Aquifer_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf){
    
    REAL flux = -0.001, S_w = 1.0, S_o = 0.;
    
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
