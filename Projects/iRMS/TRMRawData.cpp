//
//  TRMRawData.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
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
void TRMRawData::WaterReservoirBox(){
    
    std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc;
    
    // Single flow
    TPZAutoPointer<TRMPhaseProperties> water    = new TRMWaterPhase;
    TPZAutoPointer<TRMPhaseProperties> oil      = new TRMOilPhase;
    TPZAutoPointer<TRMPhaseProperties> gas      = new TRMGasPhase;
    fSystemType.Push("water");
    fPhases.Push(water);
    int n_data = fSystemType.size();
    
    // Setting up gravity
    fg.Resize(3, 0.0);
    fg[2] = 9.81;
    
    int map_model = 0; // constant
    fMap = new TRMSpatialPropertiesMap;
    fMap->SetMapModel(map_model);

    // Time control parameters
    REAL hour       = 3600.0;
    REAL day        = hour * 24.0;
    
    fn_steps  = 10;
    fdt = 1000.0*day;
    fdt_up = 1.0;
    fdt_down = 1.0;

    // Numeric controls
    fn_corrections = 2;
    fepsilon_res = 0.001;
    fepsilon_cor = 0.001;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);

//    // Singlewell mesh
//    int bc_W = 2;
//    int bc_A = 3;
//    int bc_I = 4;
//    
//    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > W(n_data);
//    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > A(n_data);
//    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > I(n_data);
//    
//    fGammaIds.Push(bc_W);
//    W[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Flux));
//    fRecurrent_bc_data.Push(W);
//    
//    fGammaIds.Push(bc_A);
//    A[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
//    fRecurrent_bc_data.Push(A);
//    
//    fGammaIds.Push(bc_I);
//    I[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Impervious));
//    fRecurrent_bc_data.Push(I);
    
    int bc_W = 11;
    int bc_E = 12;
    int bc_S = 13;
    int bc_N = 14;
    int bc_B = 15;
    int bc_T = 16;
    
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > W(n_data);
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > E(n_data);    
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > S(n_data);
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > N(n_data);
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > B(n_data);
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > T(n_data);
    
    fGammaIds.Push(bc_W);
    W[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Flux));
    fRecurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = std::make_pair(0,new TPZDummyFunction<REAL>(Pressure));
    fRecurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_N);
    N[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(N);
    
    fGammaIds.Push(bc_B);
    B[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(B);
    
    fGammaIds.Push(bc_T);
    T[0] = std::make_pair(1,new TPZDummyFunction<REAL>(Impervious));
    fRecurrent_bc_data.Push(T);
    
}

void TRMRawData::Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP)
{
    REAL p = 1.0e+7;// 1.0342e+7; // 1500 psi
    P[0] = p;
    return;
}

void TRMRawData::Flux(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    REAL f = -1.0;
    F[0] = f;
    return;
}

void TRMRawData::Impervious(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF)
{
    REAL f = 0.0;
    F[0] = f;
    return;
}

// @}