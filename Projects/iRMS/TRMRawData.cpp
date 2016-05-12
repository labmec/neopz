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
    
}

TRMRawData::~TRMRawData()
{
    
}

/** @brief Define the materials for a primitive mono-phasic example */
void TRMRawData::WaterReservoirBox(){
    
    std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > bc;
    
    // Single flow
    int n_data = 1;
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);

    int bc_W = -1;
    int bc_E = -2;
    int bc_S = -3;
    int bc_N = -4;
    int bc_B = -5;
    int bc_T = -6;
    
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > E(n_data);
    TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > W(n_data);
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
    REAL p = 1.0;//1.0342e+7; // 1500 psi
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