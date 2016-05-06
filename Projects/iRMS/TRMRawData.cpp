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
    intial_bc_data.Resize(0);
    
    /** @brief vector that stores pointers to L2 function associated with with gamma domain at given conditions */
    recurrent_bc_data.Resize(0);    
    
}

TRMRawData::~TRMRawData()
{
    
}

/** @brief Define the materials for a primitive mono-phasic example */
void TRMRawData::WaterReservoirBox(){
    
    // Single flow
    int n_data = 1;
    
    
    // Rock materials ids
    int Rock = 1;
    fOmegaIds.Push(Rock);
    
    int bc_N = -1;
    int bc_W = -2;
    int bc_E = -3;
    int bc_S = -4;
    int bc_T = -5;
    int bc_B = -6;
    
    TPZVec< TPZAutoPointer<TPZFunction<REAL> > > N(n_data);
    TPZVec< TPZAutoPointer<TPZFunction<REAL> > > W(n_data);
    TPZVec< TPZAutoPointer<TPZFunction<REAL> > > E(n_data);
    TPZVec< TPZAutoPointer<TPZFunction<REAL> > > S(n_data);
    TPZVec< TPZAutoPointer<TPZFunction<REAL> > > T(n_data);
    TPZVec< TPZAutoPointer<TPZFunction<REAL> > > B(n_data);
    
    fGammaIds.Push(bc_N);
    N[0] = new TPZDummyFunction<REAL>(Pressure);
    recurrent_bc_data.Push(N);
    
    fGammaIds.Push(bc_W);
    W[0] = new TPZDummyFunction<REAL>(Pressure);
    recurrent_bc_data.Push(W);
    
    fGammaIds.Push(bc_E);
    E[0] = new TPZDummyFunction<REAL>(Pressure);
    recurrent_bc_data.Push(E);
    
    fGammaIds.Push(bc_S);
    S[0] = new TPZDummyFunction<REAL>(Pressure);
    recurrent_bc_data.Push(S);
    
    fGammaIds.Push(bc_T);
    T[0] = new TPZDummyFunction<REAL>(Pressure);
    recurrent_bc_data.Push(T);
    
    fGammaIds.Push(bc_B);
    B[0] = new TPZDummyFunction<REAL>(Pressure);
    recurrent_bc_data.Push(B);
    
    
}

void TRMRawData::Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP)
{
    REAL p = 0.0;
    P[0] = p;
    return;
}