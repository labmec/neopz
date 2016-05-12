//
//  TRMRawData.h
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//
// This class store the computational information required for iRMS

#ifndef __PZ__TRMRawData__
#define __PZ__TRMRawData__

#include <stdio.h>
#include "pzreal.h"
#include "pzstack.h"
#include "pzfunction.h"

class TRMRawData {
    
public:
    
    /** @brief default constructor */
    TRMRawData();

    /** @brief default destructor */
    ~TRMRawData();
    
    
    /** @brief unknown part will be deprecated in time */
    REAL fLw;
    bool fHasLiner;
    bool fHasCasing;
    
    REAL fReservoirWidth;
    REAL fReservoirLength;
    REAL fReservoirHeight;
    REAL fProdVertPosition;
    REAL fWellDiam;
    
    /** @brief vector that stores all material ids associated with omega domain */
    TPZStack< int > fOmegaIds;
    
    /** @brief vector that stores all material ids associated with gamma domain */
    TPZStack< int > fGammaIds;
    
    /** @brief vector that stores all material ids associated with skeleton domain */
    TPZStack< int > fSkeletonIds;

    /** @brief vector that stores pointers to L2 function associated with with gamma domain at intial conditions */
    TPZStack< TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > >  fIntial_bc_data;
    
    /** @brief vector that stores pointers to L2 function associated with with gamma domain at given conditions */
    TPZStack< TPZVec< std::pair< int, TPZAutoPointer<TPZFunction<REAL> > > > >  fRecurrent_bc_data;
    
    /**
     * @ingroup Configuration Cases
     * @brief Define the colletion of materials ids and functions being used as boundary conditions
     * @since May 08, 2016
     */
    
    /** @brief Define the materials for a primitive mono-phasic example and their functions associated */
    void WaterReservoirBox();
    
    static void Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP);
    
    static void Flux(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    static void Impervious(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
};

#endif /* defined(__PZ__TRMRawData__) */
