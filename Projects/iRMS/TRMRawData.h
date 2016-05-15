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
     * @ingroup Required data for a simulation
     * @brief Define the colletion of materials ids and functions being used as boundary conditions
     * @since May 08, 2016
     */
    
    /** @brief Definition of the flow system one - two or three phase */
    TPZStack<std::string> fSystemType;
    
    /** @brief ntime steps */
    int fn_steps;
    
    /** @brief Time step */
    STATE fdt;
    
    /** @brief Min time step */
    STATE fdt_min;
    
    /** @brief Max time step */
    STATE fdt_max;
    
    /** @brief Increment dt factor */
    STATE fdt_up;
    
    /** @brief Decrement dt factor */
    STATE fdt_down;
    
    /** @brief number of corrections steps */
    int fn_correction;
    
    /** @brief residue overal tolerance */
    STATE fepsilon_res;
    
    /** @brief correction overal tolerance */
    STATE fepsilon_cor;
    
    // @}
    
    /**
     * @ingroup Configuration Cases :: Water flow
     * @brief Define the colletion of materials ids and functions being used as boundary conditions
     * @since May 08, 2016
     */
    
    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void WaterReservoirBox();
    
    static void Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP);
    
    static void Flux(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    static void Impervious(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    
    // @}
    
    /**
     * @ingroup Configuration Cases :: Oil flow
     * @brief Define the colletion of materials ids and functions being used as boundary conditions
     * @since May 08, 2016
     */
    
    // @}
};

#endif /* defined(__PZ__TRMRawData__) */
