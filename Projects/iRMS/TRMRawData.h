//
//  TRMRawData.h
//  PZ
//
//  Created by Philippe Devloo on 5/25/15.
//  Implemented by Omar Duran since 8/25/15.
//
// This class store the computational information required for iRMS

#ifndef __PZ__TRMRawData__
#define __PZ__TRMRawData__

#include <stdio.h>
#include "pzreal.h"
#include "pzstack.h"
#include "pzfunction.h"

#include "TRMPhaseProperties.h"
#include "TRMSpatialPropertiesMap.h"
#include "TRMWaterPhase.h"
#include "TRMOilPhase.h"
#include "TRMGasPhase.h"

class TRMRawData {
    
public:
    
    /** @brief default constructor */
    TRMRawData();
    
    /** @brief default constructor */
    TRMRawData(const TRMRawData &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor */
    TRMRawData &operator=(const TRMRawData &copy)
    {
        DebugStop();
        return *this;
    }
    
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
    
    /** @brief Porperties map */
    TPZAutoPointer<TRMSpatialPropertiesMap> fMap;
    
    /** @brief vector that stores all material ids associated with omega domain */
    TPZStack< int > fOmegaIds;
    
    /** @brief vector that stores all material ids associated with gamma domain */
    TPZStack< int > fGammaIds;
    
    /** @brief vector that stores all material ids associated with skeleton domain */
    TPZStack< int > fSkeletonIds;
    
    /** @brief vector that stores pointers to L2 function associated with with gamma domain at intial conditions */
    TPZStack< TPZVec< std::pair< int, TPZFunction<REAL> * > > >  fIntial_bc_data;
    
    /** @brief vector that stores pointers to L2 function associated with with gamma domain at given conditions */
    TPZStack< TPZVec< std::pair< int, TPZFunction<REAL> * > > >  fRecurrent_bc_data;
    
    
    /** @brief Material identifier for interfaces */
    int fInterface_mat_Id;
    
    /** @brief Definition of the flow system one - two or three phase */
    TPZStack<std::string> fSystemType;
    
    /** @brief Definition of the gravity vector flied */
    TPZVec<STATE> fg;
    
    /** @brief ntime steps */
    int fn_steps;
    
    /** @brief Store time values to be reported */
    TPZStack< std::pair< STATE , bool> , 500 > fReportingTimes;
    
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
    int fn_corrections;
    
    /** @brief residue overal tolerance */
    STATE fepsilon_res;
    
    /** @brief correction overal tolerance */
    STATE fepsilon_cor;
    
    /** @brief set the use of quasi newton method */
    bool fIsQuasiNewtonQ;
    
    /** @brief set the use p adaptation on wellbores */
    bool fIsAdataptedQ;
    
    /** @brief set the use enhanced pressure accuracy */
    bool fEnhancedPressureQ;
    
    /** @brief Use, level and resolution of MHM process */
    std::pair<bool, std::pair<int, int> > fMHMResolutionQ;
    
    /** @brief Use of increased transpor resolution transfers operators */
    std::pair<bool, int> fIncreaseTransporResolutionQ;
    
    /** @brief Gmsh grid file */
    std::string fGridName;
    
    /** @brief Set SPE10 fields file */
    std::pair< std::string , std::string > fPermPorFields;
    
    /** @brief number of blocks i, j and k  */
    TPZStack<int> fNBlocks;
    
    /** @brief size of blocks dx, dy and dz  */
    TPZStack<REAL> fBlocks_sizes;
    
    /** @brief phases = {alpha, beta, gamma} */
    TPZStack< TPZAutoPointer<TRMPhaseProperties> > fPhases;
    
    
    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void SinglePhaseReservoirHMM(bool Is3DGeometryQ);
    
    static void PressureThiem(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP);
    
    static void FluxThiem(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void SinglePhaseReservoir(bool Is3DGeometryQ);
    
    static void Pressure(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP);
    
    static void Pressure_inj(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& P, TPZFMatrix< REAL >& GradP);
    
    static void Flux(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    static void Aquifer(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    static void Impervious(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& F, TPZFMatrix< REAL >& GradF);
    
    
    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void WaterReservoirCircle(bool Is3DGeometryQ);
    
    /** @brief Define the materials for case 1 and 2: A water and gas transport tracer */
    void CaseTracerTransport(bool Is3DGeometryQ);
    
    static void CTracer_PressureOutlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void CTracer_PressureInlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void CTracer_Impervious_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    
    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void TwoPhaseWaterOilReservoir(bool Is3DGeometryQ);
    
    static void PressureOutlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void FluxInlet_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void Aquifer_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void Impervious_2p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    

    
    /** @brief Define the materials for a primitive two-phase flow example and their functions associated */
    void WaterOilReservoirVertical(bool Is3DGeometryQ);
    
    
    /** @brief Define the materials for a primitive two-phase flow example and their functions associated */
    void WaterOilReservoirCircular(bool Is3DGeometryQ);
    
    
    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void WaterOilGasReservoirBox(bool Is3DGeometryQ);
    
    static void PressureOutlet_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void FluxInlet_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void Aquifer_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    static void Impervious_3p(const TPZVec< REAL >& pt, REAL time, TPZVec< REAL >& f, TPZFMatrix< REAL >& Gradf);
    
    

    /** @brief Define the materials for a primitive one-phase flow example and their functions associated */
    void WaterOilGasReservoirCircular(bool Is3DGeometryQ);
    
    
    
    
    
    
};

#endif /* defined(__PZ__TRMRawData__) */