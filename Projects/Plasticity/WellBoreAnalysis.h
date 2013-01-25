//
//  WellBoreAnalysis.h
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef PZ_WellBoreAnalysis_h
#define PZ_WellBoreAnalysis_h

#include "pzcmesh.h"
#include "TPZSandlerDimaggio.h"
#include "TPZTensor.h"


class TPZWellBoreAnalysis
{
    
public:
    
    struct TConfig
    {
        TConfig();
        
        ~TConfig();
        
        TConfig(const TConfig &copy);
        
        TConfig &operator=(const TConfig &copy);
        
        /// Apply the deformation of the configuration to the element
        void ApplyDeformation(TPZCompEl *cel);
        
        /// Compute the maximum plastic element deformation associated with each element
        void ComputeElementDeformation();
        
        /// Delete the elements with sqj2 above the given value;
        void DeleteElementsAbove(REAL sqj2);
        
        /// Change the polynomial order of element using the plastic deformation as threshold
        void PRefineElementsAbove(REAL sqj2, int porder, set<int> &elindices);
        
        /// Divide the element using the plastic deformation as threshold
        void DivideElementsAbove(REAL sqj2, set<int> &elindices);
        
        /// Initialize the plastic history of the integration points of the element
        void InitializeElement(TConfig &from, TPZCompEl *cel);
        
        // Geometry of mesh
        /// radius of the well
        REAL fInnerRadius;
        /// radius of the computational domain
        REAL fOuterRadius;
        
        /// confinement stress
        TPZTensor<STATE> fConfinement;
        
        /// Parameters of the Sandler DiMaggio plasticity
        TPZSandlerDimaggio fSD;
        /// Fluid pressure
        REAL fFluidPressure;
        
        /// Geometric mesh
        TPZGeoMesh fGMesh;
        
        /// Computational mesh
        TPZCompMesh fCMesh;
        
        /// Matrix of incremental solutions
        TPZFMatrix<STATE> fAllSol;
        
        /// Vector containing maximum element plastic deformation
        TPZVec<REAL> fPlasticDeformSqJ2;
        
    };

    TPZWellBoreAnalysis();
    
    ~TPZWellBoreAnalysis();
    
    void ExecuteInitialSimulation();
    
    void ExecuteSimulation();
    
    void TransferSolutionTo(TConfig &config);
    
    void DeleteElementsAbove(REAL sqj2)
    {
        fCurrentConfig.DeleteElementsAbove(sqj2);
    }
    
    void PRefineElementAbove(REAL sqj2, int porder)
    {
        std::set<int> elindices;
        fCurrentConfig.PRefineElementsAbove(sqj2, porder,elindices);
        // subject to integration points to the deformation history
        ApplyHistory(elindices);
    }
    
    /// Divide the element using the plastic deformation as threshold
    void DivideElementsAbove(REAL sqj2);
    
    /// Post process the results of the current configuration
    void PostProcess(int resolution);

static void StandardConfiguration(TPZWellBoreAnalysis &obj);

    
private:
    
    /// Reset the plastic memory of the integration points of these elements
    void ApplyHistory(std::set<int> &elindices);
    
    TConfig fCurrentConfig;
    
    std::list<TConfig> fSequence;
    
    int fPostProcessNumber;
    
};

#endif
