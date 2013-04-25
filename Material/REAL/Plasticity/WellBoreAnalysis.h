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
#include "pzgeoel.h"

void CmeshWell(TPZCompMesh *CMesh, TPZMaterial * mat, TPZTensor<STATE> &Confinement, STATE pressure);


class TPZWellBoreAnalysis
{
    
public:
    
    struct TConfig
    {
        TConfig();
        
        ~TConfig();
        
        TConfig(const TConfig &copy);
        
        TConfig &operator=(const TConfig &copy);
        
        /// Write the data to the output stream
        void Write(TPZStream &out);
        
        /// Read the data from the input stream
        void Read(TPZStream &input);
        
        /// Apply the deformation of the configuration to the element
        void ApplyDeformation(TPZCompEl *cel);
        
        /// Compute the maximum plastic element deformation associated with each element
        void ComputeElementDeformation();
        
        /// Delete the elements with sqj2 above the given value;
        void DeleteElementsAbove(REAL sqj2);
        
        /// Diminish the spring which stabilizes the well by the given factor
        void RelaxWellSpring(REAL factor);
        
        /// Change the polynomial order of element using the plastic deformation as threshold
        void PRefineElementsAbove(REAL sqj2, int porder, set<int> &elindices);
        
        /// Divide the element using the plastic deformation as threshold
        void DivideElementsAbove(REAL sqj2, set<int> &elindices);
        
        /// Initialize the plastic history of the integration points of the element
        void InitializeElement(TConfig &from, TPZCompEl *cel);
        
        /// Verify the global equilibrium of the forces by boundary condition
        void VerifyGlobalEquilibrium(std::ostream &out = std::cout);
        
        /// Verify the global equilibrium of the forces by boundary condition (second version)
        void VerifyGlobalEquilibrium2(std::ostream &out = std::cout);
        
        /// Compute the Rhs for the mesh minus the elements with matid
        // this method is cumulative (sums to the rhs)
        void ComputeRhsExceptMatid(int matid, TPZFMatrix<STATE> &rhs);
        
        /// Compute the contribution of only matid
        // this method is cumulative (sums to the rhs)
        void ComputeRhsForMatid(int matid, TPZFMatrix<STATE> &rhs);
        
        /// Compute the resultant x and y force
        void ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force);
	
        /// Return gmesh
        TPZGeoMesh * GetGeoMesh();

        // Geometry of mesh
        /// radius of the well
        REAL fInnerRadius;
        /// radius of the computational domain
        REAL fOuterRadius;
        
        /// confinement stress
        TPZTensor<STATE> fConfinement;
        
        /// Parameters of the Sandler DiMaggio plasticity
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> fSD;
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
    
    /// write the object on the stream
    void Write(TPZStream &output);
    
    /// read the object from the stream
    void Read(TPZStream &input);
    
    void ExecuteInitialSimulation();
    
    void ExecuteSimulation();
    
    static void CheckDeformation(std::string filename = "deform.nb");
    
    void TransferSolutionTo(TConfig &config);
    
    void DeleteElementsAbove(REAL sqj2)
    {
        fCurrentConfig.DeleteElementsAbove(sqj2);
    }
    
    /// Diminish the spring which stabilizes the well by the given factor
    void RelaxWellSpring(REAL factor)
    {
        fCurrentConfig.RelaxWellSpring(factor);
    }
    
    void PRefineElementAbove(REAL sqj2, int porder)
    {
        std::set<int> elindices;
        fCurrentConfig.PRefineElementsAbove(sqj2, porder,elindices);
        // subject to integration points to the deformation history
        ApplyHistory(elindices);
    }
    
    void ChangeMaterialId(int idFrom, int idTo)
    {
        TPZGeoMesh *gmesh = &fCurrentConfig.fGMesh;
        int nel = gmesh->NElements();
        for (int iel=0; iel<nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel) continue;
            int matid = gel->MaterialId();
            if (matid == idFrom) {
                gel->SetMaterialId(idTo);
            }
        }
    }
    
    /// Divide the element using the plastic deformation as threshold
    void DivideElementsAbove(REAL sqj2);
    
    /// Post process the results of the current configuration
    void PostProcess(int resolution);
    
    /// verify the global equilibrium for the current configuration
    void VerifyGlobalEquilibrium(std::ostream &out = std::cout)
    {
        fCurrentConfig.VerifyGlobalEquilibrium(out);
    }
    
    TConfig * GetCurrentConfig ()
    {
        return &fCurrentConfig;
    }

static void StandardConfiguration(TPZWellBoreAnalysis &obj);

    
private:
    
    /// Reset the plastic memory of the integration points of these elements
    void ApplyHistory(std::set<int> &elindices);
    
    TConfig fCurrentConfig;
    
    std::list<TConfig> fSequence;
    
    int fPostProcessNumber;
    
};

#endif
