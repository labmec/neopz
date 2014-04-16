//
//  WellBoreAnalysis.h
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 LabMeC. All rights reserved.
//

#ifndef PZ_WellBoreAnalysis_h
#define PZ_WellBoreAnalysis_h

#define PV
#define Elastic

#include "pzcmesh.h"
#include "TPZSandlerDimaggio.h"
#include "TPZTensor.h"
#include "pzgeoel.h"
#include "pzpostprocanalysis.h"
#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "pzstack.h"
//#include "TPZMohrCoulomb.h"
//#include "TPZMohrCoulombNeto.h"
//#include "TPZMohrCoulombPV.h"

class TPZElasticityMaterial;

/// create de standard mesh
void CmeshWell(TPZCompMesh *CMesh, TPZMaterial * mat, TPZTensor<STATE> &Confinement, STATE pressure);

enum EPlasticModel  {ESandler, EMohrCoulomb, EElastic};

/// Class which simulates the stability of a wellbore
class TPZWellBoreAnalysis
{
    
public:
    
    /// this class represents a particular stage of the wellbore analysis
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
        
        /// Verify tangent elasto plastic relation
        void VerifyPlasticTangent(TPZCompEl *cel);
        
        /// Compute the maximum plastic element deformation associated with each element
        void ComputeElementDeformation();
        
        /// Compute the area of the domain at which sqJ2 is above a given value
        REAL ComputeAreaAboveSqJ2(REAL sqj2);
        
        /// Compute the area of the domain at which sqJ2 is above a given value
        REAL OpeningAngle(REAL sqj2);
        
        /// Compute the area of the domain
        REAL ComputeTotalArea();

        /// Compute the removed area of the domain
        REAL RemovedArea();
        
        /// Delete the elements with sqj2 above the given value;
        void DeleteElementsAbove(REAL sqj2);
        
        /// Diminish the spring which stabilizes the well by the given factor
        void RelaxWellSpring(REAL factor);
        
        /// Change the polynomial order of element using the plastic deformation as threshold
        void PRefineElementsAbove(REAL sqj2, int porder, set<long> &elindices);
        
        /// Divide the element using the plastic deformation as threshold
        void DivideElementsAbove(REAL sqj2, set<long> &elindices);
        
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
        
        /// Zera os componentes do rhs para connects diferentes do zero
        void FilterRhs(TPZFMatrix<STATE> &rhs);
        
        /// Compute the resultant x and y force
        void ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force);
        
        /// Add elliptic breakout
        void AddEllipticBreakout(REAL MaiorAxis, REAL MinorAxis);

        /// Create the geometric and computational mesh based on the configuration parameters
        void CreateMesh();
        
        void ModifyWellElementsToQuadratic();
        
        /// Initialize the Sandler DiMaggio object and create the computational mesh
        void CreateComputationalMesh(int porder);
        
        /// Setup post processing mesh
        void CreatePostProcessingMesh();
        

         STATE ComputeFarFieldWork();
        
        /// project the node on the boundary
        // returns true if the coordinate was changed
        bool ProjectNode(TPZVec<REAL> &co);
        
        /// return the largest y-coordinate belonging to the last ellips
        // this value will be used to identify the meaningful points to match the next ellips
        REAL MaxYfromLastBreakout();
        
        void Print(ostream &out);
        
        /// Return gmesh
        TPZGeoMesh * GetGeoMesh();

        // Geometry of mesh
        /// radius of the well
        REAL fInnerRadius;
        /// radius of the computational domain
        REAL fOuterRadius;
        
        /// number of elements in the radial and circumferential direction
        TPZManVector<int,2> fNx;
        
        /// Size of the first element in the radial direction
        REAL fDelx;
        
        /// Elliptic break out maior axis
        // Greater has to be ascending and greater than fInnerRadius
        TPZManVector<REAL,3> fGreater;
        
        /// Elliptic breakout minor axis
        // Minor has to be decreasing and smaller than fInnerRadius
        TPZManVector<REAL,3> fSmaller;
        
        /// confinement stress
        TPZTensor<STATE> fConfinement;
        
        /// Parameters
        //TPZElasticityMaterial fEl;
        
        EPlasticModel fModel;
        
#ifdef PV
        // Sandler Dimaggio PV
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> fSDPV;

        //Mohr PV
#else
        /// Parameters of the Sandler DiMaggio plasticity
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSD;
        
        
#endif
        
        
    

        /// Fluid pressure
        REAL fFluidPressure;
        
        REAL fPorePressure;
        
        /// Geometric mesh3
        TPZGeoMesh fGMesh;
        
        /// Computational mesh
        TPZCompMesh fCMesh;
        
        /// Matrix of incremental solutions
        TPZFMatrix<STATE> fAllSol;
        
        /// Vector containing maximum element plastic deformation
        TPZVec<REAL> fPlasticDeformSqJ2;
        
        /// The post processing mesh with the transferred solution
        TPZPostProcAnalysis fPostprocess;

        std::string fHistoryLog;

        
    };

    TPZWellBoreAnalysis();

    TPZWellBoreAnalysis(const TPZWellBoreAnalysis &copy);

    TPZWellBoreAnalysis &operator=(const TPZWellBoreAnalysis &copy);
    
    ~TPZWellBoreAnalysis();
    
    /// write the object on the stream
    void Write(TPZStream &output);
    
    /// read the object from the stream
    void Read(TPZStream &input);
    
    /// Computes the tension state transferring the geological stress state to the hidrostatic stress state
    /**
     * @param nsteps number of loading steps (the actual number of loading steps is one higher)
     * @param numnewton number of allowed newton iterations
     */
    void ExecuteInitialSimulation(int nsteps, int numnewton);
    
    /// Computes an equilibrium state corresponding to the current boundary conditions
    void ExecuteSimulation();
    
    /// Computes the given pressure in the specified steps
    void ExecuteSimulation(int steps,REAL pwb);
    
    /// verify the integrity of the elasto plastic material that is being used
    static void CheckDeformation(std::string filename = "deform.nb");
    
    /// verify if the stress tangent is computed correctly
    void VerifyTangentValidity();
    
    /// transfer the solution from the current configuration to the given configuration
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
    
    /// Set the polynomial order of the elements which exceed plastic deformation
    unsigned int PRefineElementAbove(REAL sqj2, int porder)
    {
        std::set<long> elindices;
        fCurrentConfig.PRefineElementsAbove(sqj2, porder,elindices);

        std::cout << "Number of elements prefined: " << elindices.size() << std::endl;
        // subject to integration points to the deformation history
        ApplyHistory(elindices);
        
        fCurrentConfig.fCMesh.Solution().Zero();
        fCurrentConfig.fAllSol = fCurrentConfig.fCMesh.Solution();

        return elindices.size();
    }
    
    /// Modify the geometry of the domain simulating an elliptic breakout
    void AddEllipticBreakout(REAL MaiorAxis, REAL MinorAxis);
    
    /// change the material id of the geometric elements of the current configuration
    void ChangeMaterialId(long idFrom, long idTo)
    {
        TPZGeoMesh *gmesh = &fCurrentConfig.fGMesh;
        long nel = gmesh->NElements();
        for (long iel=0; iel<nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel) continue;
            int matid = gel->MaterialId();
            if (matid == idFrom) {
                gel->SetMaterialId(idTo);
            }
        }
    }
    
    /// Divide the element using the plastic deformation as threshold
    unsigned int DivideElementsAbove(REAL sqj2);
    
    /// Post process the results of the current configuration
    void PostProcess(int resolution);
    
    /// Get the post processing variables
    static void PostProcessVariables(TPZStack<std::string> &scalnames, TPZStack<std::string> &vecnames);
    
    /// GetPostProcessedValues
    void PostProcessedValues(TPZVec<REAL> &x, TPZVec<std::string> &variables, TPZFMatrix<STATE> &values);
    
    /// verify the global equilibrium for the current configuration
    void VerifyGlobalEquilibrium(std::ostream &out = std::cout)
    {
        fCurrentConfig.VerifyGlobalEquilibrium(out);
    }
    
    /// Access method
    TConfig * GetCurrentConfig ()
    {
        return &fCurrentConfig;
    }

    /// Access method
    TConfig * GetConfig (int index) {
        if (index < 0 || index >= fSequence.size())
            DebugStop();

        list<TPZWellBoreAnalysis::TConfig>::iterator inte;
        int i=0;
        for (inte=fSequence.begin(); inte!=fSequence.end(); ++inte, i++)
            if (i == index)
                return &(*inte);

        DebugStop();
        return 0;
    }

    /// Return size of config list
    int GetConfigListSize () {
        return fSequence.size();
    }

    /// Initialize the object with standard parameters
    static void StandardConfiguration(TPZWellBoreAnalysis &obj);
    
    /// Configure the wellbore analysis to perform a linear analysis
    void LinearConfiguration(int porder);

    /// Define the geometry of the simulation
    void SetInnerOuterRadius(REAL inner, REAL outer)
    {
        fCurrentConfig.fInnerRadius = inner;
        fCurrentConfig.fOuterRadius = outer;
    }
    
    void SetMeshTopology(REAL delx, TPZVec<int> &nx)
    {
        fCurrentConfig.fDelx = delx;
        fCurrentConfig.fNx = nx;
    }
    
    /// Define the geological stress state and well pressure
    void SetConfinementStresses(TPZVec<STATE> &stress, STATE wellpressure)
    {
        fCurrentConfig.fConfinement.XX() = stress[0];
        fCurrentConfig.fConfinement.YY() = stress[1];
        fCurrentConfig.fConfinement.ZZ() = stress[2];
        fCurrentConfig.fFluidPressure = wellpressure;
    }
    
    void SetPorePressure(STATE porepressure)
    {
        fCurrentConfig.fPorePressure = porepressure;
    }
    
   void SetSanderDiMaggioParameters(REAL poisson, REAL Elast, REAL A, REAL B, REAL C, REAL R, REAL D, REAL W)
    {

#ifdef PV
        STATE G=Elast/(2.*(1.+poisson));
        STATE K=Elast/(3.*(1.-2*poisson));
        STATE phi=0,psi=1.,N=0.;
        fCurrentConfig.fSDPV.fYC.SetUp( A,  B, C,  D, K, G, W, R, phi, N, psi);
        fCurrentConfig.fSDPV.fER.SetUp(Elast,poisson);
#else
        fCurrentConfig.fSD.SetUp(poisson, Elast , A, B, C, R, D, W);
#endif


    }

    void Print(std::ostream &out);
    
private:
    
    /// Compute the linear elastic stiffness matrix
    void ComputeLinearMatrix();
    
    /// Set the parameters of the linear material
    void ConfigureLinearMaterial(TPZElasticityMaterial &mat);
    
    /// Recompute the plastic memory of the integration points of these elements
    void ApplyHistory(std::set<long> &elindices);
    
public:
    /** by Caju 2013 */
    /// Returns a set of points that belongs to the isoline defined by the ginen J2 value
    void GetJ2Isoline(REAL J2val, std::multimap<REAL,REAL> & polygonalChain);

    int GetPostProcessNumber () {
        return fPostProcessNumber;
    }

    void SaveConfig(stringstream &strout);
    
private:
    /// The object with the current configuration
    TConfig fCurrentConfig;
    
    /// The list of all previous configurations
    std::list<TConfig> fSequence;
    
    /// Index associated with the post processing file
    int fPostProcessNumber;
    
    /// Linear Elastic Stiffness matrix
    TPZAutoPointer<TPZMatrix<STATE> > fLinearMatrix;
    
    
};

// External variables
extern int startfrom;

#endif
