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
#include "TPZYCMohrCoulombPV.h"

class TPZElasticityMaterialSest2D;

/// create de standard mesh
void AddBoundaryConditions(TPZCompMesh *CMesh, TPZMaterial * mat, TPZTensor<STATE> &Confinement, STATE pressure);

/// Simulation models
enum EPlasticModel  {ESandler, EMohrCoulomb, EElastic};

/// Fluid model indicates whether the fluid pressure evolves around the well
enum EFluidModel  {ENonPenetrating=0,EPenetrating=1};

/// Indication of the well configuration
enum EWellConfiguration {ENoConfig, EVerticalWell, EHorizontalWellalongh, EHorizontalWellalongH};

/// Class which simulates the stability of a wellbore
class TPZWellBoreAnalysis
{
    
public:
    
    /// order in which formation stresses are stored
    enum MFormationOrder {ESh = 0, ESH = 1, ESV = 2};
    
    /// boundary condition numbers
    enum MWellBCs {EInner = -2, EBottom = -4, ELeft = -5, EOuter = -3};
    
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
        
        /// Load the solution stored in TConfig into the CompMesh and set the ZDeformation of the material
        void LoadSolution();
        
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
        
        /// Return the mesh used for computations (multiphysics mesh or fCMesh)
        TPZCompMesh *CompMeshUsed();
        
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
        
        /// Transform from physical domain to computational domain stress
        /**
         * the outcome depends on the well configuration
         */
        void FromPhysicalDomaintoComputationalDomainStress(TPZTensor<STATE> &physicalStress, TPZTensor<STATE> &computationalStress);
        
        /// print the configuration
        void Print(ostream &out);
        
        /// this method will modify the boundary condition of the computational mesh and the forcing function
        // factor is a transition parameter between the confinement tension and well pressure
        // factor = 1 corresponds to pure well pressure
        void SetWellPressure(STATE factor = 1.);
        
        /// this method will configure the forcing function and boundary condition of the computational mesh
        void ConfigureBoundaryConditions();
        
        /// Set the Z deformation (for adapting the compaction)
        void SetZDeformation(STATE epsZ);
        
        /// Compute the average vertical stress of the configuration
        STATE AverageVerticalStress();
        
        /// Return gmesh
        TPZGeoMesh * GetGeoMesh();

        // Geometry of mesh
        /// radius of the well
        REAL fInnerRadius;
        /// radius of the computational domain, whether it is vertical or horizontal
        REAL fOuterRadius;
        
        /// Variable defining the well configuration
        EWellConfiguration fWellConfig;
        
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
        
        /// confinement stress in the physical domain
        TPZTensor<STATE> fConfinementTotal;
        
        /// Parameters
        //TPZElasticityMaterial fEl;
        
        /// Plastic model
        EPlasticModel fModel;

        /// Fluid model
        EFluidModel fFluidModel;

        /// Biot coefficient
        REAL fBiotCoef;
        
#ifdef PV
        // Sandler Dimaggio PV
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> fSDPV;

        //Mohr PV
        TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> fMCPV;
#else
        /// Parameters of the Sandler DiMaggio plasticity
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSD;
#endif

        /// Wellbore effective pressure
        REAL fWellborePressure;
        
        /// Far field pore pressure
        REAL fReservoirPressure;
        
        /// Geometric mesh3
        TPZGeoMesh fGMesh;
        
        /// Computational mesh
        TPZCompMesh fCMesh;
        
        /// Matrix of incremental solutions
        TPZFMatrix<STATE> fSolution;
        
        /// Z Deformation associated with the solution
        STATE fZDeformation;
        
        /// Vector containing maximum element plastic deformation
        TPZVec<REAL> fPlasticDeformSqJ2;
        
        /// The post processing mesh with the transferred solution
        TPZPostProcAnalysis fPostprocess;

        std::string fHistoryLog;
        
        static int gNumThreads;

        
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

    
    /// Evelves the reservoir and well pressure to target pressure
    void EvolveBothPressures(int nsteps, STATE TargetWellborePressure, STATE TargetReservoirPressure);

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
#ifdef DEBUG
        std::cout << "Number of elements prefined: " << elindices.size() << std::endl;
#endif
        // subject to integration points to the deformation history
        ApplyHistory(elindices);
        
        fCurrentConfig.fCMesh.Solution().Zero();
        fCurrentConfig.fSolution = fCurrentConfig.fCMesh.Solution();

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
    
    void PopConfiguration()
    {
        if (this->fSequence.size() ==0) {
            return;
        }
        fSequence.pop_back();
        fCurrentConfig = *(this->fSequence.rbegin());
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
    void SetConfinementTotalStresses(TPZVec<STATE> &stress, STATE wellpressure)
    {
        if (fCurrentConfig.fBiotCoef < 0.) {
            DebugStop();
        }
        fCurrentConfig.fConfinementTotal.XX() = stress[ESh];
        fCurrentConfig.fConfinementTotal.YY() = stress[ESH];
        fCurrentConfig.fConfinementTotal.ZZ() = stress[ESV];
        fCurrentConfig.fWellborePressure = wellpressure;
    }
    
    
    void SetReservoirPressure(STATE reservoirpressure)
    {
        fCurrentConfig.fReservoirPressure = reservoirpressure;
    }

    static TPZVec<STATE> FromTotaltoEffective(TPZVec<STATE> &totalstress, STATE biot, STATE totalporepressure)
    {
        int nel = totalstress.size();
        TPZVec<STATE> result(nel);
        for (int i=0; i<nel; i++) {
            result[i] = totalstress[i]-biot*totalporepressure;
        }
        return result;
    }
    
    static TPZVec<STATE> FromEffectivetoTotal(TPZVec<STATE> &effectivestress, STATE biot, STATE totalporepressure)
    {
        int nel = effectivestress.size();
        TPZVec<STATE> result(nel);
        for (int i=0; i<nel; i++) {
            result[i] = effectivestress[i]+biot*totalporepressure;
        }
        return result;
        
    }
    
    static STATE FromTotalPorePressuretoEffective(STATE totalpressure, STATE biot)
    {
        return totalpressure*(1-biot);
    }
    
    static STATE FromEffectivePorePressuretoTotal(STATE effectivepressure, STATE biot)
    {
        return effectivepressure/(1.-biot);
    }
    
    void SetBiotCoefficient(STATE biot)
    {
        fCurrentConfig.fBiotCoef = biot;
    }

    void SetFluidModel(EFluidModel fmodel)
    {
        fCurrentConfig.fFluidModel = fmodel;
        fCurrentConfig.SetWellPressure();
    }

    void SetWellConfig(EWellConfiguration fwconfig)
    {
        fCurrentConfig.fWellConfig = fwconfig;
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

        fCurrentConfig.fModel = ESandler;

    }
		
	void SetMohrCoulombParameters(REAL poisson, REAL Elast, REAL c, REAL Phi, REAL Psi)
	{
#ifdef PV
		TPZElasticResponse ER;
		ER.SetUp(Elast,poisson);
		fCurrentConfig.fMCPV.fYC.SetUp(Phi, Psi, c, ER);
		fCurrentConfig.fMCPV.fER.SetUp(Elast,poisson);

    fCurrentConfig.fModel = EMohrCoulomb;
#else
		PZError << "You have to define PV to use MohrCoulombPV!" << std::endl;
		DebugStop();
#endif
	}
	
    void Print(std::ostream &out);
    
private:
    
    /// Compute the linear elastic stiffness matrix
    void ComputeLinearMatrix(TPZVec<long> &activeEquations);
    
    /// Set the parameters of the linear material
    void ConfigureLinearMaterial(TPZElasticityMaterialSest2D &mat);
    
    /// Recompute the plastic memory of the integration points of these elements
    void ApplyHistory(std::set<long> &elindices);
    
public:
    /// Test the linear matrix with vertical compaction
    void TestLinearMaterial();
    
    /** by Caju 2013 */
    /// Returns a set of points that belongs to the isoline defined by the ginen J2 value
    void GetJ2Isoline(REAL J2val, std::multimap<REAL,REAL> & polygonalChain);
    
    void ComputeAandB(REAL sqj2_refine,REAL &a,REAL &b);

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
