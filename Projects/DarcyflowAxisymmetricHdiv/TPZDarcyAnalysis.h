#ifndef TPZDarcyAnalysisH
#define TPZDarcyAnalysisH

#include "SimulationData.h"
#include "ReservoirData.h"
#include "PetroPhysicData.h"

#include "ReducedPVT.h"
#include "OilPhase.h"
#include "WaterPhase.h"
#include "GasPhase.h"

#include "pznonlinanalysis.h"

#include "pzcondensedcompel.h"

#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzstepsolver.h"

#include "pzl2projection.h"
#include "pzgradientreconstruction.h"
#include "pzbfilestream.h"


class TPZCompMesh;

/** @author Omar in 16/02/2015
 * @brief class which implements 2D analysis for monofasic axisimetric darcy flow
 */


class TPZDarcyAnalysis : public TPZNonLinearAnalysis
{
private:
    
    // Sets the material on the last state (n)
    void SetLastState();
    
    // Sets the material on the next state (n+1)
    void SetNextState();
    
    // Vector which will store tha residuum in the last state (n)
    TPZFMatrix<STATE> fResidualAtn;
    
    // Vector which will store tha residuum in the last state (n+1)
    TPZFMatrix<STATE> fResidualAtnplusOne;
    
    // Simulation data required for complete the analysis
    TPZAutoPointer<SimulationData> fSimulationData;
    
    // Reservoir datas required for material definition
    TPZVec<TPZAutoPointer<ReservoirData > > fLayers;
    
    // PetroPhysic data associated to each reservoir layer
    TPZVec<TPZAutoPointer<PetroPhysicData > > fRockPetroPhysic;
    
    // Reduced PVT data required for alpha
    TPZAutoPointer<ReducedPVT> falpha_fluid;
    
    // Reduced PVT data required for beta
    TPZAutoPointer<ReducedPVT> fbeta_fluid;
    
    // Reduced PVT data required for gamma
    TPZAutoPointer<ReducedPVT> fgamma_fluid;
   
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 4> fmeshvecini;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 4> fmeshvec;

    /** @brief Cmesh for Initial Darcy analysis */
    TPZCompMesh * fcmeshinitialdarcy;
    
    /** @brief Cmesh for Darcy analysis */
    TPZCompMesh * fcmeshdarcy;
    
    /** @brief cmesh */
    TPZCompMesh * fcmesh;
    
    /** @brief unknowns for n time step */
    TPZFMatrix<REAL> falphaAtn;
    
    /** @brief unknowns for n+1 time step */
    TPZFMatrix<REAL> falphaAtnplusOne;
    
    /** @brief Store DOF associated with Constant Saturations */    
    TPZManVector<long> fConstantSaturations;
   
    /** @brief Store DOF associated with  Saturation gradients */    
    TPZManVector<long> fGradientSaturations;    
    
    
public:
    REAL Muo;
    /// Constructor which already sets the cmesh
    TPZDarcyAnalysis(TPZAutoPointer<SimulationData> DataSimulation, TPZVec<TPZAutoPointer<ReservoirData> > Layers, TPZVec<TPZAutoPointer<PetroPhysicData> > PetroPhysic);
    
    /// Destructor
    ~TPZDarcyAnalysis();
    
    /**
     * Assemble the stiffness matrix and rhs
     **/
    void Assemble();
    
    /**
     * Assemble last step residuus
     **/
    void AssembleLastStep(TPZAnalysis *an);
    
    /**
     * Assemble the Jacobian and Residuus
     */
    void AssembleNextStep(TPZAnalysis *an);
    
    /**
     * Assemble the Residuus
     */
    void AssembleResNextStep(TPZAnalysis *an);
    
    /**
     * Computes the residuum. Used for checkconv
     */
    void Residual(TPZFMatrix<STATE> &residual, int icase);
    
    /**
     * Sets next state and computes the tangent
     */
    void ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase);
    
    /**
     * Run the simulation,
     */
    void Run();
    
    /**
     * Compute the time forward at each timestep
     */
    void TimeForward(TPZAnalysis *an);
    
    /**
     * Is is necessary to fill the vector FSolution with the correct alphaj of
     * the initial condition.
     */
    void InitializeSolution(TPZAnalysis *an);
    
    /**
     * Set the simulation data,
     */
    void SetSimulationData(TPZAutoPointer<SimulationData> SimData){fSimulationData = SimData;}
    
    /**
     * Get the simulation data,
     */
    TPZAutoPointer<SimulationData> GetSimulationData() {return fSimulationData;}
    
    
    /**
     * Set the fluids following the system type
     */
    void SetFluidData(TPZVec< TPZAutoPointer<ReducedPVT> > PVTData);
    
    /**
     * Get the alpha fluid model,
     */
    TPZAutoPointer<ReducedPVT> GetFluidAlphaData() {return falpha_fluid;}
    
    /**
     * Get the beta fluid model,
     */
    TPZAutoPointer<ReducedPVT> GetFluidBetaData() {return fbeta_fluid;}
    
    /**
     * Get the gamma fluid model,
     */
    TPZAutoPointer<ReducedPVT> GetFluidGammaData() {return fgamma_fluid;}
    
    /**
     * Rotate the geometric mesh around Z axis
     */
    void RotateGeomesh(REAL CounterClockwiseAngle);
    
    /**
     * Create geometric Mesh Based on layer average thickness and layer radius
     */
    void CreatedGeoMesh();
    
    /**
     * Create geometric Mesh for one-dimensional displacement
     */
    void GeometryLine(int nx, int ny);
    
    /**
     * Uniform Refinement
     */
    void UniformRefinement(int nh);
    
    /**
     * Uniform Refinement for specific material
     */
    void UniformRefinement(int nh, int MatId);
    
    /**
     * Uniform Refinement towards specific materials
     */
    void UniformRefinement(int nh, std::set<int> &MatToRef);
    
    /**
     * Read the geometric from dump file
     */
    void ReadGeoMesh(std::string GridFileName);
    
    /**
     * Print geometric
     */
    void PrintGeoMesh();
    
    /**
     * Create the computational mesh for Q
     */
    TPZCompMesh * CmeshFlux(int qorder);
    
    /**
     * Create the computational mesh for P
     */
    TPZCompMesh * CmeshPressure(int Porder);
    
    /**
     * Create the computational mesh for Sw
     */
    TPZCompMesh * CmeshSw(int Sworder);
    
    /**
     * Create the computational mesh for So
     */
    TPZCompMesh * CmeshSo(int Soorder);
    
    /**
     * Create the computational mixed mesh
     */
    TPZCompMesh * CmeshMixedInitial();
    
    /**
     * Create the computational mixed mesh
     */
    TPZCompMesh * CmeshMixed();
    
    /**
     * Create an analysis based on a given cmesh
     */
    TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh);
    
    /**
     * Push the initial cmesh
     */
    void PushInitialCmesh();

    /**
     * Push the current cmesh
     */
    void PushCmesh();
    
    /**
     * Print the Cmesh object
     */
    void PrintCmesh();
    
    /**
     * Create the computational continuous mesh
     */
    void CmeshH1(int porder);
    
    /**
     * Create the computational mixed
     */
    void CreateMultiphysicsMesh(int q, int p, int s);
    
    /**
     * Create vtk file
     */
    void PostProcessVTK(TPZAnalysis *an);
    
    /**
     * Create interfaces in the multiphysics mesh
     */
    void CreateInterfaces();
    
    /**
     * Apply static condensation in the multiphysics mesh
     */
    void ApplyStaticCondensation();
    
    /**
     * Print the global linear system
     */
    void PrintLS(TPZAnalysis *an);
    
    /**
     * Computes the newton Iterations
     */
    void NewtonIterations(TPZAnalysis *an);
    
    /**
     * Computes the broyden Iterations
     */
    void BroydenIterations(TPZAnalysis *an);
    
    /**
     * Computes the outer product
     */
    TPZFMatrix<STATE>  TensorProduct(TPZFMatrix<STATE> &g, TPZFMatrix<STATE> &d);
    
    /**
     * Parametric function to compute x direction
     */
    static  void Parametricfunction(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * Parametric function to compute y direction
     */
    static  void Parametricfunction2(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * Rhs function of the mass conservation equation
     */
    static  void Ffunction(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
    
    /**
     * Exact Soltuion for linear tracer
     */
    static  void LinearTracer(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Saturation, TPZFMatrix<STATE> &Grad);
    
    /**
     * Exact Soltuion for bluckley and leverett
     */
    static void BluckleyAndLeverett(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Saturation, TPZFMatrix<STATE> &Grad);
    
    /**
     * Computes the saturation at shock using the Welge method
     */
    static REAL SwatShock(REAL epsilon, REAL &ds);
    
    /**
     * Extract a value from a given list
     */
    static int Extract(REAL epsilon, TPZManVector<REAL> &list, REAL value);

    /**
     * Computes the inverse of the Global matrix
     */
    static REAL SaturationNewton( REAL x,REAL t,REAL muo, REAL muw, REAL Area,REAL q);
    
    /**
     * Computes the inverse of the Global matrix
     */
    static REAL dfdsw( REAL Sw, REAL muo,REAL muw);

    /**
     * Computes the inverse of the Global matrix
     */
    static REAL df2dsw( REAL Sw, REAL muo,REAL muw);
    
    /**
     * Computes the inverse of the Global matrix
     */
    TPZFMatrix<STATE> * ComputeInverse();
    
    /**
     * FilterEquations
     */
    void FilterSaturationGradients(TPZManVector<long> &active, TPZManVector<long> &nonactive);
    
    /**
     * Compute saturation reconstruction for Sw and So
     */
    void SaturationReconstruction(TPZAnalysis *an);
    
    /**
     * Clean up reconstructed gradient saturation for Sw and So
     */
    void CleanUpGradients(TPZAnalysis *an);    
    
    /**
     * Update state variables
     */
    void UpDateAlphaVec(TPZFMatrix<REAL> &alpha);

    /**
     * Computes the L2 projection for saturations
     */
    void SolveProjection(TPZAnalysis *an, TPZCompMesh *Cmesh);

    /**
     * Computes computational mesh for L2 projection
     */
    TPZCompMesh * L2ProjectionCmesh(TPZVec<STATE> &solini);
    
    /**
     * Computes convergence rate for an element
     */
    void CheckElementConvergence(int wichelement);
    
    /**
     * Computes convergence rate for the global problem
     */
    void CheckGlobalConvergence(TPZAnalysis * an);

    /**
     * Method to verify the jacobian
     */
    void CheckGlobalJacobian(TPZAnalysis * an);
    
    /**
     * Computes computational mesh for L2 projection
     */
    static void InitialWaterSaturation(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
    
};

#endif