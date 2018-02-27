#ifndef TPZDarcyAnalysisH
#define TPZDarcyAnalysisH

#include "SimulationData.h"
#include "ReservoirData.h"
#include "PetroPhysicData.h"

#include "Phase.h"
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
 * @brief class which implements 2D analysis for multiphasic axisimetric darcy flow
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
    TPZAutoPointer<Phase> falpha_fluid;
    
    // Reduced PVT data required for beta
    TPZAutoPointer<Phase> fbeta_fluid;
    
    // Reduced PVT data required for gamma
    TPZAutoPointer<Phase> fgamma_fluid;
    
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of initial compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2, fmeshvec[2] = SaturationL2 */
    TPZManVector<TPZCompMesh * , 4> fmeshvecini;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2, fmeshvec[2] = SaturationL2 */
    TPZManVector<TPZCompMesh * , 4> fmeshvec;
    
    /** @brief Cmesh for Initial Darcy analysis */
    TPZCompMesh * fcmeshinitialdarcy;
    
    /** @brief Cmesh for Darcy analysis */
    TPZCompMesh * fcmeshdarcy;
    
    /** @brief cmesh */
    TPZCompMesh * fcmesh;
    
    /** @brief unknowns state vars for n time step */
    TPZFMatrix<REAL> falphaAtn;
    
    /** @brief unknowns state vars for n+1 time step */
    TPZFMatrix<REAL> falphaAtnplusOne;

    /** @brief unknowns saturations for n time step */
    TPZFMatrix<REAL> fSAtn;
    
    /** @brief unknowns saturations for n+1 time step */
    TPZFMatrix<REAL> fSAtnplusOne;
    
    /** @brief Store DOF associated with active */
    TPZManVector<int64_t> fActiveEquations;
    
    /** @brief Store DOF associated with  non active */
    TPZManVector<int64_t> fNonactiveEquations;
    
    
public:
    
    /** @brief L2 norm */
    TPZManVector<REAL> fL2_norm;
    
    /** @brief L2 norm */
    TPZManVector<REAL> fL2_norm_s;
    
    /** @brief Hdiv norm */
    TPZManVector<REAL> fl2_norm_flux;

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
    void RunAnalysis();
    
    /**
     * Compute the time forward at each timestep
     */
    void TimeForward(TPZAnalysis *an);
    
    /**
     * Is necessary to fill the vector FSolution with the correct alphaj of
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
    void SetFluidData(TPZVec< TPZAutoPointer<Phase> > PVTData);
    
    /**
     * Get the alpha fluid model,
     */
    TPZAutoPointer<Phase> GetFluidAlphaData() {return falpha_fluid;}
    
    /**
     * Get the beta fluid model,
     */
    TPZAutoPointer<Phase> GetFluidBetaData() {return fbeta_fluid;}
    
    /**
     * Get the gamma fluid model,
     */
    TPZAutoPointer<Phase> GetFluidGammaData() {return fgamma_fluid;}
    
    /**
     * Rotate the geometric mesh around Z axis
     */
    void RotateGeomesh(REAL CounterClockwiseAngle);
    
    /**
     * Apply shear deformation to the geometric mesh on x-y plane
     */
    void ApplyShear(REAL CounterClockwiseAngle);
    
    /**
     * Create geometric Mesh Based on layer average thickness and layer radius
     */
    void CreatedGeoMesh();
    
    /**
     * Create geometric Mesh for one-dimensional displacement
     */
    void Geometry2D(int nx, int ny);
    
    /**
     * Apply geometric progression over the given 1-D geometric mesh -> gmesh
     */
    void ApplyPG(TPZGeoMesh * geomesh);
    
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
     * Create the computational mesh for u
     */
    TPZCompMesh * CmeshFlux(int uorder);
    
    /**
     * Create the computational mesh for P
     */
    TPZCompMesh * CmeshPressure(int Porder);
    
    /**
     * Create the computational mesh for S_alpha
     */
    TPZCompMesh * CmeshSw(int S_alpha_order);
    
    /**
     * Create the computational mesh for S_beta
     */
    TPZCompMesh * CmeshSo(int S_beta_order);
    
    /**
     * Create the initial computational mixed mesh
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
     * Create the computational mixed mesh
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
     * Computes the picard Iterations for the impes like solver
     */
    void PicardIterations(TPZAnalysis *an);
    
    /**
     * Explicit update for the impes like picard iteration
     */
    void UpdateSaturations(TPZAnalysis *an);
    
    /**
     * Computes the explicit part for the impes like picard iteration
     */
    void ComputeSaturations(TPZAnalysis *an);
    
    /**
     * Computes the outer product
     */
    TPZFMatrix<STATE>  TensorProduct(TPZFMatrix<STATE> &g, TPZFMatrix<STATE> &d);
    
    /**
     * Parametric function to compute x direction
     */
    static  void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * Parametric function to compute y direction
     */
    static  void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * Rhs function of the mass conservation equation
     */
    static  void Ffunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad);
    
    /**
     * Transient BC neumann
     */
    static  void BCNfunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad);
    
    /**
     * Transient BC dirichlet
     */
    static  void BCDfunction(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &ff, TPZFMatrix<REAL> &Grad);

    /**
     * Exact Soltuion elliptic axisymmetric darcflow
     */
    static  void Cylindrical_Elliptic(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &GradSol);
    
    /**
     * Exact Soltuion elliptic axisymmetric darcflow (Dupuit-Thiem solution)
     */
    static  void Dupuit_Thiem(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &GradSol);
    
    /**
     * Exact Soltuion parabolic axisymmetric darcflow (Morris_Muskat solution)
     */
    static  void Morris_Muskat(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &GradSol);
    
    /**
     * Exact Soltuion for linear tracer
     */
    static  void LinearTracer(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Saturation, TPZFMatrix<STATE> &Grad);
    
    /**
     * Exact Soltuion for bluckley and leverett
     */
    static void BluckleyAndLeverett(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Saturation, TPZFMatrix<STATE> &Grad);
    
    /**
     * Computes the inverse of the Global matrix
     */
    static REAL S_Newton(REAL x, REAL t, REAL u, REAL Swr, REAL Sor, REAL phi, REAL s_shok, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta, REAL epsilon);
    
    /**
     * Computes the inverse of the Global matrix
     */
    static REAL dfdsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta);
    
    /**
     * Computes the inverse of the Global matrix
     */
    static REAL df2dsw(REAL Sw, REAL Swr, REAL Sor, REAL mu_alpha, REAL mu_beta, REAL rho_alpha, REAL rho_beta);
    
    /**
     * Computes the inverse of the Global matrix
     */
    TPZFMatrix<STATE> * ComputeInverse();
    
    /**
     * FilterEquations
     */
    void FilterSaturationGradients(TPZManVector<int64_t> &active, TPZManVector<int64_t> &nonactive);
    
    /**
     * FilterEquations
     */
    void FilterSaturations(TPZManVector<int64_t> &active, TPZManVector<int64_t> &nonactive);
    
    /**
     * Compute saturation reconstruction for Sw and So
     */
    void SaturationReconstruction(TPZAnalysis *an);
    
    /**
     * Clean up reconstructed gradient saturation for Sw and So
     */
    void CleanUpGradients(TPZAnalysis *an);
    
    /**
     * Print saturations and gradients
     */
    void PrintSaturations(TPZAnalysis *an);
    
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
    static void InitialS_alpha(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

    /**
     * Computes computational mesh for L2 projection
     */
    static void P_Hydrostatic(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& P_Hydro, TPZFMatrix< STATE >& GradP_Hydro);

    /**
     * Computes the integral of the velocities for each phase
     */
    void IntegrateVelocities(TPZManVector<REAL> & velocities);
    
    /**
     * Computes the integral of the hdiv and L2 error of flux and pressure respectively
     */
    void IntegrateFluxPError(TPZManVector<REAL> & l2_norm_flux,TPZManVector<REAL> & l2_norm);
    
    /**
     * Computes the integral of the L2 error of saturation
     */
    void IntegrateL2SError(TPZManVector<REAL> & l2_norm);
    
    /**
     * Computes the integral of the production using Riemman rule
     */
    TPZFMatrix<REAL> RiemmanRule(TPZFMatrix<REAL>  current_production);
    
    /**
     * Computes the integral of the production using Trapezoidal rule
     */
    TPZFMatrix<REAL> TrapezoidalRule(TPZFMatrix<REAL>  current_production);
    
    /**
     * Comute the neigh of higher dimension
     */
    TPZGeoEl * GetVolElement(TPZGeoEl * gel);
    
    /**
     * read the rasterized list of k values
     */
    void ReadRasterized();
    
    /**
     * Compute the afin transformation from origin to destination
     */
    TPZTransform<> Transform_1D_To_2D(TPZGeoEl * gel_o, TPZGeoEl * gel_d, TPZGeoElSide & intermediate_side);
    
    
};

#endif