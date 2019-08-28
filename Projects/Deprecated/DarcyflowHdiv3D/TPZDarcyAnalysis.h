#ifndef TPZDarcyAnalysisH
#define TPZDarcyAnalysisH

#include "SimulationData.h"
#include "ReservoirData.h"
#include "PetroPhysicData.h"
#include "ReducedPVT.h"
#include "pznonlinanalysis.h"

#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzstepsolver.h"


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
    
    // Reduced PVT data required for
    TPZAutoPointer<ReducedPVT> fFluidData;
    
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> fmeshvec;
    
    /** @brief Cmesh for Darcy analysis */
    TPZCompMesh * fcmeshdarcy;
    
    /** @brief unknowns for n time step */
    TPZFMatrix<REAL> falphaAtn;
    
    /** @brief unknowns for n+1 time step */
    TPZFMatrix<REAL> falphaAtnplusOne;
    
    
public:
    
    /// Constructor which already sets the cmesh
    TPZDarcyAnalysis(TPZAutoPointer<SimulationData> DataSimulation, TPZVec<TPZAutoPointer<ReservoirData> > Layers, TPZVec<TPZAutoPointer<PetroPhysicData> > PetroPhysic, TPZAutoPointer<ReducedPVT> FluidModel);
    
    /// Destructor
    ~TPZDarcyAnalysis();
    
    /**
     * Assemble the stiffness matrix and rhs
     **/
    void Assemble();
    
    /**
     * Assemble last step residuum
     **/
    void AssembleLastStep(TPZAnalysis *an);
    
    /**
     * Assemble the Residuum
     */
    void AssembleNextStep(TPZAnalysis *an);
    
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
     * Set the fluid model,
     */
    void SetFluidData(TPZAutoPointer<ReducedPVT> PVTData){fFluidData = PVTData;}
    
    /**
     * Get the fluid model,
     */
    TPZAutoPointer<ReducedPVT> GetFluidData() {return fFluidData;}
    
    
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
    void Geometry(int nx, int ny, int nz);
    
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
    TPZCompMesh * CmeshMixed();
    
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
    static  void ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * Parametric function to compute y direction
     */
    static  void ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X);
    
    /**
     * Parametric function to compute z direction
     */
    static  void ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X);

    /**
     * Rhs function of the mass conservation equation
     */
    static  void Ffunction(const TPZVec<REAL> &pt, TPZVec<STATE> &ff);
    
    /**
     * Exact Soltuion for linear tracer
     */
    static  void AnalyticSolution(const TPZVec<REAL> &pt, REAL time, TPZVec<STATE> &Saturation, TPZFMatrix<STATE> &Grad);
    
    /**
     * Computes the inverse of the Global matrix
     */
    TPZFMatrix<STATE> * ComputeInverse();
    
    /**
     * FilterEquations
     */
    void FilterSaturations(TPZManVector<int64_t> &active, TPZManVector<int64_t> &nonactive);
    
    /**
     * Update state variables
     */
    void UpDateAlphaVec(TPZFMatrix<REAL> &alpha);
    
};

#endif