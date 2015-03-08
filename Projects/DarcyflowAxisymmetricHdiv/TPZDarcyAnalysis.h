#ifndef TPZDarcyAnalysisH
#define TPZDarcyAnalysisH

#include "SimulationData.h"
#include "ReservoirData.h"
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
	TPZFMatrix<STATE> fResidualLastState;
    
    // Simulation data required for complete the analysis
    TPZAutoPointer<SimulationData> fSimulationData;
    
    // Reservoir datas required for metrial definition
    TPZVec<TPZAutoPointer<ReservoirData > > fLayers;
    
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> fmeshvec;
    
    /** @brief Cmesh for Darcy analysis */
    TPZCompMesh * fcmeshdarcy;
    
	
public:
	
	/// Constructor which already sets the cmesh
    TPZDarcyAnalysis(TPZAutoPointer<SimulationData> DataSimulation, TPZVec<TPZAutoPointer<ReservoirData> > Layers);
	
	/// Destructor
	~TPZDarcyAnalysis();
	
	/**
	 * Assemble the stiffness matrix and rhs
	 **/
	void Assemble();
	
	/**
	 * Assemble last step residuum
	 **/
	void AssembleLastStep();
	
	/**
	 * Assemble the Residuum
	 */
	void AssembleResidual();
	
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
    void TimeForward(TPZFMatrix<STATE> &AlphasAtNplusOne, TPZFMatrix<STATE> &AlphasAtN);
    
	/**
	 * Is is necessary to fill the vector FSolution with the correct alphaj of
	 * the initial condition.
	 */
	void InitializeFirstSolution(TPZFMatrix<STATE> &AlphasAtN, REAL &ReferencePressure);
	
    /**
     * Set the simulation data,
     */
    void SetSimulationData(TPZAutoPointer<SimulationData> SimData){fSimulationData = SimData;}

    /**
     * Get the simulation data,
     */
    TPZAutoPointer<SimulationData> GetSimulationData() {return fSimulationData;}
    
    /**
     * Rotate the geometric mesh around Z axis
     */
    void RotateGeomesh(REAL CounterClockwiseAngle);

    /**
     * Create geometric Mesh Based on layer average thickness and layer radius
     */
    void CreatedGeoMesh();
    
    /**
     * Uniform Refinement
     */
    void UniformRefinement(int nh);

    /**
     * Read the geometric from dump file
     */
    void ReadGeoMesh(std::string GridFileName);   

    /**
     * Print geometric
     */
    void PrintGeoMesh();

    /**
     * Create the computational mesh for flux Q
     */
    TPZCompMesh * CmeshFlux(int qorder);
    
    /**
     * Create the computational mesh for flux P
     */
    TPZCompMesh * CmeshPressure(int Porder);
    
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
    void CreateMultiphysicsMesh(int q, int p);
    
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
    
};

#endif
