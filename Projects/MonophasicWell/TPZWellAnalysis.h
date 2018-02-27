#ifndef TPZDarcyAnalysisH
#define TPZDarcyAnalysisH


#include "pznonlinanalysis.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzstepsolver.h"


class TPZCompMesh;

/** @author Omar in 16/02/2015
 * @brief class which implements 2D analysis for monofasic axisimetric darcy flow
 */


class TPZWellAnalysis : public TPZNonLinearAnalysis
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
    
    
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 4> fmeshvec;
    
    /** @brief Cmesh for Darcy analysis */
    TPZCompMesh * fcmeshdarcy;
    
    /** @brief unknowns for n time step */
    TPZFMatrix<REAL> falphaAtn;
    
    /** @brief unknowns for n+1 time step */
    TPZFMatrix<REAL> falphaAtnplusOne;
    
    
public:
    
    /// Constructor which already sets the cmesh
    TPZWellAnalysis();
    
    /// Destructor
    ~TPZWellAnalysis();
    
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
    void GeometryLine(int ns);
    
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
     * Computes the saturation at shock using the Welge method
     */
    static REAL SwatShock(REAL epsilon, REAL ds);
    
    /**
     * Extract a value from a given list
     */
    static int Extract(REAL epsilon, TPZManVector<REAL> &list, REAL value);
    
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