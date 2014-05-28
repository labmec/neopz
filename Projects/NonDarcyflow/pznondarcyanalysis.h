#ifndef TPZNonDarcyAnalysisH
#define TPZNonDarcyAnalysisH

#include "pznonlinanalysis.h"
#include "problemdata.h"


class TPZCompMesh;

/** @author Omar Shauer in 02/04/2014
 * @brief class which implements non 2D analysis to the monofasic axisimetric 
 */
class TPZNonDarcyAnalysis : public TPZNonLinearAnalysis
{
private:
	
	// Sets the material on the last state (n)
	void SetMaterialLastState();
	
	// Sets the material on the next state (n+1)
	void SetMaterialNextState();
	
	// Vector which will store tha residuum in the last state (n)
	TPZFMatrix<STATE> fResLastState;
	
public:
	
	/// Constructor which already sets the cmesh
	TPZNonDarcyAnalysis(TPZCompMesh *cmesh, std::ostream &out = std::cout);
	
	/// Destructor
	virtual ~TPZNonDarcyAnalysis();
	
	/**
	 * Assemble the stiffness matrix and rhs
	 **/
	virtual void Assemble();
	
	/**
	 * Assemble last step residuum
	 **/
	virtual void AssembleLastStep();
	
	/**
	 * Assemble the Residuum
	 */
	virtual void AssembleResidual();
	
	/**
	 * Computes the residuum. Used for checkconv
	 */
	virtual void Residual(TPZFMatrix<STATE> &residual, int icase);
	
	/**
	 * Sets next state and computes the tangent
	 */
	virtual void ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase);
	
	/**
	 * Run all the time steps based on pini in all the fracture,
	 */
	void RunAll();
	
	/**
	 * Realizes one time step where a pressure in a step n is provided and
	 * a pressure in the step n + 1 is computed
	 * Here the Iterative process is called.
	 */
	virtual void TimeStep(TPZFMatrix<STATE> &PressureN, TPZFMatrix<STATE> &PressureNp1);
	
	/**
	 * Is is necessary to fill the vector FSolution with the correct alphaj of
	 * the initial condition. This is made here.
	 */
	void InitializeFirstSolution(TPZFMatrix<STATE> &Pressure,REAL &pini);
	
	
};

#endif
