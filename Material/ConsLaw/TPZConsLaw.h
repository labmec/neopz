/**
 * @file TPZConsLaw.h
 * @brief Contains the TPZConservationLaw class which implements the interface for conservation laws.
 */

#ifndef PZCONSLAW_H
#define PZCONSLAW_H

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatInterfaceSingleSpace.h"
#include "TPZConsLawTypes.h"
/**
 * @brief Implements the interface for conservation laws, keeping track of the timestep as well
 */
/**
 * Defines the aditional interface necessary to compute contributions over the interfaces between elements
 */
class TPZConservationLaw  :
	public TPZMatBase<STATE,
					  TPZMatSingleSpaceT<STATE>,
					  TPZMatInterfaceSingleSpace<STATE>>
{
	using TBase = TPZMatBase<STATE,
							 TPZMatSingleSpaceT<STATE>,
							 TPZMatInterfaceSingleSpace<STATE>>;
public:
	/** @brief Simple constructor with material id, time step (dt) and dimension of the spatial domain */
	TPZConservationLaw(int nummat,REAL timeStep,int dim);
	
	/** 
	 * @name Attributes and parameters
	 * @{ */
	
	/** @brief Returns the dimension of the problem */
	int Dimension() const override;
	
	/** @brief Returns the value of the time step */
	REAL TimeStep();
	
	/** 
	 * @brief Sets the time step used for time integration
	 * @param[in] timeStep Time step (dt) 
	 */
	void SetTimeStep(REAL timeStep);
	
	/**
	 * @brief Sets the time step used for time integration
	 * @return Returns the resultant time step.
	 * @param[in] maxveloc Maximal speed of flow inside the cell
	 * @param[in] deltax Greatest dimension
	 * @param[in] degree Interpolation degree
	 */
	virtual REAL SetTimeStep(REAL maxveloc,REAL deltax,int degree)=0;
	
	/** @brief Returns the CFL number */
	REAL CFL();
	
	/**
	 * @brief Sets the CFL number
	 * @param[in] CFL Value of the CFL condition
	 */
	void SetCFL(REAL CFL);
	
	/** @brief Returns the value of Gamma (constant of gas) */
	REAL Gamma();
	
	/**
	 * @brief Sets the value of Gamma (constant of gas)
	 * @param[in] gamma Gamma value to Euler equation.
	 */
	void SetGamma(int gamma);
	
	/** @brief Sets whether the contribution is advanced or referring to the last state. */
	void SetContributionTime(TPZContributeTime time);
	
	/** @brief Residual_RT for calculations and Flux_RT for convergence check. */
	void SetResidualType(TPZResidualType type);
	
	/** @brief Number of state variables according to the dimension */
	int NStateVariables() const override = 0;
	
	/**
	 * @brief Thermodynamic pressure determined by the law of an ideal gas
	 * @param[in] U vector of state variables (sol)
	 */
	virtual STATE Pressure(TPZVec<STATE> &U)=0;
	
	/**
	 * @brief Prints the state of internal variables
	 * @param out Output to print the information of the current object
	 */
	void Print(std::ostream & out) const override;
	
	/** @brief Returns the material name */
	std::string Name() const override = 0;
	
	/**
	 * @brief Returns the relative index of a variable according to its name
	 * @param[in] name Name of the variable wished.
	 */
	int VariableIndex(const std::string &name) const override = 0;
	
	int NSolutionVariables(int var) const override = 0;

	/** 
	 * @brief Returns the solution associated with the var index based on
	 * the finite element approximation
	 * @param[in] data Material data to compute the solution.
	 * @param[in] var Number of the variable wished
	 * @param[out] Solout Vector with the computed solution values
	 */
	virtual void Solution(const TPZMaterialDataT<STATE> &data, int var, TPZVec<STATE> &Solout) override = 0;
	
	/** @} */
	
	/** @name Contribute methods */
	/** @{ */
	
	/** @brief Contributes to the residual vector and tangent matrix the volume-based quantities. */
	virtual void Contribute(const TPZMaterialDataT<STATE> &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override = 0;
	/** @brief Contributes to the residual vector and tangent matrix the face-based quantities. */
	virtual void ContributeInterface(const TPZMaterialDataT<STATE> &data,
									 const TPZMaterialDataT<STATE> &dataleft,
									 const TPZMaterialDataT<STATE> &dataright,
									 REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override = 0;
	/** @brief Contributes to the residual vector the boundary conditions */
	virtual void ContributeBC(const TPZMaterialDataT<STATE> &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
							  TPZBndCondT<STATE> &bc) override = 0;
	/** @} */
	
	int ClassId() const override;
	
	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
	
protected:
	
	/** @brief Dimension of the problem */
	int fDim;
	
	/** @brief Time step used for time integration */
	REAL fTimeStep;
	
	/** @brief CFL number */
	REAL fCFL;
	
	/** 
	 * @brief Ratio between specific heat is constant and the specific heat the constant
	 * volume of a polytropic gas 
	 */
	REAL fGamma;
	
	/** @brief Variable indicating the context of the solution. */
	/** If advanced, then the implicit terms are to be
	 * contributed. If last, then the explicit.
	 */
	TPZContributeTime fContributionTime;
	
	/** @brief Variable to indicate the type of residual to be computed by Assemble. */
	/** 
	 * A Flux Evaluation type is interesting for residual evaluation, and
	 * a complete residual for the global invertion.
	 */
	TPZResidualType fResidualType;
	
};

inline int TPZConservationLaw::Dimension() const
{
	return fDim;
}

inline REAL TPZConservationLaw::CFL()
{
	return fCFL;
}

inline void TPZConservationLaw::SetCFL(REAL CFL)
{
	//if(CFL > 1e3) CFL = 1e3;
	fCFL = CFL;
}

inline void TPZConservationLaw::SetGamma(int gamma)
{
	fGamma = gamma;
}

inline REAL TPZConservationLaw::Gamma()
{
	return fGamma;
}

inline void TPZConservationLaw::SetTimeStep(REAL timeStep)
{
	fTimeStep = timeStep;
}

inline REAL TPZConservationLaw::TimeStep()
{
	if(fResidualType == Residual_RT)return fTimeStep;
	return 1.;
}

inline void TPZConservationLaw::SetContributionTime(TPZContributeTime time)
{
	fContributionTime = time;
}


inline void TPZConservationLaw::SetResidualType(TPZResidualType type)
{
	fResidualType = type;
}

#endif
