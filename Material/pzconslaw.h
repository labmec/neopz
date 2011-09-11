/**
 * @file
 * @brief Contains the TPZConservationLaw class which implements the interface for conservation laws.
 */

#ifndef PZCONSLAW_H
#define PZCONSLAW_H

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzdiscgal.h"

/**
 * @ingroup material
 * @brief Indicates the type of time discretization
 */
enum TPZTimeDiscr
{
	None_TD = -1,
	Explicit_TD = 0,
	ApproxImplicit_TD = 1,
	Implicit_TD = 2,
	Unknown_TD = 3
};

/**
 * @ingroup material
 * @brief Indicates which term is put in the right hand side and tangent matrix
 */
enum TPZContributeTime
{
	None_CT = -1,
	Last_CT = 0,
	Advanced_CT = 1
};

/**
 * @ingroup material
 * @brief Which terms are being contributed
 */
enum TPZResidualType
{
	None_RT = -1,
	Residual_RT = 0,
	Flux_RT = 1
};

/**
 * @ingroup material
 * @brief Implements the interface for conservation laws, keeping track of the timestep as well
 */
/**
 * Defines the aditional interface necessary to compute contributions over the interfaces between elements
 */
class TPZConservationLaw  : public TPZDiscontinuousGalerkin
{
public:
	
	TPZConservationLaw(int nummat,REAL timeStep,int dim);
	
	TPZConservationLaw(const TPZConservationLaw &cp) : TPZDiscontinuousGalerkin(cp),
	fDim(cp.fDim),fTimeStep(cp.fTimeStep),fCFL(cp.fCFL), fGamma(cp.fGamma),fContributionTime(cp.fContributionTime)
	,fResidualType(cp.fResidualType)
	{
	}
	
	
	virtual ~TPZConservationLaw();
	
	//------------------attributes and parameters
	
	/** @brief Returns the dimension of the problem */
	int Dimension();
	
	/** @brief Returns the value of the time step */
	REAL TimeStep();
	
	/** 
	 * @brief Sets the time step used for time integration
	 * @param timeStep [in]
	 */
	void SetTimeStep(REAL timeStep);
	
	/**
	 * @brief Sets the time step used for time integration
	 * @return Returns the resultant time step.
	 * @param maxveloc [in] maximal speed of flow inside the cell
	 * @param deltax [in] greatest dimension
	 * @param degree [in] interpolation degree
	 */
	virtual REAL SetTimeStep(REAL maxveloc,REAL deltax,int degree)=0;
	
	/** @brief Returns the CFL number */
	REAL CFL();
	
	/**
	 * @brief Sets the CFL number
	 * @param CFL [in]
	 */
	void SetCFL(REAL CFL);
	
	/* *
	 * Returns the value of delta for the artificial diffusion term
	 */
	//  REAL Delta();
	
	/**
	 * Sets the value of the artificial diffusion delta parameter
	 *
	 * @ param delta [in]
	 */
	//  void SetDelta(REAL delta);
	
	
	/** @brief Returns the value of Gamma (constant of gas) */
	REAL Gamma();
	
	/**
	 * @brief Sets the value of Gamma (constant of gas)
	 * @param gamma [in]
	 */
	void SetGamma(int gamma);
	
	/** @brief Sets whether the contribution is advanced or referring to the last state. */
	void SetContributionTime(TPZContributeTime time);
	
	/** @brief Residual_RT for calculations and Flux_RT for convergence check. */
	void SetResidualType(TPZResidualType type);
	
	/** @brief Number of state variables according to the dimension */
	virtual int NStateVariables() = 0;
	
	/**
	 * @brief Thermodynamic pressure determined by the law of an ideal gas
	 * @param U [in] vector of state variables (sol)
	 */
	virtual REAL Pressure(TPZVec<REAL> &U)=0;
	
	/**
	 * @brief Prints the state of internal variables
	 * @param out [in]
	 */
	virtual void Print(std::ostream & out);
	
	/** @brief Returns the material name */
	virtual std::string Name() = 0;
	
	/**
	 * @brief Returns the relative index of a variable according to its name
	 * @param name [in]
	 */
	virtual int VariableIndex(const std::string &name)=0;
	
	virtual int NSolutionVariables(int var)=0;
	
	/** @brief Returns the number of fluxes associated to this material */
	virtual int NFluxes();
	
	
	//------------------solutions
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,
						  TPZFMatrix &axes,int var,
						  TPZVec<REAL> &Solout)=0;
public:
	/** 
	 * @brief Returns the solution associated with the var index based on
	 * the finite element approximation 
	 */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	
	/** @name Contribute methods */
	/** @{ */
	
	/** @brief Contributes to the residual vector and tangent matrix the volume-based quantities. */
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ek,TPZFMatrix &ef)=0;
	
	/** @brief Contributes to the residual vector and tangent matrix the volume-based quantities. */
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	/** @brief Contributes to the residual vector and tangent matrix the face-based quantities. */
	virtual void ContributeInterface(TPZMaterialData &data,
									 REAL weight,
									 TPZFMatrix &ek,TPZFMatrix &ef)=0;
	/** @brief Contributes to the residual vector and tangent matrix the face-based quantities. */
	virtual void ContributeInterface(TPZMaterialData &data,
									 REAL weight,
									 TPZFMatrix &ef)
	{
		TPZDiscontinuousGalerkin::ContributeInterface(data,weight,ef);
	}
	/** @brief Contributes to the residual vector the boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,TPZFMatrix &ef,
							  TPZBndCond &bc)=0;
	/** @brief Contributes to the residual vector the boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	/** @} */
	
	//--------------------
	//virtual int IntegrationDegree() = 0;
	
	//virtual void SetIntegDegree(int degree) = 0;
	
	//---------------------Attributes
	
	
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
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

inline int TPZConservationLaw::Dimension()
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
	//   std::cout << "CFL:"<<CFL << std::endl;
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
inline int TPZConservationLaw::NFluxes()
{
	return 1;
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
