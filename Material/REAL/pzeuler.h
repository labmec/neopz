/**
 * \file
 * @brief Contains the TPZEulerEquation class which implements the weak statement of the compressible euler equations.
 */

#ifndef PZEULER_H
#define PZEULER_H

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"


#include "pzausmflux.h"
#include "pzgradientflux.h"
class TPZCompMesh;

/**
 * @ingroup material
 * @brief This material implements the weak statement of the three-dimensional compressible euler equations
 */
/**
 * It is for transient analysis, finite volume method and explicit time integrator only.
 * @note For Olivier Roussel's project.
 *
 * \par Conservation Law:
 * \f$ U_t(t, X) + \nabla \cdot F(U(t, X)) = 0 \f$
 *
 * \par Variables:
 * \f$ U = (\rho, \rho * \upsilon_x, \rho * \upsilon_y, \rho * \upsilon_z, E) \f$ \n
 * \f$ \rho = \f$ density \n
 * \f$ (\upsilon_x, \upsilon_y, \upsilon_z) = \f$ velocity \n
 * \f$ \| \mbox{velocity} \| = \| v \| = \sqrt(\upsilon_x^2 + \upsilon_y^2 + \upsilon_z^2 ) \f$ \n
 * \f$ E = \f$ energy \n
 * \f$ p = \f$ pressure \n
 * \f$ c_p = \f$ specific heat at constant pressure of the gas \n
 * \f$ c_V = \f$ specific heat at constant volume of the gas.
 *
 * \par Fluxes:
 * \f$ F_x(U) = (\rho * \upsilon_x, \rho * \upsilon_x^2 + p, \rho * \upsilon_x * \upsilon_y, \rho * \upsilon_x * \upsilon_z, \upsilon_x (E + p)) \f$ \n
 * \f$ F_y(U) = (\rho * \upsilon_y, \rho * \upsilon_x * \upsilon_y, \rho * \upsilon_y^2 + p, \rho * \upsilon_y * \upsilon_z, \upsilon_y (E + p)) \f$ \n
 * \f$ F_z(U) = (\rho * \upsilon_z, \rho * \upsilon_x * \upsilon_z, \rho * \upsilon_y * \upsilon_z, \rho * \upsilon_z^2 + p, \upsilon_z (E + p)) \f$
 *
 * \par Equation of state:
 * \f$ p = (\gamma - 1) \star (E - \frac{1}{2} \rho \| v \| ) \f$ \n
 * \f$ \gamma = \frac{c_p}{c_V} \f$ \n
 *
 * @note For tests we used \f$ \gamma = 1.4 \f$ .
 */
class TPZEulerEquation : public TPZMaterial{
	
public:
	/** @brief Type of boundary condition */
	enum BCType{EFreeSlip = 1};
	/** @brief Type of flux computing */ 
	enum CALCType{EFlux = 1, EGradient = 2};
	/** @brief Set type of flux computing */
	static void SetComputeFlux(){
		gType = EFlux;
	}
	
	static void SetComputeGradient(){
		gType = EGradient;
	}
	
#ifdef LinearConvection
	static void SetLinearConvection(TPZCompMesh * cmesh, TPZVec<REAL> &Celerity);
#endif
	
private:
	
#ifdef LinearConvection
	TPZVec<STATE> fCelerity;
#endif
	
	static CALCType gType;
	
	/** @brief Ratio between specific heat at constant pressure and the specific heat at constant volume of a polytropic gas */
	static STATE gGamma;
	
	/** @brief Convective flux object */
	TPZAUSMFlux fAUSMFlux;
	
	/** @brief Gradient flux object */
	TPZGradientFlux fGradientFlux;
	
	/** @brief Compute Euler Flux */
	void ComputeEulerFlux(TPZVec<STATE> &sol, TPZFMatrix<STATE> & F);
	
public:
	
	static REAL Gamma(){ return gGamma; }
	
	/** @brief Convert from primitive to conservative variables */
	static void FromPrimitiveToConservative(TPZVec<STATE> &sol,STATE gamma);
	
	/** @brief Convert from conservative to primitive variables */
	static void FromConservativeToPrimitive(TPZVec<STATE> &sol,STATE gamma);
	
	/** @brief Constructor with Gamma value */
	TPZEulerEquation(int nummat, STATE gamma);
	
	/** @brief Default destructor */
	~TPZEulerEquation();
	
	/** @brief Default constructor */
	TPZEulerEquation();
	
	/** @brief Copy constructor */
	TPZEulerEquation(const TPZEulerEquation &cp);
	
	/** @brief Creates a copy of this */
	TPZMaterial * NewMaterial() override;
	
	/** @brief Object-based overload */
	virtual int NStateVariables() const override;
	
	/** @brief Object-based overload */
	virtual int Dimension() const override;
	
	/** @brief Returns the pressure value */
	static REAL Pressure(TPZVec<STATE> &U, double gamma);
	
	/** @brief Computes sound speed */
	STATE cSpeed(TPZVec<STATE> & sol);
	
	/** @brief Returns \f$ u = Sqrt(u2 + v2 + w2) \f$ */
	STATE uRes(TPZVec<STATE> & sol);
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name() override {return "TPZEulerEquation";}
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;

	/** 
	 * @name Contribute methods 
	 * @brief data contains material data. data.soll and data.solr are expected in primitive variables
	 */
	/** @{ */
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override;
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override;
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<STATE> &ef,
									   TPZBndCond &bc) override;
    public:
virtual int ClassId() const override;

	/** @} */
};

#endif///PZEULER_H

