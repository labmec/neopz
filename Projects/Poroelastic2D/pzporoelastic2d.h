/*
 *  pzporoelastic2d.cpp
 *  PZ
 *
 *  Created by Agnaldo on 11/28/11.
 *	Modified and Improved by Omar Duran on 11/28/11.
 *  Copyright 2012 L@bMeC. All rights reserved.
 *
 */

#ifndef POISSONDESACOPLADOSH
#define POISSONDESACOPLADOSH

#include "pzmaterial.h"
#include "pzdiscgal.h"
#include "pzvec.h"


#include <iostream>


/**
 * @ingroup material
 * @brief Description of Biot's (1941) Linear Poroelastic system
 */
/**
 **@ingroup Linear Elastic Equation
 * \f$  div(T(u)) + b = 0  ==> Int{Grad(v).T(u)}dx - Int{v.gN}ds  = Int{b.v}dx  \f$ (Eq. 1) 
 *
 *\f$ T(u) =  lambda*Trace(E(u)I + 2*mu*(E(u)) - \f$
 *
 *\f$ E(u) =  (1/2)(Grad(u) + Transpose(Grad(u)) \f$
 * 
 *@ingroup	Diffusion equation for monophasic slightly compressible flow (e.g. oil)
 *
 *\f$ -(k/mu)*Div(Grad(p))  = d/dt{Se*p + alpha*Div(u)} (Eq. 2)  \f$ 
 *
 */

class TPZPoroElastic2d : public TPZDiscontinuousGalerkin {
	
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
	
	/** @brief Elasticity modulus */
	REAL fE;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
	/** @brief first Lame Parameter */
	REAL flambda;
	
	/** @brief Second Lame Parameter */
	REAL fmu;	
	
	/** @brief constants Biot poroelasticity */
	REAL falpha; //parAmetro poroelestico de Biot-Willis [adimensional]
	REAL fSe; //ou 1/M coeficiente poroelastico de armazenamento a volume constante [adimensional]
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Permeability of the rock and fluid viscosity*/
	REAL fk; 
	REAL fvisc;
	REAL fK;
	
	/** @brief Uses plain stress 
	* @note \f$fPlaneStress = 1\f$ => Plain stress state 
	* @note \f$fPlaneStress != 1\f$ => Plain Strain state 
	*/
	int fPlaneStress;
	
	/// timestep [s]
	REAL fTimeStep;
	
	REAL fTimeValue;
	
	REAL fCoupledfactor;
	
	REAL fmatId;
	
	REAL ftheta;
	
	/** @brief State: Stiffness or Mass Matrix Calculations */
	enum EState { ELastState = 0, ECurrentState = 1 };
	static EState gState;
	
public:
	TPZPoroElastic2d();
	
	TPZPoroElastic2d(int matid, int dim);
	
	virtual ~TPZPoroElastic2d();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPoroElastic2d"; }
	
	int Dimension() {return fDim;}
	
	virtual int NStateVariables();
	
	void SetLastState(){ gState = ELastState; }
	void SetCurrentState(){ gState = ECurrentState; }

	/** @brief Parameters of rock and fluid: */
	void SetDimension(int dimension)
	{
		fDim = dimension;

	}	
	
	
	/** @brief Parameters of rock and fluid: */
	void SetParameters(REAL perm, REAL visc)
	{
		fk = perm;
		fvisc = visc;
		fK = perm/visc;
	}
	
	/** 
	 * @brief Set parameters of elastic material:
	 * @param First  Lame Parameter Lambda
	 * @param Second Lame Parameter Mu -> G
	 * @param fx forcing function \f$ -x = fx \f$ 
	 * @param fy forcing function \f$ -y = fy \f$
	 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
	 */
	void SetParameters(REAL Lambda, REAL mu,  REAL fx, REAL fy)
	{
		fE = (mu*(3.0*Lambda+2.0*mu))/(Lambda+mu);
		fnu = (Lambda)/(2*(Lambda+mu));
	
		flambda = Lambda;
		fmu = mu;
		ff[0] = fx;
		ff[1] = fy;
	}
	
	/** @brief Set falpha parameter
	 * @param alpha : constant poroelastic Biot [dimensionless]
	 * @param Se : Coeficiente poroelastico de armazenamento a volume constante [adimensional]
	 */
	void SetBiotParameters(REAL alpha, REAL Se)
	{
		falpha = alpha;
		fSe = Se; 
	}
	
	/** @brief Set plane problem  
	 * planestress = 1 => Plain stress state 
	 * planestress != 1 => Plain Strain state 
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
	
	/// Set the timestep and Discretisation in time for more details see Ref  Finite element Method in ethe  ,.... pag 84
	void SetTimeStep(REAL Delta, REAL theta)
	{
		ftheta = theta;
		fTimeStep = Delta;
	}
	
	/// Set the timestep
	void SetTimeValue(REAL TimeValue)
	{
		fTimeValue = TimeValue;
	}	
	
	/// Set the kind of coupling: coupled -> true, uncoupled -> false
	void SetCoupled(bool coupled)
	{
		if (coupled) {
			fCoupledfactor = -1.0;
		}
		else {
			fCoupledfactor = 1.0;
		}

	}
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);

	//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
	
//	virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<REAL> &Solout);	
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, 
							 REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                       REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc);
	
	/** @name Contribute methods
	 * @{
	 */
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef) {
		DebugStop();
	}
    
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
		DebugStop();
	}
	/**
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
//    virtual void ErrorMassCal(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol,TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,TPZVec<REAL> &uexact, TPZFMatrix<REAL> &duexact,TPZVec<REAL> &val);
	
	
};
#endif
