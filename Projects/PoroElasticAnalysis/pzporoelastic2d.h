/*
 *  pzporoelastic2d.cpp
 *  PZ
 *
 *  Created by Agnaldo on 11/28/11.
 *	Modified and Improved by Omar Duran on 11/28/11.
 *  Copyright 2012 L@bMeC. All rights reserved.
 *
 */

#ifndef POISSONDESACOPLADOS_HH
#define POISSONDESACOPLADOS_HH

#include "TPZMaterial.h"

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

class TPZPoroElastic2d : public TPZMaterial {
	
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
	
	/** @brief Elasticity modulus */
	REAL fE;
	REAL fEu;	
	
	/** @brief Poison coeficient */
	REAL fnu;
	REAL fnuu;	
	
	/** @brief first Lame Parameter */
	REAL flambda;
	REAL flambdau;	
	
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
	
    int ClassId() const override;
        
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name()  override { return "TPZPoroElastic2d"; }
	
	int Dimension() const  override {return fDim;}
	
	virtual int NStateVariables() const  override {return 3;}

	
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
	void SetParameters(REAL Lambda, REAL mu, REAL Lambdau,  REAL fx, REAL fy)
	{
		fE = (mu*(3.0*Lambda+2.0*mu))/(Lambda+mu);
		fnu = (Lambda)/(2*(Lambda+mu));
		
		fEu = (mu*(3.0*Lambdau+2.0*mu))/(Lambdau+mu);
		fnuu = (Lambdau)/(2*(Lambdau+mu));
		
		flambdau = Lambdau;
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
	
	// Get Elastic Materials Parameters	
	void GetElasticParameters(REAL &Ey, REAL &EyU, REAL &nu, REAL &nuU, REAL &Lambda, REAL &LambdaU, REAL &G, REAL &Alpha)
	{
		Ey = fE;
		EyU = fEu;	
		nu =  fnu;
		nuU = fnuu;	
		Lambda =  flambda;
		LambdaU =  flambdau;	
		G = fmu;	
		Alpha = falpha;
	}		
	
	// Get Diffusion Materials Parameters	
	void GetDiffusionParameters(REAL &ALpha, REAL &Se, REAL &viscosity, REAL &Perm)
	{
		ALpha = falpha;
		Se = fSe;
		Perm = fk; 
		viscosity = fvisc;	
	}			
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;

	//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
	virtual void Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)  override {
        TPZMaterial::Solution(data,dataleftvec,datarightvec,var,Solout,Left,Right);
    }
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
							 REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec,
									 REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
        TPZMaterial::ContributeInterface(data,dataleftvec,datarightvec,weight,ek,ef);
    }
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef) override {
        TPZMaterial::ContributeInterface(data,dataleft,dataright,weight,ef);
    }
	
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
                                       REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
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
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override {
		DebugStop();
	}	
};
#endif
