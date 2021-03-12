//
//  pzporoelasticmf2d.h
//  PZ
//
//  Created by Agnaldo Farias on 7/20/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

/**
 * @file
 * @brief Contains the TPZPoroElasticMF2d class which implements the poroelastic 2d problem for simulation multi-phyisics,
 for the variables of displacement, flow and pressure.
 */

#ifndef PZ_pzporoelasticmf2d_h
#define PZ_pzporoelasticmf2d_h


#include "TPZMaterial.h"

#include "pzvec.h"

#include <iostream>


/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
/**
 **@ingroup equacao da elasticidade
 * \f$  div(T(u))  + fXf2 = 0  ==> Int{Grad(v).T(u)}dx - Int{v.gN}ds  = Int{ff.v}dx  \f$ (Eq. 1) 
 *
 *\f$ T(u) =  tr(E(u) - alpha*p*I)lambda*I + 2*nu*(E(u) - alpha*p*I)\f$
 *
 *@ingroup equacao da pressao
 * \f$ -1/visc*div(gradu k)  = 0 ==> k/visc*Int{Grad(u)Grad(v)}dx - Int{k/visv*Grad(u).n v}ds  = 0   (Eq. 2)  \f$ 
 *
 */

class TPZPoroElasticMF2d : public TPZMaterial {
	
protected:
	
	/** @brief Forcing vector */
	TPZVec<REAL>  ff;
    
    /** @brief Source term */
	STATE fSf;
	
	/** @brief Elasticity modulus */
	REAL fE;
    
    //Lame's first parameter
    REAL flambda;
    
    //Lame's second parameter
    REAL fmu;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
	/** @brief constants poroelastic Biot*/
	REAL falpha; //parameter poroelastic Biot-Willis [dimensionless]
	REAL fSe; //or 1/M poroelastic storage coefficient at constant volume [dimensionless]
    
    //pressao de referencia
    REAL fpref;
    
    //factor permeability over viscosity
    REAL fkovervisc;
	
	/** @brief Problem dimension */
	int fDim;
	
	/** @brief Permeability of the rock and fluid viscosity*/
	REAL fk;
	REAL fvisc;
    
    //comprimento de referencia
    REAL fLref;
    
    //fluid diffusivity coeffcient
    REAL fCf;
    
    //constant auxiliary adimensionalization
    REAL fSaux;
	
	/** @brief Uses plain stress
     * @note \f$fPlaneStress = 1\f$ => Plain stress state
     * @note \f$fPlaneStress != 1\f$ => Plain Strain state
     */
	int fPlaneStress;
	
	/// timestep [s]
	REAL fTimeStep;
    
    REAL fTimeValue;
    
	int fmatId;
    
    int fmatIdSourceTerm;
    
    bool fReturnSolutionDimension;
	
	/** @brief State: one ou one+1 */
	enum EState { ELastState = 0, ECurrentState = 1 };
	static EState gState;
	
public:
	TPZPoroElasticMF2d();
	
	TPZPoroElasticMF2d(int matid, int dim);
	
	virtual ~TPZPoroElasticMF2d();
    
    TPZPoroElasticMF2d(const TPZPoroElasticMF2d &copy);
    
    TPZPoroElasticMF2d &operator=(const TPZPoroElasticMF2d &copy);
    
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZPoroElasticMF2d"; }
	
	int Dimension() const {return fDim;}
	
	virtual int NStateVariables() const {return 1;}
    
    void SetIdSourceTerm(int idsourceterm){
        
        fmatIdSourceTerm = idsourceterm;
    }
	
	void SetLastState(){ gState = ELastState; }
    
	void SetCurrentState(){ gState = ECurrentState; }
    
	/** @brief Parameters of rock and fluid: */
	void SetParameters(REAL perm, REAL visc)
	{
        fk = perm;
        fvisc = visc;
//        if(fReturnSolutionDimension == true){
//            fk = 1.;
//            fvisc = 1.;
//            fLref =Lref;
//            fpref = pref;
//            fkovervisc = perm/visc;
//        }
//		else {
//            fk = perm;
//            fvisc = visc;
//        }
    }
	
	/**
	 * @brief Set parameters of elastic material:
	 * @param E elasticity modulus
	 * @param nu poisson coefficient
     * @param alpha : constant poroelastic Biot [dimensionless]
	 * @param Se : Coeficiente poroelastico de armazenamento a volume constante [adimensional]
	 * @param fx forcing function \f$ -x = fx \f$
	 * @param fy forcing function \f$ -y = fy \f$
	 */
	void SetElasticityParameters(REAL E, REAL nu,  REAL alpha, REAL Se, REAL fx, REAL fy)
	{
		fE = E;
		fnu = nu;
        //flambda = (fE*fnu)/((1.+fnu)*(1.-2.*fnu));
        //fmu = 0.5*fE/(1.+fnu);
        falpha = alpha;
		fSe = Se;
		ff[0] = fx;
		ff[1] = fy;
        
//        if(fReturnSolutionDimension == true){
//            
//            //Incompressible fluid 
//            if(Se==0){
//                
//                fSaux = falpha*falpha/(flambda+2.*fmu);
//                REAL lambdaD = flambda*fSaux;
//                REAL muD = fmu*fSaux;
//                fE = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
//                fnu = 0.5*lambdaD/(lambdaD+muD);
//                fx *=fLref/fpref;
//                fy *=fLref/fpref;
//                
//            }
//            //Compressible fluid
//            else{
//                REAL lambdau = falpha*falpha/fSe +flambda;
//                fSaux  = falpha*falpha*(lambdau+2.*fmu)/((lambdau - flambda)*(flambda+2.*fmu));
//                REAL lambdaD = flambda*fSaux;
//                REAL muD = fmu*fSaux;
//                fE = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
//                fnu = 0.5*lambdaD/(lambdaD+muD);
//                fSe = fSe/fSaux;
//                fx *=fLref/fpref;
//                fy *=fLref/fpref;
//            }
//        }
	}
    
    void SetSourceTerm(REAL Sf){
        fSf = Sf;
    }
    
    //Apos resolver o problema adimensional posso retornar as solucaoes com dimensoes
    void SetReturnSolutionDimension(REAL pref, REAL Lref, REAL Seaux)
    {
        fReturnSolutionDimension = true;
        fSaux = Seaux;
        fpref = pref;
        fLref = Lref;
    }
	
	/**
	 * @brief Set plane problem
	 * planestress = 1 => Plain stress state
	 * planestress != 1 => Plain Strain state
	 */
	void SetfPlaneProblem(int planestress)
	{
		fPlaneStress = planestress;
	}
	
	/// Set the timestep
	void SetTimeStep(REAL delt)
	{
		fTimeStep = delt;
	}
    
    /// Set the timestep
	void SetTimeValue(REAL TimeValue)
	{
		fTimeValue = TimeValue;
//        if(fReturnSolutionDimension ==true){
//            fCf = fkovervisc/fSaux;
//            if(fLref==0.){
//                
//                std::cout << " Error.The fLref value has to be set before in SetParameters()\n";
//                DebugStop();
//            }
//            fTimeStep = TimeValue*fCf/(fLref*fLref);
//        }
	}
    
    REAL GetTimeValue(){
        
        return fTimeValue;
    }
    
    int MatId()
    {
        return fmatId;
    }
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
     * @brief Applies to Dirichlet boundary condition for the elasticity equation (mechanical problem)
     */
    void ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Neumann boundary condition for elasticity equation (mechanical problem)
     */
    void ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Mixed boundary condition for elasticity equation (mechanical problem)
     */
    void ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
     *@brief Dirichlet with free boundary in the y-direction in the equation of elasticity.
     */
    void ApplyDirichletFreeY_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ApplyDirichletFreeX_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
     *@brief Neumann with free boundary in the y-direction in the equation of elasticity.
     */
    void ApplyNeumannFreeX_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    void ApplyNeumannFreeY_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    

    
    
    /*
     * @brief Applies to Dirichlet boundary condition for mixed problem (pressure and flux)
     * In the mixed formulation, the contribution of the Dirichlet boundary condition for the pressure appears in the flow equation and not in the equation of pressure
     */
    void ApplyDirichlet_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Neumann boundary condition for mixed problem (pressure and flux)
     */
    void ApplyNeumann_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to Mixed boundary condition for mixed problem (pressure and flux)
     */
    void ApplyMixed_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /*
     * @brief Applies to condition of source term in the pressure's equation to mixed problem (pressure and flux)
     */
    //void ApplySourceTerm_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
    
	//public:
	virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                       REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
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
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
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
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
		DebugStop();
	}
	
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
};

#endif
