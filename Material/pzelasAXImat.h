﻿/**
 * \file
 * @brief Contains the TPZElasticityAxiMaterial class which implements a two dimensional elastic material in plane stress or strain.
 */
#ifndef ELASAXIMATHPP
#define ELASAXIMATHPP

#include <iostream>

#include "pzmaterial.h"
#include "pzdiscgal.h"

#include <vector>
#include <math.h>

//const REAL Pi = 4.*atan(1);

/**
 * @ingroup material
 * @brief This class implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticityAxiMaterial : public TPZDiscontinuousGalerkin {
	
	public :
	
	TPZElasticityAxiMaterial();
	/** @brief Creates an elastic material with:\n
	 elasticity modulus  =   E \n
	 poisson coefficient  =   nu \n
	 forcing function -x =   fx \n
	 forcing function -y =   fy \n
	 fplainstress = 1 indicates use of plainstress
	 */
	TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy);
	
	//------------------- FEITO POR AGNALDO : 05/02/10 . Só para teste ----------
	/**
	 * @brief Creates an elastic material test with: \n
	 elasticity modulus  =   E \n
	 poisson coefficient  =   nu \n
	 forcing function -x =   fx \n
	 forcing function -y =   fy \n
	 symmetrical method = coefTheta (-1-> symmetric, 1->nonsymmetric)
	 penalty term = coefAlpha
	 */
	TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, REAL coefTheta,  REAL coefAlpha);
	
    /**
	 * @brief Copies the data of one TPZElasticityAxiMaterial object to
     another */
    TPZElasticityAxiMaterial(const TPZElasticityAxiMaterial &copy);
    
    /** @brief Creates a new material from the current object   ??*/
    virtual TPZAutoPointer<TPZMaterial> NewMaterial() { return new TPZElasticityAxiMaterial(*this);}
    
    /**Destructor*/
    virtual ~TPZElasticityAxiMaterial();
    
	
	//-----------------------------------------------------------------------------------
	
	/** @brief Set the origin of Revolution Axis (Z),
     the direction of Revolution Axis (Z), \n
     and the Radius vector (orthogonal with respect of Z axis)
	 */
	void SetOrigin(TPZManVector<REAL> &Orig, TPZManVector<REAL> &AxisZ, TPZManVector<REAL> &AxisR);
	
	REAL ComputeR(TPZVec<REAL> &x);
	
	void SetMohrCoulomb(double c, double phi)
	{
		f_c = c;
		f_phi = phi;
	}
    
    void SetThermalExpansionCoefficient(REAL alpha)
    {
        fAlpha = alpha;
    }
    
    void SetTemperature(REAL delt)
    {
        fDelTemperature = delt;
    }
	
	/** @brief Returns the dimension*/
	int Dimension() { return 2;}
	
	/** @brief Returns the number of state variables associated
     with the material
	 */
	virtual  int NStateVariables();
	
	/** @brief Prints the material data*/
	virtual void Print(std::ostream & out = std::cout);
	
	/** @brief Returns the material name*/
	std::string Name() { return "TPZElasticityAxiMaterial"; }
	
	/** @brief Return the number of components which form the flux function*/
	virtual short NumberOfFluxes(){return 3;}
	
	/** @brief Return the number of components which form the flux function*/
	virtual int NFluxes(){ return 3;}
	
	/**Cria as condicoes de contorno*/
	//virtual TPZBndCond *CreateBc(long num, int typ, TPZFMatrix &val1, TPZFMatrix &val2);
	
	/** @brief Calculates the element stiffness matrix*/
	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix &ek,TPZFMatrix &ef);
	
	/** @brief Applies the element boundary conditions*/
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
		return;
		PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	}
	
	/** @brief Returns the variable index associated with the name*/
	virtual int VariableIndex(const std::string &name);
	
	/** @brief Returns the number of variables associated with the variable
	 indexed by var. var is obtained by calling VariableIndex*/
	virtual int NSolutionVariables(int var);
	
	/** @brief returns the solution associated with the var index based on
	 * the finite element approximation*/
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
	/** @brief Compute the value of the flux function to be used
     by ZZ error estimator*/
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
	
	/** @brief Compute the error due to the difference between
     the interpolated flux and the flux computed based
     on the derivative of the solution*/
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric
	
	//virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &sol,TPZFMatrix &dsol, TPZFMatrix &axes, TPZVec<REAL> &flux,
	//		      TPZVec<REAL> &uexact,TPZFMatrix &duexact,TPZVec<REAL> &val){}
	
	/** @brief Returns the elasticity modulus E*/
	REAL E() {return fE;}
	
	/** @brief Returns the poison coefficient modulus E*/
	REAL Nu() {return fnu;}
	
	/** @brief Set PresStress Tensor*/
	void SetPreStress(REAL Sigxx, REAL Sigyy, REAL Sigxy);
    
    void SetTemperatureFunction(void (*func)(const TPZVec<REAL> &rz, REAL &temperature))
    {
        fTemperatureFunction = func;
    }
	
	virtual int ClassId() const;
	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid);
	
	TPZManVector<REAL> GetAxisR();
	TPZManVector<REAL> GetAxisZ();
	TPZManVector<REAL> GetOrigin();
	
	REAL fIntegral;
	
private:
	
	/** @brief Mohr Coulomb Plasticity Criteria Data */
	REAL f_phi;
    /** @brief Mohr Coulomb Plasticity Criteria Data */
	REAL f_c;
	
	/** @brief Elasticity modulus */
	REAL fE;
	
	/** @brief Poison coeficient */
	REAL fnu;
	
    /** @brief Thermal expansion coeficient */
    REAL fAlpha;
    
	/** @brief Forcing vector */
	REAL ff[3];
	
    /** @brief Temperature difference */
    REAL fDelTemperature;
    
	/** @brief G = E/2(1-nu)*/
	REAL fEover21PlusNu;
	
	/** @brief E/(1-nu)*/
	REAL fEover1MinNu2;
	
	/** @brief Direction of Surface*/
	TPZManVector<REAL> f_AxisR;
	
	/** @brief Revolution Axis*/
	TPZManVector<REAL> f_AxisZ;
	
	/** @brief Origin of AxisR and AxisZ*/
	TPZManVector<REAL> f_Origin;
	
	/** @brief symmetric*/ 
	REAL fSymmetric;
	
	/** @brief penalty term*/
	REAL fPenalty;
    
    /** @brief Function which defines the temperature */
    void (*fTemperatureFunction)(const TPZVec<REAL> &rz, REAL &temperature);
	
};

#endif
