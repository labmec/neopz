/**
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
 * @brief Implements a two dimensional elastic material in plane stress or strain
 */
class TPZElasticityAxiMaterial : public TPZDiscontinuousGalerkin {
	
	public :
	
	TPZElasticityAxiMaterial();
	/** 
	 * @brief Creates an elastic material with elasticity modulus (E), poisson coefficient (nu) and forcing functions (fx and fy).
	 * @param num material id
	 * @param E Elasticity modulus 
	 * @param nu \f$=\nu\f$ Poisson coefficient
	 * @param fx Forcing function \f$ -x = fx \f$
	 * @param fy Forcing function \f$ -y = fy \f$
	 * @note \f$ fplainstress = 1 \f$ indicates use of plainstress
	 */
	TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy);
	
	//------------------- FEITO POR AGNALDO : 05/02/10 . Só para teste ----------
	/**
	 * @brief Creates an elastic material test with elasticity modulus (E), poisson coefficient (nu) and forcing functions (fx, fy) and penalty term (coefAlpha).
	 * @param num material id
	 * @param E Elasticity modulus
	 * @param nu \f$=\nu\f$ Poisson coefficient
	 * @param fx Forcing function \f$ -x = fx \f$
	 * @param fy Forcing function \f$ -y = fy \f$
	 * @param coefTheta Symmetrical method (-1-> symmetric, 1->nonsymmetric)
	 * @param coefAlpha Penalty term
	 */
	TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, REAL coefTheta,  REAL coefAlpha);
	
    /** @brief Copies the data of one TPZElasticityAxiMaterial object to another */
    TPZElasticityAxiMaterial(const TPZElasticityAxiMaterial &copy);
    
    /** @brief Creates a new material from the current object   ??*/
    virtual TPZAutoPointer<TPZMaterial> NewMaterial() { return new TPZElasticityAxiMaterial(*this);}
    
    /** @brief Destructor */
    virtual ~TPZElasticityAxiMaterial();
    
	
	//-----------------------------------------------------------------------------------
	
	/** @brief Set the origin of Revolution Axis (\f$Z\f$), the direction of Revolution Axis (\f$Z\f$), \n
     * and the Radius vector (orthogonal with respect of \f$Z\f$ axis).
	 * @param Orig Origin of revolution axis
	 * @param AxisZ Direction vector of revolutin axis
	 * @param AxisR Radius vector
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
	
	/** @brief Returns the dimension */
	int Dimension() { return 2;}
	
	/** @brief Returns the number of state variables associated with the material */
	virtual  int NStateVariables();
			
	/** @brief Prints the material data */
	virtual void Print(std::ostream & out = std::cout);
	
	/** @brief Returns the material name */
	std::string Name() { return "TPZElasticityAxiMaterial"; }
	
	/** @brief Return the number of components which form the flux function */
	virtual short NumberOfFluxes(){return 3;}
	
	/** @brief Return the number of components which form the flux function */
	virtual int NFluxes(){ return 3;}
	
	/** @brief Calculates the element stiffness matrix */
	virtual void Contribute(TPZMaterialData &data, REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	/** @brief Applies the element boundary conditions */
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                     REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                       REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
		
		PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
        DebugStop();
	}
	
	/** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/** @brief Returns the number of variables associated with the variable indexed by var. \n var is obtained by calling VariableIndex */
	virtual int NSolutionVariables(int var);
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux);
	
	/**
	 * @brief Computes the error due to the difference between the interpolated flux and the flux computed based \n
     * on the derivative of the solution
	 */
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values);//Cedric
	
	/** @brief Returns the elasticity modulus E */
	REAL E() {return fE;}
	
	/** @brief Returns the poison coefficient modulus E */
	REAL Nu() {return fnu;}
	
	/** @brief Sets PresStress Tensor */
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
    
	/** @brief \f$ G = E/2(1-nu) \f$ */
	REAL fEover21PlusNu;
	
	/** @brief \f$ E/(1-nu) \f$ */
	REAL fEover1MinNu2;
	
	/** @brief Direction of Surface */
	TPZManVector<REAL> f_AxisR;
	
	/** @brief Revolution Axis */
	TPZManVector<REAL> f_AxisZ;
	
	/** @brief Origin of AxisR and AxisZ */
	TPZManVector<REAL> f_Origin;
	
	/** @brief Symmetric */ 
	REAL fSymmetric;
	
	/** @brief Penalty term */
	REAL fPenalty;
    
    /** @brief Function which defines the temperature */
    void (*fTemperatureFunction)(const TPZVec<REAL> &rz, REAL &temperature);
	
};

#endif
