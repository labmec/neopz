/**
 * \file
 * @brief Contains the TPZIncNavierStokesKEps class which implements an imcompressible Navier-Stokes formulation with modified KEpsilon turbulence model.
 */

#ifndef PZINCNSKEPS
#define PZINCNSKEPS

#include "pzmaterial.h"

#include <iostream>
#include <string>

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZBndCond;

/**
 * @ingroup material
 * @brief This class implements an imcompressible Navier-Stokes formulation with modified KEpsilon turbulence model.
 * @author Professor Paulo Vatavuk
 * @author Professor Philippe Devloo
 * @author Edimar Cesar Rylo
 * @author Luis Fernando
 * @author Roberto H. Heiderich 
 * @author Tiago Forti
 * @since June 29, 2005
 */
/** 
 * Variables are: {K, Eps, Pressure, Vx, Vy, Vz}.
 * This class is homework:
 */
class  TPZIncNavierStokesKEps : public TPZMaterial {
	
private:
	
    int fDimension;
    
    REAL fMU, fRHO, fCmu, fSigmaK, fSigmaEps, fCepsilon1, fCepsilon2;
    
    TPZVec<REAL> fBodyForce; //fc = {0,0,-fc}
    
    /** @brief Dot for matrices with same dimensions. \f$ Tr[A B]. \f$ No consistence test is made. */
    REAL Dot(TPZFMatrix<REAL> &A, TPZFMatrix<REAL> &B);
    
    /** @brief Dot of vector A with row BRow of matrix B. */
    REAL Dot(TPZVec<REAL> &A, TPZFMatrix<REAL> &B, int BRow);
	
    /** @brief Dot for vectors with same dimensions. No consistence test is made. */        
    REAL Dot(TPZVec<REAL> &A, TPZVec<REAL> &B);
	
public:
	
    enum VARIABLES {ENone = -1, EK = 0, EEpsilon, EPressure, EVx, EVy, EVz, EVvector};
	
    TPZIncNavierStokesKEps(int id, int dimension);
    
    void SetParameters(REAL MU, REAL RHO, REAL Cmu, REAL SigmaK, REAL SigmaEps, REAL Cepsilon1, REAL Cepsilon2, TPZVec<REAL> &BodyForce );
    
    void GetParameters(REAL &MU, REAL &RHO, REAL &Cmu, REAL &SigmaK, REAL &SigmaEps, REAL &Cepsilon1, REAL &Cepsilon2, TPZVec<REAL> &BodyForce );    
	
    virtual ~TPZIncNavierStokesKEps();
	
    virtual int Dimension();
	
    /** @brief Returns the number of state variables associated with the material*/
    virtual int NStateVariables();
	
    /** @brief Print out the data associated with the material*/
    virtual void Print(std::ostream &out = std::cout);
	
    /** @brief Returns the number of variables associated with the variable
     *  indexed by var. \n var is obtained by calling VariableIndex*/
    virtual int NSolutionVariables(int var);
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol,
						  TPZFMatrix<REAL> &axes, int var, TPZVec<REAL> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	
    /** @brief Computes contribution to the tangent matrix and residual at an integration point */
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<REAL> &ek,
							TPZFMatrix<REAL> &ef);
	
    /** @brief Computes contribution to the residual at an integration point */
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<REAL> &ef);
	
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at the integration point of a boundary */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ek,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc);
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at the integration point of a boundary */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc)
	{
        TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	
    /** 
	 * @brief Computes the error due to the difference between the interpolated flux and the flux computed \n
	 * based on the derivative of the solution
     */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol,
                        TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                        TPZVec<REAL> &uexact, TPZFMatrix<REAL> &duexact,
                        TPZVec<REAL> &val){
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
	
};


#endif

