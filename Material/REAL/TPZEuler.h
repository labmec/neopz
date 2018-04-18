/**
 * \file
 * @brief Contains the TPZEuler class which implements a a linear scalar convection equation.
 */

#ifndef TPZEULER_H
#define TPZEULER_H

#include "TPZMaterial.h"
#include "eulerdif.h"
#include "pzfmatrix.h"

class TPZBndCond;
class TPZEuler;


/**
 * @ingroup material
 * @brief This class implements a linear scalar convection equation with modified 
 * SUPG difusion
 */
class TPZEuler : public TPZMaterial {
public:  

virtual int ClassId() const;

	/** @brief Copy constructor */
	TPZEuler(TPZEuler & copy);
	/** @brief Simple constructor */
    TPZEuler(int id, STATE deltat) ;
	
	/** 
	 * @brief Set the state of the Euler material
	 * @param state \f$ state = 0 \f$ -> L2 projection. \f$ state = 1 \f$ -> Euler time stepping
	 */
	void SetState(int state) {
		fState = state;
	}
	
	
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const;
	
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables()  ;
	
    /** @brief Return the number of components which form the flux function */
    virtual int NFluxes() {return 2;}
	
	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
    /** @brief Computes contribution to the stiffness matrix and right hand side at an integration point */
	virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) ;
    /** @brief Computes contribution to the right hand side at an integration point */
	virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix<STATE> &ef)
	{
        TPZMaterial::Contribute(data,weight,ef);
    }

    virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ef,TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }

	/** @} */
	
    /** @brief Print out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout);
	
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
	
    virtual int NSolutionVariables(int var);
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<STATE> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
    }
	
	
    virtual void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux) {}
	
    /** @brief To create another material of the same type */
    virtual TPZMaterial * NewMaterial();
	
    /** @brief Reads data of the material from a istream (file data) */
    virtual void SetData(std::istream &data);
	
private:    
	
	static TEulerDiffusivity gEul;
	STATE fDeltaT;
	int fState;
};

#endif
