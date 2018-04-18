/**
 * \file
 * @brief Contains the TPZLinearConvection class which implements a linear scalar convection equation.
 */

#ifndef TPZLINEARCONVECTION_H
#define TPZLINEARCONVECTION_H

#include "TPZMaterial.h"
#include "pzfmatrix.h"

class TPZBndCond;
/**
 * @ingroup material
 * @brief Implements a linear scalar convection equation with modified SUPG difusion
 */
class TPZLinearConvection : public TPZMaterial {
public:
    
    virtual int ClassId() const;

    /** @brief Copy constructor */
    TPZLinearConvection(TPZLinearConvection & copy);
	/** @brief Constructor for given convection */
    TPZLinearConvection(int id,TPZVec<STATE> &conv) ;

    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const;
	
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables()  ;
	
    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {return 2;}

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
    virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
	
	
    virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
    virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix<STATE> &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	
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
    STATE fConvect[2];

};

#endif //TPZLINEARCONVECTION_H
