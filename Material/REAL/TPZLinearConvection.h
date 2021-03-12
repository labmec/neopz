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
    
    virtual int ClassId() const override;

    /** @brief Copy constructor */
    TPZLinearConvection(TPZLinearConvection & copy);
	/** @brief Constructor for given convection */
    TPZLinearConvection(int id,TPZVec<STATE> &conv) ;

    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override;
	
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override ;

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
    virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
	
	
    virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
	
    virtual void Contribute(TPZMaterialData &data, REAL weight,
							TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	
	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ef,TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }

    /** @} */
    
	/** @brief Print out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout) override;
	
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name) override;
	
    virtual int NSolutionVariables(int var) override;
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<STATE> &Solout) override
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

        Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
    }
	
    /** @brief To create another material of the same type */
    virtual TPZMaterial * NewMaterial() override;
	
    /** @brief Reads data of the material from a istream (file data) */
    virtual void SetData(std::istream &data) override;
	
private:    
    STATE fConvect[2];

};

#endif //TPZLINEARCONVECTION_H
